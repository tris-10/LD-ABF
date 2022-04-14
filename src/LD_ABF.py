# Name:        LD_ABF.py
# Purpose: Testing a vcf.gz file looking for selection signal, either balancing or positive selction, by calculating
# the linkage disequilibrium approximate Bayes factor (LD-ABF). This test statistic provides a meassure of both the
# density and strength of LD around a target variant.
#
# Author:      Tristan J. Hayeck
#
# Created:     2020
# Copyright:   (c) Tristan J. Hayeck 2020
# Licence:     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
# CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.
# -------------------------------------------------------------------------------
import os, sys, argparse, time, datetime, pstats, cProfile, scipy, cython, shutil, glob
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy as sp
import statsmodels.api as sm
from collections import OrderedDict
from itertools import cycle
from itertools import chain
from os import path
import GenicIO
import GenicUtils

#Test statistics and columns being reported
LINKAGE_SUM_COL = "linkageSum"
SUM_CORR_COL = "sumCorr"
LOG_F_ABF_COL = "LD-ABF"
DENSITY_IN_WIN_COL = "densityInWin"
#VCF like columns we'll output
CHROM_COL = "CHROM"
POS_COL  = "POS"
ID_COL = "ID"
REF_COL = "REF"
ALT_COL = "ALT"
MAF_COL="MAF"
OUTPUT_START_COLS = [CHROM_COL,POS_COL,ID_COL,REF_COL,ALT_COL]

# -------------------------------------------------------#
# --------------- R^2 based statistics ------------------#
def calcPairwiseLinkageDisequilibrium(snpA, snpB, DEBUG=False):
    """
    Calculating the pairwise LD squared correlation, r^2, for phased samples

        D = p_AB - p_A*p_B

    r^2 = D^2/ (p_A*(1-p_A) * p_B*(1-p_B) )

    :param snpA:
    :param snpB:
    :param DEBUG:
    :return:
    """
    snpACleaned, snpBCleaned = GenicUtils.removeRowsWithNulls(snpA, snpB, DEBUG)
    if snpACleaned.shape[0] == 0 or snpBCleaned.shape[0]==0:
        return 0
    snpAcopy = snpACleaned.copy()
    snpBcopy = snpBCleaned.copy()
    flipA = snpAcopy.mean() > 0.5
    if flipA:
        snpAcopy = (snpAcopy - 1) ** 2
    flipB = snpBcopy.mean() > 0.5
    if flipB:
        snpBcopy = (snpBcopy - 1) ** 2
    p_A = snpAcopy.mean()
    p_B = snpBcopy.mean()
    p_AB = ((snpAcopy == 1) & (snpBcopy == 1)).mean()
    D = p_AB - p_A * p_B
    ldCorr = D ** 2 / (p_A * (1 - p_A) * p_B * (1 - p_B))
    if DEBUG:
        print("p_AB = %s | p_A = %s (flipped_A=%s) | p_B = %s (flipped_A=%s) | D = %s |LD corr = %s" %
              (p_AB, p_A, flipA, p_B, flipB, D, ldCorr))
    return ldCorr

def runCalcPairwiseLDoverWindow(testSnp, windowSize, polymorphicGeno, DEBUG=False):
    if type(testSnp) is pd.DataFrame:
        print("\n\nError, inside runScanOverWindow mutli allelic site for this posiiton!!!")
        print(testSnp.head())
        raise ValueError
    if testSnp.mean() == 0 or testSnp.mean() == 1:
        return 0
    testSnpPos = testSnp.name
    # need to deal multi-allelic case where we have repeat positions
    positions = pd.DataFrame({POS_COL: polymorphicGeno.columns})
    posInWin = positions[positions[POS_COL].between(testSnpPos - windowSize / 2, testSnpPos + windowSize / 2)]
    # Let's ignore our test SNP position
    posInWin = posInWin[posInWin[POS_COL] != testSnpPos]
    polyGenoInWin = polymorphicGeno[posInWin[POS_COL].unique()]
    perSnpLd = polyGenoInWin.apply(lambda x: calcPairwiseLinkageDisequilibrium(testSnp, x, DEBUG=DEBUG))
    summLD = (perSnpLd.fillna(0)).sum()
    if perSnpLd.empty:
        summLD = 0.0
    return summLD

def runPersonRsqRemoveNulls(x,y,DEBUG=False):
    xACleaned, yCleaned = GenicUtils.removeRowsWithNulls(x, y, DEBUG)
    if xACleaned.shape[0]<2 or yCleaned.shape[0]<2:
        return 0
    return (stats.pearsonr(yCleaned, xACleaned)[0] ** 2)


def runCalcRsquaredOverWindow(testSnp, windowSize, polymorphicGeno, DEBUG=False):
    if type(testSnp) is pd.DataFrame:
        print("\n\nError, inside runScanOverWindow mutli allelic site for this posiiton!!!")
        print(testSnp.head())
        raise ValueError
    if testSnp.mean() == 0 or testSnp.mean() == 1:
        return 0
    testSnpPos = testSnp.name
    # need to deal multi-allelic case where we have repeat positions
    positions = pd.DataFrame({POS_COL: polymorphicGeno.columns})
    posInWin = positions[positions[POS_COL].between(testSnpPos - windowSize / 2, testSnpPos + windowSize / 2)]
    posInWin = posInWin[posInWin[POS_COL] != testSnpPos]
    polyGenoInWin = polymorphicGeno[posInWin[POS_COL].unique()]
    # Let's ignore our test SNP position
    if DEBUG:
        print("runCalcRsquaredOverWindow")
        print("polyGenoInWin")
        print(polyGenoInWin.head())
        print(polyGenoInWin.shape)
        print(polyGenoInWin.dtypes)
        print("testing")
    if polyGenoInWin.shape[1] < 1:
        return 0
    perSnpCorr = polyGenoInWin.apply(lambda x: runPersonRsqRemoveNulls(testSnp, x) )
    if DEBUG:
        print("perSnpCorr: ")
        print(perSnpCorr)
    sumRsq = perSnpCorr.sum()
    return sumRsq

def runCalcDensityOverWindow(testSnp, windowSize, polymorphicGeno, DEBUG=False):
    """
    Get the number of polymorphic SNPs in the window around the test snp
    :param testSnp:
    :param windowSize:
    :param polymorphicGeno:
    :param DEBUG:
    :return:
    """
    if type(testSnp) is pd.DataFrame:
        print("\n\nError, inside runScanOverWindow mutli allelic site for this posiiton!!!")
        print(testSnp.head())
        raise ValueError
    if testSnp.mean() == 0 or testSnp.mean() == 1:
        return -1
    testSnpPos = testSnp.name
    # need to deal multi-allelic case where we have repeat positions
    positions = pd.DataFrame({POS_COL: polymorphicGeno.columns})
    posInWin = positions[positions[POS_COL].between(testSnpPos - windowSize / 2, testSnpPos + windowSize / 2)]
    # Let's ignore our test SNP position
    posInWin = posInWin[posInWin[POS_COL] != testSnpPos]
    polyGenoInWin = polymorphicGeno[posInWin[POS_COL].unique()]
    if DEBUG:
        print("runCalcDensityOverWindow")
        print("testing % " % testSnpPos)
        print("polyGenoInWin")
        print(polyGenoInWin.head())
        print(polyGenoInWin.shape)
        print(polyGenoInWin.dtypes)
    return polyGenoInWin.shape[1]

# --------------- end of R^2 based statistics ------------------#


#--------------------------------------------------------------#
#------------ log-F Prior Approximate Bayes Factor ------------#
def runCalcLogFApproxBayesFacOverWindow(testSnp, windowSize, polymorphicGeno, DEBUG=False):
    if type(testSnp) is pd.DataFrame:
        print("\n\nError, inside runCalcLogFPriorABFOverWindow, mutli allelic site for this posiiton!!!")
        print(testSnp.head())
        raise ValueError
    if testSnp.mean() == 0 or testSnp.mean() == 1:
        return 0
    testSnpPos = testSnp.name
    # need to deal multi-allelic case where we have repeat positions
    positions = pd.DataFrame({POS_COL: polymorphicGeno.columns})
    posInWin = positions[positions[POS_COL].between(testSnpPos - windowSize / 2, testSnpPos + windowSize / 2)]
    # Let's ignore our test SNP position
    posInWin = posInWin[posInWin[POS_COL] != testSnpPos]
    polyGenoInWin = polymorphicGeno[posInWin[POS_COL].unique()]
    if DEBUG:
        print("runCalcLogFApproxBayesFacOverWindow")
    if polyGenoInWin.shape[0]==0 or polyGenoInWin.shape[1]==0 :
        return 0
    logLogfABF = polyGenoInWin.apply(lambda y: calcPairwiseLogfApproxBayesFactor(testSnp, y, DEBUG=DEBUG))
    if DEBUG:
        print(logLogfABF)
    sumLogLogfABF = (logLogfABF.fillna(0)).sum()
    return sumLogLogfABF

def testCalcPairwiseLogfApproxBayesFactor(priorOnNull):
    # Let's look at an full seperation example, should be high
    halfSamples = 5000
    data = pd.DataFrame({"y": np.concatenate((np.ones(halfSamples), np.zeros(halfSamples))),
                         "x": np.concatenate((np.ones(halfSamples), np.zeros(halfSamples)))})
    logLogfABF_fullSep = calcPairwiseLogfApproxBayesFactor(data['x'], data['y'],priorOnNull = priorOnNull, DEBUG=True)
    a=500
    b=200
    c=30
    d=0
    data = pd.DataFrame({"y": np.concatenate((np.zeros(a),np.ones(b),np.zeros(c),np.ones(d))),
                         "x": np.concatenate((np.zeros(a),np.zeros(b),np.ones(c),np.ones(d)))})
    logLogfABF_sparse = calcPairwiseLogfApproxBayesFactor(data['x'], data['y'], priorOnNull = priorOnNull,DEBUG=True)
    data = pd.DataFrame({"y": stats.binom.rvs(1, .25, 0, 200),
                         "x": stats.binom.rvs(1, .25, 0, 200)})
    logitLooRand_randomNoPrior = calcPairwiseLogfApproxBayesFactor(data['x'], data['y'], priorOnNull = False, DEBUG=True)
    logitLooRand_randomWithPrior = calcPairwiseLogfApproxBayesFactor(data['x'], data['y'], priorOnNull=True, DEBUG=True)
    print("Logf-ABF fullSep = %.2f ,sparse = %.2f , random = %.2f (no prior), random = %.2f (with prior)" %
              (logLogfABF_fullSep,logLogfABF_sparse,logitLooRand_randomNoPrior,logitLooRand_randomWithPrior))

def calcPairwiseLogfApproxBayesFactor(x, y,m=1, priorOnNull = True, DEBUG=False):
    """
    Calculating the approximate bayes factor assume log-F priors in a logistic regression setting.

    Let's say we start with some vector of values for x and y = {0,1} repressented by *
        inter x y
            1 ∗ ∗
            1 ∗ ∗
            1 ∗ ∗
            ...
            1 ∗ ∗
            1 ∗ ∗
    Then augment the data and weight it
        inter x y weight
            1 ∗ ∗ 1
            1 ∗ ∗ 1
            1 ∗ ∗ 1
            ...
            1 ∗ ∗ 1
            1 ∗ ∗ 1
            0 1 0 m/2
            0 1 0 m/2
    Then by fitting this model we get an estimate usimg log-F priors.
    We're running also adding in intercept priors, ie also add cells:
            0 0 0 m/2
            0 0 1 m/2
    :param x: the input covariate/test SNP
    :param y: the output/neighboring SNP
    :param m: hyper paramter for a symetric log-F distribution
    :param priorOnNull: include prior on the null model, in practice doesn't appear to make a big difference
    :param DEBUG:
    :return: the log of the log-F prior Approximate Bayes Factor
    """
    weightCol="weight"
    interceptCol= "intercept"
    xCol = x.name
    yCol = y.name
    xCleaned, yCleaned = GenicUtils.removeRowsWithNulls(x, y,DEBUG)
    #check if all values are null
    if xCleaned.shape[0]== 0 or yCleaned.shape[0]==0 or xCleaned.sum()==0 or yCleaned.sum()==0:
        print("Warning, setting ABF to 0, appears extensive nulls got through at: %s vs %s" % (xCol,yCol))
        print("xCleaned.shape[0]= %s or yCleaned.shape[0]= %s or xCleaned.sum()=%s or yCleaned.sum()=%s" %
              (xCleaned.shape[0], yCleaned.shape[0], xCleaned.sum(), yCleaned.sum()))
        return 0
    if xCleaned.sum()==xCleaned.shape[0] or yCleaned.sum()==xCleaned.shape[0]:
        print("Warning, setting ABF to 0, either x or y is all 1's likely related to null: %s vs %s" % (xCol,yCol))
        print("xCleaned.shape[0]==xCleaned.shape[0] %s = %s or yCleaned.shape[0]==yCleaned.shape[0] %s = %s " %
              (xCleaned.sum(), xCleaned.shape[0], yCleaned.sum(), yCleaned.shape[0]))
        return 0
    xFloat = xCleaned.copy().astype(float)
    yFloat = yCleaned.copy().astype(float)
    augmentX = pd.DataFrame(xFloat)
    augmentY = pd.DataFrame(yFloat)
    augmentX[interceptCol] = 1
    augmentX[weightCol]=1
    if augmentY[yCol].isna().sum() > 0:
        print("Error, in calcPairwiseLogfApproxBayesFactor, some y-values are null (shouldn't be possible here): ")
        print(augmentY['y'])
        raise ValueError
    #May want to use this in other settings, for now trying augmented null too assuming log-F prior on intercept only model
    if priorOnNull:
        # adding in augmentation on the null model as well
        newX = pd.DataFrame({xCol: [0, 0], interceptCol: [0, 0]})
        newX[weightCol] = m / 2.0
        augmentX = augmentX.append(newX, ignore_index=True)
        newY = pd.DataFrame({yCol: [0, 1]})
        augmentY = augmentY.append(newY, ignore_index=True)
        try:
            interceptOnlyFit = sm.Logit(augmentY, augmentX[interceptCol],
                                    freq_weights=augmentX[weightCol],
                                        disp=np.int(DEBUG)).fit()
        except Exception as e:
            print("augmentX")
            print(augmentX)
            print("augmentY")
            print(augmentY)
            print("xFloat.sum()")
            print(xFloat.sum())
            print("yFloat.sum()")
            print(yFloat.sum())
            print("Exception occurred in calcPairwiseLogfApproxBayesFactor when fitting intercept model: ")
            print(e)
            raise ValueError
        beta0inter = interceptOnlyFit.params[interceptCol]
        yPred = np.exp(beta0inter) / (1 + np.exp(beta0inter))
        logLikePerObsInter = yFloat * np.log(yPred) + (1 - yFloat) * np.log(1 - yPred)
        logLikelihoodInter = logLikePerObsInter.sum() + m / 2.0 * beta0inter - m * np.log(1 + np.exp(beta0inter))
    else:
        interceptOnlyFit = sm.Logit(yFloat, augmentX[interceptCol],disp=np.int(DEBUG)).fit()
        logLikelihoodInter=interceptOnlyFit.llf
        #now add in augmentation even though not using for this model
        newX = pd.DataFrame({xCol: [0, 0], interceptCol: [0, 0]})
        newX[weightCol] = m / 2.0
        augmentX = augmentX.append(newX, ignore_index=True)
        newY = pd.DataFrame({yCol: [0, 1]})
        augmentY = augmentY.append(newY, ignore_index=True)
    newX = pd.DataFrame({xCol: [1, 1], interceptCol: [0, 0]})
    newX[weightCol]=m/2.0
    augmentX= augmentX.append(newX,ignore_index=True)
    newY = pd.DataFrame({yCol:[0,1]})
    augmentY=augmentY.append(newY, ignore_index=True)
    try:
        logitFit = sm.Logit(augmentY.astype(float), augmentX[[interceptCol, xCol]],
                        freq_weights=augmentX[weightCol],disp=np.int(DEBUG)).fit()
    except Exception as e:
        print("augmentX")
        print(augmentX)
        print("augmentY")
        print(augmentY)
        print("xFloat.sum()")
        print(xFloat.sum())
        print("yFloat.sum()")
        print(yFloat.sum())
        print("Exception occurred in calcPairwiseLogfApproxBayesFactor when fitting full model: ")
        print(e)
        raise ValueError
    try:
        beta0 = logitFit.params[interceptCol]
        beta1 = logitFit.params[xCol]
        yPred=np.exp(beta1*xFloat + beta0) / (1+np.exp(beta1*xFloat + beta0))
        logLikePerObs = (yFloat*np.log(yPred) + (1-yFloat)*np.log(1-yPred)).sum()
        logLikelihood = logLikePerObs + m/2.0*beta0 - m*np.log(1+np.exp(beta0)) \
                        + m/2.0*beta1 - m*np.log(1+np.exp(beta1))
        logLogfABF = logLikelihood - logLikelihoodInter
    except Exception as e:
        print("Exception occurred:")
        print(e)
        print("beta0=%s beta1=%s" % (beta0, beta1))
        print("xFloat")
        print(xFloat)
        print("yFloat")
        print(yFloat)
        raise ValueError
    if DEBUG:
        print("\n\n_______________________________________")
        print("For %s with %s" % (x.name,y.name))
        a= ((xCleaned==0) & (yCleaned==0)).sum()
        b = ((xCleaned== 0) & (yCleaned ==1)).sum()
        c = ((xCleaned == 1) & (yCleaned == 0)).sum()
        d = ((xCleaned == 1) & (yCleaned == 1)).sum()
        print("| %0.2f\t%0.2f |\n| %0.2f\t%0.2f |" % (a, b, c, d))
        print("Total of %s" % (a+b+c+d))
        print("\n_______________________________________")
        print("Intercept only model:")
        print(interceptOnlyFit.summary())
        print("_______________________________________\n")
        print("Logit fit:")
        print(print(logitFit.summary()))
        print("_______________________________________\n\n")
        print("Unique y-pred: ")
        print(yPred.unique())
        diffLogLikeInter = logLikelihoodInter - interceptOnlyFit.llf
        print("Difference between null LL calc:  %.2f - %.2f = %.2f" %
              (logLikelihoodInter, interceptOnlyFit.llf, diffLogLikeInter))
        diffLogLike = logLikelihood - logitFit.llf
        print("Difference between log-likelihood calc:  %.2f - %.2f = %.2f  (without adjustment %.2f)" %
              (logLikelihood, logitFit.llf, diffLogLike,logLikePerObs))
        print("_______________________________________\n")
        print("\n_______________________________________")
        print("The approximate Bayes factor is: %.2f - %.2f = %.2f " %
              (logLikelihood, logLikelihoodInter,logLogfABF))
        pseudoR2 = 1.0 - (logLikelihood / logLikelihoodInter)
        print("Pseudo r^2 = %.4f" % pseudoR2)
        print("Pseudo r^2 = %.4f (from statsmodels output)" % logitFit.prsquared)
        print("_______________________________________\n\n")
    return logLogfABF

#------------ end of log-F Prior Approximate Bayes Factor -----#
#--------------------------------------------------------------#

def readInGenomes(windowSize, chromosome, startPos, endPos, singlePos, inputAreHaps, inputFile, outputFilePrefix, DEBUG=False):
    runningOverRegion = chromosome is not None and startPos is not None and endPos is not None
    if inputAreHaps:
        genomes = GenicIO.readHaplotypeTable(inputFile, DEBUG=DEBUG)
    else:
        # let's restrict to looking over the region plus the window (half on each end)
        if runningOverRegion:
            startPos = np.int(startPos)
            endPos = np.int(endPos)
            startPosWithWindow = np.int(np.max([0, startPos - np.ceil(windowSize / 2)]))
            endPosWithWindow = np.int(endPos + np.ceil(windowSize / 2))
        if inputFile.endswith('.gz'):
            if chromosome is None or startPos is None or endPos is None:
                print(
                    "Error, using a compressed vcf file %s so chromosome and start and end positions are required")
                print("file: %s\n region: %s:%s-%s" % (inputFile, chromosome, startPos, endPos))
                raise ValueError
            else:
                genomes = GenicIO.readCompressedVcfOverIntervalAndGetHaps(chromosome,
                                                                          startPosWithWindow, endPosWithWindow,
                                                                          inputFile, DEBUG=DEBUG)
        else:
            genomes = GenicIO.readHaplotypesFromVcf(inputFile)

    if runningOverRegion:
        genomes[POS_COL] = genomes.index
        genomes[POS_COL] = genomes[POS_COL].astype(np.int)
        genomes=genomes[genomes[POS_COL].between(startPosWithWindow, endPosWithWindow)]
        genomes = ((genomes.set_index(POS_COL)))
    return genomes

def ouputWhentUnableToRun(start,  outputFilePrefix):
    """
    Sometimes we won't be able to run the analysis, say if there's no polymorphic positions. This will still output
    the expected output files just make them blank. Also, a warning will be output.
    :param start:
    :param outputFilePrefix:
    :return:
    """
    scoresDf = pd.DataFrame(
        columns=[CHROM_COL, POS_COL, ID_COL, REF_COL, ALT_COL,
                 LOG_F_ABF_COL, LINKAGE_SUM_COL, SUM_CORR_COL,
                 DENSITY_IN_WIN_COL,MAF_COL])
    scoresDf.to_csv(outputFieP8efix + ".tsv", sep='\t', header=True, index=False)

    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print(
        "Warning, unable likely because no polymorphisms over regions, output should be blank headers for tables")
    print("The duration (seconds) run on single SNP was %s: " % elapsed.total_seconds())

# runFullyBayes,
def runScan(windowSize, chromosome, startPos, endPos, minMaf, singlePos, inputAreHaps, inputFile, outputFilePrefix,
             DEBUG=False):
    start = datetime.datetime.now()
    ableToRun = True
    # pre-specified chromosome and start and end positions to scan over, otherwise doing just 1 position if singlePos
    # is set or running on the entire input file (probably just for testing in that setting)
    runningOverRegion = chromosome is not None and startPos is not None and endPos is not None

    if path.exists(outputFilePrefix + ".tsv"):
        print("Warning!!! File already exists, for now testing if there is an issue in simulation\n"
              "when trying to run jobs in parallel there maybe an issue with multiple theano calls")
    else:
        genomes=readInGenomes(windowSize, chromosome, startPos, endPos, singlePos, inputAreHaps, inputFile, outputFilePrefix, DEBUG)
        # putting in double check removal of non-polymorphic sites for nultiallic regions
        alleleFreqs, polymorphicGenoVcf = GenicUtils.getPolyMorphicSites(genomes, minMaf, DEBUG)
        #Let's get just the genotypes but keep the
        startCols = set(OUTPUT_START_COLS) - set([POS_COL])
        hapCols = set(polymorphicGenoVcf.columns) - startCols
        polymorphicGeno= (polymorphicGenoVcf[hapCols]).T
        if polymorphicGeno.shape[1] == 0:
            print("Error, there appear to be no polymorphic sites.")
            ableToRun = False

        if singlePos is None:
            if polymorphicGeno.shape[1] == 0:
                print("Error, unable to run because there are no polymorphic regions, double check input. ")
                ableToRun = False
                scoresDf = pd.DataFrame(columns=[POS_COL, LOG_F_ABF_COL])
            if ableToRun:
                if DEBUG:
                    print("Running scan across multiple positions")

                summLD = polymorphicGeno.apply(lambda x:
                                               runCalcPairwiseLDoverWindow(testSnp=x, windowSize=windowSize,
                                                                           polymorphicGeno=polymorphicGeno,
                                                                           DEBUG=DEBUG))

                sumCorr = polymorphicGeno.apply(lambda x:
                                                runCalcRsquaredOverWindow(testSnp=x, windowSize=windowSize,
                                                     polymorphicGeno=polymorphicGeno, DEBUG=False))

                densityInWin = polymorphicGeno.apply(lambda x:
                                        runCalcDensityOverWindow(testSnp=x, windowSize=windowSize,
                                                                  polymorphicGeno=polymorphicGeno, DEBUG=False))


                logF_ABF = polymorphicGeno.apply(lambda x:
                                      runCalcLogFApproxBayesFacOverWindow(testSnp=x, windowSize=windowSize,
                                                                polymorphicGeno=polymorphicGeno, DEBUG=DEBUG))



                scoresDf = pd.DataFrame({LOG_F_ABF_COL: logF_ABF, LINKAGE_SUM_COL: summLD, SUM_CORR_COL: sumCorr,
                                             DENSITY_IN_WIN_COL: densityInWin})
            scoresDf.index.name = POS_COL

        else:
            if not np.int(singlePos) in polymorphicGeno.columns:
                print(
                    "\n\nError! Position %s being tested is not input or is not polymorphic (in MAF restriction), double check sample." % singlePos)
                print(
                    "The position could also be an indel or missing allele (like for IMGT samples). Excluding these is to be expected")
                if np.int(singlePos) in genomes.T.columns:
                    print("Unadjusted frequency: ")
                    print(genomes.T[singlePos][hapCols].mean())
                    print("Adjusted frequency (handling nulls and minor allele coded) ")
                    freqAtPos = alleleFreqs[np.int(singlePos)]
                    print("With frequency in sample %s" % freqAtPos)
                else:
                    print("not found, no allele frequency for %s" % singlePos)
                ableToRun = False
            if ableToRun:
                if DEBUG:
                    print("Running on a single position.")

                testSnp = polymorphicGeno[np.int(singlePos)]
                # if multi-allelic at input SNP, run over each one
                if type(testSnp) is pd.DataFrame:
                    print("Warning, multi-allelic site...")
                    summLD = testSnp.apply(lambda x:
                                           runCalcPairwiseLDoverWindow(testSnp=x, windowSize=windowSize,
                                                                       polymorphicGeno=polymorphicGeno, DEBUG=DEBUG))
                    sumCorr = testSnp.apply(lambda x:
                                            runCalcRsquaredOverWindow(testSnp=x, windowSize=windowSize,
                                                                      polymorphicGeno=polymorphicGeno,
                                                                      DEBUG=DEBUG))
                    densityInWin= testSnp.apply(lambda x:
                                            runCalcDensityOverWindow(testSnp=x, windowSize=windowSize,
                                                                      polymorphicGeno=polymorphicGeno,
                                                                      DEBUG=DEBUG))

                    logF_ABF = polymorphicGeno.apply(lambda x:
                                                     runCalcLogFApproxBayesFacOverWindow(testSnp=x,
                                                                                         windowSize=windowSize,
                                                                                         polymorphicGeno=polymorphicGeno,
                                                                                         DEBUG=DEBUG))
                    scoresDf = pd.DataFrame({LOG_F_ABF_COL: logF_ABF, LINKAGE_SUM_COL: summLD,
                                                  SUM_CORR_COL: sumCorr, DENSITY_IN_WIN_COL: densityInWin})
                    scoresDf.index.name = POS_COL
                else:
                    summLD = runCalcPairwiseLDoverWindow(testSnp=polymorphicGeno[np.int(singlePos)],
                                                         windowSize=windowSize,
                                                         polymorphicGeno=polymorphicGeno, DEBUG=DEBUG)

                    sumCorr = runCalcRsquaredOverWindow(testSnp=polymorphicGeno[np.int(singlePos)],
                                                        windowSize=windowSize,
                                                        polymorphicGeno=polymorphicGeno, DEBUG=DEBUG)
                    densityInWin=runCalcDensityOverWindow(testSnp=polymorphicGeno[np.int(singlePos)],
                                                        windowSize=windowSize,
                                                        polymorphicGeno=polymorphicGeno, DEBUG=DEBUG)
                    logF_ABF =runCalcLogFApproxBayesFacOverWindow(testSnp=polymorphicGeno[np.int(singlePos)],
                                                                    windowSize=windowSize,
                                                                  polymorphicGeno=polymorphicGeno, DEBUG=DEBUG)
                    if DEBUG:
                        print("Only ran ABF")
                    scoresDf = pd.DataFrame({ LOG_F_ABF_COL: logF_ABF,  LINKAGE_SUM_COL: summLD,
                                             SUM_CORR_COL: sumCorr, DENSITY_IN_WIN_COL: densityInWin},
                                            index=[np.int(singlePos)])
                    scoresDf.index.name = POS_COL

        if not ableToRun:
            ouputWhentUnableToRun(start, outputFilePrefix)
            return -1

        if runningOverRegion:
            # make sure to only include scores over the region
            if DEBUG:
                print("\n\nRestricting to test region (removing window boarders)")
                print(scoresDf.head())
            scoresDf[POS_COL] = scoresDf.index
            scoresDf[POS_COL] = scoresDf[POS_COL].astype(np.int)
            if DEBUG:
                print("Cleaning up the scores table so we're over desired region, currently %s to %s" %
                  (scoresDf[POS_COL].min(), scoresDf[POS_COL].max()))
                print(scoresDf.shape)
            scoresDf = scoresDf[scoresDf[POS_COL].between(startPos, endPos)]
            if DEBUG:
                print("then restricting min max positions with scores %s to %s for input range %s to %s" %
                  (scoresDf[POS_COL].min(), scoresDf[POS_COL].max(), startPos, endPos))
            scoresDf.set_index(POS_COL, inplace=True)
            scoresDf.index.name = POS_COL
            if DEBUG:
                print(scoresDf.shape)
                print(scoresDf.head())

        ################################
        if ableToRun:
            #Now let's merge back with vcf header information we want
            #There are better ways of doing this but would involve more of an overhaul
            #since there can be duplicate positions lets add an offset
            polymorphicGenoVcf[POS_COL] = polymorphicGenoVcf.index
            transPolymorphicGenoVcf = polymorphicGenoVcf.T
            transPolymorphicGenoVcf.columns = GenicUtils.renameDupColsIntPlusOffset(transPolymorphicGenoVcf)
            polymorphicGenoVcf = transPolymorphicGenoVcf.T
            #now for scores table so we can merge it below
            transScoresDf = scoresDf.T
            transScoresDf.columns = GenicUtils.renameDupColsIntPlusOffset(transScoresDf)
            scoresDf = transScoresDf.T

            if singlePos is None:
                scoresDf[MAF_COL]=polymorphicGenoVcf[hapCols].mean(axis=1)
                outputTable = pd.concat([polymorphicGenoVcf[OUTPUT_START_COLS], scoresDf], axis=1)
            else:
                polymorphicGenoVcfAtPos = polymorphicGenoVcf[polymorphicGenoVcf[POS_COL]== np.int(singlePos)]
                scoresDf[MAF_COL] = polymorphicGenoVcfAtPos[hapCols].mean(axis=1)
                outputTable = pd.concat([polymorphicGenoVcfAtPos[OUTPUT_START_COLS], scoresDf], axis=1)
            outputTable[LOG_F_ABF_COL]=outputTable[LOG_F_ABF_COL]/windowSize
            outputCols = OUTPUT_START_COLS+scoresDf.columns.values.tolist()
            #going back and cleaning out scores outsied of the region again
            if runningOverRegion:
                print("Trimming scores or nulls outside region (ie needed to look +/- window/2): %s to %s" % (startPos,endPos))
                print(outputTable.shape)
                outputTable=outputTable[outputTable[POS_COL].between(startPos, endPos)]
                print(outputTable.shape)
            outputTable[outputCols].to_csv(outputFilePrefix + ".tsv", sep='\t', header=True, index=False)
        else:
            print("Unable to run, probably because no polymorphic values in input to test.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputFile", type=str,default="input.vcf",
                        help="Input file with variants to be analyzed, in *.vcf.gz format. "
                             "Unless inputAreHaps set to true this should be compressed with corresponding "
                             "tabix *.vcf.gz.csi file (http://www.htslib.org/doc/tabix.html).")
    parser.add_argument("-w", "--windowSize", type=int, default=1000,
                        help="Window size for scan statistic.")
    parser.add_argument("-m", "--minMaf", type=float, default=1 / 1000000.,
                        help="Minimum minor allele frequency to consider a site polymorphic")
    parser.add_argument("-s", "--singlePos", default=None,
                        help="Running test on single position otherwise None to do full scan or scan over range")
    # scan mode terms
    parser.add_argument("--chromosome", dest='chromosome', default=None,
                        help="If running a scan over certain positions or just looking over one positions, "
                             "otherwise set to None as default")
    parser.add_argument("--start-pos", dest='startPos', default=None,type=int,
                        help="If running a scan over certain positions this is the starting position")
    parser.add_argument("--end-pos", dest='endPos', default=None,type=int,
                        help="If running a scan over certain positions this is the end position")
    parser.add_argument('--inputAreHaps', dest='inputAreHaps', action='store_true',
                        help="Flag indicating if the input file contains columns with only haplotypes instead of "
                             "biallelic phased values (ie 0 or 1 instead of {0|0,1|0,0|1,1|1}))")
    parser.add_argument('--not-inputAreHaps', dest='inputAreHaps', action='store_false')
    parser.add_argument('--debug', dest='DEBUG', action='store_true')
    parser.set_defaults(inputAreHaps=False)
    parser.set_defaults(DEBUG=False)
    parser.add_argument("-o", "--outputFilePrefix", type=str, help="Output file location prefix to be written out")
    start = datetime.datetime.now()
    args = parser.parse_args()
    print("Starting processing %s" % start)
    print(args)
    if args.DEBUG:
        print("Warning, running in debug with more details printed out...")
    runScan(args.windowSize, args.chromosome, args.startPos, args.endPos, args.minMaf, args.singlePos,
                args.inputAreHaps, args.inputFile, args.outputFilePrefix, DEBUG=args.DEBUG)
    done = datetime.datetime.now()
    elapsed = done - start
    duration = ':'.join(str(elapsed).split(':')[1:])
    print("The duration was %s: " % duration)
    print("The duration (seconds) was %s: " % elapsed.total_seconds())
    print("Finished processing %s" % done)

if __name__ == '__main__':
    main()


