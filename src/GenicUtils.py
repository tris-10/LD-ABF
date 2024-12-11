# Name:        GenicUtils.py
# Purpose: Utility functions for processing genetic data
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

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy as sp

CHROM_COL = "CHROM"
POS_COL  = "POS"
ID_COL = "ID"
REF_COL = "REF"
ALT_COL = "ALT"
OUTPUT_START_COLS = [CHROM_COL,POS_COL,ID_COL,REF_COL,ALT_COL]

def removeRowsWithNulls(x,y,DEBUG=False):
    """
    Combining x and y dropping any rows with NA values
    :param x: Series with set of values to be merged with y
    :param y: Series with set of values to be merged with x
    :param DEBUG:
    :return:
    """
    df = pd.concat([x, y], axis=1)
    dfCleaned = df.dropna(how='any',axis=0)
    if DEBUG:
        print("Inside removeRowsWithNulls started with %s rows and now %s" % (df.shape[0],dfCleaned.shape[0]))
    return dfCleaned[x.name],dfCleaned[y.name]

def renameDupColsIntPlusOffset(df, offset=0.01, DEBUG=True):
    """
    This was created to deal with multi-allelic where we have multiple allleles at the same integer position value,
    by using this we can rename the columns then later on cast them to integers again to revert to the original names
    :param df:
    :param offset:
    :param DEBUG:
    :return:
    """
    cols=pd.Series(df.columns)
    for dup in cols[cols.duplicated()].unique():
        cols[cols[cols == dup].index.values.tolist()] = [dup + i*offset if i != 0 else dup for i in range(sum(cols == dup))]
    # rename the columns with the cols list.
    return cols


def testCodeToMinorAllele():
    nrows=10
    ncols=5
    testGeno = pd.DataFrame(np.random.binomial(1, .6, ncols*nrows).reshape((nrows, ncols)))
    print("testGeno.mean()")
    print(testGeno.mean())
    testGenoMaf = codeToMinorAllele(testGeno)
    print("testGenoMaf.mean()")
    print(testGenoMaf.mean())

def codeToMinorAllele(testGeno,axis=0):
    """
    Flips any columns with minor alleles > 0.5, calling codeArrayToMinorAllele
    :param testGeno:
    :return:
    """
    return testGeno.apply(codeArrayToMinorAllele,axis=axis)

def codeArrayToMinorAllele(snp):
    snpCopy = snp.copy().astype(float)
    if np.nanmean(snpCopy.values)>0.5:
        snpCopy = (snpCopy - 1) ** 2
    return snpCopy


def getPolyMorphicSites(genomes, minMaf, DEBUG=False):
    """
    Get's the polymorphic sites for the input genomes, explicitly deals with mutli-allelic rows that
    may look like dubplicates. There are sometimes multi-allelic sites where one variant is polymorphic and the other
    isn't but, just searching and getting an index will grab both.
    This has been reset
    :param genomes: Expecting a vcf like data frame where the position is the ID (instead of a column)
                     CHROM  ID REF ALT  ...  B_58:02:01_28  B_35:01:01:02_23           ...  B_14:02:01:01_14
        POS                          ...
        31325201      6 NaN   C   A  ...            0.0               0.0              ...               0.0
        31325199      6 NaN   A   G  ...            1.0               1.0              ...               0.0
        31325193      6 NaN   T   C  ...            1.0               1.0              ...               0.0
        31325180      6 NaN   C   A  ...            0.0               0.0              ...               0.0
        31325164      6 NaN   C   T  ...            1.0               1.0              ...               0.0
    :param minMaf: minimum minor allele being allowed in, where this is considered out of the
        total number of haplotypes, even with nules
        ie genomesMafNoDupPos.sum(axis=1) / len(genomesMafNoDupPos.columns)
        This is the sum of alternate alleles / total number of samples (even if null).
        This was chosen since it's viewed as a more conservative approach.
    :return:
    """
    genomesCopy = genomes.copy()
    hapCols = list(set(genomesCopy.columns) - set(OUTPUT_START_COLS))
    genomesCopy = genomesCopy.T
    #Let's deal with multi-allelic positions adding a decimal place
    genomesCopy.columns = renameDupColsIntPlusOffset(genomesCopy)
    genomesCopy = genomesCopy.T
    # let's make sure everything is in minor allele sets
    if DEBUG:
        print(genomesCopy.head())
        incorrectlyCoded = genomesCopy[genomesCopy[hapCols].mean(1) > 0.5]
        print("Testing, before recoding to minor allele ")
        print(incorrectlyCoded.mean(1))
    genomesMafNoDupPos = codeToMinorAllele(genomesCopy[hapCols],axis=1)
    incorrectlyCoded = genomesMafNoDupPos[genomesMafNoDupPos.mean(1) > 0.5]
    if incorrectlyCoded.shape[0]>0:
        print("Error, recoding to minor allele didn't appear to work")
        print(incorrectlyCoded.mean(1))
        print("Unique:")
        for col in incorrectlyCoded.T.columns:
            print(col)
            print(incorrectlyCoded.T[col].unique())
        print(incorrectlyCoded.head())
        raise ValueError
    alleleFreqs = genomesMafNoDupPos.sum(axis=1) / len(genomesMafNoDupPos.columns)
    polymorphicSites = (alleleFreqs[alleleFreqs.between(minMaf, 1 - minMaf)]).index
    # Now let's swap out the old genomes and add the properly MAF coded ones in
    # With first step get the starting columns but, POS column is actual the index, may want to change this
    startCols = list(set(OUTPUT_START_COLS) - set([POS_COL]))
    genomesCopy=pd.concat([genomesCopy[startCols],genomesMafNoDupPos],axis=1)
    if genomes.shape[0]!=genomesCopy.shape[0] or genomes.shape[1]!=genomesCopy.shape[1]:
        print("Error, in attempting to recode alleles to MAF something went wrong oringinal")
        print("there are rows = %s vs new %s rows | then cols = %s vs new %s cols" %
              (genomes.shape[0],genomesCopy.shape[0],genomes.shape[1],genomesCopy.shape[1]))
        raise ValueError
    genomesCopy=genomesCopy.T
    polymorphicGeno = genomesCopy[polymorphicSites]
    polymorphicGeno.columns = list(map(int, polymorphicGeno.columns))
    polymorphicGeno=polymorphicGeno.T
    if DEBUG:
        print("starting genomes: ")
        print(genomes.head(5))
        print(genomes.tail(5))
        print(genomes.shape)
        print("Polymorphic sites: ")
        print(polymorphicGeno.head())
        print(polymorphicSites.shape)
    return alleleFreqs, polymorphicGeno




