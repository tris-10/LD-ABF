# Name:        GenicIO.py
# Purpose: File for processing genetic files like VCFs
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
import os, sys, argparse, time, datetime, pstats,  cProfile,  scipy, cython,pyximport
import numpy as np
import pandas as pd
import scipy.stats as stats
from itertools import cycle
from itertools import chain
import pysam
import gzip
#VCF file columns
CHROM_COL = "CHROM"
POS_COL  = "POS"
ID_COL = "ID"
REF_COL = "REF"
ALT_COL = "ALT"
QUAL_COL = "QUAL"
FILTER_COL = "FILTER"
INFO_COL = "INFO"
FORMAT_COL = "FORMAT"
VCF_START_COLS = [CHROM_COL,POS_COL,ID_COL,REF_COL,ALT_COL,QUAL_COL,FILTER_COL,INFO_COL,FORMAT_COL]
MOM_HAP_SUFFIX = "_Mom"
DAD_HAP_SUFFIX = "_Dad"
START_COL = "start"
END_COL = "end"
LABEL_COL="label"
POS_MINUS_1_COL = "posMinus1"

NULL_VALUE="."

def getFileLineCount(fileName, DEBUG=False):
    count = 0
    with open(fileName, 'r') as f:
        for line in f:
            count += 1
    if DEBUG:
        print("%s has %s lines" % (fileName,count))
    return count

def getLastLineOfFile(fileName,DEBUG=False):
    """
    Quick retrival of last line in file, based on note:
    https://stackoverflow.com/questions/46258499/read-the-last-line-of-a-file-in-python
    :param fileName:
    :param DEBUG:
    :return:
    """
    with open(fileName, 'rb') as fileOut:
        fileOut.seek(-2, os.SEEK_END)
        while fileOut.read(1) != b'\n':
            fileOut.seek(-2, os.SEEK_CUR)
        lastLine = fileOut.readline().decode()
        if DEBUG:
            print("%s has last line: " % fileName)
            print(lastLine)
        return lastLine

def findFirstOccuranceInFile(fileName, subString):
    """
    Finds the first occurance of the substring in the file, returns the index and the line
    :param fileName:
    :param subString:
    :return: tuple of index and the line where the string occured, otherwwise (-1,0)
    """
    with open(fileName,'r') as inputFile:
        for index, line in enumerate(inputFile,1):
            if subString in line:
                return index,line
    return (-1,None)

def findFirstOccuranceInCompressedFile(fileName, subString):
    """
    Finds the first occurance of the substring in the file, returns the index and the line
    :param fileName:
    :param subString:
    :return: tuple of index and the line where the string occured, otherwwise (-1,0)
    """
    with gzip.open(fileName,'rt') as inputFile:
        for index, line in enumerate(inputFile,1):
            if subString in line:
                return index,line
    return (-1,None)


def getHeaderLinesInCompressedFile(fileName, char="#"):
    """
    Gets the header lines
    :param fileName:
    :param char: Character deliniating header in formation
    :return: list of lines
    """
    headerLines=[]
    with gzip.open(fileName, 'rt') as inputFile:
        for index, line in enumerate(inputFile, 1):
            if line[0]==char:
                headerLines.append(line)
            else:
                return headerLines
    print("Error, made it to the end of file, entirely header")
    raise ValueError

def readHaplotypesFromVcf(inputFile,DEBUG=False):
    """
    Reads in a vcf file and parses to get the haplotyes
    CHROM   POS ID REF ALT  QUAL FILTER  ...  i93  i94  i95  i96  i97  i98  i99
        1     9  .   A   T  1000   PASS  ...  1|1  1|1  1|1  1|1  1|1  1|1  1|1
        1    52  .   A   T  1000   PASS  ...  0|0  0|0  0|0  0|0  0|0  0|0  0|0
        1   122  .   A   T  1000   PASS  ...  0|0  0|0  0|0  0|0  0|0  0|0  0|0

    :param inputFile:
    :return: matrix of the vcf haplotypes transposed
    """
    numRowsToSkip = findFirstOccuranceInFile(inputFile, "#CHROM")[0] - 1
    vcfTable = pd.read_csv(inputFile, sep="\t", skiprows=numRowsToSkip)
    vcfTable.rename(columns={'#CHROM': CHROM_COL}, inplace=True)
    indivColumns = set.difference(set(vcfTable.columns), set(VCF_START_COLS))
    # Creates a vcf where it splits on the haplotypes as well, assuming seperated by | which
    # for simulations should be the case.
    # This regex delimiter appears considerably more efficient than spliting after the fact.
    vcfHap = pd.read_csv(inputFile, sep="[\|\s]", skiprows=numRowsToSkip + 1, header=None,engine='python')
    hapColumns = list(chain.from_iterable((i + MOM_HAP_SUFFIX, i + DAD_HAP_SUFFIX) for i in indivColumns))
    vcfHap.columns = VCF_START_COLS + hapColumns
    columnSet = [CHROM_COL, ID_COL, REF_COL, ALT_COL] + hapColumns
    vcfTransposedGenomes = ((vcfHap.set_index(POS_COL))[columnSet])
    return vcfTransposedGenomes

def readHaplotypeTable(inputFile,removeNulls = True, removeIndels=True, DEBUG=False):
    """
    Reading in a file that is a table, likely generated from merging IMGT and IHIW resources,
    looking similar to this:

    CHROM   POS     REF     ALT     QUAL    FILTER  INFO    FORMAT  A_02:01:01:01_0 A_02:01:01:01_1
    6       0       C       *       NA      NA      NA      NA      0               0
    6       2       G       A       NA      NA      NA      NA      0               0
    6       2       G       *       NA      NA      NA      NA      0               0
    6       3       G       A       NA      NA      NA      NA      1               1

    So, this is like a vcf but we already have columns broken up by haplotype.

    :param inputFile: input table to process
    :param removeNulls: Nulls are either '*' or '.'
    :param removeIndels: any alternate allele with more than 1 character
    :param DEBUG:
    :return: matrix of the haplotypes transposed, should be same format as readHaplotypesFromVcf
    """
    hapTable = pd.read_csv(inputFile, sep="\t")
    indivColumns = set.difference(set(hapTable.columns), set(VCF_START_COLS))
    if DEBUG:
        print("Reading in haplotype table, before removing nulls/indels there %s variants for %s samples" %
              (hapTable.shape[0], len(indivColumns)))
    if removeNulls:
        hapTable = hapTable[hapTable[ALT_COL] != '*']
        hapTable = hapTable[hapTable[ALT_COL] != '.']
        if DEBUG:
            print("Removing nulls, there %s variants left" % hapTable.shape[0])
    if removeIndels:
        hapTable = hapTable[hapTable[ALT_COL].str.strip().str.len().lt(2)]
        hapTable = hapTable[hapTable[REF_COL].str.strip().str.len().lt(2)]
        if DEBUG:
            print("Removing indels, there %s variants left" % hapTable.shape[0])
    columnSet=[CHROM_COL,ID_COL,REF_COL,ALT_COL]+ list(indivColumns)
    return (hapTable.set_index(POS_COL))[columnSet]



def readVcfIntoDataFrame(vcfFile,splitPhasedHaps=True,DEBUG=False):
    """
    Read in a vcf file into a pandas data frame, ignore the header.
    from simulation results the columns look
    This is for simulated data for now
     CHROM   POS ID REF ALT  QUAL FILTER  ...  i93  i94  i95  i96  i97  i98  i99
        1     9  .   A   T  1000   PASS  ...  1|1  1|1  1|1  1|1  1|1  1|1  1|1
        1    52  .   A   T  1000   PASS  ...  0|0  0|0  0|0  0|0  0|0  0|0  0|0
        1   122  .   A   T  1000   PASS  ...  0|0  0|0  0|0  0|0  0|0  0|0  0|0

    if splitting phasing adds in the columns
    CHROM  POS ID REF ALT  ...  i98_Dad i99_Mom i99_Dad
    1    9  .   A   T  ...        1       1       1
    1   52  .   A   T  ...        0       0       0
    1  122  .   A   T  ...        0       0       0

    :param vcfFile: name of the vcf file
    :param splitPhasedHaps: flag to split phased haplotypes into seperate columns
    :param DEBUG:
    :return:
    """
    numRowsToSkip = findFirstOccuranceInFile(vcfFile, "#CHROM")[0] - 1
    vcfDf = pd.read_csv(vcfFile, sep="\t", skiprows=numRowsToSkip)
    vcfDf.rename(columns={'#CHROM': CHROM_COL}, inplace=True)
    if DEBUG:
        print(vcfDf.head())
    if splitPhasedHaps:
        indivColumns = vcfDf.columns[pd.Series(vcfDf.columns).str.startswith('i')]
        for indivId in indivColumns:
            vcfDf[[indivId + MOM_HAP_SUFFIX, indivId + DAD_HAP_SUFFIX]] = \
                vcfDf[indivId].str.split('|', expand=True)
            vcfDf[indivId + MOM_HAP_SUFFIX] = vcfDf[indivId + MOM_HAP_SUFFIX].astype(np.int64)
            vcfDf[indivId + DAD_HAP_SUFFIX] = vcfDf[indivId + MOM_HAP_SUFFIX].astype(np.int64)
        if DEBUG:
            print("Phased haplotypes split up: ")
            print(vcfDf.head())
    return vcfDf





def readCompressedVcfOverIntervalAndGetHaps(chromosome, startPos, endPos, vcfFile, removeIndels=True,
                                            KeepMonomorphicSites=False, KeepMultiAllelic=False,DEBUG=False):
    #TODO consider doing this in the future, for now let's just stop if somehow this flag accidently set
    if KeepMultiAllelic:
        print("Error, proper handling of keeping multi-allelic sites has not been implemented yet in readCompressedVcfOverIntervalAndGetHaps.")
        raise ValueError
    regionTabix = pysam.VariantFile(vcfFile, 'r')
    if endPos-startPos<1:
        print("Error, invalid intervals in readCompressedVcfOverInterval: ")
        print("%s from %s to %s  " % (chromosome, startPos, endPos))
        raise ValueError
    samples = list((regionTabix.header.samples))
    hapColumns = list(chain.from_iterable((s + MOM_HAP_SUFFIX, s + DAD_HAP_SUFFIX) for s in samples))
    headerColumns = [CHROM_COL,POS_COL, ID_COL, REF_COL, ALT_COL] + hapColumns
    outputIndex = 0
    vcfTable = pd.DataFrame(columns=headerColumns, index=range(endPos-startPos))
    if DEBUG:
        print("%s from %s to %s  " % (chromosome, startPos, endPos ))
    for record in regionTabix.fetch(chromosome, (startPos), (endPos)):
        #intervalIndicesFound.append(index)
        if DEBUG:
            print(record.pos)
            print(record.id)
            if pd.isnull(record.alts):
                print("Null at alt")
                print(record.alts)
            else:
                print(record.alts)
                print(len(record.alts))
            print(record.ref[0])
            print(record)
        # Let's remove indels
        # For now if we have a multi-allelic site and any of them are indels we remove,
        # in the future may want to handle this differently and say only remove the indel but leave the SNP
        anyAltIndel=False
        if pd.isna(record.alts):
            if DEBUG:
                print("Warning, no alternate alleles %s %s %s %s %s" %
                  (chromosome,record.pos,record.id, record.ref[0],record.alts))
        else:
            for curAlt in record.alts:
                if len(curAlt.replace(" ", "")) > 1:
                    anyAltIndel = True
        if pd.isna(record.alts):
            altSetString = "NULL"
        else:
            altSetString = ",".join(record.alts)
        if removeIndels and (len(record.ref[0].replace(" ", "")) > 1 or anyAltIndel):
                if DEBUG:
                    print("Indel at chrom %s and positions %s, ref %s and alts: %s, we're removing" % \
                        (record.chrom, record.pos, record.ref, altSetString))
        else:
            vcfTable.at[outputIndex, CHROM_COL] = record.chrom
            vcfTable.at[outputIndex, POS_COL] = record.pos
            vcfTable.at[outputIndex, ID_COL] = record.id
            vcfTable.at[outputIndex, REF_COL] = record.ref[0]
            vcfTable.at[outputIndex, ALT_COL] = NULL_VALUE if  pd.isnull(record.alts) else record.alts[0]
            curPos = record.pos
            #Introducing this term to keep track of multi-allilic sites
            alleleSet = set()
            for key, value in record.samples.iteritems():
                if DEBUG:
                    print("%s - %s" % (key,value.items()))
                    print(value.items()[0][1][0])
                    print(value.items()[0][1][0])
                    print(set(value.items()[0][1]))
                    print("test!")
                    print(len(value.items()[0][1]))
                    print(value.items()[0][1][1])
                    print("done")
                #Arbitrarily asigning mom to first and dad to second, may want to redo this so it won't cause potential
                #confusion later but, that info isn't actually used
                #allowing for nulls now
                allele1=  np.nan if value.items()[0][1][0] == None else value.items()[0][1][0]
                if len(value.items()[0][1])==1:
                    allele2 = np.nan
                else:
                    allele2 = np.nan if value.items()[0][1][1] == None else value.items()[0][1][1]
                vcfTable.at[outputIndex, key+MOM_HAP_SUFFIX] = np.double(allele1) #value.items()[0][1][0]
                vcfTable.at[outputIndex, key+DAD_HAP_SUFFIX] = np.double(allele2) #value.items()[0][1][1]
                alleleSet= alleleSet | set(value.items()[0][1]).difference(set([None]))
            #Dealing with multiallelic sites, monomorphic sites, and other odd/unexpected situations
            if len(alleleSet)>2:
                print("Warning, at chrom %s and positions %s we have a mutli-allelic site we're removing, ref %s and alts: %s" %\
                          (record.chrom,record.pos,record.ref,altSetString))
                # May want to make this more explicit, we later clean this up and remove it
                # since the chrom is set to null signaling to remove
                vcfTable.at[outputIndex, CHROM_COL] = np.nan
            #Less common situation where at a multi-allelic site in our samples it's biallelic,
            # ie only 2 alleles observed in the sample ref/alt instead of multiple alleles like ref A and alts: G,C
            if len(alleleSet)==2 and len(record.alts)>1:
                if DEBUG:
                    print("alleleSet observed:")
                    print(alleleSet)
                if alleleSet != {0,1}:
                    print("Warning, at chrom %s and positions %s we have a mutli-allelic site we're recoding "
                          "since it is biallelic in our sample,\n ref %s and alts: %s codedes as %s so %s items" % \
                        (record.chrom, record.pos, record.ref, altSetString, str(alleleSet), len(alleleSet)))
                    if DEBUG:
                        print("\n\n\n Coding BEFORE Fixed: ")
                        print(vcfTable.loc[outputIndex])
                        print("\n\n\n")
                    minAllele = min(list(alleleSet))

                    for key, value in record.samples.iteritems():
                        if vcfTable.at[outputIndex, key+MOM_HAP_SUFFIX] == minAllele:
                            vcfTable.at[outputIndex, key + MOM_HAP_SUFFIX] = 0
                        elif vcfTable.at[outputIndex, key+MOM_HAP_SUFFIX] is np.nan:
                            if DEBUG:
                                print("Warning, at chrom %s and positions %s we have a incorrectly coded site with nulls, ref %s and alts: %s" % \
                                    (record.chrom, record.pos, record.ref, altSetString))
                        else:
                            vcfTable.at[outputIndex, key + MOM_HAP_SUFFIX] = 1
                        if vcfTable.at[outputIndex, key+DAD_HAP_SUFFIX] == minAllele:
                            vcfTable.at[outputIndex, key + DAD_HAP_SUFFIX] = 0
                        elif vcfTable.at[outputIndex, key+DAD_HAP_SUFFIX] is np.nan:
                            if DEBUG:
                                print("Warning, at chrom %s and positions %s we have a incorrectly coded site with nulls, ref %s and alts: %s" % \
                                    (record.chrom, record.pos, record.ref, altSetString))
                        else:
                            vcfTable.at[outputIndex, key + DAD_HAP_SUFFIX] = 1
                    if DEBUG:
                        print("\n\n\n We filtered out mutliallelic sites thar are really biallelic in our set vcfTable.at[outputIndex, :]")
                        print(vcfTable.loc[outputIndex])
                        print("\n\n\n")
            if len(alleleSet) == 1 and not KeepMonomorphicSites:
                if DEBUG:
                    print("Warning, at chrom %s and positions %s we have a monomorphic site so we're removing, ref %s and alts: %s" % \
                        (record.chrom, record.pos, record.ref, altSetString))
                    print("alleleSet observed:")
                    print(alleleSet)
                vcfTable.at[outputIndex, CHROM_COL] = np.nan
            if len(alleleSet) == 0:
                print("Warning, no alleles found at chrom %s and positions %s we have a mutli-allelic site we're removing, ref %s and alts: %s" %\
                      (record.chrom,record.pos,record.ref,altSetString))
                print("alleleSet observed:")
                print(alleleSet)
                vcfTable.at[outputIndex, CHROM_COL] = np.nan
            if int(curPos) < startPos or int(curPos)> endPos:
                print("Warning, appears to be an issue in readCompressedVcfOverInterval where record returned values not in range")
                print("%s from %s to %s vs current pos %s\n" % (chromosome, startPos, endPos,curPos))
                print(record)
                vcfTable.at[outputIndex, CHROM_COL] = np.nan
            outputIndex+=1
    #now clean and return
    vcfTableCleaned = vcfTable[vcfTable[CHROM_COL].notna()]
    columnSet = [CHROM_COL, ID_COL, REF_COL, ALT_COL] + hapColumns
    return ((vcfTableCleaned.set_index(POS_COL))[columnSet])

