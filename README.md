# Linkage disequilibrium approximate Bayesian factor (LD-ABF) 
## Improved Characterization of Balancing Selection Genome Wide and a Detailed Look at HLA Genes

Abstract:

Balancing selection occurs when different evolutionary pressures impact the fitness of multiple alleles, resulting in increased allelic diversity in the population. A new statistical method was developed to test for selection, improving inference by using efficient Bayesian techniques to test for density and strength of linkage disequilibrium. Evolutionary simulation studies showed that the method consistently outperformed existing methods. Using this methodology, we tested for novel signals of balancing selection genome wide in 500 samples from phased trios, combining both SNP array and exome sequencing data. Several novel signals of selection appeared in CYP2A7, GPC6, and CNR2 across multiple ancestries. Additionally, tests in SIRPA demonstrate dramatically strong selection signal, significantly higher than previously observed. Well-known signals around olfactory genes and the MHC, containing HLA genes associated with the immune response, also demonstrated strong signatures of selection. So, utilizing data from the 17th IHIW, a follow up analysis was then performed by leveraging over seven thousand HLA typed samples by NGS; in contrast, the genome wide scan did not include a detailed characterization of the HLA genes.  The strongest signals observed in the IHIW samples were in DQA1 and DQB1 in or around exon 2–the portion of the gene responsible for antigen presentation and most likely to be under environmental and evolutionary pressure. Our new statistical approach and analysis suggest novel evolutionary pressure in new regions and additionally highlight the importance of improved sequencing and characterization of variation across the extended MHC and other critical regions.

Manuscript under review, further details to come.

-----------------------------------------------------
### Download LD-ABF supplemental files.

Supplemental files from the manuscript can be found online here:
 
[Genome Wide scan of CHOP Trios](https://upenn.box.com/s/yikbstl82z8g5negl7mbwd7epqapkab5), 
[HLA genes from IHIW](https://upenn.box.com/s/vqtso3lhyckwuud0qpjpkp1yd8hx636j)

This includes the statistics for each of the different populations for the two different data sets along with some additional plots and tables. 

## Running LD-ABF

LD-ABF is a test statistic that measures the amount of LD and density around a given SNP. 

![Figure 1](https://github.com/tris-10/LD-ABF/blob/main/figures/BalancingSelectionOverTime.jpg)
*Figure 1 Diagram depicting the progression of an allele under balancing selection, the green X. In the first pane the variant is introduced on a single haplotype. Then after some time has passed evolutionary pressures favoring heterozygotes at the position of focus maintaining both haplotypes with and without the polymorphism, where hitchhiking effects are observed between the variant under balancing selection–inducing LD patterns.  Recombination breaks the strong LD resulting in mosaics of the haplotypes, where strong hotspots will diffuse the effects. *

Python 3.7 and libraries required: `os, sys, argparse, time, datetime, pstats, cProfile, scipy, cython, shutil, glob, numpy, pandas, scipy, statsmodels, collections, itertools` 

Sample input files and examples coming.

To run LD-ABF: 


```
usage: LD_ABF.py [-h] [-i INPUTFILE] [-w WINDOWSIZE] [-m MINMAF]
                 [-s SINGLEPOS] [--chromosome CHROMOSOME]
                 [--start-pos STARTPOS] [--end-pos ENDPOS] [--inputAreHaps]
                 [--not-inputAreHaps] [--debug] [-o OUTPUTFILEPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        Input file with variants to be analyzed, in *.vcf.gz
                        format. Unless inputAreHaps set to true this should be
                        compressed with corresponding tabix *.vcf.gz.csi file
                        (http://www.htslib.org/doc/tabix.html).
  -w WINDOWSIZE, --windowSize WINDOWSIZE
                        Window size for scan statistic.
  -m MINMAF, --minMaf MINMAF
                        Minimum minor allele frequency to consider a site
                        polymorphic
  -s SINGLEPOS, --singlePos SINGLEPOS
                        Running test on single position otherwise None to do
                        full scan or scan over range
  --chromosome CHROMOSOME
                        If running a scan over certain positions or just
                        looking over one positions, otherwise set to None as
                        default
  --start-pos STARTPOS  If running a scan over certain positions this is the
                        starting position
  --end-pos ENDPOS      If running a scan over certain positions this is the
                        end position
  --inputAreHaps        Flag indicating if the input file contains columns
                        with only haplotypes instead of biallelic phased
                        values (ie 0 or 1 instead of {0|0,1|0,0|1,1|1}))
  --not-inputAreHaps
  --debug
  -o OUTPUTFILEPREFIX, --outputFilePrefix OUTPUTFILEPREFIX
                        Output file location prefix to be written out
```



### Running LD-ABF test scan


Let's say you want to run the LD-ABF test over a region in a bed file, say region defined by the first line of regions.bed for the compressed vcf file, input.vcf.gz with corresponding tabix input.vcf.gz.csi. Then you want to restrict to common alleles greater than 5%, ie minor allele frequency > 0.05.

```
chrom=6
startPos=149721689
endPos=149772647
outputDir="output/"
inputDir="input/"
vcfFile="${inputDir}input.vcf.gz"
python src/LD-ABF.py \
        -i ${vcfFile} --chromosome ${chrom} --start-pos ${startPos} --end-pos ${endPos} \
        -o ${outputDir}chr${chrom}_${startPos}_${endPos} --minMaf 0.05;
```

This will output a file that looks something like:
```
CHROM                 POS                   ID                    REF                   ALT                   LD-ABF                linkageSum            sumCorr               densityInWin          MAF
6                     149721690             rs237025              A                     G                     0.02258162572615963   0.9139784946236562    0.9139784946236558    1.0                   0.47058823529411764
6                     149721965             rs237024              C                     T                     0.0222778955943641    0.9139784946236562    0.9139784946236558    1.0                   0.375
...
6                     149770179             rs237005              C                     T                     0.0                   0.0                   0.0                   0.0                   0.4117647058823529
6                     149772190             rs112722576           A                     G                     0.009038690878672113  0.393939393939394     0.3939393939393939    1.0                   0.3235294117647059
6                     149772645             rs415434              G                     A                     0.008272128024190223  0.393939393939394     0.3939393939393939    1.0                   0.4583333333333333
```


Or let's say you want to run the LD-ABF test over a region in a bed file, say region defined by the first line of regions.bed for the same scenario. 
 
```
lineNum=1
bedFile="regions.bed"
curLine=$(awk "NR==${lineNum}" $bedFile)
lineAsArray=($curLine)
echo ${lineAsArray}
chrom=`echo $curLine | awk -F" " '{print $1}'`
startPos=`echo $curLine | awk -F" " '{print $2}'`
endPos=`echo $curLine | awk -F" " '{print $3}'`
outputDir="output/"
inputDir="input/"
vcfFile="${inputDir}input.vcf.gz"
python src/LD-ABF.py \
        -i ${vcfFile} --chromosome ${chrom} --start-pos ${startPos} --end-pos ${endPos} \
        -o ${outputDir}chr${chrom}_${startPos}_${endPos} --minMaf 0.05;
```

You could then modify this to run a scan across multiple regions in a bed file on a cluster just by passing in different line numbers for the bed region to look over.  