chrom=6
startPos=149721689
endPos=149772647
outputDir="output/"
inputDir="data/"
vcfFile="${inputDir}African_PanGenome_chr6_SegDupKept_filtered.vcf.gz"
python LD_ABF.py \
        -i ${vcfFile} --chromosome ${chrom} --start-pos ${startPos} --end-pos ${endPos} \
        -o ${outputDir}chr${chrom}_${startPos}_${endPos} --minMaf 0.05;
