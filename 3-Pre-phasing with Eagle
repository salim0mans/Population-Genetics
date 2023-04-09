1-Get the genetic map for Eagle from its repo:
> wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz

2-Split per Chromosome:
> for CHR in {1..23}; do
    zcat genetic_map_hg38_withX.txt.gz | \
    grep ^${CHR} | \
    sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \
    > eagle_chr${CHR}_b38.map
done

3-Make the pre-phasing

> DATASET=data

# Run phasing for each chromosome
for CHR in {20,22}; do 

    # Run phasing for each chromosome separately
    eagle \
           --vcf ${DATASET}_for_phasing.vcf.gz \
           --chrom ${CHR} \
           --geneticMapFile eagle_chr${CHR}_b38.map \
           --numThreads=8 \
           --Kpbwt=20000 \
           --outPrefix ${DATASET}_for_imputation_chr${CHR} 
done

-----------------Now we have our data ready for beagle imputation-------------------
