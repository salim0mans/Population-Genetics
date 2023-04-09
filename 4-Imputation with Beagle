1-Get the genetic map specific for beagle:

> wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

2-Extract and prepare files per Chr:

# Unzip the files
unzip plink.GRCh37.map.zip

# Rename chromosome 23
mv plink.chrX.GRCh37.map plink.chr23.GRCh37.map

###THIS COMMAND IS FOR DATA WITH hg38 BUILD... I WILL PROBABLY SKIP IT!!
# Add 'chr' tag to the beginning of the line and store the output with suitable filename
for CHR in {1..23}; do 
    cat plink.chr${CHR}.GRCh38.map | \
    sed 's/^/chr/' \
    > beagle_chr${CHR}_b38.map
done


3-Data Imputation with Beagle.22Jul22.46e: 

DATASET=data
REFERENCE_BREF=/1000GP_AF
for CHR in {20,22}; do 
    java -Xss5m -Xmx8g \
      -jar beagle.22Jul22.46e.jar \
        gt=${DATASET}_for_imputation_chr${CHR}.vcf.gz \
        ref=${REFERENCE_BREF}_chr${CHR}.vcf \
        out=${DATASET}_imputed_chr${CHR} \
        map=tbeagle_chr${CHR}_b38.map \
        nthreads=16 \
        iterations=10 \
        ne=20000 \
        impute=true \
        gp=true \
        seed=-99999 
done

#Any problems that might appear are due differences in chr tag btw panel/chip/plinkmap


-----------------------Now we Finished Imputation of each Chromosome------------------------
