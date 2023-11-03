# 5-Imputation with different Software  
  
## Introduction  
  
Imputation aims to predict genotypes based on a reference panel and calculated LD scores.  
  
There are different imputation software that rely on similar algorithms (and sometime the same base algorithm) to make the prediction. Most of these software undergo constant updates and new releases.  
  
In this guide, we will put the standard execution commands for three popular tools corrently: Beagle 5.4, Impute5 1.1.5, Minimac4.  
  
## Protocols  
  
make sure to index the vcf files upon each step with bcftools index -t \<file.vcf.gz\>  
  
### 1-Beagle  
  
-Download beagle from its official site as a jar file.  
  
-Get the genetic map for beagle (choose correct build) from <a href="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/">this site</a>.  
  
Extract and prepare files per Chr  
  
  $ unzip plink.GRCh38.map.zip  
  $ mv plink.chrX.GRCh38.map plink.chr23.GRCh38.map  
  $ for CHR in {1..23}; do                             ## This command renames chromosomes in case of using hg38 only  
    cat plink.chr${CHR}.GRCh38.map | \  
    sed 's/^/chr/' \  
    > beagle_chr${CHR}_b38.map  
  done    
  
Data Imputation with Beagle.22Jul22.46e: 
  
  $ for R in {5,10,15}; do for F in {1..10}; do for CHR in {1..22}; do        #Proceeding with our example with 3 masks and 10 folds.   
    java -Xmx16g \
      -jar beagle.22Jul22.46e.jar \
        gt=Masked_phased_chip_${R}_${F}.vcf.gz \
        ref=/path/to/final_1000G_chr${CHR}.vcf.gz \
        out=imputed_chip_${R}_${F}_chr${CHR} \
        map=/path/to/beagle_chr${CHR}_b38.map \
        nthreads=16 \
        iterations=10 \
        ne=20000 \
        impute=true \
        gp=false \
        seed=-99999 
  done;done;done  
  
  
### 2-Minimac4  
  
-Minimac4 and Minimac3 need download and installation (with make-install and make).  

-Minimac4 depends on specific formats for the reference files. First, Minimac3 is needed to convert vcf reference files into m3vcf files. Then, in the newer releases, the reference needs to be updated to msav format.  
  
  $ for CHR in {1..22}; do
    path/to/Minimac3 --refHaps final_1000G_chr${CHR}.vcf.gz \ 
    --processReference \ 
    --prefix final_1000G_chr${CHR}
  done  

  $ for CHR in {1..22}; do
  path/to/minimac4 --update-m3vcf final_1000G_chr${CHR}.m3vcf.gz > final_1000G_chr${CHR}.msav
  done

  $ for R in {5,10,15}; do for F in {1..10}; do for CHR in {1..22}; do 
    path/to/minimac4 --refHaps path/to/final_1000G_chr${CHR}.msav \
    --haps Masked_phased_chip_${R}_${F}.vcf.gz \
    --window 500000 \
    --prefix imputed_chip_${R}_${F}_chr${CHR} \
    --chr chr${CHR} \ 
    --noPhoneHome \ 
    --format GT,DS \
    --minRatio 0.00001 
  done;done;done  
  
### 3-Impute5  
  
-Download Impute5 1.1.5 from its official sites (probably in the older versions section). Genetic maps are available in the in the "misc" directory.  

-Impute5 also requires conversion of the reference files into a compatible formate (imp5).  
  
  $ for CHR in {1..22}; do
    path/to/imp5Converter_1.1.5_static \                 #Convert tool
    --h path/to/final_1000G_chr${CHR}.vcf.gz \
    --r chr${CHR} \
    --o final_1000G_chr${CHR}.imp5
  done  
  
  $ for R in {5,10,15}; do for F in {1..10}; do for CHR in {1..22}; do
    path/to/impute5_1.1.5_static \                       #Imputation tool
   --h path/to/final_1000G_chr${CHR}.imp5 \
   --m path/to/impute5_maps/chr${CHR}.b38.gmap.gz \
   --g Masked_phased_chip_${R}_${F}.vcf.gz \
   --o imputed_chip_${R}_${F}_${CHR}.vcf.gz \
   --r chr${CHR} \
   --threads 15 \
  done;done;done   


-Concatenation with beagle is also an option for split files.  

-----------------------Now we Finished Imputation of each Chromosome------------------------
