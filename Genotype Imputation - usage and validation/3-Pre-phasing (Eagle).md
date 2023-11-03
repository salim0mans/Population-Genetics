# 4-Pre-phasing (Eagle)  
    
## Introduction     
  
Phasing is a preparatory step for imputation, where genotypes undergo alignment of haplotypes according to a reference panel, or panel-free manner by a statistical model that compares alleles of samples considering the REF and ALT alleles and a genetic map.   
  
There are few well-known phasing software available publically, including Eagle, SHAPEIT, and Beagle (phasing and imputation software).  
We will perform phasing with Eagle in this guide.   
    
### 1-Download and install requirements
  
Eagle, along with the genetic maps of hg19 and hg38 assemblies are available on <a href="https://alkesgroup.broadinstitute.org/Eagle/downloads/">this site</a>.  
  
Eagle needs to be installed and configured in a linux environment.  
  
### 2-Process the genetic maps  
  
split the genetic map files by chromosome  
   
    $ for CHR in {1..22},X; do  
    zcat genetic_map_hg38_withX.txt.gz | grep ^${CHR} | sed '1ichr position COMBINED_rate(cM/Mb) Genetic_Map(cM)' \  
    > eagle_chr${CHR}_b38.map  
    done  
  
### 3-Perform a reference-based Pre-phasing with Eagle  
   
    $ for CHR in {1..22}; do    ##R represents masks and F represents folds (repeats)  
    path/to/eagle \  
       --vcfTarget /path/to/final_chip.vcf.gz \  
       --vcfRef /path/to/final_1000G_chr${CHR}.vcf.gz \
       --chrom ${CHR} \  
       --geneticMapFile path/to/eagle_chr${CHR}_b38.map \  
       --numThreads=8 \  
       --Kpbwt=20000 \  
       --outPrefix  phased_chip_${CHR}  
    done  

The output is split files (by chromosome), each file is a specific mask and a specific repeat.  
  
Warning: the reference-base procedure might be immensly long.   
  
Concatenation is possible (if preferred one file for all chromosomes per mask and repeat) using bcftools    

    $ bcftools concat \        ##an example on our chip data  
     phased_chip_1.vcf.gz phased_chip_2.vcf.gz phased_chip_3.vcf.gz phased_chip_4.vcf.gz phased_chip_5.vcf.gz phased_chip_6.vcf.gz phased_chip_7.vcf.gz phased_chip_8.vcf.gz phased_chip_9.vcf.gz phased_chip_10.vcf.gz phased_chip_11.vcf.gz phased_chip_12.vcf.gz phased_chip_13.vcf.gz phased_chip_14.vcf.gz phased_chip_15.vcf.gz phased_chip_16.vcf.gz phased_chip_17.vcf.gz phased_chip_18.vcf.gz phased_chip_19.vcf.gz phased_chip_20.vcf.gz phased_chip_21.vcf.gz phased_chip_22.vcf.gz -Oz -o phased_chip.vcf.gz  
  

-----------------Now we have a phased data, ready for the following steps-------------------
