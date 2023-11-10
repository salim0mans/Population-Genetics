# Microarray variant file (vcf) pre-processing  
  
## Introduction:  
  
Microarray data involves different early stages until it becomes annotated and formulated into a vcf file. Starting from FASTQ files through SAM/BAM files and ending with vcf files, sequences must undergo mapping and alignment to the reference genome assembly, also quality control processes, duplicate marking and recalibration (might involve bootstrapping), variant discovery, filtering and annotation. The whole pipeline is called "Variant Calling".  
  
Annotated data might contain missing genotypes, poorly genotyped sites, ambiguity, strand flips, etc,. Therefore, post calling quality control of vcf files is crucial.  
  
In addition, if the data were to be used as an input for imputation according to some reference panel, then some pre-processing involves matching and creating compatibility between the data files and the panel's files is needed and inevitable.  
  
## Pipeline:  
  
### 1-Rename Chromosomes according to genome build of choice.
  
    $ touch chr_names_chip.txt                          #Generate a chromosome renaming file  
    $ for CHR in {1..22} X ; do           
        echo ${CHR} chr${CHR}      #this direction is used for hg19 --> hg38 compatibility (otherwise, make it chr${CHR} ${CHR})  
    done >> chr_names_chip.txt  
    $ bcftools annotate --rename-chrs chr_names_chip.txt chip_data.vcf.gz -Oz -o renamed_chip.vcf.gz   
    
  
### 2-Keep only the wanted chromosomes   
(an example case is filtering 20,22)  
  
    $ chrs=$(echo chr{20,22} | tr ' ' ',')  
  
Keep only those wanted chromosomes  

    $ bcftools view -t $chrs renamed_chip.vcf.gz -Oz -o chrfiltered_chip.vcf.gz
  
### 3-QC of the chip  
3.1. Set the ID column as CHROM_POS_REF_ALT  
3.2. Align with reference genome to correct any flips  
3.3. Keep only bi-allelic records  
3.4. Remove duplicate samples, or overlapping ones with the panel  
3.5. Fill tags of AF, AC, F_MISSING, etc,.  
3.6. Exclude duplicate variants.
3.7. Set the minimum call rate (missing values), minimum AC, and minimum MAF (Recommended in general: exclude rare variant).  
  
    $ bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' chrfiltered_chip.vcf.gz -Oz -o SNPID_chip.vcf.gz #3.1    
    $ bcftools norm -f /path/to/hs37d5.fa -c ws SNPID_chip.vcf.gz -Oz | \ #3.2    
     bcftools view -m 2 -M 2 -Oz -o corrected_chip.vcf.gz    #3.3    

    $ touch duplicate_sample_IDs.txt                                                        #3.4
    $ bcftools query -l corrected_chip.vcf.gz | uniq -d >> duplicate_sample_IDs.txt    
    $ cat path/to/1000G_SampleIDs.txt >> duplicate_sample_IDs.txt   
    $ bcftools view -S ^duplicate_sample_IDs.txt --force-samples corrected_chip.vcf.gz -Oz | \  
     bcftools +fill-tags -Oz -- -t AC,AN,AF | \    #3.5    
     bcftools norm -d none -Oz | \     #3.6    
     bcftools view -e 'INFO/AC<3 | (INFO/AN-INFO/AC)<3 | INFO/MAF<0.01' -Oz -o QC_chip.vcf.gz   #3.7   

### 4-Match variants for presence and allele frequency between the panel and the chip  

    $ echo -e 'CHR\tSNP\tREF\tALT\tAF' \   #Generate a frequency table (chip_freq) for the chip     
    > chip_freq.txt     
    $ bcftools query -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' QC_chip.vcf.gz \
    >> chip_freq.txt

Compare AF in both chip and panel using an R script from <a href="https://www.protocols.io/view/genotype-imputation-workflow-v3-0-e6nvw78dlmkj/v2">this workflow by Priit Palta et al.</a>     

One can create an R file and copy the script within (the script is attached in the same "Scripts" directory under the name of "plot_AF_by_Plata.R")    

    $ touch plot_AF.R     
    $ nano plot_AF.R      #then copy the contents     

Now execute the file considering the args in the following order (chip_freq file path, prefix of output, ref_freq file path, maximum AF difference, maximum fold change):     
  
    $ Rscript --no-save /path/to/plot_AF.R chip_freq.tdt Test path/to/ref_freq_1000G.txt 0.1 5  
    
The output must include text files of excluded SNPs due to mismatch, unavailable SNPs, or very different AF, accompanied with a plot.   

We should remove every excluded SNPs    

    $ sort -V Test_exclude.txt \
    > Test_exclusion1.txt
    $ sort -V Test_nonpanel_exclude.txt \
    > Test_exclusion2.txt
    
    $ bcftools view -e ID=@Test_exclusion1.txt  QC_chip.vcf.gz | \  
    bcftools view -e ID=@Test_exclusion2.txt -Oz -o final_chip.vcf.gz
    
-----------------Now we have a cleaned vcf data file, compatible with our 1000G panel, and ready for next steps---------------------
