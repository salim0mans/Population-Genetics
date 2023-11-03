# Extraction of AFs from imputed files  

## Introduction  
  
Allele frequency comparison between samples that represent populations or groups, etc, grants insights about the effect of specific genetic variants, association between some alleles and a condition or phenotype, and so on.  
  
Starting from vcf files, we prefer imputing these files multiple times (different seeds), and taking the mean value of AF and AC.  
  
We subset the variants of interest (that we choose from literature, etc,)  
  
## Protocol  
  
### 1-Extract AF and/or AC from imputed files along with IDs (CHROM_POS_REF_ALT).  
  
    $ for F in {1..20}; do  
    echo -e "SNP,AC${F}" query_chip_${F}.csv  
    bcftools query -i ID=@SNP_IDs.txt -f "%ID,%AC\n" imputed_chip_${F}.vcf.gz  \
    >> query_chip_${F}.csv  
    done

### 2-Take mean value of all folds.  
  
The py script is available under "mean_AC.py" name, in the "Scripts" directory. Make sure to put the script in the same directory of the queries.  
  
The resulted csv file will be handled using R and/or excel.  