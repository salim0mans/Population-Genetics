# Calculate Accuracy of imputation  
  
## Introduction  
  
Since imputation is a statistical prediction method, some metrics were found to estimate its effeciency, quality, and accuracy.  
  
* R^2  
R-square represents the correlation between the imputed and actual genotype based on their concordance. It is a widely used measure to assess the quality of imputation. R^2 measure is calculated for each single variant, and can be averaged for each sample.  
* IQS (Imputation Quality Score)  
IQS is also based on concordance with some weights related to allele frequency. It is considered a strict quality measure for imputation procedure itself, and in most cases, it renders very close indicators about quality to R^2. IQS is calculated for each single variant.  
* Concordance (P0)  
Concordance is a direct measure of the percentage of matching between imputed and actual. It provides an initial insight, however, direct identical match doesn't always represents quality of prediction, since it doesn't take in consideration other aspects, such as effect of some genotype in the prediction process, the difference between homozygous and heterozygous individuals, false positives and false negatives.  
* Precision (pr)  
Precision shows consistency of imputation and how many variants were predicted correctly among the predicted variants. It is represented by the simple equation of (TP)/(TP+FP).   
* Recall (rec)  
Recall is often called "sensitivity", as it calculates how many true variants were predicted among all true variants. It is represented by the simple equation of (TP)/(TP+FN).  
* F1.Score  
A measure that considers both precision and sensitivity. It is represented by the simple equatio of (2*pr*rec)/(pr + rec).  
  
### Calculating accuracy metrics  
  
There are different tools to calculate metrics, some are built-in within imputation software (mostly R^2 and Dosage R^2), and others are standalone scripts that are publicly available.  
- Imputation Accuracy Calculator by TorkamaniLab (github), <a href="https://github.com/TorkamaniLab/imputation_accuracy_calculator">found here</a>.  
- Hap.py calculator by Illumina (github), <a href="https://github.com/Illumina/hap.py/blob/master/doc/happy.md">found here</a>.  
  
### 1-Imputation Accuracy Calculator  
This tools needs the following inputs: the imputed file to be tested, the masked input file for imputation, named genotype array (GA), and a whole genome sequencing file (WGS) to represent the ground truth reference of comparison.  
WGS files are not always available, therefore, the original phased file BEFORE MASKING can serve the same purpose.  
The script compares files, taking only the intersecting variants with the ground truth file, and excluding all variants present in the (GA) file. Therefore, it ends up only comparing the masked variants between ground truth and imputed file.  
The software calculates numerous metrics per sample and per variant.  
  
From personal experience, overwhelming the RAM might happen even in powerful servers, and the software will crash and require repetition of the whole process.  Therefore, reducing and optimizing the max_per_core and max_total_rows parameters is recommended. Better start with defaults.    
Also, files are recommended to be done with all imputed chromosomes concatenated if the purpose was to test the overall imputation and not a specific chromosome, because the tools will calculate metrics per sample based on all variants in one file.  
  
    $ for R in {5,10,15}; do for F in {1..10}; do 
        python3 Compare_imputation_to_WGS.py --ga path/to/Masked_phased_chip_${R}_${F}.vcf.gz \  
        --wgs path/to/phased_chip.vcf.gz \  
        --imputed path/to/imputed_chip_${R}_${F}.vcf.gz  
        done;done  
  
The output text files will be created by default in the same directory as the imputed ones, and named after the same prefixes.  

result tables can be directly accessed for the purpose of visualization and comparisons.  
  
### 2-Hap.py  
Hap.py is one of many scripts offered by Illumina package. It works well evaluating imputation on basis of credible estimation algorithms. It is well suited for individual sample imputations.  
  
Hap.py is built on python2, which might be problematic for most conda environments, especially finding compatible libraries. Therefore, it is better to use Docker (which is provided in the main hap.py page).  
  
Outputs of hap.py are numerous, from simple summaries (that include TP,FN, FP, etc.) and ROC coordinates. Using summaries is rather sufficient to evaluate the accuracy of imputation.  
Inputs include a reference vcf file (the original phased unmasked file), followed directly by a query file (the imputed file), in addition to a region-determining bed file (that includes only masked variants after extracting them with bcftools by comparing the masked and the phased IDs and taking only the non-interesecting ones, and them convert them into a bed file). A genome assembly fasta file is required also (must be the same build).  
  
Conversion of MaskedArea files into bed format is available in a script under the name of "bed_conversion.py" in "Scripts" directory, excuting the script must be done using python2.  
  
    $ for R in {5,10,15}; do for F in {1..10}; do 
        hap.py path/to/phased_chip.vcf.gz \  
        path/to/imputed_chip_${R}_${F}.vcf.gz \  
        -R path/to/MaskedArea_${R}_${F}.bed.gz \  
        -r path/to/hs37d5.fa \  
        -o hap.py_chip_${R}_${F} \  
        --threads 8 \  
        --write-vcf
        done;done
  
Output csv files need to be processed by a reference database. GnomAd is considered one of the largest database for AFs. We can process files in GnomAd within OpenCravat software to determine the global AF of each variant. This can be achieved by re-calculating the TPs,NPs, and FPs and further calculating the Recall, Precision, and F.score.  
  
---------------Now we get reports of our metrics, ready to be visualized after some wrangling---------------