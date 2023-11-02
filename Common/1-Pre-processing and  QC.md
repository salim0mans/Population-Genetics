## Pre-processing and Quality Control  

### Introduction  
Genotype imputation is a method for genotype prediction, that applies different statistical models.
It usually handles vcf (variant call format) files for their organized structure that facilitate the storage of information about genetic variants in numerous fields. vcf files have two main sections; the Header (all its rows start with "#") and the alignment (that contains the fields as columns). Columns include Chromosome (CHROM), ID or variant, position (POS), reference allele (REF), alternative allele (ALT), genotype quality (QUAL), FILTER field that shows the filtration status of variants against some threshold, INFO field that contains many important metrics and parameters related to describing these variants, and FORMAT fields that shows the actual haplotypes or genotypes for each sample (organism) represented by a particular system (1/0 or 0/1/2, etc.), in addition to some parameters related to individual genotypes (e.g., genotype probability, Dosage, etc.). Each column in the Format fields represents one sample.

Each row represents one variant, and these variants can be simple (SNVs, and Indels), or complex structural variants (such as CNVs).

In the general cases, we are interested only in bi-allelic simple variants (unless our study requires other aspects), which I will include.

Moreover, some casses  require special treatment, such as dealing with chromosome X in case of hemizygous individuals (e.g., human males),pedigree studies, and HLA variants.

variant calling is based on a genomic assembly (fasta file that represents multiple human genomes in one), but assemblies become updated through time, the old (still in use) assembly is hg19 (GRCh38) that was released in Feb 2009, but the present assembly is hg38 (CRCh38). Different verions differ from each other by variant positions, addition of new variants, etc.

The "liftover" process is carried out to convert vcf data from one assembly to another. 

Microarray genotype data come as two alleles per individual, with a "/" between them, where 0 is the reference allele and 1 is the alternative allele. However, after cleaning there must be a phasing algorithm to distribute alleles into two haplotypes according to a reference panel or in a reference-free statistical manner (for big samples). Then, comes the imputation algorithm that calculate the Linkage Disequilibrium (LD) score and makes prediction of missing genotypes based on it, adding all variants available in the reference panel to the imputed data and filling all genotype by prediction. 
It is proven that the larger and more population representative the reference panel, the more accurate the imputation. Also, sample size contributes a lot to the accuracy and precision.

### Pipeline  

0-Prerequisites:  

Tools required for standalone imputation pipelines include a vcf manipulation tool (plink, Bcftools, Vcftools), a phasing tool (Eagle, SHAPEIT), an imputation tool (Beagle, Impute5, Minimac4, PBWT, etc.). In addition, some tools are needed for extra tasks, such as a liftover tool (CrossMap), an accuracy calculation tools, a visualization method (R ggplot2, python matplotlib, or graphpad).  
Most tasks are performed in a linux environment (with Unix bash script).


1-Determine the genomic build (assembly) of usage:

One might need two builds: hs37d5.fa (hg19) and the Homo_sapiens.GRCh38.fa assembly (hg38), because some tasks require software and servers that currently have the hg19 only (and vice versa). However, if the task involves only quick imputation, hg38 is preferred. 

2-Obtain the reference panel files

There are only few internationally recognized open-access reference panels. 1000 Genomes is one of the most popular ones. It offers different releases and versions with size between 2504-3600 samples and 38-82 Million variants (according to the version). HGDP (Human Genome Diversity project) is another panel that contains heavy metadata and variants that are omitted in other panels, however its drawback is small size (925 samples). A similar panel is also available, SGDP (Simon's Genome Diversity Project), which contains diverse samples from a wide range of ethnicities, but also has a very small size (279 samples).
Other, more credible, reference panels are available with access restrictions, including the Haplotype Reference Consortium (HRC) which is considered one of the most useful and important references with a size of 32,470 samples gathering about 39.6 Million variants. Although, the largest known imputation reference panel is TopMed with more than 97250 samples and over 300 Million variants. These panels are available in public imputation servers (Michigan, Sanger, TopMed imputation servers).

We will take 1000G as an example
Choose the preferred release (and pay attention to the build):
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/
(downloadable with $ wget <URL>)

3-QC for the panel
  3.1-change names of chromosomes to make sure that they are consistent with the build. GRCh37 takes only 1..22,XY but GRCh38 takes chr1,chrX, etc.
  
$ touch chr_names.txt
$ for CHR in {1..22} X Y; do 
    echo -e "${CHR}\tchr${CHR}" >> chr_names.txt
done 

  3.2-Remove rare variants (with Allele Count AC<3 or Minor Allele Frequency MAF<0.01)  
  3.3-Normalize data (flips, multialleles, ambiguity, etc.)
  3.4-keep only SNPs and INDELs
  3.5-Align with reference genome to ensure that REF matches reference genome (and remove duplicate variants with -d none)
  3.6-remove multiallelic records
  3.7-remove sites with missing data
  3.8-set IDs as CHR_POS_REF_ATL

$ for CHR in {1..22}; do
#3.1
    bcftools annotate --rename-chrs chr_names.txt \
        1000G_hg38_chr${CHR}.vcf.gz -Oz | \
#3.2
    bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3 | INFO/MAF<0.01' -Oz | \
#3.3
    bcftools norm -m -any -Oz | \
#3.4   
    #(VT = variant type, it might be present int eh INFO field, but if not, maybe searching the INFO might help finding another tag for the same purpose. However if there is no such tag, then one might filter manually all complex variants by tracking a pattern in the ID, type, or INFO. For example, in some version of 1000G, there are CNVs and other variants that all have "[]" in their names or descriptive fields, hence I could filter them by knowing that)#

    bcftools view -i 'INFO/VT="SNP" | INFO/VT="INDEL"' -Oz | \ 
#3.5
    bcftools norm -f /Path/to/hs37d5.fa -d none -Oz | \
#3.6
    bcftools view -m 2 -M 2 -Oz | \
#3.7
    bcftools view -g ^miss -Oz | \
#3.8
    bcftools annotate --set-id "%CHROM\_%POS\_%REF\_%ALT" -Oz -o QC_1000G_chr${CHR}.vcf.gz
done

4-Remove duplicate IDs:
#Check and extract dublicates by a query of ID, followed by choosing only the lines that were repeated, then we put results in a text file and filter the data we have excluding the IDs in that file: 

$ for CHR in {1..22}; do
   bcftools query -f '%ID\n' QC_1000G_chr${CHR}.vcf.gz | \
    sort | uniq -d > Dup_1000G_chr${CHR}.txt

    if [[ -s Dup_1000G_chr${CHR}.txt ]]; then
    	bcftools view -e ID=@Dup_1000G_chr${CHR}.txt \
    	QC_1000G_chr${CHR}.vcf.gz \
        -Oz -o filtered_1000G_chr${CHR}.vcf.gz
    else 
    	mv QC_1000G_chr${CHR}.vcf.gz filtered_1000G_chr${CHR}.vcf.gz
    fi
done

5-Obtain Allele frequencies (AF) in the file. Then, extract all frequencies into a single file (let's call it ref_freq). It should have a column for IDs in the format CHR_POS_REF_ALT and an AF column. We concatenate all wanted chromosomes one by one in the same ref_freq file so it becomes inclusive.

$ for CHR in {1..22}; do
    bcftools +fill-tags filtered_1000G_chr${CHR}.vcf.gz -Oz -o AF_1000G_chr${CHR}.vcf.gz -- -t AF
done

$ echo -e 'CHR\tSNP\tREF\tALT\tAF' \
    > ref_freq_1000G.txt

# Query the required fields from the VCF file and append to the allele frequency file 
$ for CHR in {1..22}; do
    bcftools query -f '%CHROM\_%POS\_%REF\_%ALT\t%AF\n' AF_1000G_chr${CHR}.vcf.gz \
    >> ref_freq_1000G.txt
done

6-Generate a Panel-Sample-ID text file (to compare with any tested data for possible overlapping in samples)

$ bcftools query -l AF_1000G_chr1.vcf.gz > 1000G_SampleIDs.txt


------Now we have a reference panel that is Quality Controlled, along with its frequency table and Sample ID table------