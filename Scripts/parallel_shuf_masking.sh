#!/bin/bash

## Pay attention that upon using parallel; the first set of arguments is filename prefixes and the second set is percentages % 

total_var=$(zcat ${1}.vcf.gz | grep -cv "#")
masked_var_flt=$(echo "${2}*$total_var/100" | bc)
masked_var=$(echo ${masked_var_flt%%.*})                      # take masked portion as integer

bcftools query -f "%ID\n" ${1}.vcf.gz > snp_id_${1}_${2}_${3}.txt      ##take IDs into a text file

shuf -n ${masked_var} snp_id_${1}_${2}_${3}.txt > masked_ids_${1}_${2}_${3}.txt   ## sample a random portion (equals the portion we calculated from the percentage)

bcftools view -e ID=@masked_ids_${1}_${2}_${3}.txt ${1}.vcf.gz -Oz -o Masked_${1}_${2}_${3}.vcf.gz    #remove the sampled variant ids

rm snp_id_${1}_${2}_${3}.txt

echo file: ${1}, mask: ${2}, repeat: ${3}
echo Masked variants: ${masked_var}
echo Remaining variants: $(zcat Masked_${1}_${2}_${3}.vcf.gz |grep -vc "#")