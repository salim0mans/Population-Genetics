#!/bin/bash


echo "Determine filename (only prefix)"
read filename
echo "Determine masking percentage %"
read perc

total_var=$(zcat ${filename}.vcf.gz | grep -cv "#")
masked_var_flt=$(echo "$perc*$total_var/100" | bc)
masked_var=$(echo ${masked_var_flt%%.*})                                     # take masked portion as integer

bcftools query -f "%ID\n" ${filename}.vcf.gz > snp_id_${filename}.txt        #take IDs into a text file

shuf -n ${masked_var} snp_id_${filename}.txt > masked_ids_${filename}.txt    # sample a random portion (equals the portion we calculated from the percentage)

bcftools view -e ID==@masked_ids_${filename}.txt ${filename}.vcf.gz -Oz -o Masked_${filename}.vcf.gz    #remove the sampled variant ids

rm snp_id_${filename}.txt

echo Masked variants: ${masked_var}
echo Remaining variants: $(zcat Masked_${filename}.vcf.gz |grep -vc "#")