#!/bin/bash

echo "Determine file"
read B

D=$(zcat $B | grep -vc "#")

while read P; do
bcftools +setGT $B -- -t r:"$P" -n . --seed 100 > output.vcf
bgzip -f output.vcf
bcftools +fill-tags output.vcf.gz -Oz -o tagged.vcf.gz -- -t F_MISSING
bcftools view -e "INFO/F_MISSING>0.1" tagged.vcf.gz -Oz -o final.vcf.gz

C=$(zcat final.vcf.gz | grep -vc "#")
E=$(expr $D - $C)
echo "$P"
bc <<< "scale=4; $E*100/$D"
done < perc.txt
