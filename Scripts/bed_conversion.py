### IMPORTANT: Open this script with Python 2.7 and not 3.x

import vcf
import sys


input_vcf = f"{sys.argv[1]}.vcf"
output_bed = f"{sys.argv[1]}.bed"

with open(output_bed, "w") as bed_file:
    vcf_reader = vcf.Reader(filename=input_vcf)
    for record in vcf_reader:
        bed_file.write(
            "{}\t{}\t{}\t{}\n".format(record.CHROM, record.POS-1, record.POS, record.ID)
        )