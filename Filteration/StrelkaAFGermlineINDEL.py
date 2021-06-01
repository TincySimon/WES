#!/bin/python3
import pysam

myvcf = pysam.VariantFile(snakemake.input[0], "r")
myvcf.header.formats.add("REF_AF", ".", "Float", "AF value of REF")
myvcf.header.formats.add("ALT_AF", ".", "Float", "AF value of ALT")
vcf_out = pysam.VariantFile('snakemake.output[0]', 'w', header=myvcf.header)
with open(snakemake.output[0], "a") as out:
	for variant in myvcf:
		REF_value = float(variant['AD'][0])/(float(variant['DPI'])
		ALT_value = float(variant['AD'][1])/(float(variant['DPI'])
		variant['REF_AF'] = REF_value
		variant['ALT_AF'] = ALT_value
	out.write(str(variant))


