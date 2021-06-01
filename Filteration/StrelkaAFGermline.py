import pysam

myvcf = pysam.VariantFile(snakemake.input[0], "r")
myvcf.header.formats.add("REF_AF", ".", "Float", "AF value of REF")
myvcf.header.formats.add("ALT_AF", ".", "Float", "AF value of ALT")
with open(snakemake.output[0], "a") as out:
	for variant in myvcf:
		REF_value = ''
		ALT_value = ''
		for sample in variant.samples:
			REF_value = float(variant.samples[sample]['AD'][0])/(float(variant.samples[sample]['DP'])
				variant.samples[sample]['REF_AF'] = REF_value
			ALT_value = float(variant.samples[sample]['AD'][1])/(float(variant.samples[sample]['DP'])
				variant.samples[sample]['ALT_AF'] = ALT_value
		out.write(str(variant))


