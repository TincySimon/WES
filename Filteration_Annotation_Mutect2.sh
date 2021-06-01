#Mutect2_OneByTwo

SAMPLES = ["PNET19P","PNET21P","PNET23P","PNET24P","PNET9P","PNET10P","PNET20P","PNET22P","PNET8P"]
rule all:
	input:
          expand("Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm.vcf", sample=SAMPLES),
          expand("Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G.vcf", sample=SAMPLES),
          expand("Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G_1000GFilter.vcf", sample=SAMPLES),
          expand("Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_ExtractFields.vcf", sample=SAMPLES),
	  expand("Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields.vcf", sample=SAMPLES),
	  expand("Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields_New.vcf", sample=SAMPLES),
	  expand("Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FilterFile.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_SAFfilter.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFiltered.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFilteredIntersect.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/InitialFilter/{sample}/{sample}_FPFiltered.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf.gz", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf.gz", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_AD_filtered.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_PredAnn.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_biasRemoval.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_final_filtered.vcf", sample=SAMPLES),
	  expand("Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/CSV/{sample}/{sample}_WES.csv", sample=SAMPLES)

rule SegDup:
	input: VCF="output/Mutect2/LowStringency/{sample}/final/{sample}_sorted.vcf.gz"

        output: "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm.vcf"

        params: "/fast/projects/ngs_sers/Packages/Regions/Seg_Dup/GRCh37_SegDuplication_sort.bed"

        log: "log/SegDup/Mutect2/Normal_filtered/LowStringency/{sample}_sorted_SegDup_rm.log"

        shell:
          r"""
              	module load BEDTools/2.24.0-foss-2015a
                (bedtools intersect -v -header -a {input.VCF} -b {params} > {output}) &> {log}

          """

rule Ann1000G:
	input:
          "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm.vcf"

	output:
          "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G.vcf"

	params: "/fast/projects/ngs_sers/Packages/Regions/1000G/1000G.EUR_AF.vcf.gz"
          #"/fast/projects/ngs_sers/Packages/Regions/WES/1000G_EUR_AF_Filter.vcf.gz"

	log:
          "log/filteration/Mutect2/Normal_filtered/LowStringency/{sample}_1000G.log"

	shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        annotate {params} {input} > {output}) &> {log}

          """

rule GFilter:
	input:
          "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G.vcf"

	output:
          "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G_1000GFilter.vcf"

	log:
          "log/filteration/Mutect2/Normal_filtered/LowStringency/{sample}_1000Gfilter.log"

	shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter -n "(EUR_AF > 0.05)" {input} > {output}) &> {log}
          """

rule ProperForm:
	input: "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G_1000GFilter.vcf"

	output: "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_ExtractFields.vcf"

	log: "log/filteration/Mutect2/Normal_filtered/LowStringency/{sample}_ExtractFields.log"

	shell:
          r"""
              	grep -P '^##' {input} > {output}
                (java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        extractFields -s "\t" -e "-" {input} CHROM POS ID REF ALT QUAL FILTER \
                        ECNT HCNT MAX_ED MIN_ED NLOD TLOD EUR_AF \
                        "GEN[0].GT" "GEN[1].GT" "GEN[0].AD" "GEN[1].AD" "GEN[0].AF" "GEN[1].AF" "GEN[0].ALT_F1R2" "GEN[1].ALT_F1R2" "GEN[0].ALT_F2R1" "GEN[1].ALT_F2R1" "GEN[0].REF_F1R2" "GEN[1].REF_F2R1" >> {output}) &> {log}

	  """


rule Mutect2FPFilter:
	input: "Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_ExtractFields.vcf"

	output: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_ExtractFields.log"
	
	shell:
	  r"""
		sed 's/[, ]\+/\t/g; /^#/ d' {input} > {output}

	  """

rule RemoveHeader:
        input: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields.vcf"

        output: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields_New.vcf"

        log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_ExtractFields_new.log"

        shell:
          r"""
              	sed '1d' {input} > {output}

          """

rule CalculateSAF:
	input: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_ExtractFields_New.vcf"

	output: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FilterFile.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_FilterFile.log"
		
	shell:
	  r"""
		(awk 'BEGIN {{OFS="\t"}} {{ if ($22 > $21){{ $29=1-$21/$22 }}else{{ $29=0 }} print }}' {input} > {output}) &> {log}

	  """

rule SAFfilter:
	input: "Data/Mutect2/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FilterFile.vcf" 
	
	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_SAFfilter.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_SAFfilter.log"

	shell:
	  r"""
		(awk 'BEGIN {{OFS="\t"}} {{ if ($29 > 0.50) {{ print }} }}' {input} > {output}) &> {log}

	  """


rule AFTumorFilter:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_SAFfilter.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFiltered.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_FPFiltered.log"

	shell:
	  r"""
              	(awk 'BEGIN {{OFS="\t"}} {{ if ($22 > 0.03) {{ print }} }}' {input} > {output}) &> {log}
          """

rule addHeader:
	input:
          AVCF="Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G_1000GFilter.vcf",
	  BVCF="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFiltered.vcf"
	output:
	  Header="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_header.vcf",
	  intersectfile="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFilteredIntersect.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_Header.log"

	shell:
	  r"""
		grep -P '^#' {input.AVCF} > {output.Header}
		(sed r {output.Header} {input.BVCF} > {output.intersectfile}) &> {log}

	  """

rule IntersectVCF:
	input:
	  AVCF="Data/Mutect2/Normal_filtered/LowStringency/Filteration/{sample}/{sample}_sorted_SegDup_rm_1000G_1000GFilter.vcf",
	  BVCF="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/FPFIltering/{sample}/{sample}_FPFilteredIntersect.vcf"

	output:
	  "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/InitialFilter/{sample}/{sample}_FPFiltered.vcf"

	log: "log/FPFIltering/Mutect2/Normal_filtered/LowStringency/{sample}_FPFiltered_new.log"
	
        shell:
          r"""
              	module load BEDTools/2.24.0-foss-2015a
                (bedtools intersect -header -wa -a {input.AVCF} -b {input.BVCF} > {output}) &> {log}

          """

rule AnnotationSnpEff:
        input:
          VCF_old="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/InitialFilter/{sample}/{sample}_FPFiltered.vcf"

        output:
          VCF="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf",
          HTML="Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/InitialFilter/{sample}/{sample}_Normal_filtered_ann.html"

        params: "/fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/snpEff.config"

        log: "log/SnpEff/Mutect2_OneByTwo/Normal_filtered/LowStringency/{sample}_SegDupFiltered_Ann.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/snpEff.jar \
                        -c {params} \
                        -s {output.HTML} \
                        -canon \
                        -onlyProtein \
                        GRCh37.75 \
                        {input.VCF_old} > {output.VCF}) &> {log}
      """

rule ZipInx:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf.gz"

	log: "log/ZipInx/{sample}_WES.vcf.log"

	shell:
	  r"""
		module load SAMtools/1.3.1-foss-2015a
		(bgzip -c {input} > {output}) &> {log}

	  """

rule PASS:
	input:
          "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES.vcf"

	output:
          "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf"

	log:
          "log/filteration/Mutect2/Normal_filtered/LowStringency/{sample}_PASSfilter.log"

	shell:
          r"""
                (java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter "( na FILTER ) | (FILTER = 'PASS') | (FILTER = 'homologous_mapping_event') | (FILTER = 'clustered_events') | (FILTER = 'clustered_events;homologous_mapping_event')" -f {input} > {output}) &> {log}
          """

rule ZipInxNew:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf.gz"

	log: "log/ZipInx/{sample}_WES_PASS.vcf.log"

	shell:
          r"""
              	module load SAMtools/1.3.1-foss-2015a
		(bgzip -c {input} > {output}) &> {log}

          """
rule ADFiltered:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_AD_filtered.vcf"

	log: "log/filteration/Mutect2/Normal_filtered/{sample}_PASS_AD_filtered.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
			filter "( GEN[TUMOR].AD[1] > 2 )"  -f {input} > {output}) &> {log}
          """


rule dbNSFP:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_AD_filtered.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_PredAnn.vcf"

	log: "log/filteration/Mutect2/Normal_filtered/{sample}_PredAnn.log"

	params: "/fast/projects/ngs_sers/Packages/Regions/WES/dbNSFP2.9.txt.gz"

	shell: 
	  r"""
		(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
			dbnsfp -v -m -db {params} \
			-f Uniprot_acc,Interpro_domain,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,CADD_phred,CADD_raw_rankscore,COSMIC_ID \
			{input} > {output}) &> {log}

	  """

rule DKFZBiasFilter:
	input:
    	  bam = "output/bam/RecalibratedBam/{sample}_recalibrated.bam",
	  vcf = "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_PASS_PredAnn.vcf"

	params:
	  ref='/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa',
	  tempFolder='Data/tmp'
	output:
	  vcf = 'Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_biasRemoval.vcf'

	log: 'log/dkfzbiasfilter/{sample}.log'

	shell:
	  r"""
		(/fast/users/simont_c/scratch/miniconda3/envs/py27/bin/dkfzbiasfilter.py \
                	--tempFolder {params.tempFolder} \
                	--writeQC \
                	{input.vcf} {input.bam} {params.ref} {output.vcf}) &> {log}

      """

rule FinalFilter:
	input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_biasRemoval.vcf"

	output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_final_filtered.vcf"

	log: "log/filteration/Mutect2/Normal_filtered/{sample}_final_filtered.log"
	
        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter "( na FILTER ) | (FILTER = 'PASS') | (FILTER = 'homologous_mapping_event') | (FILTER = 'clustered_events') | (FILTER = 'clustered_events;homologous_mapping_event')" -f {input} > {output}) &> {log}
          """


rule FriendlyFormall:
        input: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/{sample}/{sample}_WES_final_filtered.vcf"

        output: "Data/Mutect2_OneByTwo/Normal_filtered/LowStringency/Results/Final/CSV/{sample}/{sample}_WES.csv"

        log: "log/filteration/Mutect2/Normal_filtered/{sample}_final_new.log"


        shell:
          r"""
              	(cat {input} \
		| /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/scripts/vcfEffOnePerLine.pl \
		| java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        extractFields - -s "," -e "-" CHROM POS ID REF ALT QUAL FILTER \
                        DB ECNT HCNT MAX_ED MIN_ED NLOD TLOD \
                        "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].FEATURE" "ANN[*].RANK" "ANN[*].HGVS_P" \
                        "GEN[*].GT" "GEN[*].AD" "GEN[*].AF" "GEN[*].FOXOG" "LOF[*].GENE" "LOF[*].NUMTR" "LOF[*].PERC" \
                        "NMD[*].GENE" "NMD[*].NUMTR" "NMD[*].PERC" \
                        dbNSFP_COSMIC_ID dbNSFP_CADD_phred dbNSFP_SIFT_pred dbNSFP_Polyphen2_HDIV_pred dbNSFP_Polyphen2_HVAR_pred dbNSFP_LRT_pred \
                        dbNSFP_MutationTaster_pred dbNSFP_MutationAssessor_pred dbNSFP_Interpro_domain > {output}) &> {log}
          """


