SAMPLES = ["PNET10P", "PNET19P", "PNET8P", "PNET20P", "PNET21P", "PNET22P", "PNET23P", "PNET24P"]

rule all:
	input:
          #expand("Data/CommonAlteration/InputFiles/{sample}/somatic.indels.vcf", sample = SAMPLES),
          #expand("output/Strelka/Results/{sample}/results/variants/somatic.indelsPASS.vcf",sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G_1000GFilter.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_w_AF.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_Indeltemp.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_biasRemovalIndel.vcf", sample = SAMPLES),
          expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf", sample = SAMPLES),
	  expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf.gz", sample = SAMPLES),
	  expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Indel_Strelka_InitialResults.tsv", sample = SAMPLES),
	  expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_WESNoHeaderIndel.tsv", sample = SAMPLES),
	  expand("Data/StrelkaVariants/Filteration/new/{sample}/{sample}_WES_FormattedIndel.tsv", sample = SAMPLES)

rule unzipIndelStrelka:
        input: "output/Strelka/Results/{sample}/results/variants/somatic.indels.vcf.gz"

        output: "Data/CommonAlteration/InputFiles/{sample}/somatic.indels.vcf"

        log: "log/unzip/{sample}indel_strelka_unzip.log"

        shell:
           r"""
               	module load SAMtools/1.3.1-foss-2015a
                (bgzip -d {input} > {output}) &> {log}
           """

rule PassStrelkaINDEL:
        input: "output/Strelka/Results/{sample}/results/variants/somatic.indels.vcf"

        output: "output/Strelka/Results/{sample}/results/variants/somatic.indelsPASS.vcf"

        log: "log/PASS/{sample}somatic.INDELPASS.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter "( na FILTER ) | (FILTER = 'PASS')" -f {input} > {output}) &> {log}
          """

rule IndelSegDup:
        input: VCF="output/Strelka/Results/{sample}/results/variants/somatic.indelsPASS.vcf"

        output: "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm.vcf"

        params: "/fast/projects/ngs_sers/Packages/Regions/Seg_Dup/GRCh37_SegDuplication_sort.bed"

        log: "log/SegDup/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm.vcf.log"

        shell:
          r"""
              	module load BEDTools/2.24.0-foss-2015a
                (bedtools intersect -v -header -a {input.VCF} -b {params} > {output}) &> {log}

          """
rule SNVAnn1000G:
        input:
          "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm.vcf"

        output:
          "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G.vcf"

        params: "/fast/projects/ngs_sers/Packages/Regions/1000G/1000G.EUR_AF.vcf.gz"
          #"/fast/projects/ngs_sers/Packages/Regions/WES/1000G_EUR_AF_Filter.vcf.gz"

        log:
          "log/StrelkaVariants/Filteration/{sample}/{sample}_Indel_1000G.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        annotate {params} {input} > {output}) &> {log}

          """
rule SNVGFilter:
        input:
          "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G.vcf"

        output:
          "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G_1000GFilter.vcf"

        log:
          "log/StrelkaVariants/Filteration/{sample}/{sample}_Indel_1000Gfilter.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter -n "(EUR_AF > 0.05)" {input} > {output}) &> {log}
          """

rule SNVAdAlleleFreqRef:
        input: "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_SegDup_rm_1000G_1000GFilter.vcf"

        output:
          "Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_w_AF.vcf"

        log:
          "log/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_w_AF.log"

        script:
          "StrelkaAFIndel.py"

rule addHeader:
        input:
          AVCF="/fast/projects/ngs_sers/2017_WES_Analysis/PNET/test.vcf",
          BVCF="Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_Strelka_w_AF.vcf"
        output:
          Header="Data/StrelkaVariants/Filteration/{sample}/{sample}_Indel_header.vcf",
          intersectfile="Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_Indeltemp.vcf"

        log: "log/Strelka/Header/{sample}_Indel_Header.log"

        shell:
          r"""
              	grep -P '^#' {input.AVCF} > {output.Header}
                (sed r {output.Header} {input.BVCF} > {output.intersectfile}) &> {log}

          """
rule SNVDKFZBiasFilter:
        input:
          bam = "output/bam/RecalibratedBamTumor/{sample}_recalibrated.bam",
          vcf = "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_Indeltemp.vcf"

        params:
          ref='/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa',
          tempFolder='Data/tmp'
        output:
          vcf = 'Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_biasRemovalIndel.vcf'

        log: 'log/Strelka/dkfzbiasfilter/{sample}.log'

        shell:
          r"""
              	(/fast/users/simont_c/scratch/miniconda3/envs/py27/bin/dkfzbiasfilter.py \
                        --tempFolder {params.tempFolder} \
                        --writeQC \
                        {input.vcf} {input.bam} {params.ref} {output.vcf}) &> {log}
          """


rule SNVfinalFilter:
        input: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Strelka_biasRemovalIndel.vcf"

        output: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf"

        log: "log/Strelka/{sample}_Strelka_InitialResults.log"

        shell:
          r"""
              	(java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        filter "( na FILTER ) | (FILTER = 'PASS')" -f {input} > {output}) &> {log}
          """

rule ZipInx:
	input: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf"

        output: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf.gz"

        log: "log/Strelka/ZipInx/{sample}_Indel_WES.vcf.log"

        shell:
          r"""
              	module load SAMtools/1.3.1-foss-2015a
                (bgzip -c {input} > {output}) &> {log}

          """
rule FriendlyFormall:
        input: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_StrelkaIndel.vcf"

        output: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Indel_Strelka_InitialResults.tsv"

        log: "log/Strelka/{sample}_Indel_csv.log"


        shell:
          r"""
              	(cat {input} \
                | /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/scripts/vcfEffOnePerLine.pl \
                | java -Xmx10g -jar /fast/projects/ngs_sers/Packages/tools/snpEff_latest-version/snpEff/SnpSift.jar \
                        extractFields - -s "," -e "-" CHROM POS ID REF ALT QUAL FILTER \
                        "GEN[*].REF_AF" "GEN[*].ALT_AF" > {output}) &> {log}
          """

rule TSVtoBed:
	input: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_Indel_Strelka_InitialResults.tsv"

        output: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_WESNoHeaderIndel.tsv"

        log: "log/filteration/Strelka/{sample}_tsvnoheaderIndel.log"

        shell:
          r"""
              	(awk -F' ' '{{if (NR!=1) {{print}}}}' {input} > {output}) &> {log}

          """

rule NewTSVtoBed:
        input: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_WESNoHeaderIndel.tsv"

        output: "Data/StrelkaVariants/Filteration/new/{sample}/{sample}_WES_FormattedIndel.tsv"

        log: "log/filteration/Strelka/{sample}_tsvformatedIndel.log"

        shell:
          r"""

              	(awk -F'\t' '{{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$4,$5,$7,$8,$9)}}' < {input} > {output}) &> {log}

          """



