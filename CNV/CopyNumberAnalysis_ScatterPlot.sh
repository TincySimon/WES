#/fast/projects/ngs_sers/2017_WES_Analysis/AC/Data/CNVkit/Results/
#import sys
#import glob

#INFILES = glob.glob('/fast/projects/ngs_sers/panel_seq/Soulafa/AC/CNVKIT/AC/*.bam')
#SAMPLES = [infile.split('/')[-1].rstrip('\.bam') for infile in INFILES]
#TYPE=["a","b"]

SAMPLES = ["PNET19P","PNET21P","PNET23P","PNET24P","PNET9P","PNET10P","PNET20P","PNET22P","PNET8P"]
rule all:
	input:
          expand("Data/CNVkit/Graph/{sample}_scatter.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/{sample}_scatter_Flasso.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_ABL1.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_ABL2.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_RICTOR.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_MEN1.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_TSC1.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_TSC2.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_ABL2.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_RICTOR.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_MEN1.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_TSC1.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_TSC2.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_ABL1.pdf", sample=SAMPLES),
	  expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_chr11gene.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_chr11gene.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_NVL.pdf", sample=SAMPLES),
          expand("Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_NVL.pdf", sample=SAMPLES)
	  

rule scatter:
	input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/{sample}_scatter.pdf"

        log: "log/CNVkit/{sample}_scatter.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -o {output}) &> {log}

           """

rule scatterFlasso:
	input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

	output: "Data/CNVkit/Graph/{sample}_scatter_Flasso.pdf"

        log: "log/CNVkit/{sample}_scatter.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -o {output}) &> {log}

           """

rule scatterABL1CBS:
	input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_ABL1.pdf"

        log: "log/CNVkit/{sample}_CBS_scatter_ABL1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 9 -g ABL1 -o {output}) &> {log}

           """

rule scatterABL1FF:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_ABL1.pdf"

        log: "log/CNVkit/{sample}_Flasso_scatter_ABL1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 9 -g ABL1 -o {output}) &> {log}

           """

rule scatterABL2CBS:
	input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_ABL2.pdf"

        log: "log/CNVkit/{sample}_CBS_scatter_ABL2.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 1 -g ABL2 -o {output}) &> {log}

           """

rule scatterABL2FF:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_ABL2.pdf"

        log: "log/CNVkit/{sample}_Flasso_scatter_ABL2.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 1 -g ABL2 -o {output}) &> {log}

           """


rule scatterRICTORCBS:
	input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_RICTOR.pdf"

        log: "log/CNVkit/{sample}_scatter_RICTOR.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 5 -g RICTOR -o {output}) &> {log}

           """

rule scatterRICTORFlasso:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_RICTOR.pdf"

        log: "log/CNVkit/{sample}_scatter_RICTOR.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 5 -g RICTOR -o {output}) &> {log}

           """

rule scatterMEN1:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_MEN1.pdf"

        log: "log/CNVkit/{sample}_scatter_MEN1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 11 -g MEN1 -o {output}) &> {log}

           """

rule scatterTSC1:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_TSC1.pdf"

        log: "log/CNVkit/{sample}_scatter_TSC1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 9 -g TSC1 -o {output}) &> {log}

           """

rule scatterTSC2:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

	output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_TSC2.pdf"

        log: "log/CNVkit/{sample}_scatter_TSC2.log"

        shell:
           r"""

                (/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 16 -g TSC2 -o {output}) &> {log}

           """

rule scatterMEN1FF:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_MEN1.pdf"

        log: "log/CNVkit/{sample}_scatter_MEN1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 11 -g MEN1 -o {output}) &> {log}

           """

rule scatterTSC1FF:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_TSC1.pdf"

        log: "log/CNVkit/{sample}_scatter_TSC1.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 9 -g TSC1 -o {output}) &> {log}

           """

rule scatterTSC2FF:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_TSC2.pdf"

        log: "log/CNVkit/{sample}_scatter_TSC2.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 16 -g TSC2 -o {output}) &> {log}

           """

rule scatterchr11geneCBS:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_chr11gene.pdf"

        log: "log/CNVkit/{sample}_scatter_chr11gene.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 11 -g TSG101,SF1,WDR74,MEN1 -o {output}) &> {log}

           """

rule scatterchr11geneFlasso:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_chr11gene.pdf"

        log: "log/CNVkit/{sample}_scatter_chr11gene.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 11 -g TSG101,SF1,WDR74,MEN1 -o {output}) &> {log}

           """

rule scatterNVLCBS:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/CBS/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/CBS/{sample}_scatter_NVL.pdf"

        log: "log/CNVkit/{sample}_scatter_NVL.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 1 -g NVL -o {output}) &> {log}

           """

rule scatterNVLFlasso:
        input:
           CNR="Data/CNVkit/Result/{sample}.cnr",
           CNS="Data/CNVkit/Result/Flasso/{sample}.cns"

        output: "Data/CNVkit/Graph/GeneChr/Flasso/{sample}_scatter_NVL.pdf"

        log: "log/CNVkit/{sample}_scatter_NVL.log"

        shell:
           r"""

               	(/fast/users/simont_c/scratch/miniconda3/bin/cnvkit.py scatter {input.CNR} -s {input.CNS} -c 1 -g NVL -o {output}) &> {log}

           """


