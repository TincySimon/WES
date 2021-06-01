#Preprocesing and generating bam and bai files from fastq files for PNET samples multiplexed onto two lanes
#installed Anaconda/miniconda

LANE = ["L1", "L2"]

SAMPLES = ["PNET19N", "PNET19P", "PNET21N", "PNET21P", "PNET23N", "PNET23P", "PNET24N", "PNET24P", "PNET9N", "PNET9P"] 

rule all:
     input: 
        expand("output/bam/{sample}/{sample}_{lane}.bam", sample=SAMPLES, lane=LANE), 
	expand("output/bam/{sample}/{sample}_merged.bam", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_sorted.bam", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_dedup.bam", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_dedup.bam.bai", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_sorted.bam.bai", sample=SAMPLES),
	expand("output/realin_inter/ReAlignerTaretCreator_{sample}.intervals", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_realigned.bam", sample=SAMPLES),
	expand("output/BQSR/{sample}/{sample}_CovarTable_Recal.table", sample=SAMPLES),
	expand("output/bam/{sample}/{sample}_recalibrated.bam", sample=SAMPLES)

rule BWA_Align:
    input: 
        READONE = "incoming/{sample}/{sample}_{lane}_R1.fastq.gz",
        READTWO = "incoming/{sample}/{sample}_{lane}_R2.fastq.gz"
    
    params:
        #RG= "@RG\tID:{lane}\tSM:{sample}\tPL:illumina\tLB:AgilentSureSelect\tPU:2000",    
        BWAREF = "/fast/projects/cubit/current/static_data/precomputed/BWA/0.7.15/GRCh37/hs37d5/hs37d5.fa" #indexed?
    
    output: "output/bam/{sample}/{sample}_{lane}.bam"

    log: "log/BWA/{sample}/{sample}_{lane}.log"

    benchmark: "benchmarks/{sample}_{lane}_BWA_Align_benchmark.txt"

    shell: 
      r"""  
      module load BWA/0.7.15-foss-2015a 
      module load SAMtools/1.3.1-foss-2015a
             
      (bwa mem -M \
	  -R "@RG\tID:{wildcards.lane}\tSM:{wildcards.sample}\tPL:illumina\tLB:AgilentSureSelect\tPU:2000" \
	  {params.BWAREF} \
	  {input.READONE} {input.READTWO} \
      | samtools view -b -S -o {output}) &> {log}
      """

def generate_mergeInput(wildcards):
    return ["output/bam/{sample}/{sample}_{lane}.bam".format(sample = wildcards.sample, lane = current_lane) for current_lane in LANE]

rule MergeBam:
   input: 
       generate_mergeInput

   output: "output/bam/{sample}/{sample}_merged.bam"
   
   log: "log/bam_merge/{sample}/{sample}_merged.log"

   benchmark: "benchmarks/{sample}_MergeBam_benchmark.txt"

   shell:
     r"""
     INPUTS=$(echo {input} | sed -e 's/output/ I=output/g')
     	
     (java -jar /fast/users/stincy/scratch/miniconda/pkgs/picard-2.11.0-py35_0/share/picard-2.11.0-0/picard.jar MergeSamFiles \
      $INPUTS \
      O={output}) &> {log}
     """

rule bam_sort:
    input: "output/bam/{sample}/{sample}_merged.bam"

    output: "output/bam/{sample}/{sample}_sorted.bam"

    log: "log/bamsort/{sample}/{sample}_sort.log"

    benchmark: "benchmarks/{sample}_bam_sort_benchmark.txt"

    shell: 
      r"""
      module load SAMtools/1.3.1-foss-2015a
      (samtools sort -o {output} {input}) &> {log}
      
      """
rule ReadDuplicates:
    input: "output/bam/{sample}/{sample}_sorted.bam"

    output: 
        RD = "output/bam/{sample}/{sample}_dedup.bam",
        MET = "log/read_dup/{sample}/metrics_{sample}.txt"

    log: "log/read_dup/{sample}/{sample}_rmdup_new.log"

    benchmark: "benchmarks/{sample}_ReadDuplicates_benchmark.txt"

    shell:
      r"""
      (java -jar /fast/users/stincy/scratch/miniconda/pkgs/picard-2.11.0-py35_0/share/picard-2.11.0-0/picard.jar MarkDuplicates \
             I={input} \
             O={output.RD} \
             M={output.MET} \
             AS=true) &> {log}
      """

rule bamindexInitial:
    input: "output/bam/{sample}/{sample}_sorted.bam"

    output: "output/bam/{sample}/{sample}_sorted.bam.bai"

    log: "log/BAI/{sample}/{sample}_baifile.log"

    benchmark: "benchmarks/{sample}_BAM_index_benchmark.txt"

    shell:
      r"""
      module load SAMtools/1.3.1-foss-2015a
      samtools index -b {input} &> {log}
      
      """
rule bamindexInitial_dedup:
    input: "output/bam/{sample}/{sample}_dedup.bam"

    output: "output/bam/{sample}/{sample}_dedup.bam.bai"

    log: "log/BAI/{sample}/{sample}_baifile_dedup.log"

    shell:
      r"""
      module load SAMtools/1.3.1-foss-2015a
      samtools index -b {input} &> {log}

      """


rule indelReAlignerTargetCreator:
    input: 
        bam = "output/bam/{sample}/{sample}_dedup.bam",
        bai = "output/bam/{sample}/{sample}_dedup.bam.bai"
    
    params: 
        REF = "/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa"

    output: "output/realin_inter/ReAlignerTaretCreator_{sample}.intervals"
  
    log: "log/ReAlignerTaretCreator/{sample}/ReAlignerTaretCreator_{sample}.log"

    benchmark: "benchmarks/{sample}_indelRealignerTargetCreator_benchmark.txt"

    shell:
      r"""
      
      (java -jar /fast/projects/ngs_sers/Packages/tools/GATK3.7/GenomeAnalysisTK.jar \
           -T RealignerTargetCreator \
           -R {params.REF} \
           -I {input.bam} \
           -o {output}) &> {log}
      """

rule indelReAligner:
    input:
        IR = "output/bam/{sample}/{sample}_dedup.bam",
        TI= "output/realin_inter/ReAlignerTaretCreator_{sample}.intervals"
    params: 
        reference= "/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa",
        #TI= "output/realin_inter/ReAlignerTaretCreator_{sample}.intervals"

    output: "output/bam/{sample}/{sample}_realigned.bam"

    log: "log/realign/{sample}/{sample}_realignment.log"

    benchmark: "benchmarks/{sample}_indel_ReAligner_benchmark.txt"

    shell:
      r"""
      (java -jar /fast/projects/ngs_sers/Packages/tools/GATK3.7/GenomeAnalysisTK.jar \
           -T IndelRealigner \
           -R {params.reference} \
           -I {input.IR} \
           -targetIntervals {input.TI} \
           -o {output}) &> {log}
      """ 

rule BQSR: 
    input: "output/bam/{sample}/{sample}_realigned.bam"

    params: 
        REF = "/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa",
	KWN = "/fast/projects/cubit/current/static_data/db/dbSNP/b147/GRCh37/All_20160408.vcf.gz"
    
    output: "output/BQSR/{sample}/{sample}_CovarTable_Recal.table"

    log: "log/BQSR/{sample}/{sample}_BQSR.log"

    benchmark: "benchmarks/{sample}_BQSR_benchmark.txt"

    shell:
      r"""
      (java -jar /fast/projects/ngs_sers/Packages/tools/GATK3.7/GenomeAnalysisTK.jar \
           -T BaseRecalibrator \
           -R {params.REF} \
           -I {input} \
           -knownSites {params.KWN} \
           -o {output}) &> {log}
      """

rule RecalibratedBam: 
    input: 
        RE = "output/bam/{sample}/{sample}_realigned.bam",
        BQSR = "output/BQSR/{sample}/{sample}_CovarTable_Recal.table"

    params:
        REF = "/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa"
	#BQSR = "output/BQSR/{sample}/{sample}_CovarTable_Recal.grp"

    output: "output/bam/{sample}/{sample}_recalibrated.bam"        

    log: "log/Recalibrate/{sample}/{sample}_Recalibrated.log"
    
    benchmark: "benchmarks/{sample}_Recalibrated_Bam_benchmark.txt"

    shell:
      r"""

      (java -jar /fast/projects/ngs_sers/Packages/tools/GATK3.7/GenomeAnalysisTK.jar \
           -T PrintReads \
           -R {params.REF} \
           -I {input.RE} \
           -BQSR {input.BQSR} \
           -o {output}) &> {log}
      """
