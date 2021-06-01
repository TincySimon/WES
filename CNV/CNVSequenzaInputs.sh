#The sequenza-utils command provides various tools; here we highlight only the basic usage:

#Process a FASTA file to produce a GC Wiggle track file:
#sequenza−utils gc_wiggle −w 50 --fasta hg19.fa -o hg19.gc50Base.wig.gz
#Process BAM and Wiggle files to produce a seqz file:
#sequenza−utils bam2seqz -n normal.bam -t tumor.bam --fasta hg19.fa \
    #-gc hg19.gc50Base.wig.gz -o out.seqz.gz
#Post-process by binning the original seqz file:
#sequenza−utils seqz_binning --seqz out.seqz.gz -w 50 -o out small.seqz.gz

SAMPLES = ["PNET19","PNET21","PNET23","PNET24","PNET10","PNET20","PNET22","PNET8"]

rule all:
	input:
	  expand('Data/Sequenza/IntermediateFIles/{sample}P_out.seqz.gz', sample=SAMPLES),
	  expand('Data/Sequenza/IntermediateFIles/{sample}P_small.seqz.gz', sample=SAMPLES),
	  expand('Data/Sequenza/IntermediateFIles/{sample}P_small.seqz', sample=SAMPLES), 
	  expand('Data/Sequenza/IntermediateFIles/{sample}P_AS_small.seqz', sample=SAMPLES),
	  expand('Data/Sequenza/IntermediateFIles/{sample}P_AS_small.seqz.gz', sample=SAMPLES) 

rule Seq:
	input:
	  T = "output/bam/{sample}P/{sample}P_recalibrated.bam",
	  N = "output/bam/{sample}N/{sample}N_recalibrated.bam"

	params:
	  reference='/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa',
	  gc='/fast/projects/ngs_sers/Packages/Regions/WES/hg19.gc50Base.wig.gz'

	output:
	  fileA = 'Data/Sequenza/IntermediateFIles/{sample}P_out.seqz.gz'

	log: 'log/Sequenza/{sample}P_out.seqz.log'

	shell:
	  r"""
		(/fast/users/simont_c/scratch/miniconda3/envs/py27/bin/sequenza-utils bam2seqz -n {input.N} -t {input.T} \
		  --fasta {params.reference} \
		  -gc {params.gc} \
		  -o {output.fileA} ) &> {log}

	  """

rule smallseq:
	input: 'Data/Sequenza/IntermediateFIles/{sample}P_out.seqz.gz'

	output: 'Data/Sequenza/IntermediateFIles/{sample}P_small.seqz.gz'

	log: 'log/Sequenza/{sample}P_small.seqz.log'

	shell:
	  r"""
		(/fast/users/simont_c/scratch/miniconda3/envs/py27/bin/sequenza-utils seqz_binning -s {input} -w 50 -o {output}) &> {log}

	  """
rule unzipsmallseq:
	input: 'Data/Sequenza/IntermediateFIles/{sample}P_small.seqz.gz'

	output: 'Data/Sequenza/IntermediateFIles/{sample}P_small.seqz'

	log: 'log/Sequenza/{sample}P_small.seqz_unzip.log'
	
	shell:
	  r"""
		(zcat {input} > {output}) &> {log}

	  """

rule ExctractAutosomeSex:
	input: 'Data/Sequenza/IntermediateFIles/{sample}P_small.seqz'

	output: 'Data/Sequenza/IntermediateFIles/{sample}P_AS_small.seqz'

	log: 'log/Sequenza/{sample}P_Extract.log'

	shell:
	  r"""
		(awk -F"\t" '{{ if (!($1 == "GL000229.1" || $1 == "GL000231.1" || $1 == "MT" || $1 == "GL000226.1")) {{ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14 }} }}' {input} > {output}) &> {log}

	  """

rule zip:
	input: 'Data/Sequenza/IntermediateFIles/{sample}P_AS_small.seqz'

	output: 'Data/Sequenza/IntermediateFIles/{sample}P_AS_small.seqz.gz'

	log: 'log/Sequenza/{sample}P_Extractzip.log'

	shell:
	  r"""
		module load SAMtools/1.3.1-foss-2015a
		(bgzip -c {input} > {output}) &> {log}

	  """

