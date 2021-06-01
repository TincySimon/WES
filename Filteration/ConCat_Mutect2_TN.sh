#expand('output/Mutect2/{sample}P/{sample}P_{panel}.vcf.gz', sample=SAMPLES, panel=PANEL)
SAMPLES = ["PNET19P","PNET21P","PNET23P","PNET24P","PNET9P","PNET10P","PNET20P","PNET22P","PNET8P"]
#python modules
import sys
import glob

INFILES = glob.glob('/fast/projects/ngs_sers/Packages/Regions/WES/scatter_region/*.bed')
PANEL = [infile.split('/')[-1].rsplit( ".", 1 )[ 0 ] for infile in INFILES]

#SAMPLES = [33]
#TYPE = ["a"]

rule all:
  input:
      expand('output/Mutect2/LowStringency/{sample}/final/{sample}.vcf.gz', sample=SAMPLES),
      expand('output/Mutect2/LowStringency/{sample}/final/{sample}.vcf.gz.done', sample=SAMPLES),
      expand('output/Mutect2/LowStringency/{sample}/final/{sample}_sorted.vcf.gz', sample=SAMPLES)

def generate_mergeInput(wildcards):
    return ['output/Mutect2/{sample}/{sample}_{panel}.vcf.gz'.format(sample = wildcards.sample, panel = currrent_panel) for currrent_panel in PANEL]


rule Mutect2_Concat:
  input: 
    generate_mergeInput

  output:
      VCF = 'output/Mutect2/LowStringency/{sample}/final/{sample}.vcf.gz',
      done = 'output/Mutect2/LowStringency/{sample}/final/{sample}.vcf.gz.done'

  log: 'log/Mutect2/LowStringency/BCFtool_Concat/{sample}_VCF_concat.log'

  shell:
      r"""

          module load BCFtools
	
          (bcftools concat -a -o {output.VCF} -O z {input} --threads 8) &> {log}

          touch {output.done}

      """

rule sort_VCF:
  input:
    'output/Mutect2/LowStringency/{sample}/final/{sample}.vcf.gz'

  output:
    sortVCF = 'output/Mutect2/LowStringency/{sample}/final/{sample}_sorted.vcf.gz'

  log: 'log/Mutect2/LowStringency/sort_VCF/{sample}_sort.log'

  shell:
      r"""
          (java -jar /fast/users/simont_c/scratch/miniconda3/pkgs/picard-2.17.6-py36_0/share/picard-2.17.6-0/picard.jar SortVcf \
              INPUT={input} \
              OUTPUT={output.sortVCF}) &> {log}

      """




