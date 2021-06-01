# WES - whole exome sequencing of Pancreatic Neuroendocrine Neoplasms (PanNENs)
Steps:
* Preprocessing
* Variant calling
* Filterations
* Annotation
* CNV analysis

Prepared in SnakeMake workflow management system

1. Preprocessing:
Files found under Preprocessing folder includes the shell script for preprocessing and json file for executing the workflow in Grid engine scheduler in HPC as well as quality check associated files

2. Variantcalling:
Files found under VariantCalling folder includes the shell script for somatic and germline calling from Strelka 

3. Filterations and Annotation:
Workflow for filteration and annotation of the mutations for identifying true pathogenic variants are found under Filteration and Annotation folder

4. CNV Analysis
CNV analysis using CNVKit and subsequent scatterplot representation of data

 
