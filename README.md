# WES - whole exome sequencing of Pancreatic Neuroendocrine Neoplasms (PanNENs)
Steps:
* Preprocessing
* Variant calling
* Filterations
* Annotation

Prepared in SnakeMake workflow management system

1. Preprocessing:
Files found under Preprocessing folder includes the shell script for preprocessing and json file for executing the workflow in Grid engine scheduler in HPC

2. Variantcalling:
Files found under VariantCalling folder includes the shell script for somatic and germline calling from Strelka 

3. Filterations:
Workflow for filteration of the mutations for identifying true pathogenic variants are found under Filteration folder

4. Annotation:
Annotation of the mutational variants are determined through the Annovar package and can be found under the folder Annotations


 
