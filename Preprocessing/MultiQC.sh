#!/bin/bash

#need to have the conda root and be in the directory to run the multiqc: so run the 2 lines below (without the hash) in the terminal before running the script:
#source /fast/users/simont_c/scratch/miniconda3/bin/activate
#cd /fast/projects/ngs_sers/2017_WES_Analysis/PNET/TechnicalAnalysis/Qualimap/ 

multiqc . --ignore *qualimapReportOutsideRegions*


