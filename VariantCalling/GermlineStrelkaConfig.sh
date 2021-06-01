#!/bin/bash

#${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
#--bam NA12878.bam \
#--referenceFasta hg19.fa \
#--runDir ${STRELKA_ANALYSIS_PATH}

#configfileworkflow: /fast/users/simont_c/scratch/miniconda3/envs/py27/bin/configureStrelkaGermlineWorkflow.py

IN='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/output/bam/RecalibratedBamNormal'
OUT='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/output/Strelka/Results'
REF='/fast/projects/cubit/current/static_data/reference/GRCh37/hs37d5/hs37d5.fa'

for file in $( find $IN/* -name '*.bam' ); do
        f=${file##*/}
        #echo $f
        F2=${f/_recalibrated.bam/}
        #echo $IN/$F2
        OUTPATH=$OUT/${f%%_*}
        #echo $OUTPATH
	/fast/users/simont_c/scratch/miniconda3/envs/py27/bin/configureStrelkaGermlineWorkflow.py --bam $file --referenceFasta $REF --runDir $OUTPATH
done

