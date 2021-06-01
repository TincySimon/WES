module load Java/1.8.0_92

IN='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/output/bam/RecalibratedBamNormal'
#OUT='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/TechnicalAnalysis/New_TechnicalAnalysis/Qualimap/
OUT='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/TechnicalAnalysis/Qualimap'

for file in $( find $IN/* -name '*.bam' ); do
	f=${file##*/}
        echo $f
        F2=${f/_recalibrated.bam/}
        echo $F2
        OUTPATH=$OUT/${f%%_*}
	echo $OUTPATH
	echo $IN/$f
	/fast/projects/ngs_sers/Packages/tools/qualimap_v2.2.1/qualimap bamqc -bam $IN/$f --java-mem-size=4G \
		-gd HUMAN \
		-gff /fast/projects/ngs_sers/Packages/Regions/WES/WES_Enriched_Region.bed \
		-sd -sdmode 2 \
		-os -outdir $OUTPATH		
done
	
