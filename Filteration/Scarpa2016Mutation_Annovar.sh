#!/bin/bash


IN='/fast/projects/ngs_sers/2017_WES_Analysis/PNET/Data/ScarpaData'
for file in $( find $IN/* -name '*.tsv' ); do
        f=${file##*/}
	#echo $f
	fT=${f%%_*}
        #echo $file
        fnew=${f%%.*}\_Annot
        #echo $IN/$fnew
	perl /fast/projects/ngs_sers/Packages/tools/annovar/table_annovar.pl \
		$file /fast/projects/ngs_sers/Packages/tools/annovar/humandb/ -buildver hg19 -out $IN/$fnew \
		-remove -protocol refGene,cytoBand,EUR.sites.2015_08,cosmic70,dbnsfp30a -operation g,r,f,f,f -nastring . --otherinfo --onetranscript -csvout
done














