## This is needed if nf-core/taxprofiler was run without bracken
## Only the different databases (kraken2, metaphlan, mOTUs etc.)
## bracken will process kraken2 output to prduce count tables (one per sample) 

echo "$HOME"
DATABASE=taxprofiler/databases/kraken2

for R in taxprofiler/kraken2/db2/*;
do
	SAMPLE=$(basename $R .report)
	echo $SAMPLE 
	bracken \
	    -d $DATABASE \
	    -i $R \
	    -o bracken_outputs/${SAMPLE}.bracken \
	    -w bracken_outputs/${SAMPLE}.bracken.report \
	    -l S
done


