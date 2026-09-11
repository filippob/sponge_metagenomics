
echo "$HOME"

for R in taxprofiler/kraken2/db2/*;
do
	SAMPLE=$(basename $R .report)
	echo $SAMPLE 
	bracken \
	    -d /path/to/kraken2_db \
	    -i $R \
	    -o bracken_outputs/${SAMPLE}.bracken \
	    -w bracken_outputs/${SAMPLE}.bracken.report \
	#   -r 150 \
	    -l S
done


