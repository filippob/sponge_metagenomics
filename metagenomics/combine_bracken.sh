## script to combine multiple bracken outputs into one single count table
## by default: this script will calculate BOTH counts and relative abundances for each sample
## hence: two solumns per sample
## I am using here a modified version of the script (provided with bracken) that only returns counts
## modifications at both header (column names) and columns (only counts are returned)
## commented out: L128, L140, L141
## modified: L129, L142
## L128-129: 	#o_file.write("\t%s_num\t%s_frac" % (name,name))
##        	o_file.write("\t%s_num" % (name))
## L140-142:
## 		#perc = float(num)/float(total_reads[sample])
##              #o_file.write("\t%i\t%0.5f" % (num,perc))
##              o_file.write("\t%i" % (num))

combine_bracken_outputs.py --files bracken_output/*.bracken  -o species_abundance.tsv

