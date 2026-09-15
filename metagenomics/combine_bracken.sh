## script to combine multiple bracken outputs into one single count table
## by default: this script will calculate BOTH counts and relative abundances for each sample
## hence: two solumns per sample
## I am using here a modified version of the script (provided with bracken) that only returns counts
## modifications at both header (column names) and columns (only counts are returned)
## commented out: L128, L140, L141
## modified: L129, L142

combine_bracken_outputs.py --files bracken_output/*.bracken  -o species_abundance.tsv

