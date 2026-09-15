
## 1) download taxdump.tar.gx from NCBI
wget http://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

## 2) extract into: /home/filippo/.taxonkit/

## 3) extract tax IDs from count table
awk 'NR>1 {print $2}' *.bracken | sort -u > taxids.txt

## 4) use taxonkit to produce the taxonomy file
taxonkit lineage taxids.txt | taxonkit reformat > taxonomy.tsv

