This repository hosts the code and workflow for the bioinformatics processing of sponge metageonmics data.
This is a collaboration between [CNR](https://www.cnr.it/en) (Italy) and [HCMR](https://www.hcmr.gr/en/) (Greece), within the "Metagenomics of Spongia officinalis from the Eastern
Mediterranean sea" [STM](https://www.cnr.it/en/short-term-mobility) (Short-term Mobility) 2026 project.

## General workflow

[bioinfo1] **Subtract host (sponge) reads**
  1. QC, adapter, trimming
  2. alignment on reference genome (which)
  3. extract unmapped reads (more details from Xenia)

[bioinfo2] **Taxonomic classification and abundance table**
  1. nf-core/taxprofiler [_taxprofiler.slurm.sh](metagenomics/_taxprofiler.slurm.sh) [stub]
  2. [optional] kraken2/bracken [bracken.sh](metagenomics/bracken.sh) [if kraken2 is run alongside other DBs in taxprofiler but bracken is not included]
  3. produce single count table [combine_bracken.sh](metagenomics/combine_bracken.sh) [!! a locally modified version of the script was used: only counts, no relative abundances]
  4. generate the taxonomy file [taxonomy.sh](metagenomics/taxonomy.sh) [steps to download taxonomy information and collate it into a file]

---

![workflow](material/workflow.png)
**Figure 1**: general workflow of the project (dashed lines refer to future lines of work)

---

[stat1] **Phyloseq**
1. import data into Phyloseq [import_to_phyloseq.r](phyloseq/import_to_phyloseq.r) [also filters for minimum counts (per taxon)]
2. CSS normalization of count (abundance) data [normalise_counts.r](phyloseq/normalise_counts.r)
3. composition of the sponge microbiome [microbiome_composition.r](phyloseq/microbiome_composition.r)


[stat1] **Clustering / PCA**
