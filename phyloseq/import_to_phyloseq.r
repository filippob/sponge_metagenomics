
library("tidyr")
library("dplyr")
library("tools")
library("phyloseq")
library("data.table")

## PARAMETERS
args = commandArgs(trailingOnly=TRUE)
if (length(args) >= 1) {
  
  #loading the parameters
  if (tools::file_ext(args[1]) %in% c("r","R")) {
    
    source(args[1])
    # source("Analysis/config.R")
  } else {
    
    load(args[1])
  }
  
} else {
  #this is the default configuration, used for development and debug
  writeLines('Using default config')
  
  #this dataframe should be always present in config files, and declared
  #as follows
  config = NULL
  config = rbind(config, data.frame(
    prjfolder = "Documents/metagenomics/sponge",
    output_folder = "analysis",
    count_table = "taxprofiler/species_abundance.tsv",
    taxonomy_file = "taxprofiler/taxonomy.tsv",
    conf_file = "data/metadata.csv",
    nfactors = 2, ## n. of design variables (e.g. treatment and timpoint --> nfactors = 2)
    min_counts = 10,
    project = "", ## USE ONLY FOR SUBSETTING !!
    treatment_column = "Type",
    sample_column = "SampleID",
    subset_variable = "",
    subset_group = "", ## subset data by sample variable (e.g. experiment, group, sex, etc.),
    force_overwrite = FALSE
  ))
}

HOME <- Sys.getenv("HOME")
prjfolder = file.path(HOME, config$prjfolder)
output_folder = file.path(prjfolder,config$output_folder)

count_table = config$count_table
taxonomy_file = config$taxonomy_file
conf_file = config$conf_file

print(paste("Abundance file is:", count_table))
print(paste("Taxonomy file is:", taxonomy_file))
print(paste("Metadata file is:", conf_file))

print(paste("Output folder is:", output_folder))

writeLines(" - creating output folder")
dir.create(output_folder, showWarnings = FALSE)
config_fname = file.path(output_folder, "import_phyloseq.config.RData")
save(config, file = config_fname)

writeLines(" - read abundance table")
fname =  file.path(prjfolder, count_table)
otu <- fread(fname)
names(otu) = gsub("_ERZ21861969_db2.kraken2.kraken2.report.txt.bracken_num", "", names(otu))

counts <- otu |> select(-c(name,taxonomy_id,taxonomy_lvl)) |> as.matrix()
rownames(counts) = otu$taxonomy_id

# taxonomy table
writeLines(" - read taxonomies")
fname =  file.path(prjfolder, taxonomy_file)
tax <- read.delim(
  fname,
  row.names = 1,
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(tax) <- c("lineage", "lineage_reduced")

tax <- tax %>%
  separate(
    lineage_reduced,
    into = c(
     # "cellular_organisms",
     # "Domain",
     # "Kingdom",
      "blank",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species",
      "str."
    ),
    sep = ";",
    fill = "right"
  ) |>
  select(-c(blank, lineage))

tax_mat <- as.matrix(tax)

# metadata
writeLines(" - read metadata ")
fname =  file.path(prjfolder, conf_file)
metadata <- read.delim(fname, row.names = 1, check.names = FALSE, sep = ",")

writeLines(" - creating phyloseq object")
physeq <- phyloseq(
  otu_table(counts, taxa_are_rows=TRUE),
  tax_table(tax_mat),
  sample_data(metadata)
)

print(paste("N. of OTUs is:", nrow(otu_table(physeq))))
print(paste("N .of samples is:", ncol(otu_table(physeq))))

writeLines(" - filtering taxa for minimum counts")
physeq_pruned <- prune_taxa(taxa_sums(physeq) > config$min_counts, physeq)

print(paste("Minimum number of counts to be retained:", config$min_counts))
print(paste("N. of OTUs after filtering is:", nrow(otu_table(physeq_pruned))))

sorted_taxa <- sort(taxa_sums(physeq_pruned))
sorted_samples <- sort(sample_sums(physeq_pruned))

fname = file.path(output_folder, "sorted_taxa.png")
png(fname)
barplot(sorted_taxa, col="red")
dev.off()

fname = file.path(output_folder, "sorted_samples.png")
png(fname)
barplot(sorted_samples, col="green")
dev.off()

## save phyloseq object
writeLines(" - saving phyloseq object with project data")
fname = file.path(output_folder, "phyloseq.RData")
save(physeq_pruned, file = fname)

print("DONE!")