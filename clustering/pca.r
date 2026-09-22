## Script that analyses the composition of a microbiome
library("tools")
library("dplyr")
library("tidyr")
library("vegan")
library("plotly")
library("ggplot2")
library("phyloseq")
library("patchwork")
library("factoextra")
library("microbiome")

################################################################################
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
    repo = "Documents/metagenomics/sponge/sponge_metagenomics",
    output_folder = "analysis",
    input_data = "analysis/phyloseq_norm.RData",
    treatment_column = "island",
    method = "bray",
    force_overwrite = FALSE
  ))
}
################################################################################

HOME <- Sys.getenv("HOME")
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$output_folder)

if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

config_fname = file.path(outdir, "pca.config.RData")
save(config, file = config_fname)

## loading data previously imported in phyloseq
writeLines(" - loading normalised phyloseq data")
fname = file.path(prjfolder, config$input_data)
load(fname)

otus <- as.data.frame(otu_table(physeq_norm))
ids <- names(otus)

otus$taxa <- row.names(otus)
temp <- otus |> 
  sample_n(2500) |>
  gather(key = "sample", value = "counts", -taxa) |>
  spread(key = taxa, value = counts)

metadata <- sample_data(physeq_norm)
temp$island <- metadata$island[match(temp$sample, row.names(metadata))]
row.names(temp) <- temp$sample
labels = select(temp, c(sample,island))
temp <- select(temp, -c(sample,island))

row.names(temp) == labels$sample
labels$island

pca_result <- prcomp(temp, 
                     scale = TRUE)

bipl1 <- fviz_pca_ind(pca_result, 
                      axes = c(1, 2),
                      label="none", 
                      habillage=labels$island,
                      addEllipses=TRUE, 
                      ellipse.level=0.95, 
                      palette = "Dark2")

bipl2 <- fviz_pca_ind(pca_result, 
                      axes = c(1, 3),
                      label="none", 
                      habillage=labels$island,
                      addEllipses=TRUE, 
                      ellipse.level=0.95, 
                      palette = "Dark2")

bipl3 <- fviz_pca_ind(pca_result, 
                      axes = c(2, 3),
                      label="none", 
                      habillage=labels$island,
                      addEllipses=TRUE, 
                      ellipse.level=0.95, 
                      palette = "Dark2")


p <- (bipl1 / bipl2 / bipl3)
fname = file.path(outdir, "figures/pca_microbiome.png")
ggsave(filename = fname, plot = p, device = "png")

