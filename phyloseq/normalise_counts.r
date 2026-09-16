## Script that applies the CSS normalisation to the microbiome count table data
library("tools")
library("phyloseq")

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
    input_data = "analysis/phyloseq.RData",
    treatment_column = "Island",
    sample_column = "SampleID",
    subset_variable = "",
    subset_group = "", ## subset data by sample variable (e.g. experiment, group, sex, etc.),
    force_overwrite = FALSE
  ))
}
################################################################################

HOME <- Sys.getenv("HOME")
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$output_folder)

config_fname = file.path(outdir, "normalise.config.RData")
save(config, file = config_fname)

source(file.path(HOME, config$repo, "support_scripts/phyloseq_transform.R")) ## from: https://github.com/vmikk/metagMisc/

## loading data previously imported in phyloseq
writeLines(" - loading phyloseq data")
fname = file.path(prjfolder, config$input_data)
load(fname)


## Preprocessing: e.g. filtering
## making normalization folder
if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

writeLines(" - CSS normalization")
physeq_norm = phyloseq_transform_css(physeq_pruned, norm = TRUE, log = FALSE)

## save phyloseq object
writeLines(" - saving phyloseq object with normalised project data")
fname = file.path(outdir, "phyloseq_norm.RData")
save(physeq_norm, file = fname)

print("DONE!")
