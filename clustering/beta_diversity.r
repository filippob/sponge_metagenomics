## Script that analyses the composition of a microbiome
library("tools")
library("dplyr")
library("vegan")
library("plotly")
library("ggplot2")
library("phyloseq")
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
    output_folder = "analysis/figures",
    input_data = "analysis/phyloseq_norm.RData",
    treatment_column = "year",
    method = "bray",
    force_overwrite = FALSE
  ))
}
################################################################################

HOME <- Sys.getenv("HOME")
prjfolder = file.path(HOME, config$prjfolder)
outdir = file.path(prjfolder,config$output_folder)

if(!file.exists(file.path(outdir))) dir.create(file.path(outdir), showWarnings = FALSE)

config_fname = file.path(outdir, "composition.config.RData")
save(config, file = config_fname)

## loading data previously imported in phyloseq
writeLines(" - loading normalised phyloseq data")
fname = file.path(prjfolder, config$input_data)
load(fname)


## beta diversity
writeLines(" - calculate microbiome-based distances")
dist.bc <- distance(physeq_norm, method = config$method)

writeLines(" - 2D-clustering")
ord.nmds <- ordinate(physeq_norm,
                     method="NMDS",
                     distance="bray")

p <- plot_ordination(physeq_norm,
                ord.nmds,
                color=config$treatment_column)

p <- p + stat_ellipse(type = "norm", alpha = 0.5, linewidth = 1.25)

fname = paste("cluster-", config$treatment_column, ".png", sep="")
fname = file.path(outdir,fname)
ggsave(filename = fname, plot = p, device = "png")

mds_D <- cmdscale(dist.bc, k = 3) |>
  as_tibble() |>
  mutate(id = labels(dist.bc)) |>
  rename(dim1 = V1, dim2 = V2, dim3 = V3)

metadata <- sample_data(physeq_norm) |> as.data.frame()
metadata$id = row.names(metadata)

mds_D <- mds_D |> inner_join(metadata, by = "id") |>
  select(-c(location, Old.New.Ranking))

writeLines(" - 3D-clustering")
group_var <- metadata[[config$treatment_column]]

p <- plot_ly(data = mds_D, 
             x = ~dim1, y = ~dim2, z = ~dim3,
             type = "scatter3d",
             color = ~ group_var,
             colors = c('#BF382A', '#0C4B8E', '#FFD700','#3CB371', '#FFB6C1')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')),
         annotations = list(
           x = 0.005,
           y = 0.01,
           text = 'Microbiome clustering',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))

fname = paste("3d-clustr-",config$treatment_column,"-", ".html", sep="")
fname = file.path(outdir, fname)
htmlwidgets::saveWidget(p, fname)

## Permanova
writeLines(" - PERMANOVA")

form <- as.formula(
  paste("dist.bc ~", config$treatment_column)
)

obj <- adonis2(form, data = data.frame(metadata))
print(obj)

pvalue = obj["Model","Pr(>F)"]
print(paste("p-value for", config$treatment_column, "is", pvalue))

print("DONE!")