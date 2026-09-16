## Script that analyses the composition of a microbiome
library("tools")
library("dplyr")
library("ggplot2")
library("phyloseq")
library("microbiome")
library("data.table")
library("RColorBrewer")

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
    treatment_column = "Island",
    topn = 15, ## max n. of phylum to display
    detection_threshold = 1e-3, ## for core microbiome
    min_prevalence = 0.8, ## for core microbiome
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

## phylum
writeLines(" - extract phylums")
ps.phylum <- tax_glom(physeq_norm, taxrank = "Phylum")

print("Calculate relative abundances")
ps.rel <- transform_sample_counts(ps.phylum, function(x) x / sum(x))
df <- psmelt(ps.rel)

# Top N phyla
writeLines(" - phylum composition")
topn <- df %>%
  group_by(Phylum) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice(1:config$topn) %>%
  pull(Phylum)

df$Phylum2 <- ifelse(df$Phylum %in% topn, df$Phylum, "Other")
df$Phylum2 = factor(df$Phylum2, levels = rev(c(topn, "Other")))

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum2)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 4))

fname = file.path(outdir, "phylum_per_sample.png")
ggsave(filename = fname, plot = p, device = "png")

## GROUP BY YEAR
df.group <- df %>%
  group_by(year, Phylum2) %>%
  summarise(Abundance = mean(Abundance),
            .groups = "drop")

cols <- c(
  "#1B9E77",
  "#D95F02",
  "#7570B3",
  "#E7298A",
  "#66A61E",
  "#E6AB02",
  "#A6761D",
  "#666666",
  "#1F78B4",
  "#33A02C",
  "#E31A1C",
  "#FF7F00",
  "#6A3D9A",
  "#B15928",
  "#17BECF",
  "#BCBD22"
)

p <- ggplot(df.group,
       aes(x = "",
           y = Abundance,
           fill = Phylum2)) + 
  geom_bar(stat="identity", width=1) + facet_wrap(~year, ncol=1) + coord_polar("y", start=0) + 
  scale_fill_manual(values = cols) +
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.y = element_blank(),
                          plot.title = element_text(hjust = 0.5, color = "#666666")) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + 
  theme(legend.title = element_text(size = 7), 
               legend.text = element_text(size = 6))

fname = file.path(outdir, "phylums_per_year.png")
ggsave(filename = fname, plot = p, device = "png")

## GROUP BY ISLAND
df.group <- df %>%
  group_by(island, Phylum2) %>%
  summarise(Abundance = mean(Abundance),
            .groups = "drop")

p <- ggplot(df.group,
            aes(x = "",
                y = Abundance,
                fill = Phylum2)) + 
  geom_bar(stat="identity", width=1) + facet_wrap(~island, ncol=2) + coord_polar("y", start=0) + 
  scale_fill_manual(values = cols) +
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.y = element_blank(),
                          plot.title = element_text(hjust = 0.5, color = "#666666")) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + 
  theme(legend.title = element_text(size = 7), 
        legend.text = element_text(size = 6))

fname = file.path(outdir, "phylums_per_island.png")
ggsave(filename = fname, plot = p, device = "png")

#############
## taxa
############
writeLines(" - Genus composition ")
ps.rel <- transform_sample_counts(physeq_norm, function(x) x / sum(x))
topn <- names(sort(taxa_sums(ps.rel), decreasing=TRUE))[1:config$topn]
ps.topn <- prune_taxa(topn, ps.rel)

df <- psmelt(ps.topn)

topn <- df %>%
  group_by(Genus) %>%
  summarise(total = sum(Abundance)) %>%
  arrange(desc(total)) %>%
  slice(1:config$topn) %>%
  pull(Genus)

df$Genus = factor(df$Genus, levels = rev(topn))

samples <- df |>
  filter(OTU == "3101277") |>
  group_by(Sample) |>
  summarise(total = sum(Abundance)) |>
  arrange(desc(total)) |>
  pull(Sample)

df$Sample = factor(df$Sample, levels = rev(samples))

p <- ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, size = 4))

fname = file.path(outdir, "genus_per_sample.png")
ggsave(filename = fname, plot = p, device = "png")

df.group <- df %>%
  group_by(island, Genus) %>%
  summarise(Abundance = mean(Abundance),
            .groups = "drop")

p <- ggplot(df.group,
       aes(x = "",
           y = Abundance,
           fill = Genus)) + 
  geom_bar(stat="identity", width=1) + facet_wrap(~island, ncol=2) + coord_polar("y", start=0) + 
  scale_fill_manual(values = cols) +
  theme_classic() + theme(axis.line = element_blank(),
                          axis.text = element_blank(),
                          axis.ticks = element_blank(),
                          axis.title.y = element_blank(),
                          plot.title = element_text(hjust = 0.5, color = "#666666")) +
  guides(color = guide_legend(override.aes = list(size = 0.5))) + 
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 7))

fname = file.path(outdir, "genus_per_island.png")
ggsave(filename = fname, plot = p, device = "png")


##########################à
## CORE MICROBIOME
## detection and prevalence must be chosen accordingly to the input data (if relative abundances or counts)

writeLines(" - analysing the core microbiome")

## get taxa names
tax <- tax_table(ps.rel)
labels <- paste(
  # tax[, "Genus"],
  tax[, "Species"],
  sep = " "
)

labels <- make.unique(labels)
taxa_names(ps.rel) <- labels


core_taxa <- core_members(ps.rel, detection = config$detection_threshold, prevalence = config$min_prevalence)
ps.core <- core(ps.rel, detection = config$detection_threshold, prevalence = config$min_prevalence)
core_taxa_names <- taxa(ps.core)
core.abundance <- sample_sums(core(ps.rel, detection = config$detection_threshold, prevalence = config$min_prevalence))
# barplot(core.abundance)

# Prevalence
prev <- apply(otu_table(ps.core), 1, function(x) mean(x > 0))
# Mean relative abundance
mean_abund <- apply(otu_table(ps.core), 1,  mean)

core_table <- data.frame(
  Taxon = names(prev),
  Prevalence = prev,
  Mean_abundance = mean_abund
)

core_table <- core_table |>
  arrange(desc(Mean_abundance))

fname = file.path(outdir, "core_microbiome.csv")
fwrite(x = core_table, file = fname, sep = ",")

#######################################
# With compositional (relative) abundances
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(ps.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")



# Core with compositionals:
# library("reshape")

prevalences <- seq(0.5, 1, .05)
detections <- round(10^seq(log10(0.001), log10(.02), length = 8), 3)

p <- plot_core(ps.rel, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE) +
  theme(axis.text.x= element_text(size=6, face="italic", hjust=1),
        axis.text.y= element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))

fname = file.path(outdir, "comre_microbiome.png")
ggsave(filename = fname, plot = p, device = "png")

print("DONE")