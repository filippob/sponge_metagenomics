
library("dplyr")
library("readxl")
library("data.table")

basefolder = "~/Documents/metagenomics/sponge"
fname = "data/Spongia officinalis 138 Samples Metadata.xlsx"
outfile = "data/metadata.csv"

writeLines(" - read metadata from Excel spreadsheet")
fname = file.path(basefolder, fname)
metadata <- readxl::read_xlsx(fname)

writeLines(" - prepare metadata file")
metadata <- metadata |> rename(`sample-id` = Sample, year = `Collection Year`, island = Island, location = Location, depth = Depth) |>
  select(-c(Longitute, Latitude))

oldc <- unique(metadata$depth)
newc <- c("4-m", "5m", "4-m", "9m", "8m", "10m", "14-15m", "13m", "25-m", "7m", "10m", "42-m", "30-m", "19m", "12m", "13m", "5m", "5m")

metadata$depth = newc[match(metadata$depth, oldc)]

writeLines(" - write out metadata to csv file")

fname = file.path(basefolder, outfile)
fwrite(x = metadata, file = fname, sep = ",")


