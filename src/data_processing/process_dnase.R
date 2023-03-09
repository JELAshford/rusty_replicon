# Load in the DNAse data, average, and bin.
# From: https://www.encodeproject.org/experiments/ENCSR000EPH/
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(tidyverse))
source("src/data_processing/bin_grange.R")

experiment_map <- list(
    "R1" = "rsc/ENCFF924FJR.bigWig",
    "R2" = "rsc/ENCFF949ANK.bigWig"
)
hg38_ranges <- readRDS("rsc/hg38_grange.RData")
hg38_binned <- bin_grange(hg38_ranges, 500)

raw_data <- lapply(experiment_map, import)

clean_data <- lapply(raw_data, function(grange) {
    filtered_data <- grange %>%
        filter(seqnames %in% paste0("chr", 1:22))
    # Assign a matching strand to contiguous values
    stranded_data <- filtered_data %>%
        mutate(strand = c("+", "-")[c(FALSE, diff(score) != 0) + 1])
    # diffs <- c(FALSE, diff(filtered_data$score) != 0)
    # strand(filtered_data) <- c("+", "-")[diffs + 1]
    combined_grange <- stranded_data %>%
        reduce_ranges_directed(score = score)
})
