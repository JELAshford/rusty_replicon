# Load in the repliseq and convert into binned profile
# from: https://www.encodeproject.org/replication-timing-series/ENCSR188HJH/
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))
source("src/data_processing/bin_grange.R")

experiment_map <- list(
    "G1" = "rsc/ENCFF001GSV.bigWig",
    "S1" = "rsc/ENCFF001GTD.bigWig",
    "S2" = "rsc/ENCFF001GTF.bigWig",
    "S3" = "rsc/ENCFF001GTH.bigWig",
    "S4" = "rsc/ENCFF001GTK.bigWig",
    "G2" = "rsc/ENCFF001GSX.bigWig"
)
grouping <- list(
    "early" = c("G1", "S1", "S2"),
    "late" = c("S3", "S4", "G2")
)
hg38_ranges <- readRDS("rsc/hg38_grange.RData")
hg38_binned <- bin_grange(hg38_ranges, 500)

# Load and bin data
raw_data <- lapply(experiment_map, import)
binned_data <- bin_granges_list(hg38_binned, raw_data, "score")
averaged_data <- lapply(grouping, function(group) {
    binned_data$average_data %>%
        select(all_of(group)) %>%
        as.matrix() %>%
        rowMeans(., na.rm = TRUE)
})

# Generate RT values and store in granges object
rt_values <- log2(averaged_data[["early"]] / averaged_data[["late"]])
rt_grange <- hg38_binned
mcols(rt_grange) <- tibble("raw_rt" = rt_values)
rt_grange <- rt_grange %>%
    plyranges::filter(!is.nan(raw_rt)) %>%
    plyranges::filter(seqnames %in% paste0("chr", 1:22))
# Smooth the data
smoothed_rt_grange <- rt_grange %>%
    as_tibble() %>%
    group_by(seqnames) %>%
    group_modify(function(chrom_data, key) {
        smooth_fit <- lowess(chrom_data$start, chrom_data$raw_rt, f = 100 / nrow(chrom_data), delta = 500)
        chrom_data %>% mutate(smooth_rt = smooth_fit$y)
    }) %>%
    ungroup() %>%
    plyranges::as_granges()
saveRDS(smoothed_rt_grange, "out/mcf7_rt_grange.RData")

# # Visualise
# smoothed_rt_grange %>%
#     as_tibble() %>%
#     filter(start > 0, end < 100000000, seqnames == "chr1", row_number() %% 40 == 0) %>%
#     pivot_longer(cols = all_of(c("raw_rt", "smooth_rt")), names_to = "smoothing", values_to = "rt") %>%
#     ggplot(aes(x = start, y = rt)) +
#     geom_point(aes(color = as.factor(smoothing)))
