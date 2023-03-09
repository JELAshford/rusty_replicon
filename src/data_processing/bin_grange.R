bin_grange <- function(full_grange, bin_size) {
    #' Create a grange of "bin_size" for full scope of the original full_grange
    #' @param grange A Granges object
    #' @param bin_size The size of the required bins
    #' @return A new Granges object with each sequence binned at bin_size
    suppressPackageStartupMessages(require("GenomicRanges"))
    binned_sequences <- endoapply(
        split(full_grange, seqnames(full_grange)),
        FUN = function(x) {
            GRanges(
                seqnames(x),
                IRanges(
                    start = seq(start(x), end(x), bin_size),
                    width = bin_size
                )
            )
        }
    )
    combined_binned_range <- Reduce(c, binned_sequences)
    return(combined_binned_range)
}


bin_granges_list <- function(target_bins, granges_list, target_col) {
    #' Bin the raw RT data into genome bins of bin_size, using GRanges objects
    #' @param target_bins Granges object to align the raw_data granges to
    #' @param granges_list List of Granges objects with data to bin
    #' @param target_col Name of the metadata column with values to aggregate
    #' @return List of two tibbles: average and count of values in each bin
    suppressPackageStartupMessages(require(GenomicRanges))
    suppressPackageStartupMessages(require(data.table))
    suppressPackageStartupMessages(require(tidyverse))
    suppressPackageStartupMessages(require(purrr))

    # Aggregate the value metadata from each set in granges_list
    cumulative_data <- map2_dfc(
        granges_list, names(granges_list),
        .f = function(set, name) {
            # Assign each bin the average of the overlapping set values
            olap <- IRanges::findOverlaps(set, target_bins)
            group <- data.table(
                to = to(olap),
                val = mcols(set)[[target_col]][from(olap)]
            )[, .(av = mean(val), n = .N), by = to]
            set_av <- rep(NA, length(target_bins))
            set_av[group$to] <- group$av
            set_n <- rep(NA, length(target_bins))
            set_n[group$to] <- group$n
            out <- list(set_av, set_n)
            names(out) <- paste0(name, c("_av", "_n"))
            return(out)
        }
    ) %>% as_tibble()

    # Extract the averages and counts separately
    average_data <- cumulative_data %>%
        dplyr::select(ends_with("_av")) %>%
        rename_with(function(x) str_replace(x, "_av", ""))
    count_data <- cumulative_data %>%
        dplyr::select(ends_with("_n")) %>%
        rename_with(function(x) str_replace(x, "_n", ""))
    return(list("average_data" = average_data, "count_data" = count_data))
}
