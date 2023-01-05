GetStrata <- function(Counts, counts_all) {
    input_counts <- Counts[, seq(1, ncol(Counts), 2)]
    mean_inputs <- rowMeans(input_counts)
    cutoffs <- as.vector(quantile(mean_inputs, probs = seq(0.25, 0.75, 0.25)))
    input_counts_all <- counts_all[, seq(1, ncol(counts_all), 2)]
    mean_inputs_all <- rowMeans(input_counts_all)
    strata_list <- list(
        strata1 = which(mean_inputs_all >= 0 & mean_inputs_all < cutoffs[1]),
        strata2 = which(mean_inputs_all >= cutoffs[1] & mean_inputs_all < cutoffs[2]),
        strata3 = which(mean_inputs_all >= cutoffs[2] & mean_inputs_all < cutoffs[3]),
        strata4 = which(mean_inputs_all >= cutoffs[3])
    )


    return(strata_list)
}
### rename this list
