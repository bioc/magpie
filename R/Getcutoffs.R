#Get quantile cutoffs based in mean input values
Getcutoffs <- function(Counts) {
    input_counts <- Counts[, seq(1, ncol(Counts), 2)]
    mean_inputs <- rowMeans(input_counts)
    cutoffs <- as.vector(quantile(mean_inputs, probs = seq(0.25, 0.75, 0.25)))
    return(cutoffs)
}
