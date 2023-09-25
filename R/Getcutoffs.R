# Get quantile cutoffs based in mean input values
Getcutoffs <- function(Counts,
                       sig_ind) {
    if (length(sig_ind) == 1){
      Counts <- matrix(Counts, nrow = 1)
      input_counts <- Counts[, seq(1, ncol(Counts), 2)]
      input_counts <- matrix(input_counts, nrow = 1)
      mean_inputs <- rowMeans(input_counts)
      cutoffs <- as.vector(quantile(mean_inputs, probs = seq(0.25, 0.75, 0.25)))
    }else{
      input_counts <- Counts[, seq(1, ncol(Counts), 2)]
      mean_inputs <- rowMeans(input_counts)
      cutoffs <- as.vector(quantile(mean_inputs, probs = seq(0.25, 0.75, 0.25)))
    }  
    return(cutoffs)
}
