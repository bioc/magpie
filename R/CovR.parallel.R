#### Calculate Cov(R) in parallel
CovR.parallel <- function(mat,
                          sf,
                          model,
                          variable,
                          Coef,
                          phi,
                          theta) {
    X <- mat[, seq(1, ncol(mat), 2)]
    sx <- sf[seq(1, ncol(mat), 2)]
    Y <- mat[, seq(2, ncol(mat), 2)]
    sy <- sf[seq(2, ncol(mat), 2)]
    D <- model.matrix(model, variable)
    #######
    res <- bplapply(seq_len(nrow(mat)), FUN = function(i) {
        iCovR(
            i, X, Y, sx, sy, D,
            Coef, phi, theta
        )
    })
    R <- matrix(NA, nrow = nrow(mat), ncol = ncol(D))
    Cov <- vector("list", length = nrow(mat))
    # for (i in 1:length(res)) {
    for (i in seq_len(length(res))) {
        R[i, ] <- res[[i]]$R
        Cov[[i]] <- res[[i]]$Cov
    }
    colnames(R) <- colnames(D)
    #########
    return(list(Coef = R, Cov = Cov))
}
