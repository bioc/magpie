SimDat <- function(nreps,
                   nsites,
                   mu,
                   phi,
                   theta,
                   sfx,
                   sfy,
                   sd_multi,
                   seed) {
    ##### use mu to generate lambda
    call.rgamma <- function(x, n) {
        rgamma(n, shape = x[seq_len(n)], scale = x[seq((n + 1), 2 * n)])
    }

    plx <- cbind((1 - mu) * (1 / phi - 1), theta$Input)
    lambda_x <- t(apply(plx, 1, call.rgamma, n = sum(nreps)))

    ply <- cbind(mu * (1 / phi - 1), theta$IP)
    lambda_y <- t(apply(ply, 1, call.rgamma, n = sum(nreps)))
    ### size factor
    sf.x <- rnorm(sum(nreps), sfx, 0.05)
    sf.y <- rnorm(sum(nreps), sfx, 0.02)
    ### sample IP and control based on lambda_x, and lambda_y
    para.x <- sweep(lambda_x, 2, sf.x, FUN = "*") * sd_multi
    para.y <- sweep(lambda_y, 2, sf.y, FUN = "*") * sd_multi
    ### input ~ poison(para.x),
    ### IP ~ poison(para.y)
    count.Input <- matrix(rpois(n = nsites * sum(nreps), as.vector(para.x)),
        nrow = nsites, byrow = FALSE
    )
    count.IP <- matrix(rpois(n = nsites * sum(nreps), as.vector(para.y)),
        nrow = nsites, byrow = FALSE
    )

    tmp <- sweep(count.IP, 2, sf.y, FUN = "/") / (sweep(count.IP, 2, sf.y, FUN = "/") +
        sweep(count.Input, 2, sf.x, FUN = "/"))

    counts <- matrix(0, nrow = nsites, ncol = 2 * sum(nreps))
    ix <- seq(1, ncol(counts), by = 2)
    counts[, ix] <- count.Input
    ix <- seq(2, ncol(counts), by = 2)
    counts[, ix] <- count.IP
    num.na <- sum(is.na(counts))
    if (length(num.na) > 0) {
        counts[is.na(counts)] <- rpois(num.na, 1e+06)
    }

    colnames(counts) <- seq_len(ncol(counts))
    colnames(counts)[seq(1, ncol(counts), 2)] <- c(
        paste0("Ctrl_Input_Rep_", seq_len(nreps[1])),
        paste0("Trt_Input_Rep_", seq_len(nreps[2]))
    )
    colnames(counts)[seq(2, ncol(counts), 2)] <- c(
        paste0("Ctrl_IP_Rep_", seq_len(nreps[1])),
        paste0("Trt_IP_Rep_", seq_len(nreps[2]))
    )
    sf <- rep(0, 2 * sum(nreps))
    sf[seq(1, ncol(counts), 2)] <- sf.x
    sf[seq(2, ncol(counts), 2)] <- sf.y
    return(list(
        counts = counts, sf = sf,
        lambda_y = lambda_y, lambda_x = lambda_x,
        para.x = para.x, para.y = para.y
    ))
}
