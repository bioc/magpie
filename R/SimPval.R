### Simulate Pvalues
SimPval <- function(N.reps.Ctrl,
                    N.reps.Trt,
                    ParaEsti,
                    idx.dmr,
                    nsim,
                    Candidates,
                    sd_multi,
                    Test_method,
                    thre,
                    model_dist) {
    PVALS <- vector("list", length = length(N.reps.Ctrl))
    for (i in seq_len(length(N.reps.Ctrl))) {
        N.reps <- c(N.reps.Ctrl[i], N.reps.Trt[i])
        res.para <- Simpara(ParaEsti,
            idx.dmr = idx.dmr,
            nreps = N.reps,
            model_dist = model_dist
        )
        message_ps <- paste0("Parameters simulation has finished for sequencing depth--", sd_multi, "x, ", N.reps[1], " Controls, ", N.reps[2], " Cases.")
        message(message_ps)
        flag <- res.para$flag
        nsites <- dim(ParaEsti$Coef)[1]
        Contrast <- c(c(0, 1), rep(0, ncol(res.para$design) - 1))

        sf.x.median <- median(Candidates$sf[seq(1, length(Candidates$sf), 2)])
        sf.y.median <- median(Candidates$sf[seq(2, length(Candidates$sf), 2)])
        # res_all = vector("list",length = 2)
        res_all <- lapply(seq_len(nsim), FUN = function(seed) {
            SimDat_apply(
                nreps = N.reps,
                nsites = nsites,
                mu = res.para$mu,
                phi = res.para$phi,
                theta = res.para$theta,
                sfx = sf.x.median,
                sfy = sf.y.median,
                res.para = res.para,
                Contrast = Contrast,
                sd_multi,
                seed,
                Test_method = Test_method,
                thre = thre
            )
        })
        PVALS[[i]] <- res_all
    }
    return(PVALS)
}
