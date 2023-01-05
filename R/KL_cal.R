KL.qp <- function(Real.dist, Sim.dist) {
    Real.range <- range(log(Real.dist[Real.dist > 0]))
    int <- seq(Real.range[1], Real.range[2], by = 0.1)
    if (max(int) > Real.range[2]) {
        xbrks <- int
    } else {
        xbrks <- c(int, max(int) + 0.1)
    }
    p.real <- hist(log(Real.dist[Real.dist > 0]),
        breaks = xbrks,
        plot = FALSE
    )$density * 0.1
    sum(p.real)
    p.real[p.real == 0] <- min(p.real[p.real > 0]) / 2
    p.real <- p.real / sum(p.real)
    min(p.real)

    #####
    Sim.dist <- Sim.dist[Sim.dist > 0]
    Sim.dist <- Sim.dist[log(Sim.dist) >= Real.range[1] & log(Sim.dist) <= Real.range[2]]
    q.sim <- hist(log(Sim.dist), breaks = xbrks, plot = FALSE)$density * 0.1

    q.sim[q.sim == 0] <- min(q.sim[q.sim > 0]) / 2
    q.sim <- q.sim / sum(q.sim)
    min(q.sim)
    KL.sTOr <- sum(p.real * log(p.real / q.sim))

    return(KL.sTOr)
}

CompareSimRealCount <- function(Sim.Candidates, Real.Candidates, PLOT = TRUE,
                                input.lim = 1000, ip.lim = 1500) {
    #### compare simulated and real candidate count density
    sim.X <- sweep(Sim.Candidates$counts[, seq(1, ncol(Sim.Candidates$counts), by = 2)], 2,
        Sim.Candidates$sf[seq(1, ncol(Sim.Candidates$counts), by = 2)],
        FUN = "/"
    )
    sim.Y <- sweep(Sim.Candidates$counts[, seq(2, ncol(Sim.Candidates$counts), by = 2)], 2,
        Sim.Candidates$sf[seq(2, ncol(Sim.Candidates$counts), by = 2)],
        FUN = "/"
    )
    Real.X <- sweep(Real.Candidates$Counts[, seq(1, ncol(Real.Candidates$Counts), 2)], 2,
        Real.Candidates$sf[seq(1, ncol(Real.Candidates$Counts), 2)],
        FUN = "/"
    )
    Real.Y <- sweep(Real.Candidates$Counts[, seq(2, ncol(Real.Candidates$Counts), 2)], 2,
        Real.Candidates$sf[seq(2, ncol(Real.Candidates$Counts), 2)],
        FUN = "/"
    )

    if (PLOT) {
        par(mfrow = c(3, 2))
        hist(log(sim.X), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.6))
        hist(log(sim.Y), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.5))
        hist(log(Real.X), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.6))
        hist(log(Real.Y), 30, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.5))

        qqplot(
            x = rowMeans(Real.X), y = rowMeans(sim.X), pch = 16, cex = 0.7,
            xlab = "Real Input", ylab = "Simulated Input", main = "Q-Q Plot",
            xlim = c(0, input.lim), ylim = c(0, input.lim)
        )
        abline(0, 1, col = "red")

        qqplot(
            x = rowMeans(Real.Y), y = rowMeans(sim.Y), pch = 16, cex = 0.7,
            xlab = "Real IP", ylab = "Simulated IP", main = "Q-Q Plot",
            xlim = c(0, ip.lim), ylim = c(0, ip.lim)
        )
        abline(0, 1, col = "red")
    }

    #### KL divergence
    KL <- (KL.qp(Real.dist = Real.X, Sim.dist = sim.X) + KL.qp(Real.dist = Real.Y, Sim.dist = sim.Y)) / 2
    # cat("The KL divergence of simulated count to real count is ", as.numeric(KL), sep = "\n")
    return(KL)
}

KL_cal <- function(N_reps,
                   ParaEsti,
                   idx.dmr,
                   Candidates,
                   sd_multi,
                   nsim = 100) {
    model_sim <- c("NB", "BB")
    KL_res <- vector("list", length = length(N_reps))
    names(KL_res) <- paste0("SampleSize_", N_reps)
    for (i in seq_len(length(N_reps))) {
        KL_val <- data.frame(NB = rep(NA, nsim), BB = rep(NA, nsim))
        for (dis in model_sim) {
            N.reps <- c(N_reps[i], N_reps[i])
            res.para <- Simpara(ParaEsti,
                idx.dmr = idx.dmr,
                nreps = N.reps,
                model_dist = dis
            )
            flag <- res.para$flag
            nsites <- dim(ParaEsti$Coef)[1]
            Contrast <- c(c(0, 1), rep(0, ncol(res.para$design) - 1))

            sf.x.median <- median(Candidates$sf[seq(1, length(Candidates$sf), 2)])
            sf.y.median <- median(Candidates$sf[seq(2, length(Candidates$sf), 2)])
            # res_all = vector("list",length = 2)
            KL_val[, dis] <- vapply(seq_len(nsim), FUN = function(seed) {
                KL_single(
                    nreps = N.reps,
                    nsites = nsites,
                    mu = res.para$mu,
                    phi = res.para$phi,
                    theta = res.para$theta,
                    sfx = sf.x.median,
                    sfy = sf.y.median,
                    sd_multi,
                    seed,
                    Candidates = Candidates
                )
            }, numeric(1))
        }
        KL_res[[i]] <- KL_val
    }

    return(KL_res)
}
