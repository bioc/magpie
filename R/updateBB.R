updateEstiPhibyBB <- function(counts, sf, design, model) {
    #### estimate dispersion using BB model
    counts.norm <- sweep(counts, 2, sf, FUN = "/")
    X.norm <- round(counts.norm[, seq(1, ncol(counts), 2)]) + 1
    Y.norm <- round(counts.norm[, seq(2, ncol(counts), 2)]) + 1

    Total <- X.norm + Y.norm
    rownames(Total) <- seq_len(nrow(Total))
    Z <- log((Y.norm + 0.5) / (X.norm + 0.5))
    Ratio <- TRESS::meRatio(counts = counts, sf = sf)
    rownames(Ratio) <- seq_len(nrow(Ratio))
    phi <- bb.pval <- rep(NA, nrow(Ratio))
    COEFS <- matrix(NA, nrow = nrow(Z), ncol = 2)
    D <- model.matrix(model, design)
    for (i in seq_len(nrow(Ratio))) {
        # cat(i, sep = "\n")
        dat <- data.frame(
            n = as.numeric(Total[i, ]),
            y = as.numeric(Y.norm[i, ]),
            Pred = factor(rep(
                c(0, 1),
                as.numeric(table(design$Trt))
            ))
        )
        resfit <- aod::betabin(cbind(y, n - y) ~ Pred, ~1, data = dat, link = "logit")
        COEFS[i, ] <- coef(resfit)
        phi[i] <- resfit@param[3]
    }


    eta <- t(model.matrix(model, design) %*% t(COEFS))
    mu <- exp(eta) / (1 + exp(eta))
    mu[mu < 1e-5] <- 0.01
    mu[mu > (1 - (1e-5))] <- 0.99

    ##### posterior of phi, in BB*logN(m,v)
    Jloglik <- function(this.phi, this.y, this.n, this.mu, m, v) {
        log.gamma <- function(x) {
            #### calculate log(gamma(x)), x can be a scaler or vector
            N <- length(x)
            log.gx <- rep(NA, length(x))
            for (ii in seq_len(length(x))) {
                if (is.infinite(log(gamma(x[ii])))) {
                    if (round(x[ii]) == x[ii]) {
                        log.gx[ii] <- sum(log(seq_len(x[ii] - 1)))
                    } else {
                        floor.x <- floor(x[ii])
                        resi.x <- x[ii] - floor(x[ii])
                        log.gx[ii] <- sum(log(seq(x[ii] - 1, resi.x, by = -1))) + log(gamma(resi.x))
                    }
                } else {
                    log.gx[ii] <- log(gamma(x[ii]))
                }
            }
            return(log.gx)
        }
        #### log(BB*logN)
        this.alpha <- this.mu * (this.phi^
            {
                -1
            } - 1)
        this.beta <- (1 - this.mu) * (this.phi^
            {
                -1
            } - 1)
        dbetabin <- sum(
            log.gamma(this.n + 1) + log.gamma(this.y + this.alpha) +
                log.gamma(this.n - this.y + this.beta) + log.gamma(this.alpha + this.beta) -
                log.gamma(this.y + 1) - log.gamma(this.n - this.y + 1) - log.gamma(this.n + this.alpha + this.beta) -
                log.gamma(this.alpha) - log.gamma(this.beta)
        ) + dlnorm(this.phi, meanlog = m, sdlog = v, log = TRUE)

        return(dbetabin)
    }

    m <- mean(log(phi)[log(phi) > -10])
    v <- sd(log(phi)[log(phi) > -10])
    postPhi <- rep(NA, length(phi))
    for (i in seq_len(nrow(mu))) {
        # cat(i, sep = "\n")
        tmp <- optimize(Jloglik,
            interval = c(0, 1),
            this.y = Y.norm[i, ],
            this.n = Total[i, ],
            this.mu = mu[i, ],
            m = m, v = v,
            maximum = TRUE
        )
        postPhi[i] <- tmp$maximum
    }

    #####
    return(postPhi)
}

updateEstiTheta <- function(counts, sf, mu0, phi0) {
    #### update theta given posterior of phi
    X <- counts[, seq(1, ncol(counts), 2)]
    Y <- counts[, seq(2, ncol(counts), 2)]
    sx <- sf[seq(1, length(sf), 2)]
    sy <- sf[seq(2, length(sf), 2)]

    log.lik.theta <- function(theta, phi, yy, xx, mmu, sx, sy) {
        loglik <- sum(dnbinom(yy,
            size = mmu * (phi^
                {
                    -1
                } - 1),
            prob = 1 / (1 + sy * theta), log = TRUE
        ) +
            dnbinom(xx,
                size = (1 - mmu) * (phi^
                    {
                        -1
                    } - 1),
                prob = 1 / (1 + sx * theta), log = TRUE
            ))
        loglik
    }
    res <- matrix(NA, nrow = nrow(mu0), ncol = 2)
    colnames(res) <- c("theta", "obj")
    for (i in seq_len(nrow(mu0))) {
        mu.i <- mu0[i, ]
        # cat(i, sep = "\n")
        if (!all(mu.i <= 0.51, na.rm = TRUE)) {
            tmp <- optimize(log.lik.theta,
                c(0, 100000),
                tol = 0.0001,
                phi = phi0[i],
                yy = Y[i, ], xx = X[i, ],
                mmu = mu.i,
                sx = sx, sy = sy,
                maximum = TRUE
            )
        }
        res[i, ] <- c(tmp$maximum, tmp$objective)
    }
    res <- as.data.frame(res)
    return(res$theta)
}
