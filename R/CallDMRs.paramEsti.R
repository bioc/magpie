### DMR modeling and estimate
CallDMRs.paramEsti <- function(counts, sf,
                               model, variable,
                               shrkPhi = TRUE) {
    #### parameter estimation based on NB model
    Ratio <- meRatio(counts = counts, sf = sf) ## methylatinon ratio
    #### 1. Preliminary estimate MLE based on NB model
    res.MLE <- MLE.parallel(
        mat = as.matrix(counts),
        sf = sf,
        D = model.matrix(model, variable)
    )
    idx <- !is.na(res.MLE$phi.nb)
    ### 2. Update phi, theta with their posterior given other parameters
    D <- model.matrix(model, variable)
    mu <- t(exp(D %*% t(res.MLE$Coef.nb)) / (1 + exp(D %*% t(res.MLE$Coef.nb))))
    if (shrkPhi) {
        PostPhi <- Posterior.phi(
            counts = counts[idx, ],
            sf = sf, D = D,
            R = res.MLE$Coef.nb[idx, ],
            phi.mom = res.MLE$phi.nb[idx],
            theta.mom = res.MLE$theta.nb[idx]
        )
        phi <- rep(NA, nrow(counts))
        phi[idx] <- as.vector(PostPhi$phi)
        theta <- rep(NA, nrow(counts))
        theta[idx] <- as.vector(PostPhi$theta)
    } else {
        phi <- res.MLE$phi.nb
        theta <- res.MLE$theta.nb
    }
    ### 3. Calculate Cov(R)

    res.Cov <- CovR.parallel(
        mat = counts, sf = sf,
        model = model, variable = variable,
        Coef = res.MLE$Coef.nb,
        phi = phi, theta = theta
    )
    ### 4. Calculate log.likelihood
    loglik <- log.lik_NB(
        x = as.matrix(counts[, seq(1, ncol(counts), 2)]),
        y = as.matrix(counts[, seq(2, ncol(counts), 2)]),
        sx = sf[seq(1, length(sf), 2)],
        sy = sf[seq(2, length(sf), 2)],
        D = D,
        R = res.MLE$Coef.nb,
        s = mylogit(phi),
        t = log(theta)
    )
    #### 5. output all estimate
    if (mean(log(phi), na.rm = TRUE) > -3) {
        phi <- phi * exp(-2.5)
        theta <- theta * exp(-2.5)
    }
    #### 6. Phi, Theta by BB
    phi.BB <- updateEstiPhibyBB(
        counts = counts,
        sf = sf,
        design = variable, model = model
    )

    theta.BB <- updateEstiTheta(
        counts = counts,
        sf = sf,
        mu0 = mu,
        phi0 = phi.BB
    )
    DMR <- list(
        Coef = res.Cov$Coef, Cov = res.Cov$Cov,
        Ratio = Ratio,
        loglik = loglik,
        Phi = round(phi, 5),
        Theta = round(theta, 5),
        Phi_BB = round(phi.BB, 5),
        Theta_BB = round(theta.BB, 5)
    )
    return(DMR)
}
