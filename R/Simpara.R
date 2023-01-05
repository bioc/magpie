Simpara <- function(ParaEsti,
                    idx.dmr,
                    nreps,
                    model_dist) {
    ### flag DMRs
    nsites <- dim(ParaEsti$Coef)[1]
    flag <- rep(FALSE, nsites)
    flag[idx.dmr] <- TRUE

    if (model_dist == "NB") {
        #### phi
        Phi0 <- ParaEsti$Phi

        #### theta

        Theta.est <- ParaEsti$Theta
    } else {
        #### phi
        Phi0 <- ParaEsti$Phi_BB

        #### theta

        Theta.est <- ParaEsti$Theta_BB
    }

    Theta0 <- list(
        Input = cbind(
            matrix(rep(Theta.est, nreps[1]),
                ncol = nreps[1], byrow = FALSE
            ),
            matrix(rep(Theta.est, nreps[2]),
                ncol = nreps[2], byrow = FALSE
            )
        ),
        IP = cbind(
            matrix(rep(Theta.est, nreps[1]),
                ncol = nreps[1], byrow = FALSE
            ),
            matrix(rep(Theta.est, nreps[2]),
                ncol = nreps[2], byrow = FALSE
            )
        )
    )
    colnames(Theta0$Input) <- c(
        paste0("Ctrl_input_rep", seq_len(nreps[1])),
        paste0("Trt_input_rep", seq_len(nreps[2]))
    )
    colnames(Theta0$IP) <- c(
        paste0("Ctrl_IP_rep", seq_len(nreps[1])),
        paste0("Trt_IP_rep", seq_len(nreps[2]))
    )

    #### mu
    alpha0 <- ParaEsti$Coef[, 1]

    beta0 <- rep(0, nsites)
    beta0[idx.dmr] <- ParaEsti$Coef[, 2][idx.dmr]


    design <- data.frame(Trt = rep(c("Ctrl", "Trt"), nreps))
    model <- ~ 1 + Trt
    model.matrix(model, design)
    eta0 <- t(model.matrix(model, design) %*% t(cbind(alpha0, beta0))) ### methylation level in both groups
    mu0 <- exp(eta0) / (1 + exp(eta0))


    colnames(mu0) <- c(paste0("Ctrl_", seq_len(nreps[1])), paste0("Trt_", seq_len(nreps[2])))


    return(list(
        mu = mu0, phi = Phi0, theta = Theta0,
        model = model, design = design, flag = flag,
        alpha = alpha0, beta = beta0
    ))
}
