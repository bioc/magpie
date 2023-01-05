FDR_t <- function(t0, padj, flag) {
    R <- sum(padj < t0, na.rm = TRUE)
    V <- sum(padj < t0 & (!flag), na.rm = TRUE)
    if (R == 0) {
        val <- NA
    } else {
        val <- V / R
    }
    return(val)
}

FDC_t <- function(t0, padj, flag) {
    TP <- sum(padj < t0 & (flag), na.rm = TRUE)
    V <- sum(padj < t0 & (!flag), na.rm = TRUE)
    if (TP == 0) {
        val <- NA
    } else {
        val <- V / TP
    }
    return(val)
}

TPR_t <- function(t0, padj, flag) {
    R <- sum(flag, na.rm = TRUE)
    V <- sum(padj < t0 & (flag), na.rm = TRUE)
    if (R == 0) {
        val <- NA
    } else {
        val <- V / R
    }
    return(val)
}

Precision_t <- function(t0, padj, flag) {
    P <- sum(padj < t0, na.rm = TRUE)
    V <- sum(padj < t0 & (flag), na.rm = TRUE)
    if (P == 0) {
        val <- NA
    } else {
        val <- V / P
    }
    return(val)
}

FDR_all <- function(padj.df, t0, flag) {
    tmp.fdr <- mean(apply(padj.df, 2, FUN = function(padj) FDR_t(t0 = t0, padj, flag = flag)), na.rm = TRUE)
    return(tmp.fdr)
}

FDC_all <- function(padj.df, t0, flag) {
    tmp.fdc <- mean(apply(padj.df, 2, FUN = function(padj) FDC_t(t0 = t0, padj, flag = flag)), na.rm = TRUE)
    return(tmp.fdc)
}

TPR_all <- function(padj.df, t0, flag) {
    tmp.tpr <- mean(apply(padj.df, 2, FUN = function(padj) TPR_t(t0 = t0, padj, flag = flag)), na.rm = TRUE)
    return(tmp.tpr)
}

Precision_all <- function(padj.df, t0, flag) {
    tmp.preci <- mean(apply(padj.df, 2, FUN = function(padj) Precision_t(t0 = t0, padj, flag = flag)), na.rm = TRUE)
    return(tmp.preci)
}

Power.cal <- function(PVALS,
                      idx.dmr,
                      N.reps,
                      thre) {
    Power.list <- list(FDR = NULL, FDC = NULL, Power = NULL, Precision = NULL)
    FDR.res <- NULL
    TPR.res <- NULL
    FDC.res <- NULL
    Precision.res <- NULL
    flag <- rep(FALSE, nrow(PVALS[[1]]))
    flag[idx.dmr] <- TRUE
    for (th in thre) {
        FDR.tmp <- data.frame(
            N.rep = N.reps,
            fdr = unlist(lapply(PVALS, FUN = function(padj.df) FDR_all(padj.df, t0 = th, flag = flag))),
            thresh = rep(th, length(N.reps))
        )
        FDR.res <- rbind(FDR.res, FDR.tmp)

        FDC.tmp <- data.frame(
            N.rep = N.reps,
            fdc = unlist(lapply(PVALS, FUN = function(padj.df) FDC_all(padj.df, t0 = th, flag = flag))),
            thresh = rep(th, length(N.reps))
        )
        FDC.res <- rbind(FDC.res, FDC.tmp)

        TPR.tmp <- data.frame(
            N.rep = N.reps,
            tpr = unlist(lapply(PVALS, FUN = function(padj.df) TPR_all(padj.df, t0 = th, flag = flag))),
            thresh = rep(th, length(N.reps))
        )
        TPR.res <- rbind(TPR.res, TPR.tmp)

        Precision.tmp <- data.frame(
            N.rep = N.reps,
            precision = unlist(lapply(PVALS, FUN = function(padj.df) Precision_all(padj.df, t0 = th, flag = flag))),
            thresh = rep(th, length(N.reps))
        )
        Precision.res <- rbind(Precision.res, Precision.tmp)
    }


    Power.list[["FDR"]] <- dcast(FDR.res, N.rep ~ thresh, value.var = "fdr")
    Power.list[["FDC"]] <- dcast(FDC.res, N.rep ~ thresh, value.var = "fdc")
    Power.list[["Power"]] <- dcast(TPR.res, N.rep ~ thresh, value.var = "tpr")
    Power.list[["Precision"]] <- dcast(Precision.res, N.rep ~ thresh, value.var = "precision")
    return(Power.list)
}


########### Strata

FDR_allstrata <- function(padj, t0, flag, strata_list) {
    val_allstr <- vapply(strata_list, FUN = function(x) {
        FDR_t(
            padj = padj[x],
            t0 = t0,
            flag = flag[x]
        )
    }, numeric(1))
    return(val_allstr)
}

FDC_allstrata <- function(padj, t0, flag, strata_list) {
    val_allstr <- vapply(strata_list, FUN = function(x) {
        FDC_t(
            padj = padj[x],
            t0 = t0,
            flag = flag[x]
        )
    }, numeric(1))
    return(val_allstr)
}

TPR_allstrata <- function(padj, t0, flag, strata_list) {
    val_allstr <- vapply(strata_list, FUN = function(x) {
        TPR_t(
            padj = padj[x],
            t0 = t0,
            flag = flag[x]
        )
    }, numeric(1))
    return(val_allstr)
}

Precision_allstrata <- function(padj, t0, flag, strata_list) {
    val_allstr <- vapply(strata_list, FUN = function(x) {
        Precision_t(
            padj = padj[x],
            t0 = t0,
            flag = flag[x]
        )
    }, numeric(1))
    return(val_allstr)
}

FDR_iter <- function(iter_list, t0, flag) {
    val_alliter <- rowMeans(vapply(iter_list, FUN = function(x) {
        FDR_allstrata(
            padj = x[[1]],
            t0 = t0,
            flag = flag,
            strata_list = x[[2]]
        )
    }, numeric(length(iter_list[[1]][[2]]))), na.rm = TRUE)


    return(val_alliter)
}

FDC_iter <- function(iter_list, t0, flag) {
    val_alliter <- rowMeans(vapply(iter_list, FUN = function(x) {
        FDC_allstrata(
            padj = x[[1]],
            t0 = t0,
            flag = flag,
            strata_list = x[[2]]
        )
    }, numeric(length(iter_list[[1]][[2]]))), na.rm = TRUE)

    return(val_alliter)
}

TPR_iter <- function(iter_list, t0, flag) {
    val_alliter <- rowMeans(vapply(iter_list, FUN = function(x) {
        TPR_allstrata(
            padj = x[[1]],
            t0 = t0,
            flag = flag,
            strata_list = x[[2]]
        )
    }, numeric(length(iter_list[[1]][[2]]))), na.rm = TRUE)

    return(val_alliter)
}

Precision_iter <- function(iter_list, t0, flag) {
    val_alliter <- rowMeans(vapply(iter_list, FUN = function(x) {
        Precision_allstrata(
            padj = x[[1]],
            t0 = t0,
            flag = flag,
            strata_list = x[[2]]
        )
    }, numeric(length(iter_list[[1]][[2]]))), na.rm = TRUE)

    return(val_alliter)
}

Power.cal_strata <- function(res_all2,
                             idx.dmr,
                             N.reps) {
    Power.list_strata <- list(FDR = NULL, FDC = NULL, Power = NULL, Precision = NULL)
    PVALS <- lapply(res_all2, FUN = function(x) vapply(x, FUN = function(y) y[[1]], numeric(length(x[[1]][[1]]))))

    flag <- rep(FALSE, nrow(PVALS[[1]]))
    flag[idx.dmr] <- TRUE
    cuts <- rowMeans(vapply(res_all2, FUN = function(x) rowMeans(vapply(x, FUN = function(y) y[[3]], numeric(length(x[[1]][[3]])))), numeric(3)))
    # print(cuts)
    # names(res_all2[[1]][[1]][[2]]) for strata names
    FDR.res <- data.frame(
        N.rep = rep(N.reps, each = 4),
        fdr = unlist(lapply(res_all2, FUN = function(x) {
            FDR_iter(
                iter_list = x,
                t0 = 0.05,
                flag = flag
            )
        })),
        strata = factor(rep(c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ), length(N.reps)), levels = c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ))
    )
    TPR.res <- data.frame(
        N.rep = rep(N.reps, each = 4),
        tpr = unlist(lapply(res_all2, FUN = function(x) {
            TPR_iter(
                iter_list = x,
                t0 = 0.05,
                flag = flag
            )
        })),
        strata = factor(rep(c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ), length(N.reps)), levels = c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ))
    )
    FDC.res <- data.frame(
        N.rep = rep(N.reps, each = 4),
        fdc = unlist(lapply(res_all2, FUN = function(x) {
            FDC_iter(
                iter_list = x,
                t0 = 0.05,
                flag = flag
            )
        })),
        strata = factor(rep(c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ), length(N.reps)), levels = c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ))
    )
    Precision.res <- data.frame(
        N.rep = rep(N.reps, each = 4),
        precision = unlist(lapply(res_all2, FUN = function(x) {
            Precision_iter(
                iter_list = x,
                t0 = 0.05,
                flag = flag
            )
        })),
        strata = factor(rep(c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ), length(N.reps)), levels = c(
            paste0("(0, ", round(cuts[1], 2), "]"),
            paste0("(", round(cuts[1], 2), ", ", round(cuts[2], 2), "]"),
            paste0("(", round(cuts[2], 2), ", ", round(cuts[3], 2), "]"),
            paste0("(", round(cuts[3], 2), ", Inf]")
        ))
    )

    Power.list_strata[["FDR"]] <- dcast(FDR.res, N.rep ~ strata, value.var = "fdr")
    Power.list_strata[["FDC"]] <- dcast(FDC.res, N.rep ~ strata, value.var = "fdc")
    Power.list_strata[["Power"]] <- dcast(TPR.res, N.rep ~ strata, value.var = "tpr")
    Power.list_strata[["Precision"]] <- dcast(Precision.res, N.rep ~ strata, value.var = "precision")
    return(Power.list_strata)
}
