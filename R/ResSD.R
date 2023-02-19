# Simulate and calculate evaluation metrics for a given sequencing depth
ResSD <- function(N.reps,
                  ParaEsti,
                  idx.dmr,
                  nsim,
                  Candidates,
                  sd_multi,
                  thre,
                  Test_method,
                  model_dist) {
    ######### 5. Simulate data and calculate pvalues
    res_all2 <- SimPval(
        N.reps.Ctrl = N.reps,
        N.reps.Trt = N.reps,
        ParaEsti = ParaEsti,
        idx.dmr = idx.dmr,
        nsim = nsim,
        Candidates = Candidates,
        sd_multi,
        Test_method = Test_method,
        model_dist = model_dist
    )
    PVALS <- lapply(res_all2, FUN = function(x) vapply(x, FUN = function(y) y[[1]], numeric(length(x[[1]][[1]]))))

    # strata_list = lapply(res_all2, FUN = function(x) x[[2]])
    ######### 6. Compute FDR, TPR, FDC
    if (sd_multi == 1) {
        Power.list <- Power.cal(
            PVALS = PVALS,
            idx.dmr = idx.dmr,
            N.reps = N.reps,
            thre = thre
        )

        Power.list_strata <- Power.cal_strata(
            res_all2 = res_all2,
            idx.dmr = idx.dmr,
            N.reps = N.reps
        )

        Power.list_all <- c(Power.list, Power.list_strata)
    } else {
        Power.list_all <- Power.cal(
            PVALS = PVALS,
            idx.dmr = idx.dmr,
            N.reps = N.reps,
            thre = thre
        )
    }


    return(Power.list_all)
}
