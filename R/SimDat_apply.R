### Simulate data wrapper
SimDat_apply <- function(nreps,
                         nsites,
                         mu,
                         phi,
                         theta,
                         sfx,
                         sfy,
                         res.para,
                         Contrast,
                         sd_multi,
                         seed,
                         Test_method) {
    res.sim <- SimDat(
        nreps = nreps,
        nsites = nsites,
        mu = mu,
        phi = phi,
        theta = theta,
        sfx = sfx,
        sfy = sfy,
        sd_multi,
        seed = seed
    )
    # input_mean = rowMeans(res.sim$counts[,seq(1, ncol(res.sim$counts), 2)])
    # ip_mean = rowMeans(res.sim$counts[,seq(2, ncol(res.sim$counts), 2)])
    # input_vars = rowVars(res.sim$counts[,seq(1, ncol(res.sim$counts), 2)])
    # ip_vars = rowVars(res.sim$counts[,seq(2, ncol(res.sim$counts), 2)])
    #
    # save(input_mean, input_vars, ip_mean, ip_vars, file = paste0("meanvars_sd", sd_multi,"_ss_",nreps[1],"_it_", seed, ".rdata"))
    # counts_test = res.sim$counts
    # strata_list = GetStrata(counts_test)
    test_res <- vector("list", length = 3)
    # print(paste0("Data simulation finished for ", N.reps[1], " Controls, ", N.reps[2], " Cases: ", "Iteration ", isim))
    if (Test_method == "TRESS") {
        ### TRESS test
        TRESS.DMR <- CallDMRs.paramEsti(
            counts = res.sim$counts,
            sf = res.sim$sf,
            variable = res.para$design,
            model = res.para$model
        )
        res.test <- TRESS_DMRtest(DMR = TRESS.DMR, contrast = Contrast)
    } else {
        ####### exomePeak2 test
        design.exome2 <- data.frame(
            Reps = c(
                rep(paste0("Rep", seq_len(nreps[1])), each = 2),
                rep(paste0("Rep", seq_len(nreps[2])), each = 2)
            ),
            IP = rep(c("Input", "IP"), sum(nreps)),
            Trt = rep(c("Ctrl", "Trt"), 2 * nreps)
        )
        model.exome2 <- ~ Reps + IP + Trt + IP * Trt
        model.matrix(model.exome2, design.exome2)
        res.test <- DiffPeak_exomePeak2(
            counts = res.sim$counts, nreps = nreps,
            model = model.exome2, design = design.exome2
        )

        #####
    }
    sig_ind <- which(res.test$padj < 0.05)

    counts_test <- res.sim$counts[sig_ind, ]
    # print(sig_ind)
    strata_list <- GetStrata(counts_test, res.sim$counts)
    cutoffs <- Getcutoffs(counts_test)
    test_res[[1]] <- res.test$padj
    test_res[[2]] <- strata_list
    test_res[[3]] <- cutoffs
    return(test_res)
}
