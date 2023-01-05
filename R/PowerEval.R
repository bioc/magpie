#### Main power function
#' Power evaluation for MeRIP-seq data under various study designs.
#'
#' This function conducts simulations with various user-defined study design parameters, including but not limited to
#' sample size, sequencing depth, and testing method. Users will need to provide either partial or whole-genome MeRIP-seq data
#' for parameter estimation purposes.
#'
#' @param Input.file A vector containing the name of BAM files of input samples.
#' @param IP.file A vector containing the name of BAM files of IP samples.
#' @param BamDir A character stating the directory path of all .BAM files.
#' @param annoDir A character stating the directory path of the ".sqlite" file for annotation.
#' @param variable A vector indicating the experimental conditions of all samples.
#' @param bam_factor A numrical value indicating the ratio of provided data to the whole genome data. Default is 0.05.
#' @param nsim An integer indicating the number of iterations to simulate under each scenario. Default is 10.
#' @param N.reps A vector of integers indicating the numbers of replicates to simulate, in both groups. Default is c(2,3,5).
#' @param depth_factor A vector of numerical values indicating how much sequencing depth of the provided data will be increased in simulations. For example, 2 means doubling the original sequencing depth. Default is c(1,2,5).
#' @param thres A vector of numerical values indicating the p-value thresholds used in power calculation. Default is c(0.01, 0.05, 0.1, 0.2).
#' @param Test_method A character indicating which DMR calling method to use. Options are "TRESS" and "exomePeak2". Default is "TRESS".
#'
#' @return A list of calculated power measurements that will be used as the input of functions \code{\link{WriteToxlsx}}, \code{\link{WriteToxlsx_strata}}, \code{\link{PlotAll}}, \code{\link{PlotRes}}, \code{\link{PlotALL_Strata}}, and \code{\link{PlotStrata}}.
#' Measurements include:
#' \item{FDR}{The ratio of number of false positives to the number of positive discoveries.}
#' \item{FDC}{The ratio of number of false positives to the number of true positives.}
#' \item{Power}{Statistical power.}
#' \item{Precision}{The ratio of number of true positives to the number of positive discoveries.}
#'
#' @import stats GenomicRanges GenomicFeatures Rsamtools BiocParallel DESeq2 TRESS
#'
#' @importFrom IRanges IRanges
#' @importFrom utils head
#' @importFrom S4Vectors subjectHits queryHits Rle
#' @importFrom methods is
#' @importFrom graphics hist abline lines title legend
#'
#' @import AnnotationDbi Matrix rtracklayer matrixStats aod
#'
#' @export
#'
#' @examples
#' \donttest{
#' power.test <- PowerEval(
#'     Input.file = c(
#'         "Ctrl1.chr1.input.bam", "Ctrl2.chr1.input.bam",
#'         "Case1.chr1.input.bam", "Case2.chr1.input.bam"
#'     ),
#'     IP.file = c(
#'         "Ctrl1.chr1.ip.bam", "Ctrl2.chr1.ip.bam",
#'         "Case1.chr1.ip.bam", "Case2.chr1.ip.bam"
#'     ),
#'     BamDir = "./data/GSE46705_split_chr",
#'     annoDir = "./data/annotation/hg18_chr1.sqlite",
#'     variable = rep(c("Ctrl", "Trt"), each = 2),
#'     bam_factor = 0.2,
#'     nsim = 2,
#'     N.reps = c(2, 3),
#'     depth_factor = c(1, 2, 5),
#'     thres = c(0.01, 0.05, 0.1),
#'     Test_method = "TRESS"
#' )
#' }
PowerEval <- function(Input.file,
                      IP.file,
                      BamDir,
                      annoDir,
                      variable,
                      bam_factor = 0.05,
                      nsim = 10,
                      N.reps = c(2, 3, 5),
                      depth_factor = c(1, 2, 5),
                      thres = c(0.01, 0.05, 0.1, 0.2),
                      Test_method = "TRESS") {
    ######## 1. Bam to bin-level data
    if (length(IP.file) == 0) {
        stop("Missing IP samples!",
            call. = TRUE, domain = NULL
        )
    }

    if (length(Input.file) == 0) {
        stop("Missing Input samples!",
            call. = TRUE, domain = NULL
        )
    }

    if (length(Input.file) != length(IP.file)) {
        stop("IP and Input samples are not paired!",
            call. = TRUE, domain = NULL
        )
    }

    if (is.na(annoDir)) {
        stop("Please provide the annotation file!",
            call. = TRUE, domain = NULL
        )
    }

    if (length(variable) == 0) {
        stop("Please provide variable for model fitting!",
            call. = TRUE, domain = NULL
        )
    }


    variable <- data.frame(Trt = variable)
    model <- ~ 1 + Trt

    allBins <- TRESS::DivideBins(
        IP.file = IP.file,
        Input.file = Input.file,
        Path_To_AnnoSqlite = annoDir,
        InputDir = BamDir
    )

    ######## 2. Expand bin-level data based on bam_factor
    # allBins = list()
    #
    # allBins$binCount = do.call("rbind",
    #                            replicate(round(1/bam_factor),
    #                                      allBins_user$binCount,
    #                                      simplify = FALSE))
    #
    # allBins$bins = do.call("rbind",
    #                        replicate(round(1/bam_factor),
    #                                  allBins_user$bins,
    #                                  simplify = FALSE))


    ######### 3. Bin-level data to region-level data
    Candidates <- TRESS::CallCandidates(
        Counts = allBins$binCount,
        bins = allBins$bins
    )
    Candidates <- TRESS::filterRegions(Candidates)

    # strata_list = GetStrata(do.call("rbind",
    #                                 replicate(round(1/bam_factor),
    #                                           Candidates$Counts,
    #                                           simplify = FALSE)))

    ######### 4. Estimation of parameters for simulation
    message("Estimating parameters...")
    ParaEsti <- CallDMRs.paramEsti(
        counts = do.call(
            "rbind",
            replicate(round(1 / bam_factor),
                Candidates$Counts,
                simplify = FALSE
            )
        ),
        sf = Candidates$sf,
        variable = variable,
        model = model,
        shrkPhi = TRUE
    )
    DMR.test <- TRESS::TRESS_DMRtest(DMR = ParaEsti, contrast = c(0, 1))
    message("Parameters estimation has finished...")
    idx.dmr <- which(DMR.test$padj < 0.05)
    ########## Extra. KL Selection
    KL.list <- KL_cal(
        N_reps = N.reps,
        ParaEsti = ParaEsti,
        idx.dmr = idx.dmr,
        Candidates = Candidates,
        sd_multi = 1,
        nsim = 100
    )

    message("KL divergence calculation has finished...")

    p.value <- wilcox.test(x = KL.list[[1]][, "NB"], y = KL.list[[1]][, "BB"], alternative = "g")$p.value

    if (p.value < 0.05) {
        dist_pref <- "BB"
    } else {
        dist_pref <- "NB"
    }

    ######### 5. Simulate data and calculate pvalues and results
    Power.list <- lapply(depth_factor, FUN = function(sd_multi) {
        ResSD(
            N.reps = N.reps,
            ParaEsti = ParaEsti,
            idx.dmr = idx.dmr,
            nsim = nsim,
            Candidates = Candidates,
            sd_multi,
            thre = thres,
            Test_method = Test_method,
            model_dist = dist_pref
        )
    })


    names(Power.list) <- paste0(depth_factor, "x")
    # save(Power.list, file = paste0(Outdir, "/", "PowerRes_",ExperimentName, ".rdata"))
    ######### 6. Plot results

    # PlotRes(Power.list = Power.list,
    #         Outdir = Outdir,
    #         ExperimentName = ExperimentName)
    # print("The .rdata file and plots are saved.")
    return(Power.list)
}
