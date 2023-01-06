#' Obtain pre-calculated results from four publicly available MeRIP-seq datasets
#'
#' This function quickly outputs pre-calculated power evaluation results from four GEO MeRIP-seq datasets:
#' (GSE46705, GSE55575, GSE115105, and GSE94613). The obtained results can be used to generate Excel files and various figures.
#'
#' {GSE46705: }{Human HeLa cell line: Two replicates of wild type (WT) and two replicates of knockdown (KD) of complex METTL3.}\cr
#' {GSE55575: }{Mouse embryonic fibroblasts: Two replicates of wild type (WT) and four replicates of knockdown (KD) of WTAP.}\cr
#' {GSE115105: }{Two sample types from WT and YTHDF1 KO mice. Each type has two replicates.}\cr
#' {GSE94613: }{Human leukemia cell line: Four replicates of wild type (WT) and eight replicates of knockdown (KD) of complex METTL3.}
#'
#'
#' @param dataset A character specifying the selected dataset. Default is 'GSE46705'. Options are
#' 'GSE46705', 'GSE55575', 'GSE115105', and 'GSE94613'.
#' @param test_method A character indicating which DMR calling method to use. Options are "TRESS" and "exomePeak2". Default is "TRESS".
#'
#' @return A list of calculated power measurements that will be used as the input of functions \code{\link{WriteToxlsx}}, \code{\link{WriteToxlsx_strata}}, \code{\link{PlotAll}}, \code{\link{PlotRes}}, \code{\link{PlotALL_Strata}}, and \code{\link{PlotStrata}}.
#' Measurements include:
#' \item{FDR}{The ratio of number of false positives to the number of positive discoveries.}
#' \item{FDC}{The ratio of number of false positives to the number of true positives.}
#' \item{Power}{Statistical power.}
#' \item{Precision}{The ratio of number of true positives to the number of positive discoveries.}
#' @importFrom utils data
#' @export
#'
#' @examples
#' library(PAM)
#' power.test <- QuickPower(dataset = "GSE55575")
#'
QuickPower <- function(dataset = "GSE46705", test_method = "TRESS") {
    power.test <- NULL
    data(list = paste0(dataset, "_", test_method, "_res"), envir = environment())
    return(power.test)
}
