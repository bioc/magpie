#### write a strata list to the same sheet
#' Write power evaulation results of four strata under all sample size scenarios to a .xlsx file.
#'
#' This function writes power evaulation results of four strata to a .xlsx file. Only results from the original sequencing depth are saved.
#'     Here, strata are determined by mean input control levels of simulated data.
#'
#' @param pl A list produced by \code{\link{PowerEval}}.
#' @param file A character indicating the name of the output .xlsx file.
#'
#' @return It outputs a .xlsx file including FDR, FDC, power, and precision under the original sequencing depth and various sample sizes
#' and input stratas.
#'
#' @import openxlsx
#'
#' @export
#'
#' @examples
#' \donttest{
#' ### Main function
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
#'     bam_factor = 0.08,
#'     nsim = 10,
#'     N.reps = c(2, 3, 5, 7),
#'     depth_factor = c(1, 2, 5),
#'     thres = c(0.01, 0.05, 0.1),
#'     Test_method = "TRESS"
#' )
#' ### write strata results
#' WriteToxlsx_strata(power.test, file = "test1_:strata.xlsx")
#' }
WriteToxlsx_strata <- function(pl, file) {
    pl <- pl[["1x"]][seq(5, 8)]
    wb <- createWorkbook() ##

    xnames <- names(pl)
    addWorksheet(wb, sheetName = paste0("Sequencing Depth--1x")) ##
    row <- 1 ##

    for (i in seq_along(pl)) {
        col <- 1
        writeData(
            wb = wb,
            sheet = paste0("Sequencing Depth--1x"),
            x = xnames[i],
            xy = c(col, row)
        )
        # cell = createCell(row, colIndex = col)##
        # setCellValue(cell[[1, 1]], xnames[i])##
        col <- col + 1
        pl[[i]][is.na(pl[[i]])] <- NA
        writeData(
            wb = wb,
            sheet = paste0("Sequencing Depth--1x"),
            x = pl[[i]],
            startRow = row, startCol = col,
            keepNA = TRUE,
            na.string = "NA"
        )
        # addDataFrame(pl[[i]], sheet,
        #              startRow = 1, startCol = col,
        #              row.names = FALSE)##
        row <- row + nrow(pl[[i]]) + 2
    }


    saveWorkbook(wb, file = file, overwrite = TRUE)
}
