#### write a list to the same sheet
#' Write power evaulation results under all scenarios to a .xlsx file.
#'
#' This function writes power evaulation results to a .xlsx file.
#'
#' @param pl A list produced by \code{\link{powerEval}}.
#' @param file A character indicating the name of the output .xlsx file.
#'
#' @return It outputs a .xlsx file including FDR, FDC, power, and precision under various sample sizes,
#' sequencing depths, and adjusted p-value thresholds.
#'
#' @import openxlsx
#'
#' @export
#'
#' @examples
#'
#' library(magpie)
#' ### Main function
#' power.test <- quickPower(dataset = "GSE46705", test_method = "TRESS")
#'
#' ### write out .xlsx
#' writeToxlsx(power.test, file = "test_TRESS.xlsx")
#'
writeToxlsx <- function(pl, file) {
    pl <- lapply(pl, FUN = function(x) x[seq_len(4)])
    wb <- createWorkbook() ##
    for (subid in seq_along(pl)) {
        xnames <- names(pl[[subid]])
        addWorksheet(wb, sheetName = paste0("Sequencing Depth--", names(pl)[subid])) ##
        row <- 1 ##

        for (i in seq_along(pl[[subid]])) {
            col <- 1
            writeData(
                wb = wb,
                sheet = paste0("Sequencing Depth--", names(pl)[subid]),
                x = xnames[i],
                xy = c(col, row)
            )
            # cell = createCell(row, colIndex = col)##
            # setCellValue(cell[[1, 1]], xnames[i])##
            col <- col + 1
            pl[[subid]][[i]][is.na(pl[[subid]][[i]])] <- NA
            writeData(
                wb = wb,
                sheet = paste0("Sequencing Depth--", names(pl)[subid]),
                x = pl[[subid]][[i]],
                startRow = row, startCol = col,
                keepNA = TRUE,
                na.string = "NA"
            )
            # addDataFrame(pl[[subid]][[i]], sheet,
            #              startRow = 1, startCol = col,
            #              row.names = FALSE)##
            row <- row + nrow(pl[[subid]][[i]]) + 2
        }
    }

    saveWorkbook(wb, file = file, overwrite = TRUE)
}
