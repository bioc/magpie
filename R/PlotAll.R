######## Plot all results
#' All plots of power evaulation results under various scenarios.
#'
#' This function plots all power measurements of a certain sequencing depth in a 2x2 panel. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @param Power.list A list produced by \code{\link{PowerEval}}.
#' @param depth_factor A numerical value indicating which sequencing depth to plot. For example, 2 means doubling the original sequencing depth. Default is 1.
#'
#' @return It plots all power measurements of a certain sequencing depth in a 2x2 panel. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @import RColorBrewer reshape2
#' @importFrom graphics axis par
#' @export
#'
#' @examples
#'
#' library(PAM)
#' ### Main function
#' power.test <- QuickPower(dataset = "GSE55575", test_method = "TRESS")
#'
#' ### plot all in a panel under sequencing depth 1x
#' PlotAll(power.test, depth_factor = 1)
#'
PlotAll <- function(Power.list,
                    depth_factor = 1) {
    SD_multiplier <- depth_factor
    Power.list <- lapply(Power.list, FUN = function(x) x[seq_len(4)])
    par(mfrow = c(2, 2))
    for (value_option in c("FDR", "FDC", "Power", "Precision")) {
        power_sub2 <- lapply(Power.list, "[", value_option)
        power_sub <- lapply(power_sub2, FUN = function(x) x[[1]])
        names(power_sub) <- names(Power.list)
        power_sub <- lapply(power_sub, FUN = function(x) {
            melt(x,
                id.vars = c("N.rep"),
                variable.name = "Threshold",
                value.name = paste0(value_option)
            )
        })

        # FDR.toP = Power.list[["FDR"]]
        # par(mfrow = c(2,2))
        df_plot <- power_sub[[paste0(SD_multiplier, "x")]]
        names(df_plot)[1] <- "Number of Replicates per Group"
        xvals <- split(df_plot[, 1], df_plot[, 2])
        yvals <- split(df_plot[, 3], df_plot[, 2])


        plot(unlist(xvals),
            xlim = c(min(unlist(xvals)), max(unlist(xvals))), ylim = c(0, max(unlist(yvals))), type = "n", xaxt = "n",
            xlab = "Number of Replicates per Group", ylab = value_option
        )

        mapply(lines, xvals, yvals,
            col = brewer.pal(n = nrow(unique(df_plot[2])), name = "Set1")[seq_len(nrow(unique(df_plot[2])))],
            pch = seq_len(nrow(unique(df_plot[2]))), type = "o"
        )
        axis(1, at = unique(df_plot[, 1]), labels = unique(df_plot[, 1]))
        title(main = paste0(value_option, " by Sample Size \n", "Sequencing Depth: ", SD_multiplier, "x"))
        legend("bottomleft",
            legend = unique(df_plot[, 2]),
            title = "P-value Threshold",
            col = brewer.pal(n = nrow(unique(df_plot[2])), name = "Set1")[seq_len(nrow(unique(df_plot[2])))], pch = seq_len(nrow(unique(df_plot[2]))), cex = 0.8, bty = "n", bg = "transparent"
        )
    }

    par(mfrow = c(1, 1))
}
