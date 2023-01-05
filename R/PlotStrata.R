#' An individual plot of power evaulation results by strata.
#'
#' This function plots a certain power measurement of the original sequencing depth by strata. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @param Power.list A list produced by \code{\link{PowerEval}}.
#' @param value_option A character indicating which measurement to plot. Options include "FDR", "FDC", "Power", and "Precision".
#'
#' @return It plots a certain power measurement of the original sequencing depth by strata. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @import RColorBrewer reshape2
#'
#'
#' @export
#' @importFrom graphics axis par
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
#' ### plot an individual strata result
#' PlotStrata(power.test, value_option = "FDR")
#' }
PlotStrata <- function(Power.list,
                       value_option = "FDR") {
    options(warn = -1)
    Power.list <- Power.list[["1x"]][5:8]
    power_sub <- Power.list[[value_option]]
    # names(power_sub) = names(Power.list)
    power_sub <- melt(power_sub,
        id.vars = c("N.rep"),
        variable.name = "Strata",
        value.name = paste0(value_option)
    )

    # FDR.toP = Power.list[["FDR"]]
    # par(mfrow = c(2,2))
    # df_plot <- power_sub[[paste0(SD_multiplier, "x")]]
    df_plot <- power_sub
    names(df_plot)[1] <- "Number of Replicates"
    xvals <- split(df_plot[, 2], df_plot[, 1])
    yvals <- split(df_plot[, 3], df_plot[, 1])

    plot(seq_along(unique(as.factor(unlist(xvals)))),
        ylim = c(0, max(unlist(yvals), na.rm = TRUE)), type = "n", xaxt = "n",
        xlab = "Average Input Strata by Percentile", ylab = value_option
    )

    mapply(lines, xvals, yvals,
        col = brewer.pal(n = nrow(unique(df_plot[1])), name = "Set1")[seq_len(nrow(unique(df_plot[1])))],
        pch = seq_len(nrow(unique(df_plot[1]))), type = "o"
    )
    axis(1,
        at = seq_along(unique(as.factor(unlist(xvals)))),
        labels = unique(as.factor(unlist(xvals)))
    )
    title(main = paste0(value_option, " by Strata \n", "Sequencing Depth: 1x"))
    # mtext(side = 3, line = 0.25, at = 1, adj = -2, "Sequencing Depth: 1x")

    legend("bottomleft",
        legend = paste0("n = ", unique(df_plot[, 1])),
        title = "",
        col = brewer.pal(n = nrow(unique(df_plot[1])), name = "Set1")[seq_len(nrow(unique(df_plot[1])))], pch = seq_len(nrow(unique(df_plot[1]))), cex = 1, bty = "n", bg = "transparent"
    )
}
