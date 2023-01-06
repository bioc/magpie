#' All plots of power evaulation results by strata.
#'
#' This function plots all power measurements of the original sequencing depth by strata in a 2x2 panel. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @param Power.list A list produced by \code{\link{PowerEval}}.
#'
#' @return It plots all power measurements of the original sequencing depth by strata in a 2x2 panel. Power measurements to plot include "FDR", "FDC", "Power", and "Precision".
#'
#' @import RColorBrewer reshape2
#' @importFrom graphics axis par
#'
#' @export
#'
#' @examples
#'
#' library(PAM)
#' ### Main function
#' power.test <- QuickPower(dataset = "GSE55575", test_method = "TRESS")
#'
#' ### plot all strata results in a panel
#' PlotALL_Strata(power.test)
#'
PlotALL_Strata <- function(Power.list) {
    options(warn = -1)
    Power.list <- Power.list[["1x"]][5:8]
    par(mfrow = c(2, 2))
    for (value_option in c("FDR", "FDC", "Power", "Precision")) {
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
            col = brewer.pal(n = nrow(unique(df_plot[1])), name = "Set1")[seq_len(nrow(unique(df_plot[1])))], pch = seq_len(nrow(unique(df_plot[1]))), cex = 0.8, bty = "n", bg = "transparent"
        )
    }

    par(mfrow = c(1, 1))
}
