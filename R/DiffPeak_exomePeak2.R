# Differential peak calling with exomePeak2
DiffPeak_exomePeak2 <- function(counts,
                                nreps,
                                model,
                                design,
                                sf = NULL) {
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = design,
        design = model
    )

    if (length(sf) == 0) {
        ## compute size factor from all counts
        sf <- colSums(counts)
        sf <- sf / median(sf)
    }
    counts.norm <- sweep(counts, 2, sf, FUN = "/")
    sizeFactors(dds) <- sf
    dds <- DESeq(dds, test = "Wald")

    ###
    # dds <- nbinomWaldTest(dds)

    ###

    resultsNames(dds) # name of each variable in the design matrix used for DESeq
    res <- DESeq2::results(dds, name = "IPIP.TrtTrt")
    dat <- data.frame(padj = res$padj, stat = res$stat)
    return(dat)
}
