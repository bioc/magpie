##Unit Test using RUnit

#test PowerEval
test_PowerEval <- function() {
  library(TBX20BamSubset)
  set.seed(123)
  power.test <- powerEval(
    Input.file = as.character(sub(".*/", "", getBamFileList()))[c(1,2,3,3)],
    IP.file = as.character(sub(".*/", "", getBamFileList()))[c(4,5,6,6)],
    BamDir = as.character(sub("/[^/]*$", "", getBamFileList()))[1],
    annoDir = file.path(system.file(package = "magpie"),
                        "extdata/mm9_chr19.sqlite"),
    variable = rep(c("Ctrl", "Trt"), each = 2),
    bam_factor = 1,
    nsim = 5,
    N.reps = c(2,3),
    depth_factor = 1,
    thres = 0.05,
    Test_method = "exomePeak2"
  )

  checkTrue(all(names(power.test) == "1x"))

  checkEquals(length(power.test), 1)
  checkEquals(length(power.test$`1x`), 8)

}


#test quickPower and plot functions
test_quickPower <- function() {
  set.seed(123)
  qp <- quickPower(dataset = "GSE94613")

  checkTrue(all(names(qp) == c("1x", "2x", "5x")))

  checkEquals(length(qp), 3)
  checkEquals(length(qp$`1x`), 8)
  checkEquals(length(qp$`2x`), 4)
  checkEquals(length(qp$`5x`), 4)

  p1 <- plotRes(qp, depth_factor = 1, value_option = "Power")

  p2 <- plotAll(qp, depth_factor = 1)

  p3 <- plotStrata(qp, value_option = "FDR")##add na

  p4 <- plotAll_Strata(qp)

  checkTrue(all(c(length(p1),
                  length(p2),
                  length(p3),
                  length(p4)) > 0))

}




