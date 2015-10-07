datapath <- system.file("extdata", package="TPP")
load(file.path(datapath, "randomESet.Rdata"))

test_that("CCR_normalization_allOK", {
  datNormalized <- tppccrNormalize(data=list("Experiment_1"=exampleSet))
  mediansNew <- apply(exprs(datNormalized[[1]]), 2, median, na.rm=T)
  expect_equal(unname(mediansNew), rep(1, ncol(exampleSet)))
})

test_that("CCR_normalization_wrongInput", {
  expect_error(tppccrNormalize(data=exampleSet))
})

test_that("CCR_checkFilters", {
  df1 <- data.frame("R_sq_T1" = seq(0.1, 1, 0.1), "R_sq_T2" = seq(0.55, 1, 0.05), 
                    "Protein_ID"=as.character(1:10))
  df2 <- data.frame("meets_FC_requirement_T1" = rep(c(FALSE, TRUE), 5),
                    "meets_FC_requirement_T2" = rep(c(TRUE, FALSE), 5))
  # NA at different positions:
  df1[1,1] <- df2[2,2] <- NA
  # NA at overlapping position:
  df1[1,2] <- df2[1,2] <- NA
  res <- checkResultCols_tppccr(curveParDF=df1, fcFilterDF=df2, minR2=0.7)
  resExp <- data.frame("Protein_ID"=as.character(1:10),
                       "passed_filter_T1" = c(F,F,F,F,F,F,F,T,F,T),
                       "passed_filter_T2" = c(F,F,F,F,T,F,T,F,T,F),
                       row.names = as.character(1:10))
  #resExp[1,3] <- NA
  print(testthat::expect_equal(res, resExp))
})
