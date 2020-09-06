# Prepare function input:
data(hdacTR_smallExample)

dat <- hdacTR_data %>% purrr::map(. %>% filter(grepl("HDAC", gene_name)))
cfg <- hdacTR_config

## ------------------------------------------------------------------------- ##
## function 'analyzeTPPTR' with option 'methods = splinefit':
## ------------------------------------------------------------------------- ##
test_that(desc="NPARC_allok", code={
  # Start analysis
  cfgIn <- cfg
  datIn <- hdacTR_data
  tpptrResults <- analyzeTPPTR(configTable = cfgIn, data = datIn, normalize = TRUE,
                               methods = "splinefit", nCores = 1, splineDF = c(3,7))
  # Are still the expected results produced?
  cols <- colnames(tpptrResults)
  colsExpected <- c("F_statistic", "F_moderated", "F_scaled", "residual_df_H1", 
                    "prior_df_H1", "df1", "df2", "df2_moderated", 
                    "posterior_var_H1", "p_NPARC", "p_adj_NPARC")
  
  expect_equal(sort(filter(tpptrResults, p_adj_NPARC <= 0.01)$Protein_ID),
               c("COPS3", "DNLZ", "HDAC1", "HDAC10", "HDAC2", "HDAC6", "HDAC8", "PSMB5", "STX4", "YBX3"))
  expect_true(all(colsExpected %in% cols))
  expect_equal(nrow(tpptrResults), 510)
})

test_that(desc="NPARC_allok_output", code={
  # Start analysis
  cfgIn <- cfg
  datIn <- dat
  dirOut <- file.path(getwd(), "TR_results")
  tpptrResults <- analyzeTPPTR(configTable = cfgIn, data = datIn, resultPath = dirOut,
                               plotCurves = FALSE, normalize = FALSE,
                               methods = "splinefit", nCores = 1, splineDF = 3)
  # Are still the expected results produced?
  cols <- colnames(tpptrResults)
  
  check1 <- all(sort(filter(tpptrResults, p_adj_NPARC <= 0.01)$Protein_ID) ==
                  c("HDAC1", "HDAC10", "HDAC2", "HDAC6", "HDAC8"))
  check2 <- !any("plot" %in% cols)
  check3 <- all.equal(round(tpptrResults$p_adj_NPARC, 10), 
                      c(0.0000000000, 0.0000000651, 0.0000000016, 0.5654779288,
                        0.2918881921, 0.9328006686, 0.0000031269, 0.1471673324,
                        0.0006474845, NA))
  check4 <- file.exists(dirOut)
  
  unlink(dirOut, recursive = TRUE)
  
  expect_true(check1 & check2 & check3 & check4)
})


test_that(desc="NPARC_allok_plot", code={
  # Start analysis
  cfgIn <- cfg
  datIn <- dat
  dirOut <- file.path(getwd(), "TR_results")
  
  expect_warning(
    tpptrResults <- analyzeTPPTR(configTable = cfgIn, data = datIn, 
                                 normalize = FALSE, resultPath = dirOut,
                                 nCores = 1, splineDF = 3)
  )
  # Are still the expected results produced?
  cols <- colnames(tpptrResults)
  
  check1 <- all(sort(filter(tpptrResults, p_adj_NPARC <= 0.01)$Protein_ID) ==
                  c("HDAC1", "HDAC10", "HDAC2", "HDAC6", "HDAC8"))
  check2 <- all(c("meltcurve_plot", "splinefit_plot") %in% cols)
  check3 <- all.equal(round(tpptrResults$p_adj_NPARC, 10), 
                      c(0.0000000000, 0.0000000651, 0.0000000016, 0.5654779288,
                        0.2918881921, 0.9328006686, 0.0000031269, 0.1471673324,
                        0.0006474845, NA))
  check4 <- c(tpptrResults$splinefit_plot, tpptrResults$meltcurve_plot) %>% file.path(dirOut, .) %>% file.exists %>% all
  
  unlink(dirOut, recursive = TRUE)
  
  expect_true(check1 & check2 & check3 & check4)
})

test_that(desc="NPARC_allok_files", code={
  
  dirIn <- system.file("example_data", package="TPP") %>% file.path("TR_example_data")
  
  files <-  c("panobinostat_1_merged_results_20150226_1013_proteins.txt",
              "panobinostat_2_merged_results_20150226_1139_proteins.txt",
              "vehicle_1_merged_results_20150227_0908_proteins.txt",
              "vehicle_2_merged_results_20150227_0847_proteins.txt")
  
  cfgIn <- openxlsx::read.xlsx(
    file.path(dirIn,"Panobinostat_TPP-TR_config.xlsx")) %>%
    mutate(Path = file.path(dirIn, files))
  
  tpptrResults <- analyzeTPPTR(configTable = cfgIn, plotCurves = FALSE,
                               methods = "splinefit", nCores = 1, splineDF = 5)
  
  # Are still the expected results produced?
  cols <- colnames(tpptrResults)
  
  colsExpected <- c("F_statistic", "F_moderated", "F_scaled", "residual_df_H1",
                    "prior_df_H1", "df1", "df2", "df2_moderated",
                    "posterior_var_H1", "p_NPARC", "p_adj_NPARC")
  
  topHits <- filter(tpptrResults, p_adj_NPARC <= 1e-5) %>%
    arrange(p_adj_NPARC) %>% extract2("Protein_ID")
  
  hitsExpected <- c("HDAC1", "CHCHD4", "TTC38", "HDAC2", "HDAC6",
                    "HDAC10", "H2AFV|H2AFZ", "GC", "STX4")
  
  check1 <- any(grepl("TPP_results", dir(dirIn)))
  check2 <- all(colsExpected %in% cols)
  check3 <- nrow(tpptrResults) == 6004
  check4 <- all(topHits == hitsExpected)
  
  unlink(file.path(dirIn, grep("TPP_results", dir(dirIn), value = TRUE)), recursive = TRUE)
  
  expect_true(check1 & check2 & check3 & check4)
})



# test_that(desc="analyzeTPPTR_9TMTlabels", code={
#   # Remove lowest temperature from config and data tables:
#   cfgIn <- hdacTR_config %>% mutate("131L" = NULL)
#   datIn <- hdacTR_data %>% purrr::map(function(d) {d$rel_fc_131L <- NULL; return(d)})
#   # Start analysis
#   tpptrResults <- analyzeTPPTR(configTable = cfgIn, data = datIn,
#                                normalize = FALSE, # under default settings, normalization requires 10 fold changes
#                                methods = c("meltcurvefit", "splinefit"), nCores = 1)
#   # View hits:
#   filter(tpptrResults, fulfills_all_4_requirements)$Protein_ID # melting curve results
#   filter(tpptrResults, p_adj_NPARC <= 0.01)$Protein_ID # smoothing spline results
#   
# })



test_that(desc="analyzeTPPTR_9TMTlabels", code={
  cfgIn <- cfg %>% mutate("131L" = NULL)
  datIn <- dat %>% purrr::map(function(d) {d$rel_fc_131L <- NULL; return(d)})
  
  # expect error because default filter criteria try to access the 10th column
  # when not adjusted by the user.
  expect_error(analyzeTPPTR(configTable = cfgIn, data = datIn,
                            methods = "meltcurve", nCores = 1))
})


test_that(desc="meltCurves_allOK_no_conditions", code={
  # Start analysis
  cfgIn <- cfg %>% select(-Condition)
  datIn <- dat
  expect_warning(
    tpptrResults <- analyzeTPPTR(configTable = cfgIn, data = datIn, normalize = FALSE,
                                 methods = "meltcurvefit", nCores = 1)
  )
  # Are still the expected results produced?
  cols <- colnames(tpptrResults)
  compCols <-  c("Panobinostat_1_vs_Vehicle_1", "Panobinostat_2_vs_Vehicle_2")
  
  check1 <- sum(c(grepl(compCols[1], cols), grepl(compCols[2], cols))) == 8
  check2 <- identical(round(tpptrResults$pVal_adj_Panobinostat_1_vs_Vehicle_1, 2),  
                      c(0.37, 0.94, 0.56, 0.56,0.56, NA, 0.94, 0.37, 0.94, NA))
  check3 <- nrow(tpptrResults) == 10
  expect_true(check1 & check2 & check3)
})