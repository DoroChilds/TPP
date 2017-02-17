test_that(desc = "2d_data_were_outcome_depends_on_OS", code = {
  y <- c(norm_rel_fc_protein_0_transformed = 0, 
         norm_rel_fc_protein_0.02_transformed = 0.4215532379877233726262,
         norm_rel_fc_protein_0.143_transformed = 1.026455783468937177361,
         norm_rel_fc_protein_1_transformed = 0.9166168400334454569034,
         norm_rel_fc_protein_5_transformed = 1)
  
  x <- c(0, 0.02, .143, 1, 5)
  xLog <- c(-15, log10(x * 10^-6)[-1])
  
  yFit <- TPP:::fitDRCurve(protID = "X020466_44.1_IPI00007956.4", expName = "dummyExp", dose = xLog, response = y, cpd_effect = "stabilized", slBds = c(1,50), verbose = FALSE)
  
  ref <- tibble(Protein_ID = "X020466_44.1_IPI00007956.4",
                pEC50 = 7.686551133741927444021,
                slope = 25.47714239657176449327,
                R_sq = 0.9905214207643284751725,
                pEC50_outside_conc_range = FALSE,
                pEC50_quality_check = "7.68655113374193",
                model_converged = TRUE,
                sufficient_data_for_fit = TRUE)
  
  check1 <- all.equal(ref$pEC50, yFit$pEC50)
  check2 <- all.equal(ref$slope, yFit$slope)
  check3 <- all.equal(ref$R_sq, yFit$R_sq)
  check4 <- all.equal(ref$pEC50_quality_check, yFit$pEC50_quality_check)
  
  expect_true(check1 & check2 & check3 & check4)
})

