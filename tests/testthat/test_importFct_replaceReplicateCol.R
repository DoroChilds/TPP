data(hdacTR_smallExample)
data(hdacCCR_smallExample)

test_that(desc="all_ok_do_replace", code={
  cfgIn <- hdacTR_config %>% mutate(Replicate = c(1,2,1,2)) %>% select(-ComparisonVT1, -ComparisonVT2)

  cfgOut <- TPP:::importFct_replaceReplicateColumn(cfg = cfgIn)
  
  check1 <- all(sort(colnames(cfgOut)) == c("126","127H","127L","128H","128L",
                                            "129H","129L","130H","130L","131L",
                                            "ComparisonVT1", "ComparisonVT2", 
                                            "Condition","Experiment"))
  
  check2 <- all(cfgOut$Comparison_1 == c("x", "", "x", ""))
  
  check3 <- all(cfgOut$Comparison_2 == c("", "x", "", "x"))
  
  expect_true(check1 & check2 & check3)
})

test_that(desc="all_ok_no_replace", code={
  cfgIn <- hdacTR_config
  
  cfgOut <- TPP:::importFct_replaceReplicateColumn(cfg = cfgIn)
  
  expect_equal(cfgIn, cfgOut)
  
})

test_that(desc="fail_because_both_column_types_defined", code={
  cfgIn <- hdacTR_config %>% mutate(Replicate = c(1,2,1,2))
  
  expect_error(TPP:::importFct_replaceReplicateColumn(cfg = cfgIn))
  
})