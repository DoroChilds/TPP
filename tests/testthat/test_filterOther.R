data(hdacTR_smallExample)

tpptrData <- suppressMessages(
  tpptrImport(configTable = hdacTR_config, data = hdacTR_data)
)

# Data columns for sorting can be inspected by:
# fData(tpptrData[[1]]) %>% head

test_that(desc="allOk_1_column", code={
  
  datIn <- tpptrData[[1]]
  colsIn <- "qssm"
  lbIn <- 4
  ubIn <- Inf
  
  out <- TPP:::filterOther(data = datIn, cols = colsIn, lb = lbIn, ub = ubIn)
  
  check1 <- nrow(out) == 315
  check2 <- all.equal(fData(out)[[colsIn]] %>% range, c(4,125))
  
  expect_true(check1 & check2)
})

test_that(desc="allOk_2_columns", code={
  
  datIn <- tpptrData[[1]]
  colsIn <- c("qssm", "qupm")
  lbIn <- c(4, 4)
  ubIn <- c(Inf, Inf)
  
  out <- TPP:::filterOther(data = datIn, cols = colsIn, lb = lbIn, ub = ubIn)
  
  check1 <- nrow(out) == 244
  check2 <- all.equal(fData(out)[,colsIn] %>% as.matrix() %>% range, c(4,125))
  
  expect_true(check1 & check2)
  
})

test_that(desc="threshold_vectors_too_short", code={
  
  datIn <- tpptrData[[1]]
  colsIn <- c("qssm", "qupm")
  lbIn <- c(4,4)
  ubIn <- Inf
  
  expect_error(
    TPP:::filterOther(data = datIn, cols = colsIn, lb = lbIn, ub = ubIn)
  )
})

test_that(desc="column_not_exists", code={
  
  datIn <- tpptrData[[1]]
  colsIn <- c("qssm", "dummy")
  lbIn <- c(4, 4)
  ubIn <- c(Inf, Inf)
  
  expect_warning(
    TPP:::filterOther(data = datIn, cols = colsIn, lb = lbIn, ub = ubIn)
  )
})

test_that(desc="no_column_exists", code={
  
  datIn <- tpptrData[[1]]
  colsIn <- c("dummy")
  lbIn <- c(4)
  ubIn <- c(Inf)
  
  expect_warning(
    TPP:::filterOther(data = datIn, cols = colsIn, lb = lbIn, ub = ubIn)
  )
})


