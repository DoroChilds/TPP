# Load function input:
data(hdacTR_smallExample)
data(hdacCCR_smallExample)

# Start tests:
test_that(desc="add_rows_from_non_unique_entries_in_single_id_column", code={
  cfgIn <- hdacTR_config
  datIn <- hdacTR_data
  #nonUnique <- datIn[[1]]$gene_name %>% grep("\\|", ., value = TRUE)
  # For test purposes only:
  singleDat <- datIn[[1]]
  ids <- c(singleDat[["gene_name"]] %>% grep("\\|", ., value = T), "HDAC1") %>% sort
  singleDat <-  singleDat %>% filter(gene_name %in% ids) %>%
   bind_rows(singleDat %>% filter(gene_name == "CALM2|CALM1") %>%
               mutate(gene_name = "CALM2")) %>%
    arrange(gene_name)
  ids <- singleDat[["gene_name"]]
  datIn[[1]] <- singleDat
  new <- tpptrImport(configTable = cfgIn, data = datIn)[[1]] %>% 
    Biobase::featureNames() %>%
    setdiff("HDAC1")
  ref <-   uniques <- unique(unlist(strsplit(grep("\\|", ids, value = TRUE), "\\|"))) %>% sort
  expect_equal(new, ref)
})
