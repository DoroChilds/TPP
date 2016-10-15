# data(panobinostat_2DTPP_smallExample)
# load(system.file("example_data/2D_example_data/fractAbundInData.RData", package="TPP"))

# test_that("evalFractAbund", code={
#   fracAbund <- tpp2dCalcFractAbundance(configTable = panobinostat_2DTPP_config, 
#                                        data = fractAbundInData,
#                                      intensityStr = "sumionarea_protein_", 
#                                      idVar = "representative")
#   expect_true(is.numeric(fracAbund$total_sumionarea_fraction))
#   expect_true(is.numeric(fracAbund$dmso1_vs_dmso2))
# })