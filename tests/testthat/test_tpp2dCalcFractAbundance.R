# filePath <- system.file("test_data", package="TPP")
# load(file.path(filePath, "panobinostat_2DTPP_smallExample.RData"))
# load(system.file("example_data/2D_example_data/fractAbundInData.RData", package="TPP"))

# test_that("evalFractAbund", code={
#   fracAbund <- tpp2dCalcFractAbundance(data = fractAbundInData)
#   expect_true(is.numeric(fracAbund$total_sumionarea_fraction))
#   expect_true(is.numeric(fracAbund$dmso1_vs_dmso2))
# })
