# dataPath <- system.file("example_data", package="TPP")
# load(file.path(dataPath, "2D_example_data/referenceImportedData.RData"))
# data(panobinostat_2DTPP_smallExample)
# 
# # test_that(desc="evalImportData", code={
# #   configTable <- panobinostat_2DTPP_config
# #   data <- panobinostat_2DTPP_data
# #   testData <- tpp2dImport(configTable=configTable, data=data, idVar="representative", 
# #                               intensityStr="sumionarea_protein_", 
# #                               addCol=c("qusm","clustername", "msexperiment_id"), 
# #                               qualColName="qupm", fcStr=NULL)
# #   expect_identical(testData, refData)
# # })
# 
# test_that(desc="evalImportDataErr1", code={
#   configTable <- panobinostat_2DTPP_config
#   data <- NULL
#   expect_error(tpp2dImport(configTable=configTable, data=data, idVar="representative", 
#                                intensityStr="sumionarea_protein_", 
#                                addCol=c("qusm","clustername", "msexperiment_id"), 
#                                qualColName="qupm", fcStr=NULL))
# })
# 
# test_that(desc="evalImportDataErr2", code={
#   configTable <- NULL
#   data <- panobinostat_2DTPP_data
#   expect_error(tpp2dImport(configTable=configTable, data=data, idVar="representative", 
#                                intensityStr="sumionarea_protein_", addCol="clustername", 
#                                qualColName="qupm"))
# })