source("legend_pipeline.R")

pathToCMData <- "artifacts/CmData_l1_t261100000_c331100000.zip" # sema vs empa example
pathToCensoringModel <- "artifacts/censoring_Cox__t261100000_c331100000_o1.rds"

ooi <- 1 # focus on outcome 1 only

Cox_censoring <- fit_censoring_model(
  pathToCMData = pathToCMData,
  outcomeOfInterestId = ooi,
  pathToCensoringModel = pathToCensoringModel,
  show_progress = TRUE
)
