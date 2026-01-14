source("legend_pipeline.R")

pathToCMData <- "artifacts/CmData_l1_t261100000_c331100000.zip" # sema vs empa example
pathToPS <- "artifacts/Ps_IPTW_t261100000_c331100000_o1.rds" # for outcome 1: 3-pt MACE
pathToCensoringModel <- "artifacts/censoring_Cox__t261100000_c331100000_o1.rds"
pathToSurvWeights <- "artifacts/surv_weights_endCohortDate_t261100000_c331100000_o1.csv"

model_results <- run_legendt2dm_pipeline(
  pathToCMData = pathToCMData,
  outcomeOfInterestId = 1,
  pathToPS = pathToPS,
  pathToCensoringModel = pathToCensoringModel,
  pathToSurvWeights = pathToSurvWeights,
  show_progress = TRUE
)

model_results$estimates
