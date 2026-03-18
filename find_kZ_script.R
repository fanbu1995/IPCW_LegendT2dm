source("legend_pipeline.R")

pathToCMData <- "artifacts/CmData_l1_t261100000_c331100000.zip" # sema vs empa example
pathToPS <- "artifacts/Ps_IPTW_t261100000_c331100000_o1.rds" # for outcome 1: 3-pt MACE
pathToCensoringModel <- "artifacts/censoring_Cox__t261100000_c331100000_o1.rds"
pathToSurvWeights <- "artifacts/surv_weights_endCohortDate_t261100000_c331100000_o1.csv"

ooi <- 1

outcomes_df.long <- compute_ipcw_weights(
  pathToCMData = pathToCMData,
  pathToPS = pathToPS,
  pathToCensoringModel = pathToCensoringModel,
  pathToSurvWeights = pathToSurvWeights,
  outcomeOfInterestId = ooi,
  show_progress = TRUE
)
