source("legend_pipeline.R")

pathToCMData <- "artifacts/CmData_l1_t261100000_c331100000.zip" # sema vs empa example
pathToPS <- "artifacts/Ps_IPTW_t261100000_c331100000_o1.rds" # for outcome 1: 3-pt MACE

ps <- run_ps_legend(
  pathToCMData = pathToCMData,
  outcomeOfInterestId = 1,
  pathToPS = pathToPS,
  show_progress = TRUE
)

