source("legend_pipeline.R")

#tcPairName = "t261100000_c331100000" # sema vs empa example
tcPairName = "t261100000_c321100000" # sema vs dapa instead

ooi = 1 # specify outcome id of interest; outcome 2: 4-pt MACE instead

db = "OptumDod" # database name 

if(!dir.exists(sprintf("artifacts_%s",db))){
  dir.create(sprintf("artifacts_%s",db))
}

# intermediate objects path -- across databases 
pathToCMData <- sprintf("../LegendT2dm_CMData/%s/CmData_l1_%s.zip", db, tcPairName) 
pathToPS <- sprintf("artifacts_%s/Ps_IPTW_%s_o%s.rds", db, tcPairName, ooi)
pathToCensoringModel <- sprintf("artifacts_%s/censoring_Cox_%s_o%s.rds", db, tcPairName, ooi)
pathToSurvWeights <- sprintf("artifacts_%s/surv_weights_endCohortDate_%s_o%s.csv", db, tcPairName, ooi)

# # intermediate objects path
# pathToCMData <- sprintf("artifacts/CmData_l1_%s.zip", tcPairName) 
# #pathToPS <- sprintf("artifacts/Ps_IPTW_%s_o1.rds", tcPairName) # for outcome 1: 3-pt MACE
# pathToPS <- sprintf("artifacts/Ps_IPTW_%s_o%s.rds", tcPairName, ooi) # outcome 2: 4-pt MACE instead
# pathToCensoringModel <- sprintf("artifacts/censoring_Cox_%s_o%s.rds", tcPairName, ooi)
# pathToSurvWeights <- sprintf("artifacts/surv_weights_endCohortDate_%s_o%s.csv", tcPairName, ooi)

# save estimate path
pathToEstimates = sprintf("artifacts_%s/estimates_%s_o%s.RDS", db, tcPairName, ooi) 

model_results <- run_legendt2dm_pipeline(
  pathToCMData = pathToCMData,
  outcomeOfInterestId = ooi,
  pathToPS = pathToPS,
  pathToCensoringModel = pathToCensoringModel,
  pathToSurvWeights = pathToSurvWeights,
  show_progress = TRUE
)

# save estimates locally as well 
saveRDS(model_results, pathToEstimates)

model_results$estimates
