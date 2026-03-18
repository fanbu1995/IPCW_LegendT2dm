source("legend_pipeline.R")

pathToSurvWeights <- "artifacts/surv_weights_endCohortDate_t261100000_c331100000_o1.csv"

if (!exists("outcomes_df.long")) {
  outcomes_df.long <- read.csv(pathToSurvWeights)
}

model_results <- fit_weighted_cox_models(outcomes_df.long)
model_results$estimates
