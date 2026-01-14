## if necessary can run in terminal

library(DatabaseConnector)
library(CohortMethod)

library(dplyr)
library(Cyclops)

# install.packages("rJava", type="source")

# install.packages("DatabaseConnector") #6.3.2 was last version
#install.packages("CohortMethod") #5.3.0 was last version
#remotes::install_github("OHDSI/Cyclops", ref = 'v3.4.0') # 3.4.1 was last version

# library(knitr)
# library(tidyr)
# library(dplyr)
# library(Cyclops)

pathToCMData = "artifacts/CmData_l1_t261100000_c331100000.zip" # sema vs empa example 
pathToPS = "artifacts/Ps_IPTW_t261100000_c331100000_o1.rds" # for outcome 1: 3-pt MACE

cohortMethodData <- CohortMethod::loadCohortMethodData(pathToCMData)

studyPop <- CohortMethod::createStudyPopulation(
  cohortMethodData = cohortMethodData,
  outcomeId = 1, # 3-pt MACE
  firstExposureOnly = FALSE,
  restrictToCommonPeriod = FALSE,
  washoutPeriod = 0,
  removeDuplicateSubjects = "keep all",
  removeSubjectsWithPriorOutcome = FALSE,
  minDaysAtRisk = 1,
  riskWindowStart = 0,
  startAnchor = "cohort start",
  riskWindowEnd = 30,
  endAnchor = "cohort end"
)

# # create propensity model
# this takes quite a while
ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop) # this takes a while 
#saveRDS(ps, "temp_ps.rds")
saveRDS(ps, pathToPS)

## load trained PS model (with IPTW)
ps <- readRDS(pathToPS)

# outcome model

# unadjusted outcome model
outcomeModel <- fitOutcomeModel(population = ps,
                                modelType = "cox",
                                inversePtWeighting = TRUE)

# outcome model but with matching
matchedPop <- matchOnPs(ps, caliper = 0.2)

outcomeModel_matching <- fitOutcomeModel(population = matchedPop,
                                         modelType = "cox")

bal_matching = computeCovariateBalance(
  population = matchedPop,
  cohortMethodData,
  subgroupCovariateId = NULL,
  maxCohortSize = 250000,
  covariateFilter = NULL
)

# outcome model, with PS stratification? 
stratPop <- stratifyByPs(ps, numberOfStrata = 5) # default 5 strata
outcomeModel_strat <- fitOutcomeModel(population = stratPop,
                                         modelType = "cox")

bal_strat = computeCovariateBalance(
  population = stratPop,
  cohortMethodData,
  subgroupCovariateId = NULL,
  maxCohortSize = 250000,
  covariateFilter = NULL
)



