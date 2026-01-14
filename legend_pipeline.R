library(CohortMethod)
library(Cyclops)
library(dplyr)
library(progress)
library(survival)
library(tibble)
library(tidyr)

validate_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path))
  }
}

progress_message <- function(message_text, show_progress = TRUE) {
  if (isTRUE(show_progress)) {
    message(message_text)
  }
}

run_ps_legend <- function(pathToCMData,
                          outcomeOfInterestId = 1,
                          pathToPS = NULL,
                          show_progress = TRUE) {
  validate_file_exists(pathToCMData, "CohortMethod data")
  if (!is.null(pathToPS) && file.exists(pathToPS)) {
    progress_message("Step 1/1: Loading existing propensity score model.", show_progress)
    return(readRDS(pathToPS))
  }
  progress_message("Step 1/4: Loading CohortMethod data.", show_progress)
  cohortMethodData <- CohortMethod::loadCohortMethodData(pathToCMData)

  progress_message("Step 2/4: Creating study population.", show_progress)
  studyPop <- CohortMethod::createStudyPopulation(
    cohortMethodData = cohortMethodData,
    outcomeId = outcomeOfInterestId,
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

  progress_message("Step 3/4: Fitting propensity score model.", show_progress)
  ps <- createPs(cohortMethodData = cohortMethodData, population = studyPop)

  if (!is.null(pathToPS)) {
    saveRDS(ps, pathToPS)
  }

  progress_message("Step 4/4: Propensity score model complete.", show_progress)
  ps
}

fit_censoring_model <- function(pathToCMData,
                                outcomeOfInterestId = 1,
                                pathToCensoringModel = NULL,
                                show_progress = TRUE) {
  validate_file_exists(pathToCMData, "CohortMethod data")
  progress_message("Step 1/3: Loading CohortMethod data for censoring model.", show_progress)
  cohortMethodData <- CohortMethod::loadCohortMethodData(pathToCMData)

  progress_message("Step 2/3: Building censoring outcomes table.", show_progress)
  cohortMethodData$outcomeOfInterest <- cohortMethodData$outcomes %>%
    filter(outcomeId == outcomeOfInterestId)
  cohortMethodData$censored_cohort <- cohortMethodData$cohorts %>%
    anti_join(cohortMethodData$outcomeOfInterest, by = "rowId")

  cohortMethodData$censored_outcomes <- cohortMethodData$censored_cohort %>%
    select(
      rowId,
      daysToEvent = daysToCohortEnd
    ) %>%
    mutate(outcomeId = 1) %>%
    select(rowId, outcomeId, daysToEvent)

  cohortMethodData$censored_outcomes <- cohortMethodData$censored_outcomes %>%
    filter(outcomeId == 1)

  cohorts <- cohortMethodData$cohorts
  censored_outcomes <- cohortMethodData$censored_outcomes

  outcomes_for_cyclops <- cohorts %>%
    left_join(censored_outcomes, by = "rowId") %>%
    mutate(y = ifelse(!is.na(outcomeId), 1, 0)) %>%
    mutate(time = ifelse(!is.na(outcomeId), daysToEvent, daysToObsEnd)) %>%
    select(rowId, y, time)

  cohortMethodData$outcomes_for_cyclops <- outcomes_for_cyclops

  censored_df <- convertToCyclopsData(
    cohortMethodData$outcomes_for_cyclops,
    cohortMethodData$covariates,
    modelType = "cox",
    addIntercept = TRUE
  )

  progress_message("Step 3/3: Fitting censoring Cox model.", show_progress)
  lassoPrior <- Cyclops::createPrior(
    priorType = "laplace",
    useCrossValidation = TRUE
  )

  Cox_censoring <- fitCyclopsModel(censored_df, prior = lassoPrior)

  if (!is.null(pathToCensoringModel)) {
    saveRDS(Cox_censoring, pathToCensoringModel)
  }

  Cox_censoring
}

breslow_est <- function(time, status, X, B, show_progress = TRUE) {
  data <- data.frame(time, status, X)
  data <- data[order(data$time), ]
  t <- unique(data$time)
  k <- length(t)
  h <- rep(0, k)
  LP_indiv <- X %*% B

  pb <- progress::progress_bar$new(
    format = "  Progress [:bar] :percent in :elapsed, ETA: :eta",
    total = length(1:k),
    clear = FALSE,
    width = 60
  )

  for (i in 1:k) {
    if (isTRUE(show_progress)) {
      pb$tick()
    }
    lp <- (LP_indiv)[data$time >= t[i]]
    risk <- exp(lp)
    h[i] <- sum(data$status[data$time == t[i]]) / sum(risk)
  }

  cumsum(h)
}

transform_to_long <- function(data, cut.times) {
  data$Tstart <- 0
  data$event <- 1 - data$y

  data.long <- survSplit(
    data = data,
    cut = cut.times,
    end = "time",
    start = "Tstart",
    event = "y"
  )
  data.long <- data.long[order(data.long$rowId, data.long$time), ]

  data.long.event <- survSplit(
    data,
    cut = cut.times,
    end = "time",
    start = "Tstart",
    event = "event"
  )
  data.long.event <- data.long.event[order(data.long.event$rowId, data.long.event$time), ]
  data.long$event <- data.long.event$event
  data.long$rowId <- as.numeric(data.long$rowId)
  data.long
}

compute_ipcw_weights <- function(pathToCMData,
                                 pathToPS,
                                 pathToCensoringModel,
                                 pathToSurvWeights = NULL,
                                 outcomeOfInterestId = 1,
                                 cut.times = NULL,
                                 show_progress = TRUE) {
  validate_file_exists(pathToCMData, "CohortMethod data")
  validate_file_exists(pathToPS, "Propensity score model")
  validate_file_exists(pathToCensoringModel, "Censoring model")

  progress_message("Step 1/6: Loading saved PS and CohortMethod data.", show_progress)
  ps <- readRDS(pathToPS)
  cohortMethodData <- CohortMethod::loadCohortMethodData(pathToCMData)
  cohortMethodData$outcomeOfInterest <- cohortMethodData$outcomes %>%
    filter(outcomeId == outcomeOfInterestId)

  progress_message("Step 2/6: Aligning cohorts to PS analysis population.", show_progress)
  cohort_ids <- ps$rowId
  cohortMethodData$cohorts <- collect(cohortMethodData$cohorts) %>%
    filter(rowId %in% cohort_ids)
  cohortMethodData$covariates <- collect(cohortMethodData$covariates) %>%
    filter(rowId %in% cohort_ids)
  cohortMethodData$outcomeOfInterest <- collect(cohortMethodData$outcomeOfInterest) %>%
    filter(rowId %in% cohort_ids)

  outcomes_for_cyclops <- ps %>%
    transmute(
      rowId = rowId,
      y = if_else(outcomeCount == 0, 1, 0),
      time = survivalTime,
      treatment = treatment,
      iptw = iptw
    )

  progress_message("Step 3/6: Preparing covariates for censoring model.", show_progress)
  Cox_censoring <- readRDS(pathToCensoringModel)
  coefs <- coef(Cox_censoring)
  non_zero_coefs <- coefs[coefs != 0]

  filtered_covariates <- collect(cohortMethodData$covariates) %>%
    filter(covariateId %in% names(non_zero_coefs))

  X_wide <- filtered_covariates %>%
    select(rowId, covariateId, covariateValue) %>%
    mutate(covariateId = as.character(covariateId)) %>%
    pivot_wider(
      names_from = covariateId,
      values_from = covariateValue,
      values_fill = 0
    )

  X_matrix <- X_wide %>%
    column_to_rownames(var = "rowId") %>%
    as.matrix()

  coeff_vector <- non_zero_coefs[colnames(X_matrix)]
  stopifnot(length(coeff_vector) == ncol(X_matrix))
  stopifnot(all(names(coeff_vector) == colnames(X_matrix)))

  progress_message("Step 4/6: Computing baseline hazards for censoring model.", show_progress)
  outcomes_df <- data.table::data.table(collect(outcomes_for_cyclops))
  outcomes_df <- outcomes_df %>%
    arrange(match(as.character(rowId), rownames(X_matrix)))

  cox_baseline <- coxph(Surv(time, y) ~ 1, data = outcomes_df)
  baseline_surv <- survfit(cox_baseline)
  baseline_surv_fun <- stepfun(baseline_surv$time, c(1, baseline_surv$surv))

  H0 <- breslow_est(
    time = outcomes_df$time,
    status = outcomes_df$y,
    X = X_matrix,
    B = coeff_vector,
    show_progress = show_progress
  )
  haz_step_fun <- stepfun(sort(unique(outcomes_df$time)), c(0, H0))

  if (is.null(cut.times)) {
    dist <- summary(outcomes_df$time)
    cut.times <- c(
      seq(from = 1, to = dist[2], by = 5),
      seq(from = dist[2] + 5, to = dist[3], by = 10),
      seq(from = dist[3] + 10, to = dist[5], by = 30),
      seq(from = dist[5] + 30, to = floor(max(outcomes_df$time) / 100) * 100, by = 100)
    )
  }

  progress_message("Step 5/6: Transforming to long format and computing IPCW.", show_progress)
  outcomes_df.long <- transform_to_long(outcomes_df, cut.times)
  eta <- X_matrix %*% coeff_vector
  outcomes_df.long$H0_Tstart <- haz_step_fun(outcomes_df.long$Tstart)

  rowid_to_eta <- data.frame(rowId = as.numeric(rownames(X_matrix)), eta = as.numeric(eta))
  outcomes_df.long <- outcomes_df.long %>%
    left_join(rowid_to_eta, by = "rowId") %>%
    mutate(KZ = exp(-H0_Tstart * exp(eta)))

  outcomes_df.long$K0_ti <- baseline_surv_fun(outcomes_df.long$Tstart)
  outcomes_df.long$Unstab_ipcw <- 1 / outcomes_df.long$KZ
  outcomes_df.long$Stab_ipcw <- outcomes_df.long$K0_ti / outcomes_df.long$KZ

  if (!is.null(pathToSurvWeights)) {
    write.csv(outcomes_df.long, pathToSurvWeights, row.names = FALSE)
  }

  progress_message("Step 6/6: IPCW weights computed.", show_progress)
  outcomes_df.long
}

prepare_weighted_dataset <- function(outcomes_df.long) {
  lo <- quantile(outcomes_df.long$Stab_ipcw, 0.01, na.rm = TRUE)
  hi <- quantile(outcomes_df.long$Stab_ipcw, 0.99, na.rm = TRUE)
  outcomes_df.long$Stab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Stab_ipcw, lo), hi)

  lo <- quantile(outcomes_df.long$Unstab_ipcw, 0.01, na.rm = TRUE)
  hi <- quantile(outcomes_df.long$Unstab_ipcw, 0.99, na.rm = TRUE)
  outcomes_df.long$Unstab_ipcw_trunc <- pmin(pmax(outcomes_df.long$Unstab_ipcw, lo), hi)

  outcomes_df.long <- outcomes_df.long %>%
    group_by(Tstart) %>%
    mutate(
      lower_iptw = quantile(iptw, 0.01, na.rm = TRUE),
      upper_iptw = quantile(iptw, 0.99, na.rm = TRUE),
      iptw_trunc = pmin(pmax(iptw, lower_iptw), upper_iptw)
    ) %>%
    ungroup()

  outcomes_df.long$comb <- outcomes_df.long$Stab_ipcw_trunc * outcomes_df.long$iptw_trunc
  outcomes_df.long
}

summarize_cox_fit <- function(fit, label) {
  fit_summary <- summary(fit)
  data.frame(
    model = label,
    coef = fit_summary$coef["treatment", "coef"],
    se_coef = fit_summary$coef["treatment", "se(coef)"],
    hr = fit_summary$conf.int["treatment", "exp(coef)"],
    hr_ci_lower = fit_summary$conf.int["treatment", "lower .95"],
    hr_ci_upper = fit_summary$conf.int["treatment", "upper .95"],
    row.names = NULL
  )
}

fit_weighted_cox_models <- function(outcomes_df.long) {
  outcomes_df.long <- prepare_weighted_dataset(outcomes_df.long)

  fit_unadjusted <- coxph(Surv(Tstart, time, event) ~ treatment, data = outcomes_df.long, id = rowId)
  fit_ipcw <- coxph(
    Surv(Tstart, time, event) ~ treatment,
    data = outcomes_df.long,
    id = rowId,
    weights = Stab_ipcw_trunc
  )
  fit_iptw <- coxph(
    Surv(Tstart, time, event) ~ treatment,
    data = outcomes_df.long,
    id = rowId,
    weights = iptw_trunc
  )
  fit_comb <- coxph(
    Surv(Tstart, time, event) ~ treatment,
    data = outcomes_df.long,
    id = rowId,
    weights = comb
  )

  estimates <- bind_rows(
    summarize_cox_fit(fit_unadjusted, "Unadjusted Cox"),
    summarize_cox_fit(fit_ipcw, "Cox with IPCW"),
    summarize_cox_fit(fit_iptw, "Cox with IPTW"),
    summarize_cox_fit(fit_comb, "Cox with IPCW + IPTW")
  )

  list(
    estimates = estimates,
    fits = list(
      unadjusted = fit_unadjusted,
      ipcw = fit_ipcw,
      iptw = fit_iptw,
      combined = fit_comb
    ),
    outcomes = outcomes_df.long
  )
}

run_legendt2dm_pipeline <- function(pathToCMData,
                                    outcomeOfInterestId = 1,
                                    pathToPS,
                                    pathToCensoringModel,
                                    pathToSurvWeights,
                                    show_progress = TRUE) {
  ps <- run_ps_legend(
    pathToCMData = pathToCMData,
    outcomeOfInterestId = outcomeOfInterestId,
    pathToPS = pathToPS,
    show_progress = show_progress
  )

  fit_censoring_model(
    pathToCMData = pathToCMData,
    outcomeOfInterestId = outcomeOfInterestId,
    pathToCensoringModel = pathToCensoringModel,
    show_progress = show_progress
  )

  outcomes_df.long <- compute_ipcw_weights(
    pathToCMData = pathToCMData,
    pathToPS = pathToPS,
    pathToCensoringModel = pathToCensoringModel,
    pathToSurvWeights = pathToSurvWeights,
    outcomeOfInterestId = outcomeOfInterestId,
    show_progress = show_progress
  )

  model_results <- fit_weighted_cox_models(outcomes_df.long)
  model_results$ps <- ps
  model_results
}
