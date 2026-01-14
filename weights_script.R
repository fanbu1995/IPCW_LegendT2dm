# run after running find_kZ_script.R
library(dplyr)
library(data.table)
library(survival)

# outcomes_df.long = read.csv("survival_weights_endObsDate.csv") %>% data.table()

# distribution of weights
# summary(outcomes_df.long$Unstab_ipcw)
# summary(outcomes_df.long$Stab_ipcw)

# Truncate stabilized IPCW weights
lo <- quantile(outcomes_df.long$Stab_ipcw, 0.01, na.rm = TRUE)
hi <- quantile(outcomes_df.long$Stab_ipcw, 0.99, na.rm = TRUE)

outcomes_df.long$Stab_ipcw_trunc <- 
  pmin(pmax(outcomes_df.long$Stab_ipcw, lo), hi)


# truncate unstabilized IPCW weights
lo <- quantile(outcomes_df.long$Unstab_ipcw, 0.01, na.rm = TRUE)
hi <- quantile(outcomes_df.long$Unstab_ipcw, 0.99, na.rm = TRUE)

outcomes_df.long$Unstab_ipcw_trunc <- 
  pmin(pmax(outcomes_df.long$Unstab_ipcw, lo), hi)


# trucate iptw weights (alr stabilized)
outcomes_df.long <- outcomes_df.long %>%
  group_by(Tstart) %>%
  mutate(
    lower_iptw = quantile(iptw, 0.01, na.rm = TRUE),
    upper_iptw = quantile(iptw, 0.99, na.rm = TRUE),
    iptw_trunc = pmin(pmax(iptw, lower_iptw), upper_iptw)
  ) %>%
  ungroup()


# combined weights
outcomes_df.long$comb = outcomes_df.long$Stab_ipcw_trunc * outcomes_df.long$iptw_trunc


#####################
# FIT MODELS
#####################

# LegendT2dm sema vs empa for 3-pt MACE (before calibration) is: 
# CCAE: 
# OptumEHR: 

# UNADJUSTED
fit_unadjusted <- coxph(Surv(Tstart, time, ami) ~ treatment, data=outcomes_df.long, 
                        id=rowId)
summary(fit_unadjusted) # coef = 0.1344  SE = 0.1347; HR = 0.8743 (0.6714, 1.138)


### IPCW STABILIZED --- censoring only 
# fit_ipcw_Stab <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
#                       id=rowId, weights = Stab_ipcw) # unstab exp(coef) = 0.51 (1/0.51 = 1.95)
# summary(fit_ipcw_Stab) 

fit_ipcw_Stab_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                             id=rowId, weights = Stab_ipcw_trunc) # unstab exp(coef) = 0.51 (1/0.51 = 1.95)
summary(fit_ipcw_Stab_trunc) # coef = -0.1756, SE = 0.1543; HR = 0.839 (0.6326, 1.113)


# IPTW weights ====== CONFOUNDING

#fit_iptw <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
#                        id=rowId, weights = iptw) 
#summary(fit_iptw) 

fit_iptw_trunc <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                        id=rowId, weights = iptw_trunc) 
summary(fit_iptw_trunc) # coef = -0.02175; SE = 0.21581; HR = 0.9785 (0.6818, 1.404)


# combined weights: IPCW & IPTW
fit_comb <- coxph(Surv(Tstart, time, ami) ~ treatment , data=outcomes_df.long, 
                  id=rowId, weights = comb) 
summary(fit_comb) # coef = -0.08727 SE =  0.24035; HR = 0.9164 (0.6199, 1.355)
