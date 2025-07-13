# ---- 0.  Packages & data ----------------------------------------------
library(haven)
library(dplyr)

dat <- read_sav("Inclusion-Exclusion.sav")

covars <- c("Age.ph1","Sex","N_Diabetes_WHO.ph1",
            "marital_status.ph1","n_education_3cat.ph1",
            "dhd_sum.ph1","nit_alcoholtot.ph1","mean_mvpa_min_wake_t.ph1",
            "dhd_sum_min_alc.ph1","N_alcohol_cat.ph1","smoking_3cat.ph1",
            "N_CVD.ph1","N_HT.ph1",
            "bmi.ph1","OSBP.ph1","ODBP.ph1")

# ---- 1.  Recreate Steps 1–4 ------------------------------------------
analysis_pool <- dat %>%
  filter(!if_any(c(mean_inbed_min_night_t.ph1,
                   mean_inbed_break_night_t.ph1,
                   mean_up_trans_night_t.ph1), is.na)) %>%
  filter(!is.na(PHQ9_score.ph1)) %>%
  filter(PHQ9_score.ph1 < 10) %>%
  filter(LD_PHQ9depr_available == 1)

# ---- 2.  Sequential covariate exclusion ------------------------------
step_log <- vector("list", length(covars))   # to store results
remaining <- analysis_pool

for (i in seq_along(covars)) {
  var <- covars[i]
  miss_idx <- is.na(remaining[[var]])
  n_excl   <- sum(miss_idx)
  n_before <- nrow(remaining)
  n_after  <- n_before - n_excl
  
  step_log[[i]] <- data.frame(
    step        = i + 4,            # continues after original 4 steps
    covariate   = var,
    excluded    = n_excl,
    remaining   = n_after
  )
  
  remaining <- remaining[!miss_idx, ]
}

covar_flow <- dplyr::bind_rows(step_log)
print(covar_flow, row.names = FALSE)
