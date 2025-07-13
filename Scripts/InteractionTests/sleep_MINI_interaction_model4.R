# --------------------------------------------------------------------
# sleep_cox_MINIinteraction_debugged_full.R
#   Fully adjusted Cox model + sleep_dur2 × MINIlifedepr.ph1 interaction
#   Ensures same sample for LR test by subsetting via complete.cases()
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist(), contrast(), Predict()
library(survival)   # Surv(), cox.zph()
library(lmtest)     # lrtest()

## 2 · Import & clean sentinel values → NA -----------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>%
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)))

## 3 · Factorise categorical variables --------------------------------
dat <- dat %>%
  mutate(
    Sex                   = factor(Sex, labels = c("Men","Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    mvpatile              = factor(mvpatile),
    MINIlifedepr.ph1      = factor(MINIlifedepr.ph1, labels = c("No","Yes")),
    sleep_dur2            = factor(sleep_dur2,
                                   levels = c(2,1,3),      # 2 = 7–9 h ref
                                   labels = c("7–9 h","<7 h","≥9 h")),
    frag_break_q4         = factor(frag_break_q4,
                                   levels = 1:4,
                                   labels = paste0("Q",1:4))
  )

cat("\nParticipants before subsetting:", nrow(dat),
    "\nIncident events:", sum(dat$LD_PHQ9depr_event, na.rm=TRUE), "\n")

## 4 · Subset to complete cases on all variables used in either model --
vars_needed <- c(
  # outcome/time
  "LD_PHQ9depr_Surv_time_standard_default", "LD_PHQ9depr_event",
  # key exposures
  "sleep_dur2", "MINIlifedepr.ph1",
  # other covariates
  "frag_break_q4", "Age.ph1", "Sex", "N_Diabetes_WHO.ph1",
  "n_education_3cat.ph1", "marital_status.ph1", "smoking_3cat.ph1",
  "N_alcohol_cat.ph1", "bmi.ph1", "N_CVD.ph1",
  "dhd_sum_min_alc.ph1", "mvpatile"
)
dat2 <- dat[ complete.cases(dat[, vars_needed]), ]
cat("Participants after subsetting:", nrow(dat2),
    "\nIncident events:", sum(dat2$LD_PHQ9depr_event), "\n")

## 5 · rms setup on cleaned data --------------------------------------
dd <- datadist(dat2)
options(datadist = "dd")

## 6 · Model without interaction (Model-4 + MINI) ---------------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 + MINIlifedepr.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2,
  x = TRUE, y = TRUE, surv = TRUE, parms = FALSE
)

## 7 · Model with sleep_dur2 × MINI interaction ----------------------
fit_int <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * MINIlifedepr.ph1 +   # main effects + interaction
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2,
  x = TRUE, y = TRUE, surv = TRUE, parms = FALSE
)

## 8 · Likelihood‐ratio test for interaction --------------------------
cat("\n--- LR test: interaction vs. no interaction ---\n")
print(lrtest(fit_base, fit_int))

## 8b · Inspect coefficient estimates ---------------------------------
cat("\n--- Coefficients (interaction model) ---\n")
print(summary(fit_int), digits = 3)   # full table with HRs & CIs

## 8c · Stratum-specific contrasts ------------------------------------
## --- 8c · Stratum-specific contrasts via coxph estimates -------------
# First refit the interaction model with survival::coxph
fit_int_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * MINIlifedepr.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2,
  ties = "efron"
)

# Extract coeffs and var–cov matrix
b <- coef(fit_int_coxph)
V <- vcov(fit_int_coxph)

# Helper to get HR, CI & p for main effects (MINI = No)
calc_main <- function(term_name) {
  β <- b[term_name]
  varβ <- V[term_name, term_name]
  se <- sqrt(varβ)
  z <- β / se
  data.frame(
    MINI  = "No",
    Sleep = sub("sleep_dur2", "", term_name),
    HR    = exp(β),
    Lower = exp(β - 1.96*se),
    Upper = exp(β + 1.96*se),
    p     = 2*(1 - pnorm(abs(z)))
  )
}

# Helper to get HR, CI & p for interaction strata (MINI = Yes)
calc_int <- function(main_term, int_term) {
  β1 <- b[main_term]
  β2 <- b[int_term]
  βsum <- β1 + β2
  var1 <- V[main_term, main_term]
  var2 <- V[int_term,   int_term]
  cov12 <- V[main_term, int_term]
  varsum <- var1 + var2 + 2*cov12
  se_sum <- sqrt(varsum)
  z <- βsum / se_sum
  data.frame(
    MINI  = "Yes",
    Sleep = sub("sleep_dur2", "", main_term),
    HR    = exp(βsum),
    Lower = exp(βsum - 1.96*se_sum),
    Upper = exp(βsum + 1.96*se_sum),
    p     = 2*(1 - pnorm(abs(z)))
  )
}

# Build the table
hr_table <- bind_rows(
  calc_main("sleep_dur2<7 h"),
  calc_main("sleep_dur2≥9 h"),
  calc_int("sleep_dur2<7 h",  "sleep_dur2<7 h:MINIlifedepr.ph1Yes"),
  calc_int("sleep_dur2≥9 h",  "sleep_dur2≥9 h:MINIlifedepr.ph1Yes")
)

# Print
cat("\n--- Stratum-specific HRs (with 95% CI) ---\n")
print(hr_table, digits = 3)

## 8d · Plot predicted survival by sleep × MINI -----------------------
plot(
  Predict(fit_int, sleep_dur2, MINIlifedepr.ph1, fun = exp),
  xlab = "Sleep duration",
  ylab = "Hazard ratio",
  main = "Sleep duration × MINI interaction"
)

## 9 · Proportional‐hazards diagnostics -------------------------------
cat("\n--- PH test (interaction model) ---\n")
print(cox.zph(fit_int))

## 10 · Clean up -------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
