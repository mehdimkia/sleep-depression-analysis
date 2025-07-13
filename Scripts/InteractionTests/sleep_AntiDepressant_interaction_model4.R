# --------------------------------------------------------------------
# sleep_cox_medInteraction.R
#   Fully adjusted Cox model + sleep_dur2 × med_depression.ph1 interaction
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
    med_depression.ph1    = factor(med_depression.ph1, labels = c("No","Yes")),
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
  "sleep_dur2", "med_depression.ph1",
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

## 6 · Model without interaction (Model-4 + med use) -----------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 + med_depression.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2,
  x = TRUE, y = TRUE, surv = TRUE, parms = FALSE
)

## 7 · Model with sleep × med use interaction ------------------------
fit_int <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * med_depression.ph1 +   # main effects + interaction
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
cat("\n--- LR test: interaction vs. no interaction (med use) ---\n")
print(lrtest(fit_base, fit_int))

## 8b · Inspect coefficient estimates ---------------------------------
cat("\n--- Interaction model coefficients (cph) ---\n")
print(summary(fit_int), digits = 3)

## 8c · Stratum-specific contrasts via coxph estimates ----------------
# Refit with survival::coxph for cross-product HRs
fit_int_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * med_depression.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, ties = "efron"
)

library(broom)
hr_med_int <- tidy(fit_int_coxph, exponentiate=TRUE, conf.int=TRUE) %>%
  filter(term %in% c("sleep_dur2<7 h:med_depression.ph1Yes",
                     "sleep_dur2≥9 h:med_depression.ph1Yes"))
cat("\n--- Cross‐product HRs (sleep × med use) ---\n")
print(hr_med_int, digits=3)

## 8d · Plot predicted HRs by sleep × med use -------------------------
plot(
  Predict(fit_int, sleep_dur2, med_depression.ph1, fun = exp),
  xlab = "Sleep duration",
  ylab = "Hazard ratio",
  main = "Sleep duration × Antidepressant use interaction"
)

## 9 · Proportional‐hazards diagnostics -------------------------------
cat("\n--- PH test (interaction model) ---\n")
print(cox.zph(fit_int))

## 10 · Clean up -------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
