# --------------------------------------------------------------------
# sleep_cox_MINI_interact.R
#   Fully adjusted model + sleep_dur2 × lifetime MINI interaction
#   – works with rms::cph for modelling
#   – refits with survival::coxph to use broom::tidy()
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist()
library(survival)   # Surv(), coxph(), cox.zph()
library(lmtest)     # lrtest()
library(broom)      # tidy()

## 2 · Import & factorise ---------------------------------------------
dat <- read_sav("Week5.sav") |> zap_labels()

dat <- dat |>
  mutate(
    Sex                   = factor(Sex, labels = c("Men","Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    MINIlifedepr.ph1      = factor(MINIlifedepr.ph1, labels = c("No","Yes")),
    sleep_dur2            = factor(sleep_dur2,
                                   levels = c(2,1,3),          # 7–9 h ref
                                   labels  = c("7–9 h","<7 h","≥9 h")),
    frag_break_q4         = factor(frag_break_q4,
                                   levels = c(1,2,3,4),
                                   labels = paste0("Q",1:4))
  )

cat("\nParticipants:", nrow(dat),
    "\nIncident PHQ-9 events:", sum(dat$LD_PHQ9depr_event), "\n")

## 3 · rms setup -------------------------------------------------------
dd <- datadist(dat); options(datadist = "dd")

## 4 · Model 4 – main effects -----------------------------------------
form_base <- Surv(LD_PHQ9depr_Surv_time_standard_default,
                  LD_PHQ9depr_event) ~
               sleep_dur2 + frag_break_q4 +
               Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
               n_education_3cat.ph1 + marital_status.ph1 +
               smoking_3cat.ph1 + N_alcohol_cat.ph1 +
               bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
               MINIlifedepr.ph1

fit_main <- cph(form_base, data = dat, x = TRUE, y = TRUE)

## 5 · Interaction model (sleep_dur2 × MINI) --------------------------
form_int <- update(form_base, . ~ . + sleep_dur2:MINIlifedepr.ph1)

fit_int <- cph(form_int, data = dat, x = TRUE, y = TRUE)

## 6 · Likelihood-ratio test ------------------------------------------
cat("\n--- LR test: sleep_dur2 × MINI interaction ---\n")
print(lrtest(fit_main, fit_int))        # df = 2

## 7 · Coefficient table (rms::cph) -----------------------------------
cat("\n--- Interaction model coefficients (cph) ---\n")
print(summary(fit_int), digits = 3)

## 8 · Refitting with coxph for tidy() --------------------------------
fit_int_coxph <- coxph(form_int, data = dat, ties = "efron")

cat("\n--- Cross-product hazard ratios (coxph + broom) ---\n")
tidy(fit_int_coxph, exponentiate = TRUE, conf.int = TRUE) |>
  filter(term %in% c("sleep_dur2<7 h:MINIlifedepr.ph1Yes",
                     "sleep_dur2≥9 h:MINIlifedepr.ph1Yes")) |>
  print(digits = 3)

## 9 · Proportional-hazards diagnostics -------------------------------
cat("\n--- Proportional-hazards test (interaction model) ---\n")
print(cox.zph(fit_int))

## 10 · Clean up -------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
