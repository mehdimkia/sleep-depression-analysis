# --------------------------------------------------------------------
# sleep_cox_sensitivity_vibration.R   
#   Cox model excluding participants with impaired vibration sense
#   Uses categorical sleep duration (sleep_dur2) + fragmentation quartiles
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist()
library(survival)   # Surv(), cox.zph()

## 2 · Import & factorise ---------------------------------------------
dat <- read_sav("Week7.sav") |> zap_labels()

dat <- dat |>
  mutate(
    Sex                   = factor(Sex, labels = c("Men", "Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    sleep_dur2            = factor(sleep_dur2,
                                   levels = c(2,1,3),      # 2 = 7–9 h ref
                                   labels = c("7–9 h","<7 h","≥9 h")),
    frag_break_q4         = factor(frag_break_q4,           # Q1 ref
                                   levels = c(1,2,3,4),
                                   labels = paste0("Q",1:4))
  )

## 3 · Exclude participants with impaired vibration sense --------------
# impaired_vibration_sense.ph1 == 1 indicates impairment in one/both toes
dat_vibfree <- dat |> filter(impaired_vibration_sense.ph1 == 0)

cat("\nParticipants retained:", nrow(dat_vibfree),
    "\nIncident events:", sum(dat_vibfree$LD_PHQ9depr_event), "\n")

## 4 · rms setup -------------------------------------------------------
dd <- datadist(dat_vibfree)
options(datadist = "dd")

## 5 · Fully adjusted Cox model ---------------------------------------
fit_vibfree <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 +                        # categorical duration
    frag_break_q4 +                     # fragmentation quartiles
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1,
  data = dat_vibfree,
  x = TRUE, y = TRUE, surv = TRUE, parms = FALSE
)

## 6 · Output HR table -------------------------------------------------
cat("\n--- Vibration-sense-free: fully adjusted model ---\n")
print(summary(fit_vibfree), digits = 3)

## 7 · PH diagnostics --------------------------------------------------
cat("\n--- Proportional-hazards test ---\n")
print(cox.zph(fit_vibfree))

## 8 · Clean up --------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
