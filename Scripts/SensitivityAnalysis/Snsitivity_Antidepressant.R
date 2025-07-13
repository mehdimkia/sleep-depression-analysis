# --------------------------------------------------------------------
# sleep_cox_medexcl.R   
#   Cox model excluding participants with baseline antidepressant use
#   Uses categorical sleep duration (sleep_dur2) + fragmentation quartiles
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist()
library(survival)   # Surv(), cox.zph()

## 2 · Import & factorise ---------------------------------------------
dat <- read_sav("Week5.sav") |> zap_labels()

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

## 3 · Exclude baseline antidepressant users --------------------------
# med_depression.ph1 == 1 indicates antidepressant use at baseline
dat_medfree <- dat |> filter(med_depression.ph1 == 0)

cat("\nParticipants retained:", nrow(dat_medfree),
    "\nIncident events:", sum(dat_medfree$LD_PHQ9depr_event), "\n")

## 4 · rms setup -------------------------------------------------------
dd <- datadist(dat_medfree)
options(datadist = "dd")

## 5 · Fully adjusted Cox model ---------------------------------------
fit_medfree <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 +                        # categorical duration
    frag_break_q4 +                     # fragmentation quartiles
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1,
  data = dat_medfree,
  x = TRUE, y = TRUE, surv = TRUE, parms = FALSE
)

## 6 · Output HR table -------------------------------------------------
cat("\n--- Antidepressant-free: fully adjusted model ---\n")
print(summary(fit_medfree), digits = 3)

## 7 · PH diagnostics --------------------------------------------------
cat("\n--- Proportional-hazards test ---\n")
print(cox.zph(fit_medfree))

## 8 · Clean up --------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
