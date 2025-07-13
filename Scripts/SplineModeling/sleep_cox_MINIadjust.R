# --------------------------------------------------------------------
# sleep_cox_MINIadjust.R
#   Cox model *adjusting* for lifetime MINI depression (Model-4)
#   Categorical sleep duration (sleep_dur2) + fragmentation quartiles
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph()
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
    mvpatile              = factor(mvpatile),
    MINIlifedepr.ph1      = factor(MINIlifedepr.ph1, labels = c("No","Yes")),
    sleep_dur2            = factor(sleep_dur2,
                                   levels = c(2,1,3),      # 2 = 7–9 h ref
                                   labels  = c("7–9 h","<7 h","≥9 h")),
    frag_break_q4         = factor(frag_break_q4,
                                   levels = c(1,2,3,4),
                                   labels = paste0("Q",1:4))
  )

cat("\nParticipants analysed:", nrow(dat),
    "\nIncident events:", sum(dat$LD_PHQ9depr_event), "\n")

## 3 · rms setup -------------------------------------------------------
dd <- datadist(dat); options(datadist = "dd")

## 4 · Fully adjusted Cox model (Model-4 + MINI) -----------------------
fit_MINIadj <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 +                         # categorical duration
    frag_break_q4 +                      # fragmentation quartiles
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile + MINIlifedepr.ph1,                    # NEW covariate
  data = dat,
  x = TRUE, y = TRUE
)

## 5 · Output HR table -------------------------------------------------
cat("\n--- Model 4 + lifetime MINI adjustment ---\n")
print(summary(fit_MINIadj), digits = 3)

## 6 · PH diagnostics --------------------------------------------------
cat("\n--- Proportional-hazards test ---\n")
print(cox.zph(fit_MINIadj))

## 7 · Clean up --------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
