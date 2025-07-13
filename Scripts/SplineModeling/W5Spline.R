#=== 1) Load packages ============================================
library(haven)      # for read_sav() + zap_labels()
library(survival)   # for Surv() + cph()
library(rms)        # for datadist(), rcs(), Predict()

#=== 2) Import your SPSS data ===================================
dat <- read_sav("Week5.sav") |> 
  zap_labels()    # drop SPSS value‐labels so R sees pure numerics

#=== 3) Tell rms where to place knots ===========================
dd <- datadist(dat)
options(datadist="dd")

#=== 4) Fit the fully adjusted Cox model with a 4‐knot spline =====
fit_full <- cph(
  formula = Surv(
    LD_PHQ9depr_Surv_time_standard_default,
    LD_PHQ9depr_event
  ) ~
    rcs(mean_inbed_min_night_t.ph1, 4)     +   # 4‐knot spline on sleep
    frag_break_q4                          +   # fragmentation quartiles
    Age.ph1 + Sex + N_Diabetes_WHO.ph1     +   # Model 1 covariates
    n_education_3cat.ph1 + marital_status.ph1 +   # Model 2 additions
    smoking_3cat.ph1 + N_alcohol_cat.ph1   +   # Model 3 additions
    bmi.ph1 + N_CVD.ph1                    +   # clinical covariates
    dhd_sum_min_alc.ph1 + mvpatile,           # diet & physical activity
  data = dat,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

#=== 5) Plot HR vs sleep, anchored at 420 min (7 h) ==============
#=== 5) Plot HR vs sleep, anchored at 420 min (7 h) ==============
plot(
  Predict(
    fit_full,
    mean_inbed_min_night_t.ph1,
    ref.zero = 420    # anchor HR=1 at 7h (420 min)
  ),
  xlab = "Minutes in bed per night",
  ylab = "Hazard ratio (ref = 7 h)",
  main = "Fully adjusted Cox (RCS) – HR anchored at 7 h"
)

