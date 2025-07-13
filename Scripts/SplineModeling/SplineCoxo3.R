#--------------------------------------------------------------------
#  full_sleep_spline.R  –  Model-3 spline Cox curve (ref = 7 h)
#--------------------------------------------------------------------

## 1.  Load libraries  -------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(rms)        # cph(), rcs(), Predict()
library(survival)   # Surv()
library(dplyr)      # data wrangling
library(ggplot2)    # plotting

## 2.  Read and prepare data  -----------------------------------------
dat <- read_sav("Week5.sav") |> zap_labels()

dat <- dat |>
  mutate(
    Sex                   = factor(Sex),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    mvpatile              = factor(mvpatile),
    frag_break_q4         = factor(frag_break_q4)      # NEW: treat as categorical
  )

dd <- datadist(dat);  options(datadist = "dd")         # tells rms where to find percentiles

## 3.  Fit fully adjusted Cox model with 4-knot restricted cubic spline
fit_full <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,
       LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat, x = TRUE, y = TRUE
)

## 4.  Predict HRs every 10 min from 300 to 720 min -------------------
newgrid <- data.frame(mean_inbed_min_night_t.ph1 = seq(300, 720, by = 10))
pred    <- Predict(fit_full, mean_inbed_min_night_t.ph1 = newgrid$mean_inbed_min_night_t.ph1,
                   fun = exp) |> as.data.frame()

# Re-centre so HR = 1 at 7 h (420 min)
ref   <- pred$yhat[pred$mean_inbed_min_night_t.ph1 == 420]
pred  <- pred |>
  mutate(HR    = yhat  / ref,
         lower = lower / ref,
         upper = upper / ref)

## 5.  Plot ------------------------------------------------------------
ggplot(pred, aes(mean_inbed_min_night_t.ph1, HR)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 420, linetype = "dashed") +
  scale_x_continuous(
    "Minutes in bed per night",
    breaks = seq(300, 720, 60),
    sec.axis = dup_axis(~ ./60, name = "Hours")
  ) +
  coord_cartesian(ylim = c(0.5, 3)) +
  labs(
    y = "Hazard ratio (95% CI, ref = 7 h/night)",
    title = "Restricted cubic-spline Cox model – fully adjusted"
  ) +
  theme_minimal(base_size = 14)

#--------------------------------------------------------------------
# End of script
#--------------------------------------------------------------------
