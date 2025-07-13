# --------------------------------------------------------------------
# sleep_cox_medfree.R
#   Spline–Cox analysis *excluding* baseline psychotropic users
#   – prints HR table, spline χ² tests, PH tests, event count,
#     and draws an HR curve (ref = 7 h)
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # data wrangling
library(rms)        # cph(), rcs(), Predict(), datadist()
library(survival)   # Surv(), cox.zph()
library(ggplot2)    # plotting

## 2 · Import & basic cleaning ----------------------------------------
dat <- read_sav("AV744_Seyedmahdi_Mirkialangaroodi_20250326 copy 2.sav") |>
  zap_labels()

dat <- dat |>
  mutate(
    Sex                   = factor(Sex),          # 1 = Men (ref), 2 = Women
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    sleepfrag_group.      = factor(sleepfrag_group),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    med_depression.ph1    = factor(med_depression.ph1),   # 0 = no, 1 = yes
    med_sleep.ph1         = factor(med_sleep.ph1)         # 0 = no, 1 = yes
  )

## 3 · Exclude baseline psychotropic users ----------------------------
dat_medfree <- dat |>
  filter(impaired_vibration_sense.ph1 == 0)

cat("\nParticipants retained:", nrow(dat_medfree),
    "\nIncident events:", sum(dat_medfree$LD_PHQ9depr_event), "\n")

## 4 · Tell rms about percentiles -------------------------------------
dd <- datadist(dat_medfree); options(datadist = "dd")

## 5 · 4-knot spline Cox model (Model-3 covariate set) -----------------
fit_medfree <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +    # minutes in bed spline
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + sleepfrag_group + dhd_sum_min_alc.ph1,
  data = dat_medfree,
  x = TRUE, y = TRUE
)

## 6 · Print coefficient table (HR, 95 % CI, Wald p) ------------------
cat("\n--- Model summary (psychotropic-free cohort) ---\n")
print(summary(fit_medfree), digits = 3)

## 7 · Spline block & non-linearity tests -----------------------------
cat("\n--- Likelihood-ratio χ² tests (overall & non-linear spline) ---\n")
print(anova(fit_medfree, test = "LR"))

## 8 · Proportional-hazards check -------------------------------------
cat("\n--- Proportional-hazards test (cox.zph) ---\n")
print(cox.zph(fit_medfree))

## 9 · Predicted HR curve (ref = 420 min = 7 h) -----------------------
pred_medfree <- Predict(fit_medfree,
                        mean_inbed_min_night_t.ph1 = seq(300, 720, 10),
                        fun = exp,
                        ref.zero = 420) |>
  as.data.frame()

ggplot(pred_medfree,
       aes(x = mean_inbed_min_night_t.ph1, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.20) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 420, linetype = "dashed") +
  scale_x_continuous("Minutes in bed per night",
                     breaks = seq(300, 720, 60),
                     sec.axis = dup_axis(~ ./60, name = "Hours")) +
  coord_cartesian(ylim = c(0.5, 3)) +
  labs(y = "Hazard ratio (95% CI, ref = 7 h/night)",
       title = "Restricted cubic-spline Cox model (psychotropic users excluded)") +
  theme_minimal(base_size = 14)

## 10 · Reset rms options ---------------------------------------------
options(datadist = NULL)

# --------------------------------------------------------------------
# End of script
# --------------------------------------------------------------------
