#------------------------------------------------------------------------------#
# sleep_cox_full_CI.R — Full spline–Cox analysis with 95% CI ribbon
#------------------------------------------------------------------------------#

# 1) (Once) install required packages — then you can comment this out:
# install.packages(c("haven","rms","survival","splines","dplyr","ggplot2"))

# 2) Load libraries
library(haven)      # read_sav(), zap_labels()
library(rms)        # cph(), rcs(), Predict()
library(survival)   # Surv()
library(splines)    # ns() if you ever want base‐R splines
library(dplyr)      # data‐wrangling
library(ggplot2)    # plotting

# 3) (Optional) set working directory to where Week5.sav lives
# setwd("/Volumes/education/ThesisFeb/DMS_744_S_Mirkialangarooodi")

# 4) Import the SPSS dataset
dat <- read_sav("Week5.sav")

# 5) Strip SPSS labels so continuous vars stay numeric
dat <- zap_labels(dat)

# 6) Convert categorical codes to factors
dat <- dat %>% mutate(
  Sex                   = factor(Sex),                    # 1=Men,2=Women
  N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),     # 0=No,1=Type2
  n_education_3cat.ph1  = factor(n_education_3cat.ph1),   # 1=Low,2=Mid,3=High
  marital_status.ph1    = factor(marital_status.ph1),     # e.g. 1=Single,2=Partner…
  smoking_3cat.ph1      = factor(smoking_3cat.ph1),       # 0=Never,1=Former,2=Current
  N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),      # 0=None,1=Moderate,2=High
  N_CVD.ph1             = factor(N_CVD.ph1),              # 0=None,1=Yes
  mvpatile              = factor(mvpatile)                # 1=Low,2=Med,3=High
)

# 7) Tell rms about your data distribution (for knot placement)
dd <- datadist(dat)
options(datadist = "dd")

# 8) Fit the full restricted-cubic-spline Cox model
fit_full <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,
       LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +  # spline on minutes in bed
    frag_break_q4 +                       # fragmentation quartile
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +  # core covariates
    n_education_3cat.ph1 +                # education
    marital_status.ph1 +                  # marital status
    smoking_3cat.ph1 +                    # smoking
    N_alcohol_cat.ph1 +                   # alcohol use
    bmi.ph1 +                             # BMI
    N_CVD.ph1 +                           # CVD history
    dhd_sum_min_alc.ph1 +                 # diet score
    mvpatile,                             # MVPA tertile
  data = dat,
  x = TRUE, y = TRUE
)
print(summary(fit_full))

# 9) Use rms::Predict() to get HR + 95% CI, re-centered so HR=1 at 480 min
p_full <- Predict(
  fit_full,
  mean_inbed_min_night_t.ph1,
  ref.zero = TRUE,  # re-center so spline = 0 → HR = 1
  fun      = exp    # exponentiate to get HR
)
df_pred <- as.data.frame(p_full)
# df_pred has columns: mean_inbed_min_night_t.ph1, yhat, lower, upper

# 10) Plot with 95% CI ribbon
ggplot(df_pred, aes(x = mean_inbed_min_night_t.ph1, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x     = "Minutes in bed per night",
    y     = "Hazard ratio (95% CI)",
    title = "Full spline–Cox model with 95% confidence bands"
  ) +
  theme_minimal(base_size = 14)

# 11) OPTIONAL: If you’d rather use 7 h (420 min) as your reference—re-center manually:
p2    <- Predict(fit_full, mean_inbed_min_night_t.ph1, fun = exp)
df2   <- as.data.frame(p2)
ref420 <- df2$yhat[df2$mean_inbed_min_night_t.ph1 == 420]

df2 <- df2 %>% mutate(
  HR     = yhat / ref420,
  lower2 = lower / ref420,
  upper2 = upper / ref420
)

ggplot(df2, aes(x = mean_inbed_min_night_t.ph1, y = HR)) +
  geom_ribbon(aes(ymin = lower2, ymax = upper2), alpha = 0.2) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x        = "Minutes in bed per night",
    y        = "Hazard ratio (95% CI)",
    subtitle = "Re-centered at 7 h (420 min)",
    title    = "Full spline–Cox model"
  ) +
  theme_minimal(base_size = 14)

#------------------------------------------------------------------------------#
# End of script
#------------------------------------------------------------------------------#
