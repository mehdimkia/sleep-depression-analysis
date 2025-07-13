#------------------------------------------------------------------------------#
# sleep_cox_full_with_CI.R — Full spline–Cox analysis of Week5.sav (Model 3)
#                              with 95% confidence intervals
#------------------------------------------------------------------------------#

# 1) Install (once) –– then you can comment this out:
# install.packages(c("haven","rms","survival","splines","dplyr","ggplot2"))

# 2) Load libraries
library(haven)      # read_sav(), zap_labels()
library(rms)        # cph(), rcs()
library(survival)   # Surv()
library(splines)    # ns()
library(dplyr)      # data‐wrangling
library(ggplot2)    # plotting

# 3) (Optional) set your working directory
# setwd("/Volumes/education/ThesisFeb/DMS_744_S_Mirkialangarooodi")

# 4) Import the SPSS dataset
dat <- read_sav("Week5.sav")

# 5) Strip SPSS labels (so numeric stays numeric)
dat <- zap_labels(dat)

# 6) Convert all categorical codes into R factors
dat <- dat %>% mutate(
  Sex                  = factor(Sex),                    # 1=Men, 2=Women
  N_Diabetes_WHO.ph1   = factor(N_Diabetes_WHO.ph1),     # 0=No, 1=Type2
  n_education_3cat.ph1 = factor(n_education_3cat.ph1),   # 1=Low,2=Mid,3=High
  marital_status.ph1   = factor(marital_status.ph1),     # e.g. 1=Single,2=Partner…
  smoking_3cat.ph1     = factor(smoking_3cat.ph1),       # 0=Never,1=Former,2=Current
  N_alcohol_cat.ph1    = factor(N_alcohol_cat.ph1),      # 0=None,1=Moderate,2=High
  N_CVD.ph1            = factor(N_CVD.ph1),              # 0=None,1=Yes
  mvpatile             = factor(mvpatile)                # 1=Low,2=Med,3=High
)

# 7) Inform rms about your data (for knot placement)
dd <- datadist(dat)
options(datadist = "dd")

# 8) Fit the full spline–Cox model
fit_full <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,
       LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +    # spline on minutes in bed
    frag_break_q4 +                         # sleep fragmentation quartile
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +    # core covariates
    n_education_3cat.ph1 +                  # education level
    marital_status.ph1 +                    # marital status
    smoking_3cat.ph1 +                      # smoking status
    N_alcohol_cat.ph1 +                     # alcohol use category
    bmi.ph1 +                               # body‐mass index
    N_CVD.ph1 +                             # history of CVD
    dhd_sum_min_alc.ph1 +                   # diet score
    mvpatile,                               # MVPA tertile
  data = dat,
  x = TRUE, y = TRUE
)
print(summary(fit_full))

# 9) Build data.frame for prediction at typical covariate values
newdat_full <- data.frame(
  mean_inbed_min_night_t.ph1 = seq(300, 720, by = 10),
  frag_break_q4             = median(dat$frag_break_q4, na.rm = TRUE),
  Age.ph1                   = mean(dat$Age.ph1, na.rm = TRUE),
  Sex                       = factor(1, levels = levels(dat$Sex)),               # Men
  N_Diabetes_WHO.ph1        = factor(0, levels = levels(dat$N_Diabetes_WHO.ph1)),# No diabetes
  n_education_3cat.ph1      = factor(1, levels = levels(dat$n_education_3cat.ph1)),  # Low
  marital_status.ph1        = factor(2, levels = levels(dat$marital_status.ph1)),   # Partner
  smoking_3cat.ph1          = factor(0, levels = levels(dat$smoking_3cat.ph1)),     # Never
  N_alcohol_cat.ph1         = factor(0, levels = levels(dat$N_alcohol_cat.ph1)),    # None
  bmi.ph1                   = mean(dat$bmi.ph1, na.rm = TRUE),
  N_CVD.ph1                 = factor(0, levels = levels(dat$N_CVD.ph1)),             # No CVD
  dhd_sum_min_alc.ph1       = mean(dat$dhd_sum_min_alc.ph1, na.rm = TRUE),
  mvpatile                  = factor(1, levels = levels(dat$mvpatile))              # Low MVPA
)

# 10) Predict linear predictor + SE, then re-center & compute 95% CI
pred <- predict(
  fit_full,
  newdata = newdat_full,
  type    = "lp",
  se.fit  = TRUE
)

newdat_full <- newdat_full %>%
  mutate(
    lp      = pred$fit,
    se      = pred$se.fit,
    # reference at 480 min (8 h):
    logHR   = lp - lp[which(newdat_full$mean_inbed_min_night_t.ph1 == 480)],
    HR      = exp(logHR),
    lowerCI = exp(logHR - 1.96 * se),
    upperCI = exp(logHR + 1.96 * se)
  )

# 11) Plot U‐shaped curve with 95% CI ribbon
ggplot(newdat_full, aes(x = mean_inbed_min_night_t.ph1, y = HR)) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = 0.2) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x     = "Minutes in bed per night",
    y     = "Hazard ratio of depression",
    title = "Full spline–Cox Model (with 95% CI)"
  ) +
  theme_minimal(base_size = 14)

#------------------------------------------------------------------------------#
# End of script
#------------------------------------------------------------------------------#
