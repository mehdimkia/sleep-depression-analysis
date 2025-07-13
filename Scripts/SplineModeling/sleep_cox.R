#------------------------------------------------------------------------------#
# sleep_cox.R — U-shaped spline–Cox analysis of Week5.sav
#------------------------------------------------------------------------------#

# 1) (Once) install required packages, then you can comment this out:
# install.packages(c("haven","rms","survival","splines","dplyr","ggplot2"))

# 2) Load libraries
library(haven)      # read_sav()
library(rms)        # cph(), rcs()
library(survival)   # Surv()
library(splines)    # ns() if you prefer
library(dplyr)      # data prep
library(ggplot2)    # plotting

# 3) Set working directory to where Week5.sav lives
#setwd("/Volumes/education/ThesisFeb/DMS_744_S_Mirkialangarooodi")

# 4) Import the SPSS file
dat <- read_sav("Week5.sav")

# 5) Strip away SPSS labels (so continuous vars stay numeric)
dat <- zap_labels(dat)

# 6) Convert categorical codes to R factors
dat <- dat %>% mutate(
  Sex               = factor(Sex),               # 1=Men, 2=Women
  N_Diabetes_WHO.ph1 = factor(N_Diabetes_WHO.ph1) # 0=No, 1=Type2
)

# 7) Tell rms about your data distribution (for knot placement)
dd <- datadist(dat)
options(datadist = "dd")

# 8) Fit the restricted-cubic-spline Cox model
fit <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,  # follow-up time (years)
       LD_PHQ9depr_event)                       # event flag 0/1
  ~ rcs(mean_inbed_min_night_t.ph1, 4)        # spline on minutes in bed
  + Age.ph1                                  # age at visit1
  + Sex                                      # sex
  + N_Diabetes_WHO.ph1                       # diabetes status
  + frag_break_q4,                           # quartile of night breaks
  data = dat, x = TRUE, y = TRUE
)
print(summary(fit))

# 9) Build a plain data.frame for prediction
newdat <- data.frame(
  mean_inbed_min_night_t.ph1 = seq(300, 720, by = 10),
  Age.ph1                   = mean(dat$Age.ph1, na.rm = TRUE),
  Sex                       = factor(1, levels = levels(dat$Sex)),               # reference = 1
  N_Diabetes_WHO.ph1        = factor(0, levels = levels(dat$N_Diabetes_WHO.ph1)),# reference = 0
  frag_break_q4             = median(dat$frag_break_q4, na.rm = TRUE)
)

# Inspect it
str(newdat)
#> 'data.frame': 43 obs. of  5 variables:
#>  $ mean_inbed_min_night_t.ph1: num  300 310 320 ... 720
#>  $ Age.ph1                   : num  57 57 57 ... 57
#>  $ Sex                       : Factor w/ 2 levels "1","2": 1 1 1 ... 1
#>  $ N_Diabetes_WHO.ph1        : Factor w/ 2 levels "0","1": 1 1 1 ... 1
#>  $ frag_break_q4             : num  2 2 2 ... 2

# 10) Predict log‐hazard and re‐center at 480 min
lp     <- predict(fit, newdata = newdat, type = "lp")
ref_lp <- lp[newdat$mean_inbed_min_night_t.ph1 == 480]
newdat$HR <- exp(lp - ref_lp)

# 11) Plot
library(ggplot2)
ggplot(newdat, aes(x = mean_inbed_min_night_t.ph1, y = HR)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Minutes in bed per night",
    y = "Hazard ratio (depression)",
    title = "U-shaped spline Cox model"
  ) +
  theme_minimal(base_size = 14)
