# --------------------------------------------------------------------
# sleep_spline_interactions.R
#   Test effect modification of the sleep–depression spline by sex
#   and by type-2 diabetes status
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>%, mutate()
library(rms)        # cph(), datadist(), rcs()
library(survival)   # Surv()
library(lmtest)     # lrtest()

## 2 · Import & clean sentinel → NA -----------------------------------
dat <- read_sav("Week7.sav") |> 
  zap_labels() %>%
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)
  ))

## 3 · Derive continuous sleep in hours or minutes --------------------
# here we use minutes in bed directly:
# mean_inbed_min_night_t.ph1 must already exist
# If you prefer hours: uncomment the next line
# dat <- dat %>% mutate(sleep_hours = mean_inbed_min_night_t.ph1 / 60)

## 4 · Factorise covariates --------------------------------------------
dat <- dat %>%
  mutate(
    Sex     = factor(Sex,         levels=1:2, labels=c("Men","Women")),
    # recode type-2 diabetes: WHO code 2 = T2D
    T2diab  = factor(N_Diabetes_WHO.ph1 == 2,
                     levels=c(FALSE,TRUE),
                     labels=c("NoT2D","T2D")),
    frag_q4 = factor(frag_break_q4, levels=1:4, labels=paste0("Q",1:4)),
    # other covariates as before…
    n_edu   = factor(n_education_3cat.ph1),
    marst   = factor(marital_status.ph1),
    smoke   = factor(smoking_3cat.ph1),
    alc     = factor(N_alcohol_cat.ph1),
    CVD     = factor(N_CVD.ph1),
    mvpat   = factor(mvpatile)
  )

## 5 · Subset to complete cases on all model vars ---------------------
vars_needed <- c(
  "LD_PHQ9depr_Surv_time_standard_default","LD_PHQ9depr_event",
  "mean_inbed_min_night_t.ph1",  # continuous exposure
  "Sex", "T2diab",               # modifiers
  "frag_q4","Age.ph1","n_edu","marst",
  "smoke","alc","bmi.ph1","CVD","dhd_sum_min_alc.ph1","mvpat"
)
dat2 <- dat[complete.cases(dat[, vars_needed]), ]

## 6 · rms setup -------------------------------------------------------
dd <- datadist(dat2)
options(datadist="dd")

## 7 · Base spline model (no interactions) ----------------------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +  # 4-knot spline
    frag_q4 + Age.ph1 + Sex + T2diab +
    n_edu + marst + smoke + alc +
    bmi.ph1 + CVD + dhd_sum_min_alc.ph1 + mvpat,
  data=dat2, x=TRUE, y=TRUE, surv=TRUE
)

## 8 · Interaction with **Sex** ----------------------------------------
# 8a · Add spline × Sex
fit_int_sex <- update(fit_base,
                      . ~ rcs(mean_inbed_min_night_t.ph1, 4) * Sex +
                        frag_q4 + Age.ph1 + T2diab +
                        n_edu + marst + smoke + alc +
                        bmi.ph1 + CVD + dhd_sum_min_alc.ph1 + mvpat
)

# 8b · LRT
cat("\n--- LR test: spline × Sex interaction ---\n")
print(lrtest(fit_base, fit_int_sex))

## 9 · Interaction with **T2D** ----------------------------------------
# 9a · Add spline × T2diab
fit_int_diab <- update(fit_base,
                       . ~ rcs(mean_inbed_min_night_t.ph1, 4) * T2diab +
                         frag_q4 + Age.ph1 + Sex +
                         n_edu + marst + smoke + alc +
                         bmi.ph1 + CVD + dhd_sum_min_alc.ph1 + mvpat
)

# 9b · LRT
cat("\n--- LR test: spline × T2diab interaction ---\n")
print(lrtest(fit_base, fit_int_diab))

## 10 · (Optional) Plotting predicted curves --------------------------
# to visualize separately by Sex or T2diab, e.g.:
# library(ggplot2); library(ggpubr)
# pred_sex  <- Predict(fit_int_sex,  mean_inbed_min_night_t.ph1, Sex, fun=exp)
# autoplot(pred_sex) + labs(x="Minutes in bed", y="Hazard ratio", color="Sex")
#
# pred_diab <- Predict(fit_int_diab, mean_inbed_min_night_t.ph1, T2diab, fun=exp)
# autoplot(pred_diab) + labs(x="Minutes in bed", y="Hazard ratio", color="T2D status")

## 11 · Reset -----------------------------------------------------------
options(datadist=NULL)
# --------------------------------------------------------------------
