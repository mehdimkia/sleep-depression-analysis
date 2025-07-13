# --------------------------------------------------------------------
# sleep_cox_anyVulnInteraction.R
#   Fully adjusted Cox model + sleep_dur2 × any‐vulnerability interaction
#   (any MINI history, antidepressant use, or neuropathy)
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist(), Predict()
library(survival)   # Surv(), coxph(), cox.zph()
library(lmtest)     # lrtest()
library(broom)      # tidy()

## 2 · Import & clean sentinel → NA -----------------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>% 
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)))

## 3 · Factorise & composite modifier ---------------------------------
dat <- dat %>%
  mutate(
    Sex                         = factor(Sex,                   labels=c("Men","Women")),
    N_Diabetes_WHO.ph1          = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1        = factor(n_education_3cat.ph1),
    marital_status.ph1          = factor(marital_status.ph1),
    smoking_3cat.ph1            = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1           = factor(N_alcohol_cat.ph1),
    N_CVD.ph1                   = factor(N_CVD.ph1),
    mvpatile                    = factor(mvpatile),
    MINIlifedepr.ph1            = factor(MINIlifedepr.ph1,      labels=c("No","Yes")),
    med_depression.ph1          = factor(med_depression.ph1,    labels=c("No","Yes")),
    impaired_vibration_sense.ph1= factor(impaired_vibration_sense.ph1, labels=c("No","Yes")),
    any_vuln                    = factor(
      (MINIlifedepr.ph1=="Yes") |
        (med_depression.ph1=="Yes") |
        (impaired_vibration_sense.ph1=="Yes"),
      labels=c("No","Yes")
    ),
    sleep_dur2                  = factor(sleep_dur2,
                                         levels=c(2,1,3),
                                         labels=c("7–9 h","<7 h","≥9 h")),
    frag_break_q4               = factor(frag_break_q4,
                                         levels=1:4,
                                         labels=paste0("Q",1:4))
  )

cat("\nTotal participants:", nrow(dat),
    "\nTotal events:", sum(dat$LD_PHQ9depr_event, na.rm=TRUE), "\n")

## 4 · Subset complete cases ------------------------------------------
vars_needed <- c(
  "LD_PHQ9depr_Surv_time_standard_default","LD_PHQ9depr_event",
  "sleep_dur2","any_vuln",
  "frag_break_q4","Age.ph1","Sex","N_Diabetes_WHO.ph1",
  "n_education_3cat.ph1","marital_status.ph1","smoking_3cat.ph1",
  "N_alcohol_cat.ph1","bmi.ph1","N_CVD.ph1","dhd_sum_min_alc.ph1","mvpatile"
)
dat2 <- dat[ complete.cases(dat[, vars_needed]), ]
cat("After complete-case subset:", nrow(dat2),
    "\nEvents:", sum(dat2$LD_PHQ9depr_event), "\n")

## 5 · rms setup -------------------------------------------------------
dd <- datadist(dat2); options(datadist="dd")

## 6 · Model without interaction --------------------------------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 + any_vuln +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, x=TRUE, y=TRUE, surv=TRUE, parms=FALSE
)

## 7 · Model with sleep × any_vuln interaction ------------------------
fit_int <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * any_vuln +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, x=TRUE, y=TRUE, surv=TRUE, parms=FALSE
)

## 8 · LR test omnibus interaction ------------------------------------
cat("\n--- LR test: sleep_dur2 × any_vulnerability ---\n")
print(lrtest(fit_base, fit_int))

## 9 · Interaction coefficients & contrasts ---------------------------
cat("\n--- Interaction model coefficients (cph) ---\n")
print(summary(fit_int), digits=3)

# Cross‐product HRs via survival::coxph + broom::tidy
fit_int_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * any_vuln +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, ties="efron"
)

cp <- tidy(fit_int_coxph, exponentiate=TRUE, conf.int=TRUE) %>%
  filter(grepl("sleep_dur2.*:any_vulnYes", term))
cat("\n--- Stratum‐specific HRs for sleep by any_vuln ---\n")
print(cp, digits=3)

## 10 · Plot predicted HR curves --------------------------------------
plot(
  Predict(fit_int, sleep_dur2, any_vuln, fun=exp),
  xlab="Sleep duration", ylab="Hazard ratio",
  main="Sleep × Any-Vulnerability interaction"
)

## 11 · PH diagnostics -----------------------------------------------
cat("\n--- PH test (interaction model) ---\n")
print(cox.zph(fit_int))

## 12 · Clean up -------------------------------------------------------
options(datadist=NULL)
# --------------------------------------------------------------------
