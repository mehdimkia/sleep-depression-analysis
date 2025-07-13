# --------------------------------------------------------------------
# sleep_cox_sexDiabInteractions.R
#   Fully adjusted Cox models + sleep_dur2 × Sex and sleep_dur2 × T2diab interactions
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% 
library(rms)        # cph(), datadist()
library(survival)   # Surv(), coxph(), cox.zph()
library(lmtest)     # lrtest()
library(broom)      # tidy()

## 2 · Import & clean sentinel values → NA -----------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>%
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)))

## 3 · Factorise categorical variables --------------------------------
dat <- dat %>%
  mutate(
    # Sex: 1 = Men, 2 = Women
    Sex = factor(
      Sex,
      levels = c(1, 2),
      labels = c("Men","Women")
    ),
    
    # Diabetes status: 0 = No diabetes, 1 = Prediabetes, 2 = Type 2 diabetes
    N_Diabetes_WHO.ph1 = factor(
      N_Diabetes_WHO.ph1,
      levels = c(0, 1, 2),
      labels = c("No diabetes","Prediabetes","Type 2 diabetes")
    ),
    
    # Binary T2-diabetes flag (for the interaction test)
    T2diab = factor(
      N_Diabetes_WHO.ph1 == "Type 2 diabetes",
      levels = c(FALSE, TRUE),
      labels = c("No","Yes")
    ),
    
    # Sleep duration categories
    sleep_dur2 = factor(
      sleep_dur2,
      levels = c(2,1,3), 
      labels = c("7–9 h","<7 h","≥9 h")
    ),
    
    # Fragmentation quartiles
    frag_break_q4 = factor(
      frag_break_q4,
      levels = 1:4,
      labels = paste0("Q",1:4)
    ),
    
    # (and the rest as before…)
    n_education_3cat.ph1 = factor(n_education_3cat.ph1),
    marital_status.ph1   = factor(marital_status.ph1),
    smoking_3cat.ph1     = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1    = factor(N_alcohol_cat.ph1),
    N_CVD.ph1            = factor(N_CVD.ph1),
    mvpatile             = factor(mvpatile)
  )

cat("\nParticipants before subsetting:", nrow(dat),
    "\nIncident events:", sum(dat$LD_PHQ9depr_event, na.rm=TRUE), "\n")

## 4 · Subset to complete cases on key variables -----------------------
vars_needed <- c(
  "LD_PHQ9depr_Surv_time_standard_default","LD_PHQ9depr_event",
  "sleep_dur2","Sex","T2diab",
  "frag_break_q4","Age.ph1",
  "n_education_3cat.ph1","marital_status.ph1",
  "smoking_3cat.ph1","N_alcohol_cat.ph1",
  "bmi.ph1","N_CVD.ph1","dhd_sum_min_alc.ph1","mvpatile"
)
dat2 <- dat[ complete.cases(dat[, vars_needed]), ]
cat("Participants after subsetting:", nrow(dat2),
    "\nIncident events:", sum(dat2$LD_PHQ9depr_event), "\n")

## 5 · rms setup on cleaned data --------------------------------------
dd <- datadist(dat2)
options(datadist="dd")

## 6 · Base Cox model (no interaction) -------------------------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,LD_PHQ9depr_event) ~
    sleep_dur2 + Sex + T2diab +
    frag_break_q4 + Age.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data=dat2, x=TRUE, y=TRUE, surv=TRUE, parms=FALSE
)

# ---- Interaction with Sex ----
## 7a · Fit sleep_dur2 × Sex
fit_sex_int <- update(fit_base,
                      . ~ . + sleep_dur2:Sex
)

## 7b · LRT for Sex interaction
cat("\n--- LR test: sleep_dur2 × Sex ---\n")
print(lrtest(fit_base, fit_sex_int))

## 7c · Stratum-specific HRs for sleep by Sex
fit_sex_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,LD_PHQ9depr_event) ~
    sleep_dur2 * Sex +
    frag_break_q4 + Age.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, ties="efron"
)
cat("\n--- Stratum-specific HRs (sleep × Sex) ---\n")
print(
  tidy(fit_sex_coxph, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(grepl("sleep_dur2.*:SexWomen", term)),
  digits=3
)

# ---- Interaction with T2diab ----
## 8a · Fit sleep_dur2 × T2diab
fit_diab_int <- update(fit_base,
                       . ~ . + sleep_dur2:T2diab
)

## 8b · LRT for T2diab interaction
cat("\n--- LR test: sleep_dur2 × T2diab ---\n")
print(lrtest(fit_base, fit_diab_int))

## 8c · Stratum-specific HRs for sleep by T2diab
fit_diab_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default,LD_PHQ9depr_event) ~
    sleep_dur2 * T2diab +
    frag_break_q4 + Age.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat2, ties="efron"
)
cat("\n--- Stratum-specific HRs (sleep × T2diab) ---\n")
print(
  tidy(fit_diab_coxph, exponentiate=TRUE, conf.int=TRUE) %>%
    filter(grepl("sleep_dur2.*:T2diabYes", term)),
  digits=3
)

## 9 · PH diagnostics --------------------------------------------------
cat("\n--- PH test (Sex‐interaction model) ---\n")
print(cox.zph(fit_sex_int))

cat("\n--- PH test (T2diab‐interaction model) ---\n")
print(cox.zph(fit_diab_int))

## 10 · Clean up --------------------------------------------------------
options(datadist=NULL)
# --------------------------------------------------------------------
