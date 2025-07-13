# --------------------------------------------------------------------
# sleep_cox_MINIinteraction_diag.R
#   Fully adjusted Cox model + sleep_dur2 × MINIlifedepr.ph1 interaction
#   + diagnostics: cell counts, event counts, Firth-penalised model
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(rms)        # cph(), datadist(), contrast(), Predict()
library(survival)   # Surv(), coxph(), cox.zph()
library(lmtest)     # lrtest()
library(coxphf)     # Firth-penalised Cox (diagnostic block C)

## 2 · Import & clean sentinel values → NA -----------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>%
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)))

## 3 · Factorise categorical variables --------------------------------
dat <- dat %>%
  mutate(
    Sex                   = factor(Sex, labels = c("Men","Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    mvpatile              = factor(mvpatile),
    MINIlifedepr.ph1      = factor(MINIlifedepr.ph1, labels = c("No","Yes")),
    sleep_dur2            = factor(
      sleep_dur2,
      levels = c(2,1,3),                  # 2 = 7–9 h reference
      labels = c("7–9 h","<7 h","≥9 h")
    ),
    frag_break_q4         = factor(
      frag_break_q4,
      levels = 1:4,
      labels = paste0("Q",1:4)
    )
  )

cat("\nParticipants before subsetting:", nrow(dat),
    "\nIncident events:", sum(dat$LD_PHQ9depr_event, na.rm = TRUE), "\n")

## 4 · Subset to complete cases on all variables used -----------------
vars_needed <- c(
  "LD_PHQ9depr_Surv_time_standard_default", "LD_PHQ9depr_event",
  "sleep_dur2", "MINIlifedepr.ph1",
  "frag_break_q4", "Age.ph1", "Sex", "N_Diabetes_WHO.ph1",
  "n_education_3cat.ph1", "marital_status.ph1", "smoking_3cat.ph1",
  "N_alcohol_cat.ph1", "bmi.ph1", "N_CVD.ph1",
  "dhd_sum_min_alc.ph1", "mvpatile"
)
dat2 <- dat[complete.cases(dat[, vars_needed]), ]

cat("Participants after subsetting:", nrow(dat2),
    "\nIncident events:", sum(dat2$LD_PHQ9depr_event), "\n")

## ---- Block A · Raw counts by MINI × sleep --------------------------
cat("\n--- Raw participant counts (MINI × sleep) ---\n")
print(
  dat2 %>% count(MINIlifedepr.ph1, sleep_dur2),
  n = Inf
)

## ---- Block B · Events & participants per stratum -------------------
cat("\n--- Events & participants per stratum ---\n")
print(
  dat2 %>% 
    group_by(MINIlifedepr.ph1, sleep_dur2) %>% 
    summarise(events = sum(LD_PHQ9depr_event),
              participants = n(),
              .groups = "drop"),
  n = Inf
)

## 5 · rms setup on cleaned data --------------------------------------
dd <- datadist(dat2)
options(datadist = "dd")

## 6 · Model without interaction (Model-4 + MINI) ---------------------
fit_base <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 + MINIlifedepr.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data  = dat2,
  x     = TRUE, y = TRUE, surv = TRUE
)

## 7 · Model with sleep_dur2 × MINI interaction -----------------------
fit_int <- cph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * MINIlifedepr.ph1 +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data  = dat2,
  x     = TRUE, y = TRUE, surv = TRUE
)

## 8 · Likelihood-ratio test for interaction --------------------------
cat("\n--- LR test: interaction vs. no interaction ---\n")
print(lrtest(fit_base, fit_int))

## 8b · Coefficient table ---------------------------------------------
cat("\n--- Coefficients (interaction model) ---\n")
print(summary(fit_int), digits = 3)

## 8c · Stratum-specific contrasts (unchanged code) -------------------
# ... (keep your existing contrast-building section here) ...

## 8d · Plot predicted survival by sleep × MINI -----------------------
plot(
  Predict(fit_int, sleep_dur2, MINIlifedepr.ph1, fun = exp),
  xlab = "Sleep duration",
  ylab = "Hazard ratio",
  main = "Sleep duration × MINI interaction"
)

## 9 · Proportional-hazards diagnostics -------------------------------
cat("\n--- PH test (interaction model) ---\n")
print(cox.zph(fit_int))

## ---- Block C · Firth-penalised sensitivity model -------------------
cat("\n--- Firth-penalised sensitivity model (sleep × MINI) ---\n")
fit_firth <- coxphf(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    sleep_dur2 * MINIlifedepr.ph1 + Age.ph1 + Sex + bmi.ph1,
  data = dat2
)
print(summary(fit_firth), digits = 3)

## 10 · Clean up -------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------------
