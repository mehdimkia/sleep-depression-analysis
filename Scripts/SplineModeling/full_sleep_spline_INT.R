# --------------------------------------------------------------
# sleep_interaction_models_with_export.R  
#   –  Cox spline + interaction tests + export to HTML/Word
# --------------------------------------------------------------

## 1 ·  Libraries ------------------------------------------------------
library(haven)        # read_sav(), zap_labels()
library(dplyr)        # data wrangling
library(rms)          # cph(), rcs()
library(survival)     # Surv()
library(lmtest)       # lrtest()
library(broom)        # tidy()
library(stargazer)    # export HTML tables
library(flextable)    # pretty tables
library(officer)      # write Word docs

## 2 ·  Read and prepare data -----------------------------------------
dat <- read_sav("Week5.sav") |>
  zap_labels()

# Factorise categorical covariates
dat <- dat |>
  mutate(
    Sex                   = factor(Sex, levels = c(1, 2), labels = c("Men", "Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1,
                                   levels = c(0, 1, 2),
                                   labels = c("Normo", "PreDM", "T2DM")),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    mvpatile              = factor(mvpatile),
    frag_break_q4         = factor(frag_break_q4)
  )

## 2b ·  rms setup -----------------------------------------------------
dd <- datadist(dat)
options(datadist = "dd")

## 3 ·  Base Cox model with 4-knot restricted cubic spline -------------
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

cat("\n=== Base model (fit_full):\n")
s_full <- summary(fit_full)
print(s_full, digits = 3)

## 4 ·  Interaction with SEX ------------------------------------------
fit_sexInt <- update(fit_full,
                     . ~ rcs(mean_inbed_min_night_t.ph1, 4) * Sex +
                       frag_break_q4 +
                       Age.ph1 + N_Diabetes_WHO.ph1 +
                       n_education_3cat.ph1 + marital_status.ph1 +
                       smoking_3cat.ph1 + N_alcohol_cat.ph1 +
                       bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
                       mvpatile
)

cat("\n=== Interaction model: sleep spline × Sex (fit_sexInt):\n")
s_sexInt <- summary(fit_sexInt)
print(s_sexInt, digits = 3)

lr_sex <- lrtest(fit_full, fit_sexInt)
cat("\nLR test for interaction with sex:\n")
print(lr_sex)

## 5 ·  Interaction with DIABETES STATUS ------------------------------
fit_dmInt <- update(fit_full,
                    . ~ rcs(mean_inbed_min_night_t.ph1, 4) * N_Diabetes_WHO.ph1 +
                      frag_break_q4 +
                      Age.ph1 + Sex +
                      n_education_3cat.ph1 + marital_status.ph1 +
                      smoking_3cat.ph1 + N_alcohol_cat.ph1 +
                      bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
                      mvpatile
)

cat("\n=== Interaction model: sleep spline × Diabetes (fit_dmInt):\n")
s_dmInt <- summary(fit_dmInt)
print(s_dmInt, digits = 3)

lr_dm <- lrtest(fit_full, fit_dmInt)
cat("\nLR test for interaction with diabetes status:\n")
print(lr_dm)

# --------------------------------------------------------------
# Step 6 ·  Refit with coxph() so broom::tidy() works
# --------------------------------------------------------------

# 6a) Base model with coxph (using the same spline & covariates)
fit_full_coxph <- coxph(
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event) ~
    rcs(mean_inbed_min_night_t.ph1, 4) +
    frag_break_q4 +
    Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
    n_education_3cat.ph1 + marital_status.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
    mvpatile,
  data = dat,
  ties = "efron"
)

# 6b) Interaction models with coxph
fit_sexInt_coxph <- update(fit_full_coxph,
                           . ~ rcs(mean_inbed_min_night_t.ph1, 4)*Sex +
                             frag_break_q4 + Age.ph1 + N_Diabetes_WHO.ph1 +
                             n_education_3cat.ph1 + marital_status.ph1 +
                             smoking_3cat.ph1 + N_alcohol_cat.ph1 +
                             bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
                             mvpatile
)

fit_dmInt_coxph  <- update(fit_full_coxph,
                           . ~ rcs(mean_inbed_min_night_t.ph1, 4)*N_Diabetes_WHO.ph1 +
                             frag_break_q4 + Age.ph1 + Sex +
                             n_education_3cat.ph1 + marital_status.ph1 +
                             smoking_3cat.ph1 + N_alcohol_cat.ph1 +
                             bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
                             mvpatile
)

# 6c) Tidy into data‐frames
library(broom)
coefs_full   <- tidy(fit_full_coxph,   exponentiate = TRUE, conf.int = TRUE)
coefs_sexInt <- tidy(fit_sexInt_coxph, exponentiate = TRUE, conf.int = TRUE)
coefs_dmInt  <- tidy(fit_dmInt_coxph,  exponentiate = TRUE, conf.int = TRUE)

# --------------------------------------------------------------
# Step 7 ·  Export with stargazer & flextable/officer
# --------------------------------------------------------------

# 7a) HTML table of the three models
library(stargazer)
stargazer(
  fit_full_coxph, fit_sexInt_coxph, fit_dmInt_coxph,
  type          = "html",
  out           = "cox_models.html",
  column.labels = c("Base", "×Sex", "×Diabetes"),
  ci            = TRUE,
  single.row    = TRUE,
  digits        = 3
)

# 7b) Word doc of the base‐model coefficients
library(flextable)
library(officer)
ft <- flextable(coefs_full)
doc <- read_docx() |>
  body_add_par("Full CoxPH Model Coefficients", style = "heading 1") |>
  body_add_flextable(ft)
print(doc, target = "coefs_full.docx")

# --------------------------------------------------------------
# Clean up
# --------------------------------------------------------------
options(datadist = NULL)
# --------------------------------------------------------------

# --------------------------------------------------------------
# End of script
# --------------------------------------------------------------
