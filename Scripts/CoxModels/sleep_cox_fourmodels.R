# --------------------------------------------------------------------
# sleep_cox_fourmodels.R  –  Recreate Table of HRs for Models 1–4
# --------------------------------------------------------------------
#  * Week7.sav assumed in working directory
#  * Uses survival::coxph (simpler extraction than rms::cph for tables)
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # data wrangling pipes
library(survival)   # coxph(), Surv()
library(broom)      # tidy() for model output

## 2 · Import & factorise ---------------------------------------------
dat <- read_sav("Week7.sav") %>% zap_labels()

dat <- dat %>%
  mutate(
    # Key exposures ────────────────────────────────────────────────
    sleep_dur2 = factor(
      sleep_dur2,
      levels = c(2, 1, 3),                # 2 = “Normal (7–9 h)” reference
      labels = c("Normal", "Short", "Long")
    ),
    frag_break_q4 = factor(
      frag_break_q4,
      levels = 1:4,
      labels = paste0("Q", 1:4)
    ),
    # Covariates ───────────────────────────────────────────────────
    Sex                  = factor(Sex),                     # 1=Men, 2=Women
    N_Diabetes_WHO.ph1   = factor(N_Diabetes_WHO.ph1),      # 0=No, 1=Type-2
    n_education_3cat.ph1 = factor(n_education_3cat.ph1),    # 1=Low,2=Mid,3=High
    marital_status.ph1   = factor(marital_status.ph1),      # 1=Single,2=Partner,…
    smoking_3cat.ph1     = factor(smoking_3cat.ph1),        # 0=Never,1=Former,2=Curr
    N_alcohol_cat.ph1    = factor(N_alcohol_cat.ph1),       # 0=None,1=Mod,2=High
    N_CVD.ph1            = factor(N_CVD.ph1),               # 0=No CVD, 1=Yes
    mvpatile             = factor(mvpatile)                 # 1=Low,2=Med,3=High
  )

## 3 · Survival object -------------------------------------------------
surv_obj <- with(
  dat,
  Surv(LD_PHQ9depr_Surv_time_standard_default, LD_PHQ9depr_event)
)

## 4 · Model formulas --------------------------------------------------
form_m1 <- surv_obj ~ sleep_dur2 + frag_break_q4

form_m2 <- update(
  form_m1,
  . ~ . + Age.ph1 + Sex + N_Diabetes_WHO.ph1
)

form_m3 <- update(
  form_m2,
  . ~ . + marital_status.ph1 + n_education_3cat.ph1
)

form_m4 <- update(
  form_m3,
  . ~ . + bmi.ph1 + N_CVD.ph1 +
    smoking_3cat.ph1 + N_alcohol_cat.ph1 +
    dhd_sum_min_alc.ph1 + mvpatile
)

models <- list(
  "Model 1 (Crude)"                       = form_m1,
  "Model 2 (+ Age, Sex, T2D)"             = form_m2,
  "Model 3 (+ Marital, Education)"        = form_m3,
  "Model 4 (+ BMI, CVD, lifestyle)"       = form_m4
)

## 5 · Fit models & collect results -----------------------------------
extract_rows <- c(
  "sleep_dur2Short",
  "sleep_dur2Long",
  "frag_break_q4Q2",
  "frag_break_q4Q3",
  "frag_break_q4Q4"
)

tbl_list <- lapply(names(models), function(mn) {
  fit <- coxph(models[[mn]], data = dat)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    filter(term %in% extract_rows) |>
    mutate(Model = mn,
           term  = recode(term,
                          sleep_dur2Short   = "Short (<7 h / day)",
                          sleep_dur2Long    = "Long (>9 h / day)",
                          frag_break_q4Q2   = "Q2",
                          frag_break_q4Q3   = "Q3",
                          frag_break_q4Q4   = "Q4 (highest)")) |>
    select(Model, `Sleep / Fragmentation level` = term,
           HR = estimate, `Lower 95%` = conf.low, `Upper 95%` = conf.high,
           `p value` = p.value)
})

results_tbl <- bind_rows(tbl_list) |>
  mutate(across(c(HR, `Lower 95%`, `Upper 95%`),
                ~ sprintf("%.3f", .x)),
         CI = paste0(HR, " (", `Lower 95%`, "–", `Upper 95%`, ")")) |>
  select(Model, `Sleep / Fragmentation level`, CI, `p value`)

## 6 · Print nicely ----------------------------------------------------
print(results_tbl, n = Inf, row.names = FALSE)

## 7 · (Optional) save to CSV -----------------------------------------
# write.csv(results_tbl, "Table_sleep_models.csv", row.names = FALSE)

# --------------------------------------------------------------------
# End of script
# --------------------------------------------------------------------
