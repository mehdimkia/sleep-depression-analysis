# --------------------------------------------------------------------
# sleep_cutoff_sensitivity_debugged.R
#   Re‐estimate the fully‐adjusted Cox model under alternative
#   sleep‐duration cut‐offs, using survival::coxph() for tidy()
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>% |>
library(purrr)      # pmap_dfr()
library(broom)      # tidy()
library(survival)   # Surv(), coxph()

## 2 · Import & clean -------------------------------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>%
  mutate(across(everything(),
                ~ replace(., . %in% c(-777, -888, -999), NA_real_)
  )) %>%
  mutate(
    # continuous sleep in hours
    sleep_hours = mean_inbed_min_night_t.ph1 / 60,
    # factorize all covariates exactly as your main model
    Sex                   = factor(Sex, labels = c("Men","Women")),
    N_Diabetes_WHO.ph1    = factor(N_Diabetes_WHO.ph1),
    n_education_3cat.ph1  = factor(n_education_3cat.ph1),
    marital_status.ph1    = factor(marital_status.ph1),
    smoking_3cat.ph1      = factor(smoking_3cat.ph1),
    N_alcohol_cat.ph1     = factor(N_alcohol_cat.ph1),
    N_CVD.ph1             = factor(N_CVD.ph1),
    mvpatile              = factor(mvpatile),
    frag_break_q4         = factor(frag_break_q4,
                                   levels = 1:4,
                                   labels = paste0("Q",1:4))
  )

## 3 · Helper: build sleep_cat, fit coxph, extract sleep HRs ---------
run_cutoff_model <- function(short_thr, long_thr) {
  d2 <- dat %>%
    mutate(
      sleep_cat = factor(
        case_when(
          sleep_hours <  short_thr ~ sprintf("<%.1f h", short_thr),
          sleep_hours >= long_thr  ~ sprintf("≥%.1f h", long_thr),
          TRUE                     ~ sprintf("%.1f–%.1f h", short_thr, long_thr)
        ),
        levels = c(
          sprintf("%.1f–%.1f h", short_thr, long_thr),
          sprintf("<%.1f h", short_thr),
          sprintf("≥%.1f h", long_thr)
        )
      )
    ) %>%
    # drop any missing in outcome, sleep_cat, or covariates
    filter(!is.na(sleep_cat)) %>%
    drop_na(
      LD_PHQ9depr_Surv_time_standard_default,
      LD_PHQ9depr_event,
      frag_break_q4, Age.ph1, Sex, N_Diabetes_WHO.ph1,
      n_education_3cat.ph1, marital_status.ph1, smoking_3cat.ph1,
      N_alcohol_cat.ph1, bmi.ph1, N_CVD.ph1, dhd_sum_min_alc.ph1,
      mvpatile
    )
  
  # fit with coxph()
  fit_r <- coxph(
    Surv(LD_PHQ9depr_Surv_time_standard_default,
         LD_PHQ9depr_event) ~
      sleep_cat + frag_break_q4 +
      Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
      n_education_3cat.ph1 + marital_status.ph1 +
      smoking_3cat.ph1 + N_alcohol_cat.ph1 +
      bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
      mvpatile,
    data = d2,
    ties = "efron"
  )
  
  # extract just the sleep_cat rows
  broom::tidy(fit_r, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(grepl("^sleep_cat", term)) %>%
    mutate(
      short_thr = short_thr,
      long_thr  = long_thr,
      N         = nrow(d2),
      events    = sum(d2$LD_PHQ9depr_event)
    ) %>%
    select(short_thr, long_thr, N, events,
           term, HR = estimate, CI_low = conf.low,
           CI_high = conf.high, p.value)
}

## 4 · Build the grid & run --------------------------------------------
cutoff_grid <- tibble::tribble(
  ~short_thr, ~long_thr,
  6.0,         9.5,
  6.0,        10.0,
  6.5,         9.5,
  6.5,        10.0
)

sensitivity_HRs <- pmap_dfr(cutoff_grid, run_cutoff_model)

## 5 · Display ---------------------------------------------------------
cat("\n--- Sleep-duration sensitivity analyses ---\n")
print(sensitivity_HRs, digits = 3)
