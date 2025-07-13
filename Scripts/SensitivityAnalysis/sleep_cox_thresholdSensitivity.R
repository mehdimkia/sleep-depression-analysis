# --------------------------------------------------------------------
# sleep_cox_thresholdSensitivity_debugged2.R
#   Fully adjusted Cox model under alternative cutoffs
#   (automatically matches the actual cph() term names)
# --------------------------------------------------------------------

## 1 · Libraries -------------------------------------------------------
library(haven)      # read_sav(), zap_labels()
library(dplyr)      # %>%, mutate(), case_when()
library(tidyr)      # drop_na()
library(rms)        # cph()
library(survival)   # Surv()
library(purrr)      # map_dfr()
library(tibble)     # tibble()

## 2 · Load & sentinel → NA -------------------------------------------
dat <- read_sav("Week7.sav") |> zap_labels()
dat <- dat %>% mutate(across(everything(),
                             ~ replace(., . %in% c(-777, -888, -999), NA_real_)))

## 3 · Helper: fit + extract “long-sleep” HR robustly ------------------
make_model <- function(short_cut, long_cut) {
  df <- dat %>%
    mutate(
      bed_h = mean_inbed_min_night_t.ph1 / 60,
      sleep_cat = case_when(
        bed_h <  short_cut ~ paste0("<", short_cut, " h"),
        bed_h >= long_cut  ~ paste0("≥", long_cut, " h"),
        TRUE               ~ paste0(short_cut, "–", long_cut, " h")
      ),
      sleep_cat = factor(
        sleep_cat,
        levels = c(
          paste0(short_cut, "–", long_cut, " h"),
          paste0("<", short_cut, " h"),
          paste0("≥", long_cut, " h")
        )
      )
    ) %>%
    drop_na(
      sleep_cat, frag_break_q4, Age.ph1, Sex, N_Diabetes_WHO.ph1,
      n_education_3cat.ph1, marital_status.ph1, smoking_3cat.ph1,
      N_alcohol_cat.ph1, bmi.ph1, N_CVD.ph1, dhd_sum_min_alc.ph1,
      mvpatile, LD_PHQ9depr_Surv_time_standard_default,
      LD_PHQ9depr_event
    )
  
  # fit the model
  fit <- cph(
    Surv(LD_PHQ9depr_Surv_time_standard_default,
         LD_PHQ9depr_event) ~ sleep_cat + frag_break_q4 +
      Age.ph1 + Sex + N_Diabetes_WHO.ph1 +
      n_education_3cat.ph1 + marital_status.ph1 +
      smoking_3cat.ph1 + N_alcohol_cat.ph1 +
      bmi.ph1 + N_CVD.ph1 + dhd_sum_min_alc.ph1 +
      mvpatile,
    data = df, x=TRUE, y=TRUE, surv=TRUE, parms=FALSE
  )
  
  # pull coef names and find the one that ends with our long_cut
  cs   <- coef(fit)
  sc_terms <- grep("^sleep_cat", names(cs), value=TRUE)
  # look for the numeric long_cut (e.g. "9.5") in those names
  match_term <- sc_terms[grep(as.character(long_cut), sc_terms)]
  
  if (length(match_term) != 1) {
    # no exact match → return an NA‐row
    return(tibble(
      short_cut = short_cut,
      long_cut  = long_cut,
      N         = nrow(df),
      events    = sum(df$LD_PHQ9depr_event),
      HR        = NA_real_,
      CI_low    = NA_real_,
      CI_high   = NA_real_,
      p_value   = NA_real_
    ))
  }
  
  # compute HR, CI, p
  vc   <- vcov(fit)
  β    <- cs[match_term]
  s_e  <- sqrt(vc[match_term, match_term])
  HR   <- exp(β)
  low  <- exp(β - 1.96 * s_e)
  high <- exp(β + 1.96 * s_e)
  pval <- 2*(1 - pnorm(abs(β / s_e)))
  
  tibble(
    short_cut = short_cut,
    long_cut  = long_cut,
    N         = nrow(df),
    events    = sum(df$LD_PHQ9depr_event),
    HR        = HR,
    CI_low    = low,
    CI_high   = high,
    p_value   = pval
  )
}

## 4 · Apply over all four specs --------------------------------------
shorts <- c(6.0, 6.5)
longs  <- c(9.5, 10.0)

results <- map_dfr(shorts, function(s) {
  map_dfr(longs, function(l) {
    make_model(s, l)
  })
})

## 5 · View results ---------------------------------------------------
print(results, digits = 3)
