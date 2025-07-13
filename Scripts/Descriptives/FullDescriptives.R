############################################################
##  Table S1  ·  Baseline characteristics by inclusion status
############################################################

## 1 ·  Packages -------------------------------------------------------
libs <- c("tidyverse", "haven", "rstatix", "effectsize")
invisible(lapply(libs, \(p) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p)))
lapply(libs, library, character.only = TRUE)

## 2 ·  Read data ------------------------------------------------------
dat <- read_sav("Inclusion-Exclusion.sav")     # adjust path

## 3 ·  Inclusion cascade ---------------------------------------------
## 3a ·  List every covariate that must be complete
covars <- c(
  "Age.ph1", "Sex", "N_Diabetes_WHO.ph1",
  "marital_status.ph1", "n_education_3cat.ph1",
  "dhd_sum.ph1", "nit_alcoholtot.ph1",        # alcohol (continuous)
  "mean_mvpa_min_wake_t.ph1",                 # mean daily MVPA
  "dhd_sum_min_alc.ph1", "N_alcohol_cat.ph1", # diet & alcohol category
  "smoking_3cat.ph1",
  "N_CVD.ph1", "N_HT.ph1",
  "bmi.ph1", "OSBP.ph1", "ODBP.ph1"           # BP
)

## 3b ·  Row-wise flag: TRUE only if *all* inclusion criteria are met
analysis_flag <- with(dat,
                      !is.na(mean_inbed_min_night_t.ph1)  &          # accel. data complete
                        !is.na(mean_inbed_break_night_t.ph1) &
                        !is.na(mean_up_trans_night_t.ph1)    &
                        !is.na(PHQ9_score.ph1)               &          # baseline PHQ-9 present
                        PHQ9_score.ph1 < 10                  &          # no baseline depression
                        LD_PHQ9depr_available == 1                       # follow-up PHQ-9 present
)

## 3c ·  Add covariate completeness requirement (drops the extra 538)
analysis_flag <- analysis_flag &
  rowSums(is.na(dat[ , covars])) == 0   # no NA in covars

## 3d ·  Mark participants
dat <- dat %>%
  mutate(included = factor(
    if_else(analysis_flag, "Included", "Excluded"),
    levels = c("Excluded", "Included")
  ))

## 4 ·  Variable lists -------------------------------------------------
cat_vars <- c(
  "Sex", "n_education_3cat.ph1", "marital_status.ph1",
  "sleep_dur2", "frag_break_q4", "med_sleep.ph1",
  "smoking_3cat.ph1", "N_alcohol_cat.ph1",
  "N_Diabetes_WHO.ph1", "N_CVD.ph1",
  "med_depression.ph1", "impaired_vibration_sense.ph1",
  "mvpatile"
)

cont_mean   <- c("mean_inbed_min_night_t.ph1", "dhd_sum_min_alc.ph1")   # mean ± SD
cont_median <- c("Age.ph1", "mean_up_trans_night_t.ph1",
                 "mean_mvpa_min_wake_t.ph1", "bmi.ph1")                 # median [IQR]

all_vars <- unique(c(cat_vars, cont_mean, cont_median))

## 5 ·  Clean SPSS sentinel missings & convert labels -----------------
to_na <- function(x) if (is.numeric(x)) replace(x, x <= -700, NA) else x

dat_cln <- dat |>
  mutate(across(everything(), to_na)) |>
  mutate(across(all_of(cat_vars), \(x) haven::as_factor(x, levels = "labels"))) |>
  mutate(across(all_of(c(cont_mean, cont_median)), as.numeric))

## 6 ·  Helpers --------------------------------------------------------
fmt_pct <- function(n, tot) sprintf("%d (%.1f%%)", n, 100*n/tot)

fmt_p <- function(p) {
  if (is.na(p))       return("")
  if (p < 0.001)      return("<.001")
  sprintf("%.3f", round(p, 3))
}

summarise_cat <- function(df, var) {
  if (!var %in% names(df)) return(NULL)
  tmp <- df |> filter(!is.na(.data[[var]]))
  if (nrow(tmp) == 0L) return(NULL)
  
  tmp[[var]] <- droplevels(tmp[[var]])
  
  counts <- tmp |>
    count(included, !!sym(var), name = "n") |>
    group_by(included) |>
    mutate(pct = fmt_pct(n, sum(n))) |>
    ungroup() |>
    select(included, level = !!sym(var), value = pct) |>
    pivot_wider(names_from = included, values_from = value)
  
  totals <- tmp |>
    count(!!sym(var), name = "n") |>
    mutate(Total = fmt_pct(n, sum(n))) |>
    select(level = !!sym(var), Total)
  
  pval <- NA_real_
  tbl  <- table(tmp$included, tmp[[var]])
  if (nrow(tbl) >= 2 && ncol(tbl) >= 2)
    pval <- tryCatch(chisq.test(tbl)$p.value, error = \(e) NA_real_)
  
  counts |>
    left_join(totals, by = "level") |>
    mutate(Variable = var,
           `p value` = fmt_p(pval),
           .before = level)
}

summarise_cont <- function(df, var, style = c("mean", "median")) {
  if (!var %in% names(df)) return(NULL)
  tmp <- df |> filter(!is.na(.data[[var]]))
  if (nrow(tmp) == 0L) return(NULL)
  
  style <- match.arg(style)
  
  grp <- tmp |>
    group_by(included) |>
    summarise(value = if (style == "mean")
      sprintf("%.2f ± %.2f",
              mean(.data[[var]], na.rm = TRUE),
              sd  (.data[[var]], na.rm = TRUE))
      else
        sprintf("%.2f [%.2f–%.2f]",
                median(.data[[var]], na.rm = TRUE),
                quantile(.data[[var]], .25, na.rm = TRUE),
                quantile(.data[[var]], .75, na.rm = TRUE)),
      .groups = "drop") |>
    pivot_wider(names_from = included, values_from = value)
  
  overall <- if (style == "mean")
    sprintf("%.2f ± %.2f",
            mean(tmp[[var]], na.rm = TRUE),
            sd  (tmp[[var]], na.rm = TRUE))
  else
    sprintf("%.2f [%.2f–%.2f]",
            median(tmp[[var]], na.rm = TRUE),
            quantile(tmp[[var]], .25, na.rm = TRUE),
            quantile(tmp[[var]], .75, na.rm = TRUE))
  
  pval <- NA_real_
  if (style == "mean")
    pval <- tryCatch(t.test(tmp[[var]] ~ tmp$included)$p.value,
                     error = \(e) NA_real_)
  else
    pval <- tryCatch(wilcox.test(tmp[[var]] ~ tmp$included)$p.value,
                     error = \(e) NA_real_)
  
  grp |>
    mutate(Total     = overall,
           level     = var,
           Variable  = var,
           `p value` = fmt_p(pval),
           .before = level)
}

## 7 ·  Build table ----------------------------------------------------
table_S1 <- bind_rows(
  map(cat_vars,    summarise_cat,    df = dat_cln),
  map(cont_mean,   summarise_cont,   df = dat_cln, style = "mean"),
  map(cont_median, summarise_cont,   df = dat_cln, style = "median")
)

## 8 ·  Add n-missing column ------------------------------------------
miss_overall <- dat_cln |>
  summarise(across(all_of(all_vars), ~ sum(is.na(.)))) |>
  pivot_longer(everything(),
               names_to  = "Variable",
               values_to = "n missing")

table_S1 <- table_S1 |>
  left_join(miss_overall, by = "Variable") |>
  select(Variable, level,
         Excluded = `Excluded`,
         Included = `Included`,
         Total, `p value`, `n missing`)

## 9 ·  Inspect & export ----------------------------------------------
print(table_S1, n = Inf, width = Inf)
write.csv(table_S1, "Table_S1_values.csv", row.names = FALSE)
