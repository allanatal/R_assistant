# Median follow-up via reverse Kaplan-Meier
#
# Inputs:   data frame with time-to-event and event (1/0).
# Outputs:  Overall median follow-up + 95% CI; stratified by arm if requested.
# Why:      median(time) and max(time) are NOT follow-up estimators.
#           Reverse-KM swaps event/censoring indicators so that "events" are
#           administrative censoring — what we actually want to estimate.

# ---- packages ---------------------------------------------------------------
library(prodlim)
library(survival)
library(dplyr)
library(here)

# ---- data -------------------------------------------------------------------
df_analysis <- readr::read_csv(here::here("data", "trial.csv"))

stopifnot(all(df_analysis$os_event %in% c(0, 1)))   # event coding sanity

# ---- overall median follow-up (prodlim) ------------------------------------
fu_fit <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ 1,
  data    = df_analysis,
  reverse = TRUE
)
fu_overall <- quantile(fu_fit, q = 0.5)
print(fu_overall)

# ---- overall via survival package (cross-check) ----------------------------
fu_fit2 <- survival::survfit(
  survival::Surv(os_months, 1 - os_event) ~ 1,
  data = df_analysis
)
print(fu_fit2)   # `median` row in the output = median follow-up + 95% CI

# ---- stratified by arm ------------------------------------------------------
# >>> EDIT: change `arm` to whichever stratification variable you want
fu_by_arm <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ arm,
  data    = df_analysis,
  reverse = TRUE
)
fu_strat <- quantile(fu_by_arm, q = 0.5)
print(fu_strat)

# ---- formatted summary for the methods section -----------------------------
get_median_ci <- function(qres) {
  # qres is the data frame returned by quantile(prodlim_fit, q = 0.5)
  paste0(
    formatC(qres$quantile, format = "f", digits = 1), " ",
    "(95% CI ",
    formatC(qres$lower, format = "f", digits = 1), "-",
    formatC(qres$upper, format = "f", digits = 1), ")"
  )
}

cat(sprintf(
  "Median follow-up (reverse KM): %s\n",
  get_median_ci(fu_overall)
))

# Stratified (one line per arm)
cat("\nMedian follow-up by arm:\n")
print(fu_strat)

# ---- save tidy summary ------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(
  tibble::as_tibble(fu_overall) |> dplyr::mutate(group = "Overall"),
  file.path(out_dir, "median_followup_overall.csv")
)
readr::write_csv(
  tibble::as_tibble(fu_strat),
  file.path(out_dir, "median_followup_by_arm.csv")
)

cat("Saved follow-up summaries to ", out_dir, "\n")

# ---- caveat -----------------------------------------------------------------
# If the median is reported as NA, less than half the cohort would have been
# censored by that time — the median is not estimable. Report the
# 25th/75th percentile or the maximum observed follow-up instead.
