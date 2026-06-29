# Univariable Cox proportional hazards — single covariate
#
# Inputs:   data frame with time/event/group; reference level set explicitly.
# Outputs:  HR (95% CI), p-value, PH check, tidy table; printed and saved.
# Assumes:  event coding verified (1 = event, 0 = censored); reference level set.

# ---- packages ---------------------------------------------------------------
library(survival)
library(broom)
library(gtsummary)
library(dplyr)
library(forcats)
library(here)

# ---- data -------------------------------------------------------------------
# >>> EDIT: load data + set reference level
df_analysis <- readr::read_csv(here::here("data", "trial.csv")) |>
  dplyr::mutate(arm = forcats::fct_relevel(arm, "Control"))

# ---- fit --------------------------------------------------------------------
# >>> EDIT: replace `arm` with the covariate of interest
cox_uni <- survival::coxph(
  survival::Surv(os_months, os_event) ~ arm,
  data = df_analysis
)
summary(cox_uni)

# ---- tidy output (HR with 95% CI) -------------------------------------------
hr_tab <- broom::tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE) |>
  dplyr::transmute(
    Covariate = term,
    HR        = round(estimate, 2),
    CI_low    = round(conf.low, 2),
    CI_high   = round(conf.high, 2),
    p         = signif(p.value, 3)
  )
print(hr_tab)

# ---- proportional hazards check ---------------------------------------------
zph <- survival::cox.zph(cox_uni)
print(zph)

# Schoenfeld residual plot — flag any non-flat smooth
# pdf(here::here("output", "fig_phcheck_uni.pdf"), width = 5, height = 4)
# plot(zph); abline(h = 0, lty = 2)
# dev.off()

if (any(zph$table[, "p"] < 0.05, na.rm = TRUE)) {
  warning(
    "PH assumption likely violated. ",
    "Consider RMST (scripts/rmst.R), stratification, or time-varying coefficient."
  )
}

# ---- manuscript table -------------------------------------------------------
tab <- cox_uni |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3),
    label = list(arm ~ "Treatment arm")     # >>> EDIT label
  ) |>
  gtsummary::bold_labels() |>
  gtsummary::bold_p(t = 0.05)
tab

# ---- export -----------------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tab |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = file.path(out_dir, "table_cox_univariable.docx"))

cat("Saved: ", file.path(out_dir, "table_cox_univariable.docx"), "\n")
