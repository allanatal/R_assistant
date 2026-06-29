# Multivariable Cox proportional hazards — clinically meaningful adjustment
#
# Inputs:   data frame with time/event/treatment + adjustment covariates;
#           covariates chosen by clinical reasoning, NOT stepwise selection.
# Outputs:  Adjusted HR (95% CI) for each covariate, PH diagnostics, forest plot,
#           univariable + multivariable side-by-side gtsummary table.
# Assumes:  event coding verified; factor levels set; >= ~10 events per covariate.

# ---- packages ---------------------------------------------------------------
library(survival)
library(gtsummary)
library(broom.helpers)
library(forestmodel)
library(dplyr)
library(forcats)
library(here)

# ---- data + factor setup ----------------------------------------------------
df_analysis <- readr::read_csv(here::here("data", "trial.csv")) |>
  dplyr::mutate(
    arm   = forcats::fct_relevel(arm,   "Control"),         # >>> EDIT
    ecog  = forcats::fct_relevel(factor(ecog), "0"),        # >>> EDIT
    stage = forcats::fct_relevel(factor(stage), "I"),       # >>> EDIT
    sex   = forcats::fct_relevel(sex,   "Female")           # >>> EDIT
  )

# Event-per-variable sanity check
n_events <- sum(df_analysis$os_event == 1, na.rm = TRUE)
n_covs   <- 5   # >>> EDIT: number of covariates in the model
if (n_events < 10 * n_covs) {
  warning(sprintf(
    "Only %d events for %d covariates (< 10 events/covariate). HR estimates may be unstable.",
    n_events, n_covs
  ))
}

# ---- fit --------------------------------------------------------------------
# >>> EDIT: pick covariates based on clinical reasoning (confounders, effect modifiers)
cox_mv <- survival::coxph(
  survival::Surv(os_months, os_event) ~ arm + age + sex + ecog + stage,
  data = df_analysis
)
summary(cox_mv)

# ---- PH diagnostics ---------------------------------------------------------
zph <- survival::cox.zph(cox_mv)
print(zph)

violations <- rownames(zph$table)[zph$table[, "p"] < 0.05]
if (length(violations) > 0) {
  warning(
    "PH violated for: ", paste(violations, collapse = ", "), ". ",
    "Consider stratifying these variables (`+ strata(var)`), ",
    "fitting a time-varying effect (`tt()`), or reporting RMST."
  )
}

# Schoenfeld plots (one panel per covariate)
# pdf(here::here("output", "fig_phcheck_mv.pdf"), width = 7, height = 5)
# par(mfrow = c(2, 3)); plot(zph)
# dev.off()

# ---- manuscript table: univariable + multivariable -------------------------
uni <- df_analysis |>
  gtsummary::tbl_uvregression(
    method  = survival::coxph,
    y       = survival::Surv(os_months, os_event),
    include = c(arm, age, sex, ecog, stage),
    exponentiate = TRUE,
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3)
  )

mv <- cox_mv |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3),
    label = list(                                            # >>> EDIT labels
      arm   ~ "Treatment arm",
      age   ~ "Age (per year)",
      sex   ~ "Sex",
      ecog  ~ "ECOG performance status",
      stage ~ "Disease stage"
    )
  )

cox_table <- gtsummary::tbl_merge(
  list(uni, mv),
  tab_spanner = c("**Univariable**", "**Multivariable**")
) |>
  gtsummary::bold_labels()

cox_table

# ---- forest plot ------------------------------------------------------------
# forestmodel produces a one-call manuscript-ready forest
fm <- forestmodel::forest_model(cox_mv)
print(fm)

# ---- export -----------------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Table → Word
cox_table |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = file.path(out_dir, "table_cox_multivariable.docx"))

# Forest → PDF
ggplot2::ggsave(
  filename = file.path(out_dir, "fig_forest_mv.pdf"),
  plot     = fm,
  width    = 7, height = 5, units = "in",
  device   = grDevices::cairo_pdf
)

cat("Saved table_cox_multivariable.docx and fig_forest_mv.pdf to ", out_dir, "\n")
