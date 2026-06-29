# Reproducible analysis skeleton — clinical oncology project
#
# Use this as the starting point for a new analysis. Adapt section-by-section.
# Project structure assumed:
#   project_root/
#     ├── data/                   ← raw user data
#     ├── R/                      ← reusable functions
#     ├── output/                 ← exported tables and figures
#     ├── reports/                ← Rmd/Quarto manuscripts
#     └── this script             ← end-to-end analysis driver

# ---- packages ---------------------------------------------------------------
library(here)
library(readxl)
library(readr)
library(janitor)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(survival)
library(survminer)
library(gtsummary)
library(broom)
library(flextable)
library(ggplot2)

# ---- session reproducibility ------------------------------------------------
set.seed(42)
sessionInfo() |> capture.output(file = here::here("output", "sessionInfo.txt"))

# ---- 1. LOAD ----------------------------------------------------------------
# >>> EDIT: path + sheet
data_path <- here::here("data", "trial.xlsx")
readxl::excel_sheets(data_path)            # confirm which sheet to load

df_raw <- readxl::read_excel(
  path  = data_path,
  sheet = "patients",                       # >>> EDIT
  na    = c("", "NA", "N/A", "NULL", ".", "?")
)

# ---- 2. INSPECT (see references/07-data-inspection.md) ---------------------
df_clean <- df_raw |> janitor::clean_names()

dplyr::glimpse(df_clean)

# Missingness summary
df_clean |>
  dplyr::summarise(dplyr::across(everything(), ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  dplyr::mutate(pct_missing = round(100 * n_missing / nrow(df_clean), 1)) |>
  dplyr::arrange(dplyr::desc(n_missing)) |>
  print()

# Event coding check
print(table(df_clean$os_event, useNA = "always"))

# Time range check
print(summary(df_clean$os_months))

# ---- 3. PREPARE -------------------------------------------------------------
# >>> EDIT: recode events if needed; set factor levels with clinical reference
df_analysis <- df_clean |>
  dplyr::mutate(
    # Example recodes (uncomment as needed):
    # os_event = dplyr::case_when(
    #   os_event %in% c(1, "Dead", "Yes")  ~ 1L,
    #   os_event %in% c(0, "Alive", "No")  ~ 0L,
    #   TRUE ~ NA_integer_
    # ),
    arm   = forcats::fct_relevel(arm,   "Control"),
    ecog  = forcats::fct_relevel(factor(ecog), "0"),
    stage = forcats::fct_relevel(factor(stage), "I")
  ) |>
  dplyr::filter(!is.na(os_months), !is.na(os_event))

# Cohort flow
cat(sprintf(
  "Raw n = %d -> analysis n = %d (excluded %d for missing time/event)\n",
  nrow(df_raw), nrow(df_analysis), nrow(df_raw) - nrow(df_analysis)
))

# ---- 4. DESCRIBE — Table 1 -------------------------------------------------
# (see scripts/baseline_table1.R for a full version)
tab1 <- df_analysis |>
  dplyr::select(age, sex, ecog, stage, arm) |>
  gtsummary::tbl_summary(
    by = arm,
    missing = "ifany",
    statistic = list(
      gtsummary::all_continuous()  ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    )
  ) |>
  gtsummary::add_overall() |>
  gtsummary::bold_labels()
tab1

# ---- 5. ANALYZE — survival -------------------------------------------------
# Median follow-up via reverse KM (see scripts/reverse_km_followup.R)
fu_fit <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ 1,
  data = df_analysis, reverse = TRUE
)
print(quantile(fu_fit, q = 0.5))

# Kaplan-Meier
os_fit <- survival::survfit(
  survival::Surv(os_months, os_event) ~ arm,
  data = df_analysis
)
print(os_fit)
print(summary(os_fit, times = c(12, 24, 36, 60)))

# Cox model (univariable, see scripts/cox_univariable.R for full version)
cox_uni <- survival::coxph(
  survival::Surv(os_months, os_event) ~ arm,
  data = df_analysis
)
broom::tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE) |> print()

# PH check
zph <- survival::cox.zph(cox_uni)
print(zph)
if (any(zph$table[, "p"] < 0.05)) {
  warning("PH violated. Consider RMST (scripts/rmst.R) or stratification.")
}

# ---- 6. EXPORT --------------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Table 1 -> Word
tab1 |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = file.path(out_dir, "table1.docx"))

# KM figure -> PDF
km <- survminer::ggsurvplot(
  os_fit, data = df_analysis,
  risk.table = TRUE, conf.int = FALSE, pval = TRUE,
  xlab = "Time (months)", ylab = "Overall survival",
  ggtheme = ggplot2::theme_classic(base_size = 12),
  palette = c("#1F77B4", "#D62728")
)
pdf(file.path(out_dir, "fig_km_os.pdf"), width = 6.5, height = 5)
print(km)
dev.off()

cat("All outputs saved to ", out_dir, "\n")

# ---- 7. NEXT STEPS ----------------------------------------------------------
# - Add multivariable Cox: see scripts/cox_multivariable.R
# - If PH violated, add RMST: see scripts/rmst.R
# - For meta-analysis or indirect comparison: see scripts/ipd_from_km.R
# - For a forest plot: see scripts/forest_plot.R
