# Baseline characteristics — manuscript-ready Table 1
#
# Inputs:   patient-level data frame; variables to summarize defined below.
# Outputs:  gtsummary Table 1 with appropriate summary stats and tests,
#           exported to .docx for journal submission.
# Default:  continuous = median (IQR), categorical = n (%); appropriate tests.
# Edit:     `selected_vars` and `label` list to match your dataset.

# ---- packages ---------------------------------------------------------------
library(gtsummary)
library(flextable)
library(dplyr)
library(forcats)
library(labelled)
library(here)

# ---- data + factor setup ----------------------------------------------------
df_analysis <- readr::read_csv(here::here("data", "trial.csv")) |>
  dplyr::mutate(
    arm   = forcats::fct_relevel(arm,   "Control"),                # >>> EDIT
    ecog  = forcats::fct_relevel(factor(ecog), "0"),
    stage = forcats::fct_relevel(factor(stage), "I"),
    sex   = forcats::fct_relevel(sex,   "Female")
  ) |>
  labelled::set_variable_labels(                                   # >>> EDIT labels
    age     = "Age (years)",
    sex     = "Sex",
    ecog    = "ECOG performance status",
    stage   = "Disease stage (AJCC v8)",
    smoking = "Smoking history",
    ldh     = "LDH (U/L)"
  )

# ---- journal theme ----------------------------------------------------------
gtsummary::theme_gtsummary_journal(journal = "jama")               # >>> EDIT journal

# ---- build Table 1 ----------------------------------------------------------
# >>> EDIT: which variables to include
selected_vars <- c("age", "sex", "ecog", "stage", "smoking", "ldh")

tab1 <- df_analysis |>
  dplyr::select(dplyr::all_of(selected_vars), arm) |>
  gtsummary::tbl_summary(
    by = arm,
    missing = "ifany",
    statistic = list(
      gtsummary::all_continuous()  ~ "{median} ({p25}, {p75})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(gtsummary::all_continuous() ~ 1)
  ) |>
  gtsummary::add_n() |>
  gtsummary::add_overall() |>
  # Optional: add p-values. For randomized trials, consider OMITTING — see references/04-baseline-tables.md.
  gtsummary::add_p(
    test = list(
      gtsummary::all_continuous()  ~ "wilcox.test",
      gtsummary::all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3)
  ) |>
  gtsummary::modify_header(label = "**Characteristic**") |>
  gtsummary::bold_labels()

tab1

# ---- export to Word ---------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tab1 |>
  gtsummary::as_flex_table() |>
  flextable::set_table_properties(layout = "autofit", width = 1) |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::font(fontname = "Times New Roman", part = "all") |>
  flextable::save_as_docx(path = file.path(out_dir, "table1_baseline.docx"))

# HTML option:
# tab1 |> gtsummary::as_gt() |> gt::gtsave(file.path(out_dir, "table1_baseline.html"))

cat("Saved: ", file.path(out_dir, "table1_baseline.docx"), "\n")

# ---- when to skip p-values --------------------------------------------------
# In a randomized trial, baseline imbalances are by chance and p-values are not
# informative. Many journals (NEJM, JAMA) discourage them in Table 1. To remove:
#
#   tab1_no_p <- df_analysis |>
#     dplyr::select(all_of(selected_vars), arm) |>
#     gtsummary::tbl_summary(by = arm, missing = "ifany") |>
#     gtsummary::add_overall() |>
#     gtsummary::add_n()
