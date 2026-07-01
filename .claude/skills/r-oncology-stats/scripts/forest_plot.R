# Forest plot -- three flavors
#
# A. Cox MV coefficients        (custom oncology layout + forestmodel QC)
# B. Subgroup forest             (treatment HR per stratum + interaction p)
# C. Meta-analysis forest        (pooled HRs across studies)
#
# Sections A + B: see references/11-forest-plots.md §3-§6 for methodology
# (precision-scaled squares, "HR = 1.00 (1.00-1.00)" trap, forestmodel crash).
# Section C: see references/05-manuscript-figures.md Meta-analysis section.
#
# Pick the section that matches the use case.

# ---- packages ---------------------------------------------------------------
library(survival)
library(broom.helpers)
library(ggplot2)
library(patchwork)
library(scales)
library(gtsummary)
library(flextable)
library(forestmodel)
library(forcats)
library(purrr)
library(dplyr)
library(here)
library(tibble)
# library(meta)         # uncomment for section C

# Load the 7 helpers (plot_forest_oncology, plot_forest_subgroup,
# nonest_terms, mark_nonest_tbl, save_forestmodel, build_subgroup_tidy,
# interaction_p). Mirrored from references/11-forest-plots.md.
source(here::here("scripts", "forest_helpers.R"))

out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# A. Cox MV coefficients -- custom oncology layout + forestmodel QC
# =============================================================================
# >>> EDIT: assume you already have a fitted multivariable Cox model `cox_mv`
# cox_mv <- coxph(Surv(os_months, os_event) ~ arm + age + sex + ecog + stage,
#                 data = df_analysis)

# Print summary to the run log so the plot can be cross-checked against
# the raw Wald HRs. Also flag any non-estimable terms upfront.
summary(cox_mv)
bad <- nonest_terms(cox_mv)
if (length(bad) > 0) {
  message("Non-estimable Cox terms detected (rendered as 'Not estimable' ",
          "in the forest and blanked in the .docx table): ",
          paste(bad, collapse = ", "))
}

# Manuscript figure -- precision-scaled two-panel layout
plot_forest_oncology(
  model          = cox_mv,
  file_base      = "fig_forest_cox_mv",
  endpoint_label = "Overall survival",           # >>> EDIT endpoint
  out_dir        = out_dir
)

# QC sanity check -- forestmodel default (uniform squares), with degenerate
# variables auto-excluded so the base-R log-axis tick generator can't crash
# with "_LARGE_ range: invalid {xy}axp or par".
save_forestmodel(cox_mv,
                 file_base = "fig_forest_cox_mv_qc",
                 out_dir   = out_dir)

# =============================================================================
# B. Subgroup forest -- treatment HR within each subgroup
# =============================================================================
# Fits ONE Cox model per subgroup level (e.g., one per age band). Reports
# ONE interaction p-value per subgroup variable (not per row). Use ONLY for
# pre-specified subgroups -- see reference doc §5.1 for the fishing warning.

df_analysis <- readr::read_csv(here::here("data", "trial.csv"))

# >>> EDIT: subgroup variables, time / event / treatment columns
subgroup_vars <- c("age_group", "sex", "ecog", "stage")
time_var  <- "os_months"
event_var <- "os_event"
treat_var <- "arm"

# One forest per subgroup variable
for (sv in subgroup_vars) {
  sub_tidy <- build_subgroup_tidy(
    data         = df_analysis,
    time_var     = time_var,
    event_var    = event_var,
    treat_var    = treat_var,
    subgroup_var = sv,
    min_events   = 10                # strata below this are "Not estimable"
  )
  p_int <- interaction_p(df_analysis, time_var, event_var, treat_var, sv)
  message(sprintf("Interaction p (%s * %s): %.4f", treat_var, sv, p_int))

  plot_forest_subgroup(
    tidy_df        = sub_tidy,
    file_base      = paste0("fig_forest_subgroup_", sv),
    endpoint_label = "Overall survival",       # >>> EDIT endpoint
    p_interaction  = p_int,
    out_dir        = out_dir
  )
}

# =============================================================================
# C. Meta-analysis forest -- pooled HRs across studies
# =============================================================================
# Suppose you have a tibble: study, log_hr, se_log_hr (one row per study).
# Typical source: reconstructed IPD per published curve (see ipd_from_km.R)
# or extracted directly from each paper's reported HR + CI.
#
# library(meta)
#
# pooled <- tibble::tibble(
#   study     = c("Trial A", "Trial B", "Trial C"),
#   log_hr    = c(log(0.72), log(0.81), log(0.65)),
#   se_log_hr = c(0.12,      0.18,      0.10)
# )
#
# m <- meta::metagen(
#   TE = log_hr, seTE = se_log_hr, studlab = study,
#   data = pooled, sm = "HR",
#   common = FALSE, random = TRUE      # random-effects model
# )
# summary(m)
#
# pdf(file.path(out_dir, "fig_forest_metaanalysis.pdf"), width = 8, height = 4)
# meta::forest(m, leftcols = c("studlab"), rightcols = c("effect", "ci"),
#              xlab = "Hazard ratio (95% CI)", smlab = "Pooled HR (random effects)")
# dev.off()

cat("Forest plots saved to ", out_dir, "\n")
