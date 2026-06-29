# Forest plot — three flavors
#
# A. Model-based forest from a Cox model      (forestmodel::forest_model)
# B. Subgroup forest with treatment HR per stratum  (custom ggplot)
# C. Meta-analysis forest of pooled HRs       (meta::metagen + meta::forest)
#
# Pick the section that matches the use case.

# ---- packages ---------------------------------------------------------------
library(forestmodel)
library(ggplot2)
library(dplyr)
library(forcats)
library(here)
library(survival)
# library(meta)         # uncomment for section C

# =============================================================================
# A. Model-based forest from a Cox model
# =============================================================================
# >>> EDIT: assume you already have a fitted multivariable Cox model `cox_mv`
# cox_mv <- coxph(Surv(os_months, os_event) ~ arm + age + sex + ecog + stage,
#                 data = df_analysis)

fm <- forestmodel::forest_model(
  cox_mv,
  recalculate_width = TRUE
)
print(fm)

out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  filename = file.path(out_dir, "fig_forest_cox.pdf"),
  plot     = fm,
  width    = 7, height = 5, units = "in",
  device   = grDevices::cairo_pdf
)

# =============================================================================
# B. Subgroup forest — treatment HR within each subgroup
# =============================================================================
# Computes a Cox HR for treatment within each level of each subgroup variable.
# Reports interaction p (treatment × subgroup) for each variable.

df_analysis <- readr::read_csv(here::here("data", "trial.csv"))

# >>> EDIT: subgroup variables
subgroup_vars <- c("age_group", "sex", "ecog", "stage")

# Helper: compute treatment HR within one subgroup level
subgroup_hr <- function(data, sub_var, sub_level) {
  d <- dplyr::filter(data, .data[[sub_var]] == sub_level)
  if (nrow(d) < 10 || sum(d$os_event) < 5) {
    return(tibble::tibble(subgroup = sub_var, level = sub_level,
                          n = nrow(d), events = sum(d$os_event),
                          hr = NA, lower = NA, upper = NA))
  }
  f <- survival::coxph(survival::Surv(os_months, os_event) ~ arm, data = d)
  ci <- exp(confint(f))
  tibble::tibble(
    subgroup = sub_var,
    level    = as.character(sub_level),
    n        = nrow(d),
    events   = sum(d$os_event),
    hr       = exp(coef(f))[1],
    lower    = ci[1, 1],
    upper    = ci[1, 2]
  )
}

# Build long table over all subgroup × level combos
sub_tab <- purrr::map_dfr(subgroup_vars, function(sv) {
  levs <- unique(df_analysis[[sv]])
  purrr::map_dfr(levs, ~ subgroup_hr(df_analysis, sv, .x))
}) |>
  dplyr::mutate(label = paste0(subgroup, ": ", level))

# Interaction p-values (LRT) per subgroup variable
int_p <- purrr::map_dfr(subgroup_vars, function(sv) {
  f0 <- survival::coxph(
    as.formula(paste0("survival::Surv(os_months, os_event) ~ arm + ", sv)),
    data = df_analysis
  )
  f1 <- survival::coxph(
    as.formula(paste0("survival::Surv(os_months, os_event) ~ arm * ", sv)),
    data = df_analysis
  )
  tibble::tibble(
    subgroup = sv,
    p_interaction = anova(f0, f1)$"P(>|Chi|)"[2]
  )
})
print(int_p)

# Forest
sub_plot <- ggplot(sub_tab, aes(x = hr, y = forcats::fct_rev(label))) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.35) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_log10() +
  labs(
    x = "Hazard ratio (treatment vs control), 95% CI",
    y = NULL,
    caption = "Log-scale x-axis. HR < 1 favors treatment."
  ) +
  theme_classic(base_size = 11)

print(sub_plot)

ggplot2::ggsave(
  filename = file.path(out_dir, "fig_forest_subgroup.pdf"),
  plot     = sub_plot,
  width    = 7, height = 5, units = "in",
  device   = grDevices::cairo_pdf
)

# =============================================================================
# C. Meta-analysis forest — pooled HRs across studies
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
