# Kaplan-Meier curve template — manuscript-ready
#
# Inputs:   a data frame with one row per patient containing:
#             - time-to-event variable (numeric, in months unless otherwise stated)
#             - event/status variable  (1 = event, 0 = censored)
#             - group variable         (factor; reference level set explicitly)
# Outputs:  KM plot with risk table, log-rank p-value, manuscript theme
#           saved as PDF to output/
# Assumes:  event coding has been verified; time unit is months; specific endpoint name known.
# Edit:     >>> EDIT markers indicate places to swap for the user's data.

# ---- packages ---------------------------------------------------------------
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(forcats)
library(here)

# ---- data -------------------------------------------------------------------
# >>> EDIT: load your data
df_analysis <- readr::read_csv(here::here("data", "trial.csv"))

# >>> EDIT: set reference level for group variable
df_analysis <- df_analysis |>
  dplyr::mutate(arm = forcats::fct_relevel(arm, "Control"))

# ---- sanity checks ----------------------------------------------------------
stopifnot(
  all(df_analysis$os_event %in% c(0, 1)),                # event coding
  all(df_analysis$os_months >= 0, na.rm = TRUE),         # non-negative time
  is.factor(df_analysis$arm)
)
cat("Events per arm:\n")
print(
  df_analysis |>
    dplyr::group_by(arm) |>
    dplyr::summarise(n = dplyr::n(),
                     events = sum(os_event),
                     median_t = median(os_months, na.rm = TRUE))
)

# ---- fit --------------------------------------------------------------------
# >>> EDIT: swap variable names if needed
os_fit <- survival::survfit(
  survival::Surv(time = os_months, event = os_event) ~ arm,
  data = df_analysis
)

# Median survival + 95% CI
print(os_fit)

# Landmark survival probabilities (1y / 2y / 3y / 5y)
print(summary(os_fit, times = c(12, 24, 36, 60)))

# Log-rank test
print(
  survival::survdiff(survival::Surv(os_months, os_event) ~ arm, data = df_analysis)
)

# ---- plot -------------------------------------------------------------------
# >>> EDIT: axis labels, legend labels, palette, dimensions
km <- survminer::ggsurvplot(
  os_fit, data = df_analysis,

  risk.table        = TRUE,
  risk.table.height = 0.22,
  conf.int          = FALSE,
  censor            = TRUE,
  censor.shape      = "|",
  censor.size       = 2.5,

  pval              = TRUE,
  pval.method       = TRUE,
  pval.coord        = c(3, 0.05),
  pval.size         = 4,

  xlab              = "Time since randomization (months)",
  ylab              = "Overall survival",            # >>> EDIT: specific endpoint
  xlim              = c(0, 60),
  break.x.by        = 12,
  ylim              = c(0, 1),
  break.y.by        = 0.2,

  legend.title      = "",
  legend.labs       = c("Standard of care", "Experimental"),   # >>> EDIT
  legend            = c(0.8, 0.85),

  ggtheme           = ggplot2::theme_classic(base_size = 12),
  palette           = c("#1F77B4", "#D62728"),

  risk.table.y.text = FALSE,
  tables.theme      = survminer::theme_cleantable()
)

print(km)

# ---- export -----------------------------------------------------------------
# Manuscript-quality vector PDF; survminer composite needs base device + dev.off()
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(out_dir, "fig_km_os.pdf"),
    width = 6.5, height = 5, useDingbats = FALSE)
print(km)
dev.off()

# Alternative: TIFF 600 DPI if a journal requires raster
# tiff(file.path(out_dir, "fig_km_os.tiff"),
#      width = 6.5, height = 5, units = "in", res = 600, compression = "lzw")
# print(km)
# dev.off()

cat("Saved: ", file.path(out_dir, "fig_km_os.pdf"), "\n")
