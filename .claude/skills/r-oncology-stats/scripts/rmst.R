# Restricted Mean Survival Time (RMST) — survRM2
#
# Inputs:   data frame with time, status (0/1), arm coded 0/1.
# Outputs:  RMST per arm + 95% CI, RMST difference + 95% CI + p-value,
#           RMST ratio + 95% CI; KM overlay annotated with tau.
# Assumes:  event coding verified; tau chosen with explicit clinical or
#           at-risk-based rationale; both arms have meaningful data up to tau.
#
# When to prefer RMST over Cox HR:
#   - Proportional hazards assumption violated (cox.zph p < 0.05).
#   - Curves cross visibly.
#   - The clinical message is "average months alive over X months", not a constant ratio.

# ---- packages ---------------------------------------------------------------
library(survRM2)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(here)

# ---- data + binary arm coding ----------------------------------------------
df_analysis <- readr::read_csv(here::here("data", "trial.csv")) |>
  dplyr::mutate(
    arm_binary = dplyr::if_else(arm == "Experimental", 1L, 0L)   # >>> EDIT levels
  )

# ---- choose tau -------------------------------------------------------------
# Tau MUST be:
#   (a) <= last observed time in each arm
#   (b) at a point where both arms have meaningful at-risk count (>= ~10% of original n)
# Document the rationale.

# >>> EDIT: pick tau (in same time unit as `time`)
tau_months <- 36     # e.g., 3-year RMST, pre-specified per protocol

# Diagnostic: check at-risk count at tau per arm
at_risk_at_tau <- df_analysis |>
  dplyr::group_by(arm_binary) |>
  dplyr::summarise(
    n_total    = dplyr::n(),
    n_at_risk  = sum(os_months >= tau_months),
    pct_at_risk = round(100 * sum(os_months >= tau_months) / dplyr::n(), 1)
  )
print(at_risk_at_tau)

if (any(at_risk_at_tau$pct_at_risk < 10)) {
  warning(sprintf(
    "At tau = %s, at-risk fraction < 10%% in at least one arm. ",
    tau_months
  ),
  "RMST estimates near tau will be unstable. Consider a smaller tau.")
}

# ---- fit RMST ---------------------------------------------------------------
rmst_res <- survRM2::rmst2(
  time   = df_analysis$os_months,
  status = df_analysis$os_event,
  arm    = df_analysis$arm_binary,
  tau    = tau_months
)
print(rmst_res)

# Key outputs:
#   - rmst_res$RMST.arm0$result : RMST in arm 0 + 95% CI
#   - rmst_res$RMST.arm1$result : RMST in arm 1 + 95% CI
#   - rmst_res$unadjusted.result : differences and ratios (Est, CI, p)

# ---- formatted summary ------------------------------------------------------
fmt <- function(x) formatC(x, format = "f", digits = 2)
diff_row  <- rmst_res$unadjusted.result["RMST (arm=1)-(arm=0)", ]
ratio_row <- rmst_res$unadjusted.result["RMST (arm=1)/(arm=0)", ]

cat(sprintf(
  "RMST (Arm 0 = Control): %s months (95%% CI %s-%s)\n",
  fmt(rmst_res$RMST.arm0$result["Est.", "Est."]),
  fmt(rmst_res$RMST.arm0$result["Est.", "lower .95"]),
  fmt(rmst_res$RMST.arm0$result["Est.", "upper .95"])
))
cat(sprintf(
  "RMST (Arm 1 = Experimental): %s months (95%% CI %s-%s)\n",
  fmt(rmst_res$RMST.arm1$result["Est.", "Est."]),
  fmt(rmst_res$RMST.arm1$result["Est.", "lower .95"]),
  fmt(rmst_res$RMST.arm1$result["Est.", "upper .95"])
))
cat(sprintf(
  "RMST difference (Exp - Ctrl): %s months (95%% CI %s-%s), p = %s\n",
  fmt(diff_row["Est."]),
  fmt(diff_row["lower .95"]),
  fmt(diff_row["upper .95"]),
  signif(diff_row["p"], 3)
))
cat(sprintf(
  "RMST ratio (Exp / Ctrl): %s (95%% CI %s-%s), p = %s\n",
  fmt(ratio_row["Est."]),
  fmt(ratio_row["lower .95"]),
  fmt(ratio_row["upper .95"]),
  signif(ratio_row["p"], 3)
))
cat(sprintf("Tau used: %s months\n", tau_months))

# ---- KM overlay annotated with tau -----------------------------------------
os_fit <- survival::survfit(
  survival::Surv(os_months, os_event) ~ arm_binary,
  data = df_analysis
)
km <- survminer::ggsurvplot(
  os_fit, data = df_analysis,
  risk.table = TRUE, conf.int = FALSE, censor = TRUE,
  xlab = "Time (months)", ylab = "Overall survival",
  legend.labs = c("Control", "Experimental"),
  ggtheme = ggplot2::theme_classic(base_size = 12),
  palette = c("#1F77B4", "#D62728")
)
km$plot <- km$plot +
  ggplot2::geom_vline(xintercept = tau_months, linetype = "dashed", color = "grey40") +
  ggplot2::annotate("text", x = tau_months, y = 0.05,
                    label = paste0("tau = ", tau_months, " mo"),
                    hjust = -0.1, size = 3.2, color = "grey30")

# ---- export -----------------------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(out_dir, "fig_km_rmst.pdf"),
    width = 6.5, height = 5, useDingbats = FALSE)
print(km)
dev.off()

# Save RMST results as tidy data frame
rmst_summary <- tibble::tibble(
  comparison = c("RMST_arm0", "RMST_arm1", "RMST_diff", "RMST_ratio"),
  estimate   = c(rmst_res$RMST.arm0$result["Est.", "Est."],
                 rmst_res$RMST.arm1$result["Est.", "Est."],
                 diff_row["Est."],
                 ratio_row["Est."]),
  ci_low     = c(rmst_res$RMST.arm0$result["Est.", "lower .95"],
                 rmst_res$RMST.arm1$result["Est.", "lower .95"],
                 diff_row["lower .95"],
                 ratio_row["lower .95"]),
  ci_high    = c(rmst_res$RMST.arm0$result["Est.", "upper .95"],
                 rmst_res$RMST.arm1$result["Est.", "upper .95"],
                 diff_row["upper .95"],
                 ratio_row["upper .95"]),
  p_value    = c(NA, NA, diff_row["p"], ratio_row["p"]),
  tau        = tau_months
)
readr::write_csv(rmst_summary, file.path(out_dir, "rmst_summary.csv"))

cat("Saved RMST KM figure + summary CSV to ", out_dir, "\n")
