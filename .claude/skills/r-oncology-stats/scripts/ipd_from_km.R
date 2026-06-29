# Reconstruct pseudo-IPD from a published KM curve
#
# Workflow:
#   Step 1 (outside R): digitize the curve coordinates with WebPlotDigitizer.
#                        Save one CSV per arm with columns (time, survival).
#                        Manually transcribe the numbers-at-risk table.
#   Step 2 (this script): reconstruct pseudo-IPD with IPDfromKM.
#   Step 3 (this script): validate against the original (KM overlay + numeric checks).
#
# Reference: Guyot P et al. BMC Med Res Methodol 2012;12:9.

# ---- packages ---------------------------------------------------------------
library(IPDfromKM)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(here)

# ---- inputs you must provide ------------------------------------------------
# >>> EDIT all four pieces:
#  1. Digitized (time, survival) CSV per arm.
#  2. Time points of the numbers-at-risk table.
#  3. Numbers at risk at those time points, per arm.
#  4. (Optional but recommended) Total events per arm from the paper.

# Arm A
digi_arm_a   <- readr::read_csv(here::here("data", "digi_arm_a.csv"),   # cols: time, survival
                                col_types = readr::cols(time = "n", survival = "n"))
trisk_arm_a  <- c(0, 6, 12, 18, 24, 36, 48)
nrisk_arm_a  <- c(120, 95, 71, 52, 38, 22, 9)
events_arm_a <- 78

# Arm B
digi_arm_b   <- readr::read_csv(here::here("data", "digi_arm_b.csv"),
                                col_types = readr::cols(time = "n", survival = "n"))
trisk_arm_b  <- c(0, 6, 12, 18, 24, 36, 48)
nrisk_arm_b  <- c(122, 110, 92, 78, 65, 48, 30)
events_arm_b <- 52

# Y-axis convention (0-1 or 0-100)
maxy_val <- 1   # use 100 if your digitized survival is on 0%-100% scale

# ---- preprocess + reconstruct ----------------------------------------------
pre_a <- IPDfromKM::preprocess(
  dat    = as.data.frame(digi_arm_a),
  trisk  = trisk_arm_a,
  nrisk  = nrisk_arm_a,
  maxy   = maxy_val
)
ipd_a <- IPDfromKM::getIPD(prep = pre_a, armID = 1, tot.events = events_arm_a)$IPD

pre_b <- IPDfromKM::preprocess(
  dat    = as.data.frame(digi_arm_b),
  trisk  = trisk_arm_b,
  nrisk  = nrisk_arm_b,
  maxy   = maxy_val
)
ipd_b <- IPDfromKM::getIPD(prep = pre_b, armID = 2, tot.events = events_arm_b)$IPD

# Combine
ipd_full <- dplyr::bind_rows(ipd_a, ipd_b) |>
  dplyr::rename(time = time, event = status, arm = treat) |>
  dplyr::mutate(arm = factor(arm, levels = c(1, 2), labels = c("Arm A", "Arm B")))

# ---- validation: visual overlay ---------------------------------------------
reco_fit <- survival::survfit(
  survival::Surv(time, event) ~ arm,
  data = ipd_full
)

km_overlay <- survminer::ggsurvplot(
  reco_fit, data = ipd_full,
  risk.table = TRUE, conf.int = FALSE, censor = FALSE,
  xlab = "Time (months)", ylab = "Survival probability (reconstructed)",
  legend.title = "",
  ggtheme = ggplot2::theme_classic(base_size = 12),
  palette = c("#1F77B4", "#D62728")
)
print(km_overlay)

# ---- validation: numeric checks ---------------------------------------------
landmarks <- c(12, 24, 36)
reco_summary <- summary(reco_fit, times = landmarks)
reco_tab <- tibble::tibble(
  arm   = rep(c("Arm A", "Arm B"), each = length(landmarks)),
  time  = rep(landmarks, 2),
  surv  = reco_summary$surv,
  lower = reco_summary$lower,
  upper = reco_summary$upper
)
print(reco_tab)

cat("\nReconstructed median survival (with 95% CI):\n")
print(reco_fit)

cat("\nReconstructed total events:\n")
print(
  ipd_full |>
    dplyr::group_by(arm) |>
    dplyr::summarise(n = dplyr::n(),
                     events = sum(event),
                     median_t = median(time))
)

# Compare reco_tab and reco_fit medians against the original paper's reported values.
# Tolerance: ~5% on survival probabilities, ±1 unit on median, ±1 on events.
# If outside tolerance, REPEAT digitization (Step 1) with more points or better calibration.

# ---- save reconstructed IPD -------------------------------------------------
out_dir <- here::here("output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(ipd_full, file.path(out_dir, "ipd_reconstructed.csv"))

pdf(file.path(out_dir, "fig_km_reconstructed.pdf"),
    width = 6.5, height = 5, useDingbats = FALSE)
print(km_overlay)
dev.off()

cat("Saved reconstructed IPD + KM overlay to ", out_dir, "\n")
cat("\nNEXT: visually compare fig_km_reconstructed.pdf against the published figure.\n",
    "If they don't match within tolerance, repeat digitization, not the reconstruction step.\n",
    sep = "")
