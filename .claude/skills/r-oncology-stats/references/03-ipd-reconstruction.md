# IPD reconstruction from published Kaplan–Meier curves

The workflow: extract (x, y) coordinates from a published KM figure → combine with numbers-at-risk → reconstruct pseudo-individual-patient data → validate against the original.

This is a two-step process. **Do not conflate them.** Digitization (Step 1) is image processing. IPD reconstruction (Step 2) is statistical and uses the Guyot et al. (2012) algorithm.

---

## Workflow overview

```
Published KM figure (PNG/PDF)
        │
        ▼
[Step 1] Digitize coordinates with WebPlotDigitizer (or similar)
        │   → CSV of (time, survival_probability) pairs per arm
        │
        ▼
[Step 2] Reconstruct pseudo-IPD with IPDfromKM
        │   → preprocess() → getIPD()
        │   → data.frame(time, event, arm)
        │
        ▼
[Step 3] Validate: overlay reconstructed KM on the original figure;
         check median, landmark survival, total events
        │
        ▼
[Step 4] Use the pseudo-IPD for downstream analyses
         (HR re-estimation, meta-analysis, indirect comparison)
```

---

## Step 1 — Digitization with WebPlotDigitizer

WebPlotDigitizer (https://apps.automeris.io/wpd/) runs in-browser. The user does this outside R; you advise on the workflow.

**What the user should extract:**

1. **Axis calibration** — define 4 reference points on the figure: (x_min, y_min), (x_max, y_min), (x_min, y_max), (x_max, y_max). For KM curves: x is time, y is survival probability (0 to 1, or 0% to 100%).
2. **One curve per arm.** Trace each arm separately. Save each as a separate CSV (or as separate datasets within WebPlotDigitizer).
3. **Sample density.** Sample heavily where the curve changes (early steep drops); sparser where the curve plateaus. Aim for 30–80 points per arm depending on the curve.
4. **Numbers at risk table.** Manually transcribe the numbers-at-risk table beneath the figure into a small CSV: rows = arms, columns = time points. **This is NOT optional** — the reconstruction algorithm needs it to estimate event/censoring distributions.
5. **Document everything**: figure source (paper, figure number, page), DPI used, any visual quirks (overlapping arms, low-resolution print).

**What to ask the user before reconstruction:**

- Did you extract numbers at risk? (If no, reconstruction quality will be poor.)
- Total events per arm (if reported in the paper's text/table)? — used as a quality check.
- Median survival per arm (if reported)? — used as a quality check.

---

## Step 2 — Reconstruct pseudo-IPD with `IPDfromKM`

```r
library(IPDfromKM)

# Per arm: digitized coordinates (time, survival)
# Plus numbers-at-risk table for that arm.

arm_a_pre <- preprocess(
  dat       = digi_arm_a,           # data.frame(time, survival) from WPD
  trisk     = c(0, 6, 12, 18, 24, 36, 48),   # time points where N-at-risk is reported
  nrisk     = c(120, 95, 71, 52, 38, 22, 9), # N at risk at those time points
  maxy      = 1                     # if survival on 0–1; use 100 if 0%–100%
)

arm_a_ipd <- getIPD(
  prep      = arm_a_pre,
  armID     = 1,                    # arbitrary numeric label
  tot.events = 78                   # total events reported in the paper (optional, improves fit)
)$IPD                               # → data.frame(time, status, treat)

# Repeat for arm B with armID = 2, then bind:
ipd_full <- dplyr::bind_rows(arm_a_ipd, arm_b_ipd)
```

Notes:

- The `IPDfromKM` package implements the Guyot et al. (2012) algorithm — see the citation footer for the paper.
- `trisk` and `nrisk` must align: same length, same arm.
- `maxy` flips the y-scale convention: `1` for 0–1 survival, `100` for 0–100%.
- Provide `tot.events` when known — it constrains the algorithm and improves accuracy substantially.

---

## Step 3 — Validate

This step is **mandatory**. If reconstruction is bad, downstream analyses are garbage.

### Visual overlay

```r
library(ggplot2)
library(survival)

# Reconstructed KM
reco_fit <- survfit(Surv(time, status) ~ treat, data = ipd_full)

# Plot reconstructed curve over the original
survminer::ggsurvplot(reco_fit, data = ipd_full,
  conf.int = FALSE, risk.table = TRUE,
  palette = c("#1F77B4", "#D62728"),
  ggtheme = theme_classic()
)
```

Compare side-by-side with the published figure. The reconstruction should be visually indistinguishable from the original at the resolution of the plot.

### Numeric checks

| Metric | Reconstructed | Reported in paper | Acceptable difference |
|---|---|---|---|
| Median survival, arm A | ? | ? | ± 5% or ± 1 unit (mo/yr) |
| Median survival, arm B | ? | ? | ± 5% or ± 1 unit |
| Landmark 12-mo survival | ? | ? | ± 2 percentage points |
| Landmark 24-mo survival | ? | ? | ± 2 percentage points |
| Total events, arm A | ? | ? | ± 1 |
| Total events, arm B | ? | ? | ± 1 |
| Hazard ratio (Cox) | ? | ? | ± 0.05 absolute |

If any of these are off by more than the tolerance, repeat digitization (Step 1) with more points or better axis calibration — do not "fix" reconstruction by tweaking the algorithm parameters.

---

## Step 4 — Downstream uses

Pseudo-IPD is suitable for:

- Re-estimating HR with 95% CI (Cox) when the original paper only reported a KM curve.
- Calculating RMST (`survRM2::rmst2`) for indirect treatment comparison.
- Pooling across trials (meta-analysis of HRs after refit, or pooled IPD meta-analysis).
- Computing landmark survival differences with CIs.
- Comparing different sub-arms when the published Cox model is not what you want.

Pseudo-IPD is NOT suitable for:

- Individual-level subgroup analyses (the algorithm assumes uniform censoring; subgroups are not preserved).
- Multivariable adjustment (no baseline covariates to adjust on).
- Any analysis where the published numbers-at-risk are missing.

---

## What to document in the methods section

For any paper using reconstructed IPD:

1. Source figure (paper citation + figure number + page).
2. Digitization tool (WebPlotDigitizer version) and operator.
3. Reconstruction tool (`IPDfromKM` version) and algorithm (cite Guyot 2012).
4. Validation results (visual overlay + numeric agreement vs published).
5. Sensitivity analyses if reconstruction was repeated.
6. Limitations: pseudo-IPD inherits uncertainty from digitization and from the published curve's resolution.

---

## Common failure modes

- **No numbers at risk in the figure.** Reconstruction quality drops sharply. If the user really needs IPD, ask whether the paper's supplement or original trial registry has them.
- **Curves overlap heavily on a low-resolution PNG.** Hand-tracing is noisy; sample more points and validate aggressively.
- **Stepped vs smooth digitization.** KM curves are step functions. When digitizing, click on the corners of the steps, not the diagonals — otherwise the reconstruction infers events at wrong times.
- **Y-axis as percentage but `maxy = 1`** (or vice versa) — the survival probabilities become wildly wrong. Always confirm the y-axis scale matches `maxy`.
- **Mismatched `trisk` and `nrisk` lengths** — `IPDfromKM::preprocess` will throw a cryptic error. Re-check the at-risk table.

---

## Reference

Guyot, P., Ades, A. E., Ouwens, M. J. N. M., & Welton, N. J. (2012). Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan–Meier survival curves. *BMC Medical Research Methodology*, 12, 9.
