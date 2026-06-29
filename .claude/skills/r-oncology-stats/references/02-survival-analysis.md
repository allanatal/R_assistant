# Survival analysis methodology

Covers: Kaplan–Meier estimation, Cox proportional hazards, PH diagnostics, RMST, reverse-KM follow-up, landmark analysis. The corresponding paste-and-adapt code lives in `scripts/`.

---

## 1. Kaplan–Meier

### Fit

```r
os_fit <- survival::survfit(
  survival::Surv(time = os_months, event = os_event) ~ arm,
  data = df_analysis
)
```

- `Surv()` arguments: `time` = follow-up time, `event` = 1/0 (event/censored). If event is coded differently, recode first; never paper over miscoding with `event == "Dead"` inside `Surv()`.
- For one-arm (whole-cohort) estimation: `~ 1`.
- For stratified curves: `~ strata(stratum) + arm` — note this is for stratified Cox; for KM curves stratification means separate curves per stratum level.

### Read off

- `summary(os_fit)` — survival probabilities at each event time.
- `summary(os_fit, times = c(12, 24, 60))` — landmark survival at chosen time points.
- `survival::survdiff(Surv(time, event) ~ arm, data = df)` — log-rank test.
- Median survival: `os_fit` print method shows medians + 95% CI. If `NA`, the curve has not reached 0.5 — report "median not reached (NR)" and quote landmark probabilities instead.

### Visualize

Use `survminer::ggsurvplot()` for most cases (more layout knobs, risk table is mature). Use `ggsurvfit::ggsurvfit()` for tidier ggplot composability when you need to layer or facet heavily.

Minimal manuscript-quality call:

```r
survminer::ggsurvplot(
  os_fit, data = df_analysis,
  risk.table = TRUE, risk.table.height = 0.22,
  conf.int = FALSE, censor = TRUE,
  pval = TRUE, pval.method = TRUE,
  xlab = "Time (months)",
  ylab = "Overall survival",          # ⚠ specific endpoint, not "Cumulative survival"
  legend.title = "", legend.labs = c("Arm A", "Arm B"),
  ggtheme = ggplot2::theme_classic(base_size = 12),
  palette = c("#1F77B4", "#D62728")
)
```

See `references/05-manuscript-figures.md` for export, theme, and dimension guidance.

---

## 2. Cox proportional hazards

### Univariable

```r
cox_uni <- survival::coxph(
  survival::Surv(os_months, os_event) ~ arm,
  data = df_analysis
)
broom::tidy(cox_uni, exponentiate = TRUE, conf.int = TRUE)
```

- Output the HR (`exponentiate = TRUE`), 95% CI, p-value.
- For categorical covariates with > 2 levels, present every non-reference level; do NOT collapse silently.

### Multivariable

```r
cox_mv <- survival::coxph(
  survival::Surv(os_months, os_event) ~ arm + age + sex + ecog + stage,
  data = df_analysis
)
```

Principles:

- **Clinically meaningful covariate selection.** Discourage stepwise (AIC/BIC) selection unless the user has a real reason.
- **One event per ~10 covariates** is the classic rule of thumb (Peduzzi 1995). Warn if violated.
- **Categorical reference levels** must make clinical sense (ECOG 0, stage I, wild-type, control arm).
- **Continuous covariates**: don't dichotomize unnecessarily. If you must categorize, use clinically established thresholds, not data-driven median splits.
- **Interactions**: only fit pre-specified ones, or treatment × subgroup interactions for subgroup forest plots.

### PH diagnostics

Always run after any Cox fit:

```r
zph <- survival::cox.zph(cox_mv)
print(zph)            # per-covariate and global Schoenfeld test
plot(zph)             # Schoenfeld residuals over time per covariate
```

Interpretation:

- Global p > 0.05 and per-covariate p > 0.05 → PH assumption not contradicted; proceed.
- Any covariate p < 0.05 OR visibly non-flat smooth → PH likely violated for that covariate.

If violated, options (let the user choose):

- **Stratify** the offending covariate (`+ strata(var)`) — removes its effect from estimation while controlling for it.
- **Time-varying coefficient** with `tt()` in `coxph` — preserves the effect estimate but allows it to vary.
- **Switch to RMST** for the primary treatment effect — model-free, no PH assumption (see section 4).
- **Restricted-time Cox** (`tstart`/`tstop` setup with piecewise effects).
- **Milestone survival** — report 1y/2y/5y survival differences with CIs as a non-HR effect measure.

Document the choice and rationale in the analysis output.

---

## 3. Reverse-KM follow-up

Median follow-up = the median time patients would have been followed if no events had occurred. The standard estimator is reverse KM (Schemper & Smith 1996).

### Overall

```r
# Option A: prodlim
fu_fit <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ 1,
  data = df_analysis, reverse = TRUE
)
quantile(fu_fit, q = 0.5)   # median follow-up with 95% CI

# Option B: survival package
fu_fit2 <- survival::survfit(
  survival::Surv(os_months, 1 - os_event) ~ 1,
  data = df_analysis
)
print(fu_fit2)              # median + 95% CI
```

Both give the same answer. `prodlim` integrates nicely with multi-state / competing-risks setups.

### Stratified

```r
fu_by_arm <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ arm,
  data = df_analysis, reverse = TRUE
)
quantile(fu_by_arm, q = 0.5)
```

Report follow-up by treatment arm if the arms have meaningfully different accrual or censoring patterns (e.g., asymmetric administrative censoring).

### Anti-patterns

- ❌ `median(df$os_months)` — that's median observed time, not median follow-up.
- ❌ `max(df$os_months)` — that's maximum observed time, also not follow-up.
- ❌ `median(df$os_months[df$os_event == 0])` — median among censored only; biased toward early censoring.

---

## 4. Restricted Mean Survival Time (RMST)

RMST = expected survival time up to a fixed truncation time τ. Model-free, no PH assumption, interpretable as "average months alive in the first τ months".

### When to use

- PH assumption is violated.
- Curves cross or the late tail is informative (immunotherapy, long-term survivorship).
- Communicating effect size to clinicians ("the new treatment buys 3.2 months over 24 months") is more useful than an HR.

### Fit

```r
rmst_res <- survRM2::rmst2(
  time   = df_analysis$os_months,
  status = df_analysis$os_event,
  arm    = df_analysis$arm_binary,   # MUST be 0/1
  tau    = 36                        # truncation in same unit as `time`
)
print(rmst_res)
```

### Choosing τ

τ must be ≤ the smaller of: (a) the last observed event time in each arm, (b) a time at which both arms have meaningful numbers at risk (rule of thumb: ≥ 10% of the original cohort). Document the rationale.

Common choices:

- Pre-specified clinical milestone (24 mo, 36 mo, 60 mo).
- Time at which the smaller arm's risk set drops below a threshold (e.g., 5–10%).
- Maximum observation time that is reasonably common to both arms.

### Report

- RMST per arm + 95% CI.
- RMST difference + 95% CI + p-value.
- RMST ratio + 95% CI (alternative).
- The chosen τ and its justification.

### Caveats

- RMST requires both arms to have data extending to τ — extrapolation is not valid.
- For > 2 groups, use pairwise RMST differences (no native multi-group test in `survRM2`).
- `survRM2::rmst2` requires arm coded 0/1. Recode before calling.

---

## 5. Landmark analysis

Useful when:

- Comparing groups defined by an event that happens during follow-up (e.g., responders vs non-responders), to avoid immortal-time bias.
- Reporting survival probabilities at fixed time points (1y, 2y, 5y).

### Landmark for response (immortal-time)

```r
landmark_t <- 3   # months
df_landmark <- df_analysis |>
  dplyr::filter(os_months >= landmark_t) |>
  dplyr::mutate(
    os_months_lm = os_months - landmark_t,
    response_by_lm = ifelse(time_to_response <= landmark_t & response == 1, 1, 0)
  )
```

Then run KM / Cox on `os_months_lm` against `response_by_lm`.

### Landmark survival probabilities

```r
summary(os_fit, times = c(12, 24, 36, 60))
# extract surv, lower, upper, n.risk per arm at each time
```

Report as % (3-sig-fig CIs): "12-month OS was 84.2% (95% CI 78.1–90.7) in arm A vs 71.5% (95% CI 64.8–78.9) in arm B."

---

## 6. Common pitfalls

- **Time-zero ambiguity.** Confirm what time 0 means — randomization, diagnosis, treatment start, surgery. Different choices give different KM curves.
- **Left truncation / delayed entry.** If patients enter after time 0 (e.g., registry data with diagnosis-to-enrollment lag), use `Surv(tstart, tstop, event)` syntax — failing to do so biases the estimate.
- **Recurrent events.** A standard Cox treats only the first event. If recurrence matters, consider Andersen-Gill, PWP, or frailty models — flag and discuss; out of scope for v0.1.0 templates.
- **Competing risks.** If non-cancer death is non-trivial (older / comorbid cohorts), a cause-specific Cox biases the cause-specific effect; consider Fine-Gray subdistribution hazard (`cmprsk` or `tidycmprsk`). Flag as a v0.2.0 addition.
- **Crossing curves.** If KM curves cross, the log-rank test loses power and a single HR is misleading. Report RMST and milestone survival, not just an HR.
- **Tied event times.** Default Cox tie handler is Efron — usually fine. Switch to `method = "exact"` only when ties are extreme and the sample is small.

---

## References (further reading, not for citation insertion)

- Schemper & Smith. Stat Med 1996 — reverse-KM follow-up.
- Royston & Parmar. BMC Med Res Methodol 2013 — RMST methodology.
- Uno et al. JCO 2014 — RMST for clinical communication.
- Therneau & Grambsch — "Modeling Survival Data: Extending the Cox Model" (canonical Cox reference).
