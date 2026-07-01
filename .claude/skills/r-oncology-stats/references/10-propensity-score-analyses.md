# Propensity-score analyses for retrospective oncology cohorts

Covers: when (and when NOT) to use propensity-score (PS) methods, estimand selection, NCDB-specific guardrails, covariate timing, balance assessment, and outcome modeling for survival endpoints. Includes paste-and-adapt code for PS-adjusted Cox, matching (`MatchIt`), IPTW (`WeightIt`), overlap weighting, RMST after PS adjustment, and manuscript output.

PS methods are easy to misuse. **Do not generate code first.** Clarify the causal question, pick an estimand, check overlap and balance, then fit the outcome model.

---

## 1. When to use this reference

Consult this reference whenever the user asks about any of the following:

- Propensity-score matching
- Propensity-score adjustment (PS as covariate in a Cox model)
- PS-adjusted Cox models
- Inverse probability of treatment weighting (IPTW)
- ATT, ATE, or ATC weighting
- Overlap weighting (ATO)
- Matching with `MatchIt`
- Weighting with `WeightIt` or `PSweight`
- Balance assessment, Love plots, standardized mean differences (SMDs)
- Reviewer-requested confounding adjustment or sensitivity analyses
- NCDB matched or weighted survival analyses
- Any retrospective oncology cohort study where treatment/exposure assignment is non-random

For mechanics of KM, Cox, PH diagnostics, reverse-KM follow-up, and RMST, see `references/02-survival-analysis.md`. For data inspection see `references/07-data-inspection.md`. For Table 1 see `references/04-baseline-tables.md`. For figure/table export see `references/05-manuscript-figures.md` and `references/06-manuscript-tables.md`.

---

## 2. Decide whether PS analysis is appropriate

PS methods address **measured confounding** in observational/retrospective cohorts where treatment assignment is non-random. They are **not** a magic bridge to randomized-trial-quality causal inference.

PS methods do NOT fix:

- Unmeasured confounding (the variable you wish you had)
- Selection bias (who entered the cohort)
- Immortal time bias (treatment defined by surviving long enough to receive it)
- Informative censoring
- Incorrect endpoint definitions
- Poor data quality / coding errors
- Missing key confounders
- Poor overlap / positivity violations
- Residual confounding from registry limitations (granularity, missing fields)

Never frame a PS analysis as making a retrospective study "equivalent to a randomized trial." Frame it as **addressing measured confounding only**, with explicit limitations.

When ordinary multivariable Cox already adjusts for the same covariates, PS methods provide a complementary view (different functional form for the confounding adjustment, plus explicit balance diagnostics). They do not provide a stronger causal claim from the same data.

---

## 3. Define the clinical question and estimand BEFORE coding

Clarify these before estimating any propensity score. If unclear, ask the user — do not assume.

### Causal-question checklist

- **Treatment / exposure**: what defines the exposed group? (drug, regimen, surgery, RT, sequence, timing)
- **Comparator / reference**: what is "control"? (no treatment, alternative treatment, standard of care)
- **Primary endpoint**: OS, DFS, PFS, EFS, RFS, cancer-specific survival, recurrence — and is it actually available in the dataset?
- **Time zero / index date**: diagnosis, randomization-analog, treatment start, surgery date?
- **Eligibility window**: who is in the analytic cohort and why?
- **Censoring definition**: what censors a patient?
- **Was treatment baseline, time-fixed, time-dependent, or duration-based?**
- **Immortal-time risk**: must the patient survive (event-free) to "receive" the exposure? If yes, address it before any PS work (landmark, time-dependent treatment, eligibility restriction).
- **Target population**: who do we want to generalize to?

### Estimands — pick one and report it explicitly

| Estimand | Interpretation | Typical method |
|---|---|---|
| **ATE** (average treatment effect) | Effect if the whole eligible population were treated vs untreated | IPTW with `estimand = "ATE"`, PS covariate adjustment |
| **ATT** (avg. treatment effect on the treated) | Effect among patients who actually received treatment | 1:1 / k:1 matching, IPTW with `estimand = "ATT"` |
| **ATC** (avg. treatment effect on the controls) | Effect among untreated patients (rare in oncology) | IPTW with `estimand = "ATC"` |
| **ATO / overlap** | Effect among patients with clinical equipoise (good overlap) | Overlap weighting (`estimand = "ATO"`) |
| **Matched estimand** | Effect in the retained matched population | Depends on matching algorithm; often approximates ATT |

Start with the clinical question and target population, then choose the estimand, then choose the method. **Do not start with the statistical method.**

---

## 4. NCDB-specific considerations

NCDB (National Cancer Database) is a **hospital registry** of Commission on Cancer (CoC)-accredited facilities — not a population-based registry like SEER, and not a randomized trial. NCDB analyses are observational registry analyses requiring extra care.

### NCDB cohort-definition checklist

- **Diagnosis year range** — confirm the user's PUF version supports the years they want.
- **Cancer site / histology / stage** — use ICD-O-3 site/histology codes and the AJCC edition appropriate for the diagnosis years (editions change cutoff thresholds).
- **Endpoint availability** — confirm the endpoint is actually in NCDB before promising it.
- **Survival time** — measured from diagnosis to last contact / death (months).
- **First-course treatment** — most NCDB treatment variables describe the initial treatment plan; salvage / second-line therapy is generally absent.
- **Treatment timing** — surgery date, days from diagnosis to chemo/RT/immunotherapy start; check whether the variables are populated and credible.
- **Zero-month survival** — patients with `dx_lastcontact_death_months == 0` need careful handling (data-entry artifact vs early death); decide and document.
- **Vital status** — exclude unknown vital status from OS analyses.
- **Facility-level variables** — facility type, location, volume; consider clustering by facility ID or robust SEs.
- **Missingness** — substantial for many fields; assess whether it is informative.
- **Generalizability** — NCDB covers CoC-accredited facilities (~70% of incident US cancer cases); do not overgeneralize to the entire US population.

### Endpoints in NCDB — what is and isn't supported

Standard NCDB PUF data **commonly supports**:

- **Overall survival** (vital status + months from diagnosis to last contact / death)
- **30-day or 90-day post-operative mortality** when surgical date and survival time allow
- **Treatment receipt** (which patients got which first-course treatment)
- **Pathologic outcomes** when relevant pathology variables are present (e.g., margin status, pCR for some sites)

Standard NCDB PUF data **generally does NOT support**:

- Disease-free survival (DFS)
- Progression-free survival (PFS)
- Recurrence-free survival (RFS)
- Cancer-specific survival (CSS) — cause of death is not in the standard PUF

If the user asks for DFS / PFS / RFS / CSS in NCDB, **warn explicitly** and pivot to OS unless they confirm they have a non-standard dataset that includes the required variables.

### Immortal-time bias in NCDB

NCDB survival time starts at diagnosis. If the exposure occurs after diagnosis (surgery, adjuvant chemo, immunotherapy, radiotherapy), patients must survive long enough to "receive" the exposure — this creates immortal time bias. Address it before any PS analysis:

- **Landmark analysis** — restrict to patients alive and event-free at the landmark time (e.g., 90 or 180 days post-diagnosis), shift time zero to the landmark, exclude pre-landmark events.
- **Time-dependent treatment** — code treatment as a time-varying covariate (more complex; rarely needed for NCDB sensitivity analyses).
- **Eligibility restriction** — require receipt of treatment by a specified time post-diagnosis.

Apply the landmark BEFORE estimating propensity scores.

### Statistical significance vs clinical relevance in NCDB

NCDB sample sizes are often very large. Trivial absolute differences can be highly statistically significant. Emphasize:

- Standardized mean differences (SMDs) for balance, not p-values
- Absolute risk differences and effect sizes in outcomes, not just HRs and p-values
- Clinical interpretability of the magnitude of effect

---

## 5. Choose only pre-treatment covariates

The propensity-score model includes variables measured **before** treatment / exposure assignment and **before** the index date. This is non-negotiable.

### Do NOT include

- ❌ Post-treatment variables (response, recurrence, ypT/ypN if exposure is neoadjuvant therapy)
- ❌ Mediators of the treatment effect (variables on the causal pathway between treatment and outcome)
- ❌ Colliders (variables caused by both treatment and outcome)
- ❌ Variables measured after follow-up begins
- ❌ Outcome-related variables not available at baseline
- ❌ Treatment-response variables
- ❌ Variables that themselves define eligibility (unless clearly justified)

### Clinically justified baseline covariate set (typical for retrospective oncology)

Adapt to the dataset and question:

- Age (continuous or clinically grouped)
- Sex
- Race / ethnicity
- Insurance status
- Area-level income / education (when available)
- Facility type, facility location, distance to facility (NCDB)
- Year of diagnosis
- Charlson-Deyo comorbidity score
- ECOG performance status (if available; rarely in NCDB)
- Clinical stage (cT, cN, cM) — at diagnosis, not pathologic
- Tumor grade
- Histology
- Tumor location / site / laterality
- Biomarker status (when available and measured at baseline)
- Baseline tumor size
- Treatment intent (curative vs palliative) if recorded
- Prior therapy (when relevant)

For NCDB specifically, do not assume ECOG, recurrence, progression, cause of death, or detailed systemic-therapy regimens are available — they generally are not in the standard PUF.

Use a clinically motivated, pre-specified covariate set. Avoid automated stepwise selection on the PS model.

---

## 6. Check data structure and missingness BEFORE modeling

Before fitting a PS model, run the inspection from `references/07-data-inspection.md` plus a treatment-group breakdown.

Inspect:

- Sample size by treatment group
- Number of events by treatment group
- Median follow-up by group (reverse KM — see `references/02-survival-analysis.md`)
- Missingness per variable, overall and by treatment group
- Variable types and factor levels
- Rare factor levels that may need collapsing
- Complete-case sample size with the planned covariate set
- Whether missingness is likely informative
- Whether "Unknown" is a meaningful category or a missing-data sentinel
- Whether any covariates are structurally missing for an entire treatment group

❌ Do not silently drop hundreds of patients to complete-case analysis without telling the user. Report the analytic-cohort N and the reason for exclusion.

### NCDB variables with frequently meaningful missingness

Pay particular attention to:

- Tumor grade
- Tumor size
- Lymphovascular invasion (LVI)
- Biomarker fields (ER/PR/HER2/EGFR/KRAS/MMR/etc.)
- Treatment timing variables (days to chemo, days to RT)
- Radiation dose
- Surgical margins
- Facility-level variables when restricted to a subset

Consider multiple imputation (`mice`) when missingness is substantial and not plausibly MCAR — but only after discussing with the user and confirming the analytic plan supports it.

---

## 7. Estimate propensity scores or balancing weights — packages

Use the following packages (see `references/09-package-quickref.md` for one-line role summaries):

| Task | Package | Notes |
|---|---|---|
| PS matching | `MatchIt` | Returns a matched dataset; supports nearest-neighbor, optimal, full, exact matching |
| PS weighting (IPTW / ATT / ATE / ATO) | `WeightIt` | Unified interface for many weighting estimands |
| PS weighting + augmented estimators | `PSweight` | Strong overlap-weight support; built-in variance estimation |
| Balance diagnostics | `cobalt` | SMDs, Love plots, balance tables; integrates with `MatchIt`/`WeightIt`/`PSweight` |
| Survival outcome modeling | `survival` | `coxph`, `cox.zph`, `Surv` |
| Weighted Cox / robust SEs | `survey` | `svydesign` + `svycoxph` for IPTW outcome models |
| Manuscript regression tables | `gtsummary` (+ `broom`, `broom.helpers`) | Tidy HR tables with labels |
| Table export | `flextable`, `officer`, `gt` | Word / HTML / multi-table docs |
| Balance SMD table fallback | `tableone` | `CreateTableOne(..., smd = TRUE)` is a quick PS-balance table |
| Multiple imputation | `mice` | When missingness is substantial; combine with PS workflow carefully |

Roles in one line:

- `MatchIt` → matching
- `WeightIt` → weighting
- `cobalt` → balance diagnostics
- `PSweight` → weighting with strong overlap and variance support
- `survival` / `survey` → outcome modeling
- `gtsummary` / `broom` / `gt` / `flextable` / `officer` → manuscript output

---

## 8. Assess overlap, positivity, and extreme weights

Before trusting any weighted or matched estimate, check:

- Propensity-score distributions by treatment group (histograms / density plots)
- Common support / overlap region
- Positivity: are there regions of covariate space with no controls (or no treated)?
- Weight summary statistics (min, median, max, IQR)
- Extreme weights (the upper tail of IPTW weights)
- Effective sample size (ESS) after weighting
- Whether trimming is needed (and the rule used)
- Whether stabilized weights would help
- Whether overlap weights are a better choice
- Whether the estimand should be reframed because of poor overlap

❌ Do not push through with IPTW when many propensity scores are near 0 or 1 — weights become unstable, the ESS collapses, and the point estimate is dominated by a handful of patients.

When overlap is limited and the scientific question can be reframed to the clinical-equipoise population, **overlap weighting (ATO)** is often the better choice. Explicitly report when trimming, overlap weighting, or matching changes the target population — the estimand changes too.

---

## 9. Assess covariate balance BEFORE outcome modeling

Balance assessment is a **gate**. Do not interpret outcome models until balance is acceptable.

Use:

- `cobalt::bal.tab()` — standardized mean differences (SMDs) per covariate, before vs after adjustment
- `cobalt::love.plot()` — visual balance summary
- `cobalt::bal.plot()` — variable-level distribution before vs after
- Effective sample size after weighting / matched N after matching
- Inspection of propensity-score overlap before and after

Typical convention: **absolute SMD < 0.1 is acceptable**, < 0.05 is excellent. These are conventions, not guarantees — clinical importance of imbalance varies by covariate.

❌ Do not use p-values (e.g., t-test, χ²) to judge covariate balance after matching or weighting. P-values conflate balance with sample size — large matched samples will flag trivial imbalances as significant; small samples will miss meaningful ones.

If balance is inadequate:

- Revise the PS model (add interactions, splines, missed covariates)
- Change the estimand (ATE → ATT or ATO)
- Trim non-overlap regions
- Switch to overlap weighting
- Collapse rare factor levels with `forcats::fct_lump()`
- Reconsider whether PS analysis is appropriate for the data

Report balance diagnostics in the manuscript (balance table + Love plot, or both, depending on journal).

---

## 10. Decision tree for reviewer-requested PS analyses

When a reviewer asks for "a PS-adjusted Cox model or IPTW" (the most common phrasing):

1. **Start with diagnostics, not the model.** Inspect sample size by group, event counts, missingness, and PS overlap before choosing the method.
2. **Small sample / few events / poor overlap** → PS-adjusted Cox model (PS as covariate) or overlap weights; explain the choice and limitations.
3. **Good overlap and stable weights** → IPTW is reasonable.
4. **Goal is comparability of treated vs control** → ATT matching or ATT weighting.
5. **Target is the full eligible cohort with adequate positivity** → ATE weighting.
6. **Question is about patients who could plausibly have received either treatment** → overlap weighting (ATO).
7. **Matching discards too many patients or leaves too few events** → switch to weighting or PS-adjusted Cox; explain why matching was not ideal.
8. **IPTW produces extreme weights or low ESS** → stabilized weights, trimming, overlap weights, or ATT weights.
9. **PH assumption violated** → RMST after the PS adjustment / matching / weighting; define τ explicitly.
10. **Exposure defined by receiving treatment after diagnosis or completing N cycles** → assess immortal time bias first. Pragmatic default: landmark analysis (§14.7). When treatment timing is well recorded and landmark misclassification is likely to matter, prefer time-dependent (start-stop) Cox (§11.F, code in §13.H). Apply either BEFORE the PS step.
11. **NCDB analysis requesting DFS / PFS / RFS / CSS** → warn the user that standard NCDB PUF does not support those endpoints; pivot to OS unless they have additional variables.

### Preferred default for retrospective oncology survival papers

1. **Primary analysis**: multivariable Cox with clinically selected covariates (see `references/02-survival-analysis.md`).
2. **Sensitivity analysis**: PS-based approach chosen from the decision tree above.
3. **Report balance diagnostics BEFORE the outcome estimate.**
4. **Do not over-interpret a PS-adjusted HR when balance is poor.**
5. **State explicitly that the analysis addresses measured confounding only.**

Reviewer-requested PS analyses are almost always framed as **sensitivity analyses**, not primary analyses, unless pre-specified in the protocol.

---

## 11. Outcome modeling for survival endpoints

For OS, DFS, PFS, EFS, RFS, etc. (when the endpoint is supported by the data).

### A. PS covariate adjustment

- Estimate the propensity score on baseline covariates (logistic regression by default).
- Fit a Cox model including treatment plus the propensity score as a covariate.
- Optionally include clinically essential covariates that are not fully captured by the PS (e.g., for double adjustment).
- Run `cox.zph()` to check PH.
- Report HR, 95% CI, p-value, N, events.
- Simpler than matching/weighting but more dependent on the PS model being correct.

### B. Matching

- Build the matched dataset with `MatchIt`.
- Report retained N and events; quantify discarded patients and how they differ.
- Assess balance with `cobalt`.
- Fit the Cox model in the matched dataset.
- Use robust SEs or cluster by matched-pair / subclass when appropriate.
- Do not interpret the treatment effect if post-matching balance is poor — revisit the match specification.

### C. IPTW / weighting

- Estimate weights with `WeightIt` (or `PSweight`) for the chosen estimand (ATE / ATT / ATC).
- Assess balance and ESS.
- Inspect the weight distribution; consider stabilization or trimming when extreme.
- Fit a weighted Cox model.
- Use robust / sandwich SEs (`coxph(..., weights = w, robust = TRUE)` or `survey::svycoxph`).
- Report the estimand explicitly.
- Report any stabilization or trimming and how it affects the estimand.

### D. Overlap weighting (ATO)

- Estimate overlap weights with `WeightIt(estimand = "ATO")` or `PSweight`.
- Assess overlap and ESS.
- Report that the estimand applies to the clinical-equipoise / overlap population.
- Fit a weighted Cox model with robust SEs.
- Present as a clinically interpretable alternative when ATE / IPTW positivity is poor.

### E. Non-proportional hazards

If `cox.zph()` flags PH violations, or curves visibly cross:

- Inspect Schoenfeld residual plots.
- Consider RMST as the primary effect measure (see `references/02-survival-analysis.md` section 4).
- Define τ explicitly and justify the choice.
- Do not rely only on a Cox HR when PH is strongly violated.

### F. Time-dependent treatment (start-stop Cox) — alternative to landmark

Landmark analysis is the pragmatic first-line correction for immortal-time bias, but it has a known residual attenuation: patients whose treatment starts **after** the landmark are classified as "untreated by landmark" while still enjoying post-treatment protection during their follow-up. This biases the HR toward the null.

Time-dependent treatment (start-stop Cox / counting-process form) eliminates that misclassification by letting each patient contribute untreated person-time before their treatment date and treated person-time after it. Use it when:

- The exposure is a post-diagnosis event (surgery, chemotherapy, immunotherapy start, RT start).
- The timing of exposure is recorded per patient (`days_to_treatment` or equivalent).
- Landmark misclassification is clinically meaningful (e.g., treatment routinely initiated over a wide window).
- Sample size allows (start-stop expands each treated patient into two rows).

**Do NOT use time-dependent treatment when:**

- Exposure timing is unknown or unreliable in the dataset.
- The exposure is truly baseline (e.g., a diagnostic biomarker, not a post-diagnosis treatment).
- The `treat × time` interaction is the scientific question — that is non-proportional hazards, not time-dependent covariates; use §E instead.

The estimand is the same as under landmark (effect of exposure vs no exposure on the hazard), but the estimator is less attenuated. Report time-dependent Cox as either:
- **Primary analysis** when timing data are clean and clinically important.
- **Sensitivity analysis** alongside landmark, to demonstrate robustness.

Weighted start-stop Cox (PS-based) is straightforward: estimate baseline-covariate PS/weights **once** on the row-per-patient dataset (weights are baseline properties, not time-varying), then attach the weight to every counting-process row for that patient. Code template in §13.H.

---

## 12. Recommended analysis workflow

Follow these steps in order:

1. Clarify clinical question, exposure, comparator, endpoint, time zero, follow-up, eligibility, and estimand (section 3).
2. Confirm endpoint availability — especially for NCDB (section 4).
3. Inspect data: group sizes, event counts, follow-up, missingness (section 6).
4. Define a clinically justified pre-treatment covariate set (section 5).
5. Choose the method:
   - PS-adjusted Cox for a simple reviewer-requested sensitivity analysis
   - Matching when comparability of treated/control is the goal and N allows
   - IPTW for ATE / ATT / ATC when overlap is adequate and weights are stable
   - Overlap weighting when common support is limited and the overlap population is clinically meaningful
   - RMST after PS adjustment when PH is violated or curves cross
6. Estimate the propensity score or weights.
7. Diagnose overlap, positivity, weight distribution, ESS (section 8).
8. Assess balance with SMDs and Love plots (section 9).
9. Modify the design if balance is inadequate (revise PS model, change estimand, trim, switch method).
10. Fit the survival outcome model only after acceptable balance (section 11).
11. Export manuscript-ready tables and figures (section 13G).
12. Clearly document limitations (section 15).

---

## 13. R code templates

Realistic object names: `df_analysis` (analytic cohort), `treat` (0/1 treatment indicator), `ps_fit` (PS model), `ps` (predicted PS), `match_obj` / `df_match` (matched objects), `w_obj` / `wt` (weights), `cox_ps` / `cox_match` / `cox_iptw` / `cox_ow` (outcome models).

### A. Data inspection and missingness

```r
# Group sizes and events
df_analysis |>
  dplyr::group_by(treat) |>
  dplyr::summarise(
    n         = dplyr::n(),
    events    = sum(os_event, na.rm = TRUE),
    censored  = sum(os_event == 0, na.rm = TRUE),
    median_t  = median(os_months, na.rm = TRUE)
  )

# Reverse-KM median follow-up by group (see references/02-survival-analysis.md)
fu_by_treat <- prodlim::prodlim(
  prodlim::Hist(os_months, os_event) ~ treat,
  data = df_analysis, reverse = TRUE
)
quantile(fu_by_treat, q = 0.5)

# Missingness by variable, by group
df_analysis |>
  dplyr::group_by(treat) |>
  dplyr::summarise(dplyr::across(
    dplyr::all_of(c("age", "sex", "ecog", "cstage", "grade", "histology",
                    "comorbidity", "facility_type", "insurance")),
    ~ round(100 * mean(is.na(.x)), 1)
  ))

# Overall structure and per-variable summary
dplyr::glimpse(df_analysis)
skimr::skim(df_analysis)
```

### B. Propensity-score covariate adjustment

```r
# 1. Estimate PS via logistic regression on baseline covariates
ps_fit <- glm(
  treat ~ age + sex + cstage + grade + histology + comorbidity +
          facility_type + insurance + dx_year,
  family = binomial(link = "logit"),
  data   = df_analysis
)

df_analysis$ps <- predict(ps_fit, type = "response")

# 2. Quick overlap check
ggplot2::ggplot(df_analysis, ggplot2::aes(x = ps, fill = factor(treat))) +
  ggplot2::geom_density(alpha = 0.4) +
  ggplot2::labs(x = "Estimated propensity score", fill = "Treatment") +
  ggplot2::theme_classic(base_size = 12)

# 3. Cox model with treatment + PS
cox_ps <- survival::coxph(
  survival::Surv(os_months, os_event) ~ treat + ps,
  data = df_analysis
)

# 4. PH check and tidy output
survival::cox.zph(cox_ps)
cox_ps |>
  gtsummary::tbl_regression(exponentiate = TRUE,
                            pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3))
```

### C. Propensity-score matching with `MatchIt`

```r
library(MatchIt)
library(cobalt)

# 1. 1:1 nearest-neighbor matching on the logit of the PS, with a caliper
match_obj <- MatchIt::matchit(
  treat ~ age + sex + cstage + grade + histology + comorbidity +
          facility_type + insurance + dx_year,
  data     = df_analysis,
  method   = "nearest",
  distance = "glm",
  link     = "logit",
  caliper  = 0.2,            # in SD units of the logit-PS
  ratio    = 1,
  replace  = FALSE
)

summary(match_obj)             # retained N, discarded N, balance summary
df_match <- MatchIt::match.data(match_obj)

# 2. Balance diagnostics
cobalt::bal.tab(match_obj, m.threshold = 0.1, un = TRUE)
love_p <- cobalt::love.plot(
  match_obj,
  thresholds = c(m = 0.1),
  abs        = TRUE,
  stars      = "raw",                            # disambiguates SMD vs raw-diff axis
  var.order  = "unadjusted",
  colors     = c("#1F77B4", "#D62728")
)
love_p

# 3. Retained N and events
df_match |>
  dplyr::group_by(treat) |>
  dplyr::summarise(n = dplyr::n(), events = sum(os_event))

# 4. Cox model on matched data; cluster on matched-pair subclass for SE
cox_match <- survival::coxph(
  survival::Surv(os_months, os_event) ~ treat,
  data    = df_match,
  cluster = subclass,
  robust  = TRUE
)
summary(cox_match)
survival::cox.zph(cox_match)
```

### D. IPTW with `WeightIt`

```r
library(WeightIt)
library(cobalt)

# 1. Estimate weights (ATE here; switch estimand = "ATT" / "ATC" as needed)
w_obj <- WeightIt::weightit(
  treat ~ age + sex + cstage + grade + histology + comorbidity +
          facility_type + insurance + dx_year,
  data     = df_analysis,
  method   = "glm",            # logistic regression PS
  estimand = "ATE",
  stabilize = TRUE
)

summary(w_obj)                 # weight distribution + ESS

# 2. Balance diagnostics
cobalt::bal.tab(w_obj, m.threshold = 0.1, un = TRUE)
cobalt::love.plot(w_obj, thresholds = c(m = 0.1), abs = TRUE, stars = "raw")

# 3. Inspect extreme weights
df_analysis$wt <- w_obj$weights
summary(df_analysis$wt)
quantile(df_analysis$wt, probs = c(0.5, 0.9, 0.95, 0.99, 1))

# Optional: truncate at the 99th percentile if extreme (report the rule)
# w99 <- quantile(df_analysis$wt, 0.99)
# df_analysis$wt_trunc <- pmin(df_analysis$wt, w99)

# 4. Weighted Cox with robust SE (cluster on patient id)
cox_iptw <- survival::coxph(
  survival::Surv(os_months, os_event) ~ treat,
  data    = df_analysis,
  weights = wt,
  robust  = TRUE,
  cluster = patient_id
)
summary(cox_iptw)
survival::cox.zph(cox_iptw)

# Alternative: survey-design weighted Cox
# des <- survey::svydesign(ids = ~ patient_id, weights = ~ wt, data = df_analysis)
# cox_iptw_svy <- survey::svycoxph(survival::Surv(os_months, os_event) ~ treat, design = des)
```

Report the estimand (ATE / ATT / ATC), whether weights were stabilized, and any trimming or truncation rule applied.

### E. Overlap weighting

```r
library(WeightIt)

w_overlap <- WeightIt::weightit(
  treat ~ age + sex + cstage + grade + histology + comorbidity +
          facility_type + insurance + dx_year,
  data     = df_analysis,
  method   = "glm",
  estimand = "ATO"             # overlap (Li, Morgan, Zaslavsky 2018)
)

summary(w_overlap)
cobalt::bal.tab(w_overlap, m.threshold = 0.1)
cobalt::love.plot(w_overlap, thresholds = c(m = 0.1), abs = TRUE, stars = "raw")

df_analysis$wt_ow <- w_overlap$weights

cox_ow <- survival::coxph(
  survival::Surv(os_months, os_event) ~ treat,
  data    = df_analysis,
  weights = wt_ow,
  robust  = TRUE,
  cluster = patient_id
)
summary(cox_ow)
```

State explicitly that the ATO estimand applies to the **clinical-equipoise / overlap population**, not the full eligible cohort. This is a feature, not a limitation, when the question is "for patients who could plausibly have received either treatment, what is the effect?"

`PSweight` is an alternative with strong built-in variance estimation:

```r
# library(PSweight)
# ow_fit <- PSweight::PSweight(
#   ps.formula = treat ~ age + sex + cstage + grade + histology + comorbidity,
#   yname      = "os_event",            # for binary outcomes; survival via separate workflow
#   data       = df_analysis,
#   weight     = "overlap"
# )
# summary(ow_fit)
```

### F. RMST after PS adjustment / matching / weighting

When `cox.zph()` flags PH violations or KM curves visibly cross, present RMST as the primary effect measure or as a sensitivity analysis.

```r
# 1. PH check first
survival::cox.zph(cox_iptw)

# 2. Choose tau (see references/02-survival-analysis.md section 4)
tau_months <- 36   # justify clinically + by risk-set size in each arm

# 3. RMST in the matched dataset (simplest, well-supported)
rmst_match <- survRM2::rmst2(
  time   = df_match$os_months,
  status = df_match$os_event,
  arm    = df_match$treat,         # must be 0/1
  tau    = tau_months
)
print(rmst_match)
```

**Caveat for weighted RMST:** `survRM2::rmst2` does not natively accept weights. For weighted RMST after IPTW or overlap weighting, options include (a) pseudo-observations with weighted GEE, (b) `adjustedCurves` / `RMST` packages that support weighting, or (c) bootstrap of the weighted KM-based RMST. This is non-trivial — flag for statistician review when high-stakes.

### G. Manuscript tables and figures

```r
# G.1 Baseline Table 1 by treatment, BEFORE adjustment (see references/04-baseline-tables.md)
tab1_before <- df_analysis |>
  dplyr::select(age, sex, ecog, cstage, grade, histology, comorbidity,
                facility_type, insurance, dx_year, treat) |>
  gtsummary::tbl_summary(
    by      = treat,
    missing = "ifany",
    statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")
  ) |>
  gtsummary::add_overall() |>
  gtsummary::bold_labels()

# G.2 Balance table after matching/weighting (cobalt → data.frame → flextable)
bal <- cobalt::bal.tab(w_obj, m.threshold = 0.1, un = TRUE)$Balance
bal_tbl <- bal |>
  tibble::rownames_to_column("variable") |>
  dplyr::select(variable,
                smd_unadj = Diff.Un,
                smd_adj   = Diff.Adj)

# G.3 Love plot export (see references/05-manuscript-figures.md — device = "pdf" is
# the portable default; cairo_pdf requires Cairo/XQuartz and fails silently otherwise)
ggplot2::ggsave(
  filename = here::here("output", "fig_loveplot.pdf"),
  plot     = love_p,
  width    = 6.5, height = 5, units = "in",
  device   = "pdf"
)

# G.4 Cox HR table (use the model from the chosen method)
tab_cox <- cox_iptw |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    label = list(treat ~ "Treatment (vs control)"),
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3)
  ) |>
  gtsummary::modify_header(estimate = "**Weighted HR**", ci = "**95% CI**") |>
  gtsummary::bold_labels()

# G.5 KM curves before adjustment
km_before <- survminer::ggsurvplot(
  survival::survfit(survival::Surv(os_months, os_event) ~ treat, data = df_analysis),
  data = df_analysis,
  risk.table = TRUE, conf.int = FALSE, censor = TRUE, pval = TRUE,
  xlab = "Time since diagnosis (months)", ylab = "Overall survival",
  legend.title = "", legend.labs = c("Control", "Treated"),
  ggtheme = ggplot2::theme_classic(base_size = 12)
)

# G.6 KM curves AFTER matching (use df_match)
km_after_match <- survminer::ggsurvplot(
  survival::survfit(survival::Surv(os_months, os_event) ~ treat, data = df_match),
  data = df_match,
  risk.table = TRUE, conf.int = FALSE, censor = TRUE, pval = TRUE,
  xlab = "Time since diagnosis (months)", ylab = "Overall survival (matched cohort)",
  ggtheme = ggplot2::theme_classic(base_size = 12)
)

# G.7 Adjusted survival curve from a weighted Cox (one option for IPTW visualization)
# adjustedCurves::adjustedsurv(...)  # see package vignette; non-trivial

# G.8 Export to Word (see references/06-manuscript-tables.md)
tab_cox |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = here::here("output", "table_cox_iptw.docx"))
```

### H. Time-dependent treatment (start-stop Cox) with PS weights

Alternative to landmark when the exposure is post-baseline and timing is recorded per patient. Removes the residual attenuation that landmark introduces by misclassifying late-treated patients (see §11.F).

**Object conventions.** `df_analysis` is one row per patient with baseline covariates, `os_months`, `os_event`, `days_to_treatment` (NA for never-treated). `wt` is a baseline PS-based weight already estimated on the row-per-patient dataset. Weights are baseline properties; they attach to *every* counting-process row for a given patient.

```r
# 1. Prepare per-patient row-per-patient data with the exposure "start" time
df_pp <- df_analysis |>
  dplyr::mutate(
    t_rx_months = days_to_treatment / 30.4375,   # NA if never treated
    # Ceiling to os_months so patients treated after their last-contact date
    # are handled as never-treated (this can happen with recording lag)
    t_rx_months = ifelse(!is.na(t_rx_months) & t_rx_months < os_months,
                         t_rx_months, NA_real_)
  )

# 2. Split each patient's follow-up on t_rx_months using survSplit()
#    - Never-treated patients get a single row (0, os_months]
#    - Treated patients get two rows: (0, t_rx_months] treated=0 and
#                                     (t_rx_months, os_months] treated=1
df_ss <- df_pp |>
  dplyr::mutate(row_id = dplyr::row_number()) |>
  survival::survSplit(
    formula = survival::Surv(os_months, os_event) ~ .,
    cut     = df_pp$t_rx_months,          # patient-specific cuts (ignore NAs)
    end     = "os_stop",
    start   = "os_start",
    event   = "os_event_ss",
    episode = "episode"
  ) |>
  dplyr::group_by(row_id) |>
  dplyr::mutate(
    treated_now = dplyr::case_when(
      is.na(t_rx_months)               ~ 0L,                   # never treated
      os_start >= t_rx_months          ~ 1L,                   # post-treatment
      TRUE                             ~ 0L                    # pre-treatment
    )
  ) |>
  dplyr::ungroup()

# NB: survSplit's `cut = df_pp$t_rx_months` splits at every unique non-NA cut
# time for every patient — that is NOT what we want. The correct pattern is to
# split each patient at their OWN cut. Do it in a loop or with `tmerge()`:
df_ss <- survival::tmerge(
  data1 = df_pp,
  data2 = df_pp,
  id    = patient_id,
  os_event_ss = event(os_months, os_event),
  treated_now = tdc(t_rx_months)   # time-dependent covariate: 0 before t_rx, 1 after
)

# 3. Fit the weighted start-stop Cox
cox_tdep <- survival::coxph(
  survival::Surv(tstart, tstop, os_event_ss) ~ treated_now,
  data    = df_ss,
  weights = wt,                        # baseline PS weight, replicated per row
  robust  = TRUE,
  cluster = patient_id                 # cluster on patient (multiple rows)
)
summary(cox_tdep)

# 4. PH check — cox.zph works on start-stop objects
survival::cox.zph(cox_tdep)

# 5. Baseline covariate adjustment (double robustness): add PS-selected
#    covariates. They are time-fixed, so `tmerge` does not need to touch them.
cox_tdep_adj <- survival::coxph(
  survival::Surv(tstart, tstop, os_event_ss) ~ treated_now +
    age + sex + cstage + grade + histology + comorbidity +
    facility_type + insurance + dx_year,
  data = df_ss, weights = wt, robust = TRUE, cluster = patient_id
)
summary(cox_tdep_adj)
```

**Doubly-robust interpretation.** In practice, the covariate-adjusted, PS-weighted start-stop Cox (`cox_tdep_adj` above) is the estimator you want. PS weighting alone can leave residual bias if the baseline hazard depends on covariates that are not fully captured by the weight. Combining PS weights with outcome-model covariate adjustment gives a **doubly-robust** estimator: it is consistent if EITHER the PS model OR the outcome model is correctly specified. Empirically this behaves much closer to the true HR than either estimator alone; landmark, by contrast, cannot be made doubly robust in the same way because it discards information at the design stage.

**Guardrails specific to start-stop:**

- `tmerge` is the ergonomic way to build the counting-process dataset. `survSplit` also works, but requires a per-patient loop or careful cut construction. Prefer `tmerge`.
- Always `cluster = patient_id` in the Cox call. Otherwise the sandwich variance treats the two rows of a treated patient as independent, undercounting variance.
- The weight `wt` is a **baseline** property; do not re-estimate it inside each interval. `tmerge` copies baseline columns to every row for a patient, which is what you want.
- After `tmerge`, verify `table(df_ss$treated_now)` gives the number of rows-under-treatment vs rows-untreated, and that the total person-time (`sum(tstop - tstart)`) equals the sum of `os_months` in the original per-patient dataset. Mismatches indicate a split error.
- Do **not** re-weight or re-match on the split dataset — matching / IPTW estimand is defined on the baseline population, not on person-time.
- Report the same estimand language as under landmark (ATE / ATT / ATO) — the estimator changed, the estimand did not.

---

## 14. NCDB-focused code templates and pseudocode

Do not assume specific NCDB variable names. Ask the user for the data dictionary or inspect `names(dat)` first. The code below uses placeholder names — replace with the user's actual variable names.

```r
# 14.1 Filter diagnosis years, site, histology, stage
df_ncdb <- df_raw |>
  janitor::clean_names() |>
  dplyr::filter(
    year_of_diagnosis  %in% 2010:2017,
    primary_site       %in% c("C50.1", "C50.2", "C50.3", "C50.4",
                              "C50.5", "C50.8", "C50.9"),     # example: breast
    histology %in% c(8500, 8501, 8520),                       # ductal/lobular
    analytic_stage_group %in% c(1, 2, 3)                      # AJCC I–III
  )

# 14.2 Define treatment groups from first-course treatment
df_ncdb <- df_ncdb |>
  dplyr::mutate(
    treat = dplyr::case_when(
      rx_summ_surg_prim_site %in% c(20:80) & rx_summ_radiation > 0 ~ 1L,  # surgery + RT
      rx_summ_surg_prim_site %in% c(20:80) & rx_summ_radiation == 0 ~ 0L, # surgery only
      TRUE ~ NA_integer_
    )
  ) |>
  dplyr::filter(!is.na(treat))

# 14.3 Treatment-timing sanity check
df_ncdb |>
  dplyr::group_by(treat) |>
  dplyr::summarise(
    median_days_to_rt = median(dx_rad_started_days, na.rm = TRUE),
    pct_missing       = round(100 * mean(is.na(dx_rad_started_days)), 1)
  )

# 14.4 OS time and event derivation
df_ncdb <- df_ncdb |>
  dplyr::mutate(
    os_months = dx_lastcontact_death_months,
    os_event  = dplyr::case_when(
      puf_vital_status == 0 ~ 1L,    # 0 = dead in standard PUF coding (CONFIRM with codebook)
      puf_vital_status == 1 ~ 0L,    # 1 = alive
      TRUE                  ~ NA_integer_
    )
  )

# 14.5 Exclude unknown vital status and handle zero-month survival
df_ncdb <- df_ncdb |>
  dplyr::filter(!is.na(os_months), !is.na(os_event)) |>
  dplyr::filter(os_months > 0)        # decide and document; alternative: keep + sensitivity

# 14.6 Factor labels and rare-category collapse
df_ncdb <- df_ncdb |>
  dplyr::mutate(
    facility_type = factor(facility_type,
                           levels = c(1, 2, 3, 4),
                           labels = c("Community", "Comprehensive",
                                      "Academic/Research", "Integrated Network")),
    insurance     = forcats::fct_lump_min(factor(insurance), min = 100,
                                          other_level = "Other/Unknown"),
    grade         = forcats::fct_explicit_na(factor(grade), na_level = "Unknown")
  )

# 14.7 Landmark structure for post-diagnosis exposure (immortal time)
# Pragmatic first-line correction. Trades away some information: patients
# treated after the landmark are misclassified as "untreated by landmark",
# which attenuates the HR toward the null. See §11.F / §13.H for time-
# dependent (start-stop) Cox, which does not have this attenuation when
# treatment timing is well recorded.
landmark_days <- 180
df_ncdb_lm <- df_ncdb |>
  dplyr::mutate(days_to_chemo = dx_chemo_started_days) |>
  dplyr::filter(os_months * 30.44 >= landmark_days) |>     # alive at landmark
  dplyr::mutate(
    treat_lm = dplyr::case_when(
      !is.na(days_to_chemo) & days_to_chemo <= landmark_days ~ 1L,
      is.na(days_to_chemo) | days_to_chemo > landmark_days   ~ 0L
    ),
    os_months_lm = os_months - (landmark_days / 30.44)
  )

# 14.8 Cluster on facility for robust SE (illustrative — confirm cluster variable)
# cox_facility <- survival::coxph(
#   survival::Surv(os_months, os_event) ~ treat + age + sex + cstage + grade,
#   data    = df_ncdb,
#   cluster = facility_id,
#   robust  = TRUE
# )
```

When using NCDB, always inspect the data dictionary (PUF codebook for the relevant diagnosis year) and confirm coding of vital status, treatment fields, and survival time **before** writing analysis code. Vital-status coding has changed across PUF versions.

---

## 15. Manuscript language (examples)

Keep wording cautious, oncology-focused, and explicit about limitations.

### Methods — propensity-score model

> "Propensity scores were estimated using logistic regression on baseline pre-treatment covariates selected a priori based on clinical relevance, including age, sex, race/ethnicity, insurance status, Charlson-Deyo comorbidity score, year of diagnosis, clinical stage, tumor grade, histology, and facility type."

### Methods — matching

> "Patients in the [treated] group were matched 1:1 to patients in the [control] group on the logit of the propensity score using nearest-neighbor matching with a caliper of 0.2 standard deviations. Matched pairs were retained only when balance, assessed by standardized mean differences (SMDs), was below 0.1 across all covariates."

### Methods — IPTW

> "Inverse probability of treatment weights were estimated for the average treatment effect (ATE). Weights were [stabilized / truncated at the 99th percentile if needed]. Weighted Cox proportional hazards models for overall survival used robust sandwich variance estimators."

### Methods — overlap weighting

> "Given limited overlap of propensity scores between groups, overlap weights (Li et al. 2018) were used to estimate the average treatment effect in the population with clinical equipoise. This estimand applies to patients who could plausibly have received either treatment, rather than the entire eligible cohort."

### Methods — balance assessment

> "Covariate balance was assessed using absolute standardized mean differences (SMDs), with values below 0.1 considered acceptable. Outcome models were interpreted only after balance diagnostics were reviewed."

### Methods — weighted/matched Cox

> "Hazard ratios and 95% confidence intervals were estimated using [weighted / matched-cohort] Cox proportional hazards models. The proportional hazards assumption was assessed using scaled Schoenfeld residuals."

### Methods — RMST sensitivity

> "When the proportional hazards assumption was violated, restricted mean survival time (RMST) differences were estimated at τ = [X] months as a sensitivity analysis."

### Methods — NCDB cohort

> "Patients were identified from the National Cancer Database (NCDB), a hospital-based registry sponsored by the American College of Surgeons Commission on Cancer and the American Cancer Society, capturing approximately 70% of newly diagnosed cancer cases in the United States from Commission on Cancer–accredited facilities. The analysis used NCDB Participant Use File data for diagnosis years [YYYY–YYYY]."

### Limitations

> "This analysis is subject to the limitations of observational registry data. Although propensity-score methods address measured confounding, residual confounding from unmeasured factors — including performance status, molecular biomarkers not consistently captured, and detailed treatment regimens — cannot be excluded. Because cause of death and recurrence are not captured in the standard NCDB PUF, survival analyses were restricted to overall survival. Generalizability is limited to patients treated at Commission on Cancer–accredited facilities."

---

## 16. Guardrails and warnings

Surface these explicitly to the user when relevant:

- ❌ Propensity-score methods do not fix unmeasured confounding.
- ❌ Poor overlap / positivity violations make IPTW unstable; consider overlap weighting or ATT instead.
- ❌ Post-treatment covariates, mediators, and colliders must not enter the PS model.
- ❌ Immortal time bias must be addressed (landmark, time-dependent treatment, eligibility restriction) BEFORE PS estimation when the exposure is defined by post-diagnosis treatment receipt or duration.
- ❌ Do not interpret outcome models before balance diagnostics are reviewed.
- ❌ Do not use p-values to judge covariate balance after matching or weighting.
- ❌ Small sample sizes or few events make PS methods unstable; consider PS-adjusted Cox or alternative designs.
- ❌ Matching reduces sample size and changes the target population; report N retained and N discarded.
- ❌ Trimming and overlap weights change the estimand — report this explicitly.
- ❌ Reviewer-requested PS analyses are sensitivity analyses unless pre-specified.
- ❌ NCDB is a hospital registry of CoC-accredited facilities, not a population registry, and not a randomized trial.
- ❌ Standard NCDB PUF data generally does not support DFS, PFS, RFS, or cancer-specific survival — confirm endpoint availability before promising.
- ❌ NCDB vital status and follow-up have limitations; verify coding against the codebook for the diagnosis-year version in use.
- ❌ Treatment timing in NCDB can create immortal time bias if not handled.
- ❌ NCDB sample sizes are large — emphasize SMDs, absolute differences, and clinical relevance over p-values.
- ❌ Weighted RMST is non-trivial to implement correctly; flag for statistician review.
- ⚠ Complex or high-stakes PS analyses (regulatory submissions, primary endpoints, treatment-decision-informing) should be reviewed with a statistician — frame this as a partnership, not a deflection.

---

## 17. Style and integration

This reference is intentionally **not a duplicate** of:

- `references/02-survival-analysis.md` — KM, Cox, PH diagnostics, reverse-KM follow-up, RMST mechanics live there
- `references/04-baseline-tables.md` — Table 1 construction and test selection live there
- `references/05-manuscript-figures.md` — KM/forest plot themes, axis labels, export formats live there
- `references/06-manuscript-tables.md` — Word/HTML/CSV export, number formatting, labels live there
- `references/07-data-inspection.md` — first-5-minutes data inspection lives there
- `references/08-code-review-checklist.md` — code-review checklist lives there
- `references/09-package-quickref.md` — one-line package roles live there (and PS packages are listed in its "Propensity-score / causal inference" section)

When a PS analysis touches survival mechanics, Table 1, figure export, or table export, cross-reference the relevant file rather than restating.
