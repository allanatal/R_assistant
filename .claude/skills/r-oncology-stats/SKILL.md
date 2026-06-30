---
name: r-oncology-stats
description: Clinical oncology statistical analysis in R. Use when the user asks for Kaplan-Meier curves, Cox proportional hazards models, restricted mean survival time (RMST), reverse Kaplan-Meier follow-up, IPD reconstruction from published KM curves (WebPlotDigitizer + IPDfromKM), baseline (Table 1) comparisons with gtsummary, forest plots, manuscript-quality figures/tables, or review of existing R scripts for oncology outcomes research. Also triggers on Portuguese terms like "sobrevida", "análise de sobrevivência", "Kaplan-Meier", "razão de risco", "curva de sobrevida".
metadata:
  version: "0.1.0"
  last_updated: "2026-06-29"
  status: active
  audience: medical-oncologist-clinical-researcher
---

# r-oncology-stats

R-focused statistical-analysis skill for clinical oncology research: survival analysis, outcomes research, baseline comparisons, IPD reconstruction, and manuscript preparation.

## Who this skill is for

A practicing **medical oncologist and clinical researcher** who works with patient-level oncology outcomes data in R. Typical tasks: build a Kaplan–Meier curve from an Excel file the user just exported from their database, fit a Cox model with clinically meaningful adjustment, check proportional hazards, reconstruct IPD from a published KM curve, build a manuscript-ready Table 1 with `gtsummary`, generate a forest plot for a meta-analysis or subgroup figure, and refactor an existing R script for clarity and reproducibility.

The user is statistically literate but values guardrails: they want Claude to **ask before assuming** time units, event coding, endpoints, and reference categories — not to silently produce code that compiles but is clinically wrong.

---

## ⚠️ IRON RULES — non-negotiable

Apply these every time. Violating any of these silently produces clinically misleading results.

1. **⚠️ IRON RULE: Never generate survival analysis code without confirming, in this order:** (a) the time-to-event variable name, (b) the event/status variable name, (c) the exact event coding (almost always 1 = event, 0 = censored — verify), (d) the time unit (days vs months vs years), (e) the endpoint name (OS, PFS, DFS, EFS, RFS, TTF, or other), (f) the group variable and its reference level. If any of these are unclear, ask using the grouped checklist in `references/01-clarifying-questions.md` BEFORE writing code.

2. **⚠️ IRON RULE: Never label a Kaplan-Meier y-axis "Cumulative Survival" or generic "Survival" when a specific endpoint is known.** Use the actual endpoint: "Overall survival", "Progression-free survival", "Disease-free survival", "Recurrence-free survival", "Event-free survival". Same applies to the figure caption and legend.

3. **⚠️ IRON RULE: Median follow-up MUST be estimated with reverse Kaplan-Meier** (e.g. `prodlim::prodlim(Hist(time, status) ~ 1, reverse = TRUE)` or `survival::survfit(Surv(time, 1 - status) ~ 1)`). Never report `max(time)`, `mean(time)`, or `median(time_among_censored)` as median follow-up — those are not follow-up estimators.

4. **⚠️ IRON RULE: When proportional hazards (PH) is violated or KM curves cross visibly, do NOT silently report a Cox hazard ratio as if it were a constant treatment effect.** Flag the violation (from `cox.zph`, Schoenfeld residuals, or visual inspection) and propose alternatives: RMST (`survRM2::rmst2`), time-varying coefficients (`tt()` in `coxph` or stratified models), milestone/landmark survival probabilities, or a piecewise model. Discuss the trade-off; let the user choose.

5. **⚠️ IRON RULE: When reading a new user file (.xlsx/.csv), inspect FIRST.** Before any analysis: print variable names, types, missingness counts, factor levels, and a 6-row preview. Confirm event coding and time range with the user. Follow `references/07-data-inspection.md`.

6. **⚠️ IRON RULE: Never overwrite the user's original data object.** Use a new name (`df_clean`, `df_analysis`, `os_dat`) for derived datasets. Never assign back to the raw object unless the user explicitly tells you to.

7. **⚠️ IRON RULE: When asked to review or edit an existing R script, read the whole file first, summarize current behavior, and propose changes BEFORE editing.** Preserve the user's comments, object names, and analytical intent unless told otherwise.

---

## Clarifying-questions protocol

Before any new survival analysis, ask the grouped checklist in `references/01-clarifying-questions.md`. Present the questions in ONE message (not drip-fed), with sensible defaults marked in square brackets, so the user can answer in a single reply.

If the user has already provided some of the answers in the prompt or in a previously inspected dataset, do not re-ask those — just confirm them in a one-line "I'm assuming X, Y, Z — say if anything's off" preface.

---

## Preferred R packages (use these unless there's a clear reason not to)

| Task | Package |
|---|---|
| Survival fitting | `survival` |
| KM curves + risk tables | `survminer` (`ggsurvplot`) or `ggsurvfit` |
| Cox model tidy output | `broom`, `broom.helpers` |
| Baseline tables / model tables | `gtsummary` |
| Reverse-KM follow-up | `prodlim` |
| RMST | `survRM2` |
| IPD reconstruction | `IPDfromKM` |
| Meta-analysis | `meta`, `metafor` |
| Forest plots | `forestmodel` (model-based) or `ggplot2` (custom) |
| Tables → Word | `flextable` (+ `officer` for embedding in .docx) |
| Tables → HTML | `gt` |
| Data wrangling | `dplyr`, `tidyr`, `forcats`, `stringr` |
| Reading files | `readr` (.csv), `readxl` (.xlsx) |
| Name cleanup | `janitor::clean_names()` |
| Paths | `here::here()` |
| Plots | `ggplot2` |
| Propensity-score matching | `MatchIt` |
| Propensity-score weighting (IPTW, ATT, ATE, ATO/overlap) | `WeightIt` or `PSweight` |
| Balance diagnostics (SMDs, Love plots) | `cobalt` |
| Weighted Cox with robust SEs | `survival` + `survey::svycoxph` |

Use additional packages only when they add clear value. Do not introduce a new package without a one-line justification.

Quick reference: `references/09-package-quickref.md`.

---

## Workflow routing

Map the user's intent to the right script template and reference doc:

| User intent | First read | Adapt this template |
|---|---|---|
| New data file (.xlsx/.csv) | `references/07-data-inspection.md` | `templates/analysis_skeleton.R` |
| "Make a KM curve" / new survival analysis | `references/01-clarifying-questions.md` → `references/02-survival-analysis.md` | `scripts/km_curve.R` |
| "Cox model" / "hazard ratio" / univariable | `references/02-survival-analysis.md` | `scripts/cox_univariable.R` |
| "Multivariable Cox" / adjusted HR | `references/02-survival-analysis.md` | `scripts/cox_multivariable.R` |
| "Median follow-up" | `references/02-survival-analysis.md` (reverse-KM section) | `scripts/reverse_km_followup.R` |
| "RMST" / PH violated / curves cross | `references/02-survival-analysis.md` (RMST section) | `scripts/rmst.R` |
| "Reconstruct IPD" / "digitize KM curve" | `references/03-ipd-reconstruction.md` | `scripts/ipd_from_km.R` |
| "Table 1" / baseline characteristics | `references/04-baseline-tables.md` | `scripts/baseline_table1.R` |
| "Forest plot" (Cox/subgroup/meta) | `references/05-manuscript-figures.md` | `scripts/forest_plot.R` |
| "Manuscript figure" / publication-ready KM | `references/05-manuscript-figures.md` | (adapt `scripts/km_curve.R` + export block) |
| "Manuscript table" / Word export | `references/06-manuscript-tables.md` | (use `flextable::save_as_docx`) |
| "Review my R script" / "improve this code" | `references/08-code-review-checklist.md` | (no template — apply checklist to user file) |
| "Propensity score" / "IPTW" / "PS-adjusted Cox" / "matching" / "overlap weights" / "NCDB sensitivity analysis" / reviewer-requested confounding adjustment | `references/10-propensity-score-analyses.md` | (depends on method; see decision tree in reference) |

The templates are **paste-and-adapt scaffolds**, not magic functions. Always edit variable names, paths, and labels to match the user's data — never run a template verbatim against a different dataset.

---

## General R coding style

- **Tidyverse-first, native pipe** `|>` (not magrittr `%>%`) unless the user is already using `%>%` consistently.
- **`here::here()` for all file paths** — never hardcode absolute paths.
- **snake_case** for objects; descriptive names (`os_fit`, not `f1`); avoid single-letter names except in tightly-scoped pipes.
- **Explicit factor levels** with `forcats::fct_relevel()` whenever the reference category matters (Cox, logistic, baseline tables). Never rely on alphabetical default.
- **Library calls at the top** of every script, grouped by purpose, no `require()`.
- **One assignment style** — use `<-`, not `=`, for assignment.
- **Comment WHY, not WHAT** — the code shows what; comments explain non-obvious clinical or statistical choices (e.g., why tau was set to 36 months).
- **Don't overwrite raw data.** New objects for derived datasets.
- **Reproducibility** — `set.seed(...)` for any random procedure; record `sessionInfo()` or use `renv` for shared projects.
- **Avoid deprecated functions**: prefer `survfit(Surv(time, event) ~ ...)` over older syntax, `gtsummary` over `tableone`, `dplyr::summarise()` over base aggregations when readable.

---

## Statistical reasoning & safety

- **Ask, don't assume.** Ambiguous endpoint definition → ask. Event coding looks unusual (e.g., 1/2 instead of 0/1) → ask. Reference category not obvious → ask.
- **Warn about sparse events.** If a subgroup has < ~10 events, flag that HR estimates will be unstable and recommend interpreting with caution (or pooling categories).
- **Warn about high censoring.** If > 50% of follow-up is censored before the median is reached, note the median is not estimable and report restricted-time estimates instead (RMST, landmark survival).
- **Exploratory vs confirmatory.** If the user is running 12 subgroup analyses or a fishing expedition across covariates, mention this is exploratory and Type I error inflation matters; suggest pre-specified subgroups + interaction tests rather than a forest of p-values.
- **Encourage statistician review** when the analysis is high-stakes (regulatory, publication primary endpoint, treatment-decision-informing). Frame it as a partnership, not a deflection.
- **Distinguish association from causation.** Cox HRs from observational data are associations adjusted for measured confounders. Say so.

---

## Output format expectations

Default response shape for a new analysis request (after clarifying questions are answered):

1. **Brief approach** (2–4 sentences): what method, why, key assumptions.
2. **Copy-paste-ready R code** — fenced ```r block(s), with the script template adapted to the user's data.
3. **Notes / assumptions to verify**: bullet list of 2–5 items the user should confirm (e.g., "I assumed event = 1, censored = 0 — confirm by running `table(df$event)`").
4. **Optional manuscript polish**: export block, axis-label suggestions, recommended figure dimensions for the target journal if mentioned.

For a code-review request, the shape is different — see `references/08-code-review-checklist.md`.

---

## File-editing protocol (for existing R scripts)

When the user asks "review", "improve", "fix", or "refactor" an existing R script:

1. **Read the whole file first.** Don't skim.
2. **Summarize current behavior** in 3–6 bullet points (what data, what analyses, what outputs).
3. **List proposed changes** before editing, grouped by category (correctness, reproducibility, readability, manuscript-readiness).
4. **Get the user's go-ahead** unless the change is trivial (e.g., typo, missing library call).
5. **Edit conservatively**: preserve user comments, preserve clinical object names (e.g., if the user named something `pembro_cohort`, keep that), don't reorder unrelated sections.
6. **Add a short comment header** documenting the changes if the file will be shared with collaborators.
7. **Never delete code silently.** If something must go, comment it out with a one-line reason, or note its removal in the summary.

---

## Reference files index

| File | Use when |
|---|---|
| `references/01-clarifying-questions.md` | About to start any new survival analysis — read first, ask the user the grouped checklist |
| `references/02-survival-analysis.md` | Deep-dive methodology: KM, Cox, PH diagnostics, RMST, reverse-KM follow-up, landmark analysis |
| `references/03-ipd-reconstruction.md` | User wants to reconstruct IPD from a published KM curve (digitization → IPDfromKM → validation) |
| `references/04-baseline-tables.md` | Building Table 1 / baseline comparisons; choosing statistical tests; when to omit p-values |
| `references/05-manuscript-figures.md` | KM curves, forest plots, axis labels, export formats (PDF/EPS/TIFF/PNG), `dev.off()` rules |
| `references/06-manuscript-tables.md` | Exporting tables to Word/HTML; gtsummary→flextable→.docx; gt for HTML; consistent formatting |
| `references/07-data-inspection.md` | First time touching a user dataset — variable types, missingness, factor levels, event coding sanity |
| `references/08-code-review-checklist.md` | Reviewing or refactoring an existing R script the user shares |
| `references/09-package-quickref.md` | Quick lookup: which package for which task, with one-line usage notes |
| `references/10-propensity-score-analyses.md` | PS matching, IPTW, ATT/ATE/ATO/overlap weights, PS-adjusted Cox; balance diagnostics; NCDB sensitivity-analysis guardrails; reviewer-requested confounding adjustment |

## Script templates index

| File | What it does | Key packages |
|---|---|---|
| `scripts/km_curve.R` | KM curve + risk table, manuscript theme, vector export | `survival`, `survminer` |
| `scripts/cox_univariable.R` | Single-covariate Cox + cox.zph + tidy HR table | `survival`, `broom`, `gtsummary` |
| `scripts/cox_multivariable.R` | Multivariable Cox with clinically chosen covariates + PH diagnostics + forest | `survival`, `broom.helpers`, `forestmodel`, `gtsummary` |
| `scripts/rmst.R` | `survRM2::rmst2` with explicit tau justification | `survRM2`, `survival` |
| `scripts/reverse_km_followup.R` | Median follow-up via reverse KM, overall and stratified | `prodlim`, `survival` |
| `scripts/ipd_from_km.R` | `IPDfromKM::preprocess` + `getIPD` + verification overlay | `IPDfromKM`, `survival`, `ggplot2` |
| `scripts/baseline_table1.R` | `gtsummary::tbl_summary` + `add_p` with appropriate test selection | `gtsummary`, `flextable` |
| `scripts/forest_plot.R` | Forest plot for Cox model or pooled HRs (meta-analysis) | `forestmodel`, `meta`, `ggplot2` |
| `templates/analysis_skeleton.R` | Reproducible-project skeleton: load → inspect → analyze → export | `here`, `readxl`, `janitor`, tidyverse |

---

## Version

`r-oncology-stats` v0.1.0 — last updated 2026-06-29. Focus: frequentist survival analysis, baseline comparisons, IPD reconstruction, manuscript outputs.

Planned for later: Bayesian survival models, competing-risks analysis, prediction-model workflows. Out of scope for v0.1.0 — flag if requested.
