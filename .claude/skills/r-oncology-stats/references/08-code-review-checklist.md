# Code review and refactoring checklist for existing R scripts

When the user asks you to "review", "improve", "fix", or "refactor" an R script, follow this protocol.

---

## Protocol

1. **Read the WHOLE file first** — never skim. If it's long, read in chunks and note structure.
2. **Summarize current behavior** in 4–8 bullets: data sources, analyses performed, outputs produced.
3. **Identify issues** by category (see checklist below).
4. **Propose changes** before editing. Get the user's go-ahead unless changes are trivial (typos, missing libraries).
5. **Edit conservatively**: preserve comments, object names, and analytical intent.
6. **Document changes** with a short comment header if the script will be shared.

---

## The review categories

### A. Correctness (highest priority)

- **Event coding**: is `Surv(time, event)` getting 1=event/0=censored, or something else? Check `table(data$event)` against script assumptions.
- **Reference level**: are factor levels in a clinically meaningful order, or alphabetical default? Look for `factor()`, `forcats::fct_relevel`, or absence thereof.
- **Time unit**: are axis labels, median survival reports, and tau values all in the same unit? Mixed units silently break interpretation.
- **PH assumption**: is `cox.zph()` run after every Cox model? If not, flag it.
- **Median follow-up**: is it computed via reverse KM or via `median(time)`? The latter is wrong.
- **Censoring vs competing events**: does the analysis treat non-cancer death as censoring when competing risks should apply? (Often a flag, not a fix.)
- **Multiple testing**: if there are subgroup analyses or many univariable screens, does the script adjust or note exploratory status?
- **Missing-value handling**: is `na.action` set? Is `na.omit` silently dropping rows? Are missing categorical levels collapsed to "Unknown" silently?
- **Reference-data leak**: are diagnostic plots or summaries computed on the modeling data (vs a validation subset, when intended)?

### B. Reproducibility

- **`set.seed()`** before any random procedure (bootstrap CI, multiple imputation, random forest splits).
- **`here::here()`** for file paths, not absolute paths like `"/Users/.../data.xlsx"`.
- **Library calls at top** — not scattered through the script.
- **`sessionInfo()`** or `renv` lockfile referenced if the script will be re-run later.
- **Hard-coded values** (cutoffs, thresholds, dates) — should be named constants at the top.
- **Output paths** — should write to a clear `output/` directory using `here::here("output", "...")`.

### C. Readability / maintainability

- **Variable names**: `f1`, `mod`, `m` → unclear; `os_cox_mv`, `pfs_km_fit` → clear.
- **Magic numbers**: `tau <- 36` is fine if commented (`# 3-year RMST per protocol`); unexplained `0.7` deep in the script is suspect.
- **Comment hygiene**: comments should explain WHY (clinical or statistical rationale), not WHAT (the code shows what).
- **Repeated blocks**: if the same analysis runs for OS and PFS with copy-pasted code, propose a small function.
- **Long script**: > 300–400 lines suggests modular organization (one file per analytical step, sourced from a master script).
- **Mixed pipe styles**: scattered `%>%` and `|>` — unify on one (prefer `|>` for new code).

### D. Manuscript-readiness

- **KM y-axis label**: is it specific (e.g., "Overall survival") or generic ("Cumulative survival")?
- **Risk table**: present on every KM figure?
- **Number formatting**: HRs to 2 decimals, p-values 3 sig figs with `< 0.001` floor, CIs with en-dash?
- **Table labels**: variable names cleaned for the manuscript (`os_months` → "OS, months")?
- **Export format**: PDF/EPS for vectors, TIFF 600 DPI when required, `dev.off()` after base devices?
- **Footnotes**: tests and stats explained on Table 1?

### E. Outdated / deprecated patterns

- **`tableone`** → migrate to `gtsummary`.
- **Base R `survfit` formula with `~ as.factor(arm)`** when `arm` is already a factor → unnecessary.
- **`as.numeric(as.character(factor))`** → use `as.numeric(levels(factor))[factor]` or stay as factor.
- **`dplyr::summarise_at` / `_if` / `_all`** → migrate to `across()`.
- **`reshape2::melt` / `dcast`** → migrate to `tidyr::pivot_longer` / `pivot_wider`.
- **`plyr::ddply`** → migrate to `dplyr::group_by + summarise`.
- **`ggsurv` from old GGally** → migrate to `survminer::ggsurvplot` or `ggsurvfit`.

---

## Summary report template

When you report findings, structure them like this:

```
## Review of `analysis.R`

### What the script does
- Reads `data/trial.xlsx`, sheet "patients".
- Cleans names, recodes ECOG to 3 levels.
- Fits univariable + multivariable Cox for OS and PFS.
- Generates Table 1 (gtsummary) and KM curves (survminer).
- Exports tables to .docx, figures to .pdf.

### Issues found (in priority order)

**Correctness**
1. Event coding not verified: `Surv(os_months, status)` assumes status = 1/0 but
   `table(df$status)` shows 1/2. → Recode before fitting.
2. No `cox.zph()` after Cox models. → Add PH check; flag any violations.
3. Median follow-up computed as `median(df$os_months)`. → Replace with reverse KM
   (`prodlim::prodlim(Hist(os_months, status) ~ 1, reverse = TRUE)`).

**Reproducibility**
4. Hard-coded path `/Users/allan/data/trial.xlsx`. → Use `here::here("data", "trial.xlsx")`.
5. No `set.seed()` before bootstrap CI (line 142).

**Readability**
6. Variables `m1`, `m2`, `m3` for three Cox models — rename to `cox_os_uni`,
   `cox_os_mv`, `cox_pfs_mv`.
7. Same KM-plot block copy-pasted four times — extract to a helper function.

**Manuscript-ready**
8. KM y-axis is "Survival" — should be "Overall survival" (or specific endpoint).
9. p-values formatted as `0.000` in the output — switch to `style_pvalue(digits = 3)`.

### Suggested next step
I can apply fixes 1–5 (correctness + reproducibility) now if you want, or you can
review and tell me which to change first. Fixes 6–9 (readability + manuscript polish)
are less urgent — happy to do them in a second pass.
```

---

## When NOT to "improve"

- The user wrote the script in a particular style intentionally — don't impose tidyverse on a base-R purist.
- The script works and is clear — don't refactor for refactoring's sake.
- The script uses an older package for a specific reason (e.g., reviewer requested it) — flag the modern alternative as a comment but don't switch.
- The user is co-authoring with someone who only uses Stata-style — keep the script approachable.

Ask if you're unsure about scope. "I see these patterns I'd usually update — want me to leave them as-is?" is better than over-editing.

---

## Edit-mode rules

- Edit one logical chunk at a time, not the whole file at once.
- Show the diff (or the proposed `Edit` block) before applying to anything non-trivial.
- Preserve original comments unless they're factually wrong; even then, update rather than delete.
- Never delete code silently — comment out with a one-liner reason, or note the removal in your summary.
- After editing, run `lintr::lint()` or just re-read the file to confirm nothing's broken.

---

## Anti-patterns in code review

- ❌ Saying "looks good, just modernize the syntax" without checking event coding.
- ❌ Rewriting the script in your preferred style without asking.
- ❌ Reordering function definitions for "logical flow" — breaks user's mental map.
- ❌ Removing seemingly-dead code without confirming it's truly unused.
- ❌ Auto-running and applying `styler::style_file()` without warning — overwrites user formatting choices.
- ❌ Pasting back the entire script with one tiny change — show only what changed.
