# Baseline (Table 1) and group-comparison tables

The canonical tool for Table 1 in oncology manuscripts is `gtsummary::tbl_summary()`. This file covers default choices, when to add p-values, which tests to use, and how to export.

---

## Minimal Table 1

```r
library(gtsummary)
library(dplyr)

tab1 <- df_analysis |>
  dplyr::select(age, sex, ecog, stage, smoking, arm) |>
  gtsummary::tbl_summary(
    by = arm,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",   # median (IQR) by default
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1)
  ) |>
  gtsummary::add_overall() |>
  gtsummary::add_n() |>
  gtsummary::modify_header(label = "**Characteristic**") |>
  gtsummary::bold_labels()

tab1
```

---

## Continuous variables — median (IQR) vs mean (SD)

**Default to median (IQR).** Most clinical variables (age in elderly cohorts excepted) are skewed, and median + IQR is robust.

Use mean (SD) only when:

- The variable is genuinely normal-looking (visually or by Shapiro on n < 50).
- The reader expects mean (e.g., baseline lab values commonly reported as mean in oncology trials).

Override per-variable:

```r
tbl_summary(
  statistic = list(
    age ~ "{mean} ({sd})",                              # mean (SD) for age
    all_continuous() ~ "{median} ({p25}, {p75})"        # median (IQR) for everything else
  )
)
```

---

## Categorical variables

Default `{n} ({p}%)` — count and percent. `gtsummary` handles missing-as-category if you set `missing = "always"` (forces a row), `missing = "ifany"` (only when present), or `missing = "no"` (drops missing).

Manuscript norm: `missing = "ifany"` with a separate "Missing" row when relevant.

---

## When to add p-values

`gtsummary::add_p()` adds significance tests per row, by-group.

**Add p-values when:**

- The table is a clinical-trial style baseline characteristic table where the convention is to show them (even though it's debated).
- The user is comparing groups defined post hoc (e.g., responders vs non-responders) and the comparison IS the point.

**Skip p-values when:**

- The table is purely descriptive (e.g., a registry cohort with no group comparison intended).
- Group assignment is randomized — by design, p-values are uninformative and journals like NEJM/JAMA no longer require them.
- The reader will misinterpret n.s. p-values as "the groups are similar" (they don't; a baseline imbalance is meaningful regardless of p).

---

## Choosing the right test (when adding p-values)

`add_p()` defaults are often fine, but be explicit when sample size is small or distributions are skewed:

```r
tbl_summary(by = arm) |>
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",           # non-parametric default
      age               ~ "t.test",                # override for one variable
      all_categorical() ~ "fisher.test"            # Fisher when expected counts < 5
    ),
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3)
  )
```

### Decision matrix

| Variable type | Default test | When to override |
|---|---|---|
| Continuous, ≥ 2 groups | Wilcoxon (2 grp) / Kruskal-Wallis (≥3) | t-test / ANOVA only if normality + equal variance hold and large n |
| Continuous, paired | Wilcoxon signed-rank | Paired t-test if normal |
| Categorical 2×2 | Fisher's exact | χ² only if all expected cell counts ≥ 5 |
| Categorical r×c | Fisher's exact (small) / χ² (large) | Same expected-count rule (Cochran) |
| Ordered categorical | Kruskal-Wallis on ranks | Or trend test (Mantel–Haenszel) |

`gtsummary::add_p()` auto-picks based on sample size and cell counts, but spell it out for the manuscript methods section regardless.

---

## Footnote conventions

Most journals expect:

- Continuous: "Median (IQR)" or "Mean ± SD" footnoted.
- Categorical: "n (%)" footnoted.
- Tests used (Fisher exact, Wilcoxon, etc.) footnoted.
- Missing data row labeled.

`gtsummary` handles all of this automatically; verify by printing the table and reading the footer.

---

## Export to Word / HTML

```r
# Word (preferred for journal submission)
library(flextable)
tab1 |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = here::here("output", "table1.docx"))

# HTML (preferred for sharing in markdown or notebooks)
tab1 |>
  gtsummary::as_gt() |>
  gt::gtsave(filename = here::here("output", "table1.html"))

# CSV (for collaborators using Excel/Stata/etc.)
tab1 |>
  gtsummary::as_tibble() |>
  readr::write_csv(here::here("output", "table1.csv"))
```

See `references/06-manuscript-tables.md` for advanced formatting (merging columns, custom footnotes, journal-specific styling).

---

## Stratified Table 1 (3+ groups + overall)

```r
df_analysis |>
  dplyr::select(age, sex, ecog, stage, arm) |>
  tbl_summary(by = arm) |>
  add_overall(last = TRUE) |>
  add_p(test = list(
    all_continuous()  ~ "kruskal.test",
    all_categorical() ~ "fisher.test"
  )) |>
  bold_labels()
```

---

## Anti-patterns

- ❌ Reporting mean ± SD for every variable regardless of distribution.
- ❌ Defaulting to χ² without checking expected cell counts.
- ❌ Hiding the "Missing" row to make the table look clean.
- ❌ Stripping out p-values to "let the data speak" — be explicit either way.
- ❌ Reporting `p = 0.000` — use `p < 0.001`.
- ❌ Showing 4-decimal-place percentages (53.7392%) — round to 1 decimal or whole numbers.
- ❌ Using `tableone::CreateTableOne` for new code — `gtsummary` is the modern standard with better export options.

---

## Outcomes / results tables (not baseline)

For results tables (treatment effects, HRs, ORs):

```r
# Single multivariable Cox model → manuscript-ready table
cox_mv |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3),
    label = list(
      arm     ~ "Treatment arm",
      age     ~ "Age (years)",
      ecog    ~ "ECOG performance status",
      stage   ~ "Disease stage"
    )
  ) |>
  bold_p() |>
  bold_labels()
```

Combine univariable + multivariable side-by-side:

```r
uni  <- tbl_uvregression(df_analysis, method = coxph,
          y = Surv(os_months, os_event),
          include = c(arm, age, sex, ecog, stage),
          exponentiate = TRUE)
mv   <- tbl_regression(cox_mv, exponentiate = TRUE)
tbl_merge(list(uni, mv), tab_spanner = c("**Univariable**", "**Multivariable**"))
```
