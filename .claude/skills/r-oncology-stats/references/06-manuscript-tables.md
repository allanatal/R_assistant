# Manuscript-quality tables (Word / HTML / CSV export)

When tables move from "analytical output" to "manuscript object", four things matter: consistent number formatting, clean labels, footnotes that mean something, and an export format that survives copy-paste into the final document.

Use `gtsummary` to build, then convert downstream:

- `gtsummary` → `flextable` → `.docx` (journal submission)
- `gtsummary` → `gt` → `.html` (web / notebook)
- `gtsummary` → tibble → `.csv` (collaborator handoff)

---

## Number formatting standards

| Quantity | Format | Example |
|---|---|---|
| p-values | 3 sig figs, `p < 0.001` floor | `0.034`, `< 0.001` |
| HR / OR / RR | 2 decimals | `0.72` |
| 95% CI | 2 decimals, en-dash | `(0.58–0.89)` |
| Median survival | 1 decimal in the units used | `18.4 months` |
| Percent | 1 decimal or whole | `62.4%` / `62%` |
| Count | integer | `127` |

Set globally for `gtsummary`:

```r
gtsummary::theme_gtsummary_journal(journal = "jama")   # or "nejm", "lancet"
# OR custom:
gtsummary::set_gtsummary_theme(
  list(
    "style_number-arg:big.mark" = "",
    "tbl_regression-arg:pvalue_fun" = function(x) gtsummary::style_pvalue(x, digits = 3),
    "tbl_regression-arg:estimate_fun" = function(x) gtsummary::style_number(x, digits = 2)
  )
)
```

`theme_gtsummary_journal()` applies sensible defaults per journal — set it at the top of your analysis script.

---

## Labels — clean, manuscript-ready

Variable names in code are not labels in a manuscript. Set them explicitly:

```r
tbl_summary(
  data = df_analysis,
  by   = arm,
  label = list(
    age      ~ "Age, years",
    sex      ~ "Sex",
    ecog     ~ "ECOG performance status",
    stage    ~ "Disease stage (AJCC v8)",
    smoking  ~ "Smoking history",
    ldh      ~ "Lactate dehydrogenase, U/L"
  )
)
```

Or label once on the dataset (preferred — labels propagate):

```r
df_analysis <- df_analysis |>
  labelled::set_variable_labels(
    age   = "Age, years",
    ecog  = "ECOG performance status",
    stage = "Disease stage (AJCC v8)"
  )
```

Categorical level labels (e.g., `0 = "PS 0"`, `1 = "PS 1"`) should be set on the factor:

```r
df_analysis$ecog <- factor(df_analysis$ecog,
  levels = c(0, 1, 2),
  labels = c("ECOG 0", "ECOG 1", "ECOG 2"))
```

---

## Export to Word (`flextable` → `.docx`)

```r
library(flextable)
library(officer)

doc_path <- here::here("output", "table1.docx")

tab1 |>
  gtsummary::as_flex_table() |>
  flextable::set_table_properties(layout = "autofit", width = 1) |>
  flextable::fontsize(size = 10, part = "all") |>
  flextable::font(fontname = "Times New Roman", part = "all") |>
  flextable::save_as_docx(path = doc_path)
```

For embedding multiple tables in one .docx with text in between, use the `officer` package:

```r
doc <- officer::read_docx() |>
  officer::body_add_par("Table 1. Baseline characteristics by treatment arm.", style = "heading 2") |>
  flextable::body_add_flextable(value = gtsummary::as_flex_table(tab1)) |>
  officer::body_add_break() |>
  officer::body_add_par("Table 2. Cox model for overall survival.", style = "heading 2") |>
  flextable::body_add_flextable(value = gtsummary::as_flex_table(tab_cox))

print(doc, target = here::here("output", "all_tables.docx"))
```

---

## Export to HTML (`gt`)

```r
tab1 |>
  gtsummary::as_gt() |>
  gt::tab_header(
    title    = "Table 1. Baseline characteristics",
    subtitle = "By treatment arm"
  ) |>
  gt::gtsave(filename = here::here("output", "table1.html"))
```

HTML is best for sharing in Quarto/RMarkdown documents and for online supplementary material.

---

## Export to CSV (collaborator handoff)

```r
tab1 |>
  gtsummary::as_tibble() |>
  readr::write_csv(here::here("output", "table1.csv"))
```

CSV loses formatting and footnotes — only use when the recipient needs to manipulate the data in Excel/Stata.

---

## Cox / regression results tables

For an analysis-results table (HRs from a Cox model), use `tbl_regression`:

```r
tab_cox <- cox_mv |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    label = list(
      arm     ~ "Treatment arm",
      age     ~ "Age (per year)",
      ecog    ~ "ECOG (per unit increase)",
      stage   ~ "Stage (vs. I)"
    ),
    pvalue_fun = function(x) gtsummary::style_pvalue(x, digits = 3)
  ) |>
  gtsummary::bold_p(t = 0.05) |>
  gtsummary::bold_labels() |>
  gtsummary::modify_header(estimate = "**HR**", ci = "**95% CI**")

tab_cox |>
  gtsummary::as_flex_table() |>
  flextable::save_as_docx(path = here::here("output", "table2_cox.docx"))
```

---

## Combining tables (univariable + multivariable)

```r
uni <- gtsummary::tbl_uvregression(
  data    = df_analysis,
  method  = survival::coxph,
  y       = survival::Surv(os_months, os_event),
  include = c(arm, age, sex, ecog, stage),
  exponentiate = TRUE
)

mv <- gtsummary::tbl_regression(cox_mv, exponentiate = TRUE)

gtsummary::tbl_merge(
  list(uni, mv),
  tab_spanner = c("**Univariable**", "**Multivariable**")
)
```

---

## Footnotes that mean something

Defaults from `gtsummary` include:

- Statistic explanation ("Median (IQR); n (%)").
- Test explanation ("Wilcoxon rank-sum test; Pearson's χ²").

Add domain-specific footnotes manually:

```r
tab_cox |>
  gtsummary::modify_footnote(
    estimate = "Adjusted hazard ratio from Cox proportional hazards model.",
    ci = "95% Wald confidence interval."
  )
```

Footnote anti-patterns:

- Footnoting *every* cell — clutters the table.
- Footnoting obvious things (e.g., "% = percent").
- Footnoting statistical tests that aren't in the table.

---

## Anti-patterns specific to manuscript tables

- ❌ Different number-of-decimals across rows of the same column.
- ❌ Mixing comma and period as thousands separator (be consistent: `12,000` OR `12.000`, follow the journal).
- ❌ Reporting `p = 0.0000` — use `p < 0.001`.
- ❌ Bolding everything to make the table look important.
- ❌ Including 4-level categorical with single observations per cell (collapse first or note sparsity).
- ❌ Forgetting to convert variable names to clean labels before export — `os_months` instead of `OS (months)` looks unfinished.
- ❌ Saving a .docx with no caption — captions live in the figure/table itself or in the manuscript text, but Table 1 / Table 2 should be labeled.

---

## Quick-reference table

| Goal | Build with | Convert with | Export to |
|---|---|---|---|
| Baseline / Table 1 | `tbl_summary` | `as_flex_table` | `save_as_docx` |
| Cox / regression table | `tbl_regression` | `as_flex_table` | `save_as_docx` |
| Univariable screen | `tbl_uvregression` | `as_flex_table` | `save_as_docx` |
| Multiple tables in one doc | `tbl_*` × N | `as_flex_table` | `officer::body_add_flextable` |
| HTML report | `tbl_*` | `as_gt` | `gtsave(... .html)` |
| CSV handoff | `tbl_*` | `as_tibble` | `write_csv` |
