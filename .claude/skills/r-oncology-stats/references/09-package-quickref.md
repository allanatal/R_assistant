# R package quick reference

One-line reminders for the preferred packages, grouped by task. When in doubt, use these. Add others only with a clear justification.

---

## Survival core

| Package | Use for | Key functions |
|---|---|---|
| **survival** | KM fitting, Cox models, PH diagnostics; the foundation. | `Surv()`, `survfit()`, `coxph()`, `cox.zph()`, `survdiff()` |
| **survminer** | Manuscript-ready KM curves with risk tables. | `ggsurvplot()`, `ggcoxzph()`, `surv_pvalue()` |
| **ggsurvfit** | Modern, ggplot-composable survival plots; alternative to survminer. | `ggsurvfit()`, `add_risktable()`, `add_pvalue()` |
| **prodlim** | Reverse-KM follow-up (median follow-up time), product-limit estimation. | `prodlim()` (with `reverse = TRUE`), `Hist()` |
| **survRM2** | Restricted Mean Survival Time. | `rmst2()` |

## IPD reconstruction

| Package | Use for | Key functions |
|---|---|---|
| **IPDfromKM** | Reconstruct pseudo-IPD from digitized KM curves (Guyot algorithm). | `preprocess()`, `getIPD()` |

## Meta-analysis

| Package | Use for | Key functions |
|---|---|---|
| **meta** | General-purpose meta-analysis (HRs, ORs, means). | `metagen()`, `metabin()`, `forest()` |
| **metafor** | Advanced meta-analysis: meta-regression, multi-level, network. | `rma()`, `rma.mv()`, `forest()` |

## Tables

| Package | Use for | Key functions |
|---|---|---|
| **gtsummary** | Build Table 1, regression tables, table merges. | `tbl_summary()`, `tbl_regression()`, `tbl_uvregression()`, `add_p()`, `tbl_merge()` |
| **flextable** | Export tables to Word (.docx) and PowerPoint. | `save_as_docx()`, `set_table_properties()` |
| **gt** | Export tables to HTML, also LaTeX/RTF. | `gt()`, `gtsave()`, `tab_header()` |
| **officer** | Build multi-table .docx with text between tables. | `read_docx()`, `body_add_par()`, `body_add_flextable()` |

## Tidy model output

| Package | Use for | Key functions |
|---|---|---|
| **broom** | Convert model objects to tidy data frames. | `tidy()`, `glance()`, `augment()` |
| **broom.helpers** | Tidy bridge from `coxph` / `glm` to a manuscript-ready data frame with `reference_row` / `header_row` flags. Foundation for custom forest plots. | `tidy_plus_plus(model, exponentiate = TRUE, conf.int = TRUE, add_reference_rows = TRUE, add_header_rows = TRUE)` |
| **forestmodel** | One-call forest plots from a `coxph` / `lm` / `glm` model. Use as QC sanity check next to the custom oncology layout; see `references/11-forest-plots.md` §6. | `forest_model()` |

## Propensity-score / causal inference

For full guidance on when and how to use these (estimand selection, balance assessment, NCDB-specific guardrails), see `references/10-propensity-score-analyses.md`.

| Package | Use for | Key functions |
|---|---|---|
| **MatchIt** | Propensity-score matching; returns a matched dataset. | `matchit()`, `match.data()`, `summary()` |
| **WeightIt** | Balancing weights: IPTW (ATE/ATT/ATC), overlap (ATO), multi-category and continuous treatments. | `weightit()`, `summary()`, `trim()` |
| **PSweight** | Weighting (incl. overlap), augmented estimators, variance estimation, diagnostics. | `SumStat()`, `PSweight()`, `summary()` |
| **cobalt** | Balance diagnostics: standardized mean differences, Love plots, balance tables. | `bal.tab()`, `love.plot()`, `bal.plot()` |
| **survey** | Weighted Cox / robust-SE workflows for IPTW outcome models. | `svydesign()`, `svycoxph()` |
| **tableone** | Quick balance-table SMD display (alternative to `cobalt::bal.tab` for a printable summary). | `CreateTableOne()`, `print(..., smd = TRUE)` |

Preferred manuscript path: `gtsummary` + `cobalt::bal.tab()` for the balance table; `tableone` is a fallback when an SMD-focused printed table is needed.

## Data wrangling

| Package | Use for | Key functions |
|---|---|---|
| **dplyr** | Filter, mutate, group, summarise; verb-based wrangling. | `filter()`, `mutate()`, `group_by()`, `summarise()`, `across()` |
| **tidyr** | Reshape (long ↔ wide), separate, unite, nest. | `pivot_longer()`, `pivot_wider()`, `drop_na()`, `nest()` |
| **stringr** | Consistent string manipulation. | `str_detect()`, `str_replace()`, `str_extract()` |
| **forcats** | Factor manipulation: reorder levels, lump rare categories. | `fct_relevel()`, `fct_recode()`, `fct_lump()`, `fct_infreq()` |
| **janitor** | Name cleanup, tabulation, Excel date handling. | `clean_names()`, `tabyl()`, `excel_numeric_to_date()` |

## Reading files

| Package | Use for | Key functions |
|---|---|---|
| **readr** | Read CSV, TSV, fixed-width; fast and consistent. | `read_csv()`, `read_csv2()` (Euro), `read_tsv()` |
| **readxl** | Read Excel `.xls` / `.xlsx`. | `read_excel()`, `excel_sheets()` |

## Project structure

| Package | Use for | Key functions |
|---|---|---|
| **here** | Sane file paths relative to the project root. | `here()` |

## Plotting

| Package | Use for | Key functions |
|---|---|---|
| **ggplot2** | All custom plots; foundation for survminer, ggsurvfit, forestmodel. | `ggplot()`, `geom_*()`, `theme_classic()`, `ggsave()` |
| **patchwork** | Combine multiple ggplots into multi-panel layouts (e.g., the two-panel forest in `references/11-forest-plots.md`). | `plot_layout()`, `wrap_plots()`, `+` / `/` operators |
| **scales** | Axis transforms, log-tick break helpers, precision-based size rescaling for forest squares. | `rescale()`, `log_breaks()`, `pretty_breaks()` |

---

## When to reach beyond this list

These are not in the default skill list — use them when the situation calls for it, and note the reason:

| Package | Use case |
|---|---|
| **tidycmprsk** | Competing-risks analysis (Fine-Gray, cumulative incidence). Out of v0.1.0 scope but worth knowing. |
| **cmprsk** | Lower-level competing-risks (older API; tidycmprsk wraps it). |
| **rms** | Frank Harrell's regression modeling strategies — flexible parametric, restricted cubic splines, validation. |
| **mice** | Multiple imputation for missing data — including for PS-model covariates when missingness is substantial. |
| **boot** | Bootstrap CIs when you need them (most survival tools have built-in CIs). |
| **labelled** | Variable labels that propagate through gtsummary. |
| **skimr** | First-look data summary; complements `glimpse`. |
| **naniar** | Missingness visualization. |
| **lubridate** | Date/time parsing and arithmetic. |
| **renv** | Per-project package versions for reproducibility. |

---

## Anti-list — avoid in new code

| Old | Use instead |
|---|---|
| `tableone::CreateTableOne` | `gtsummary::tbl_summary` |
| `reshape2::melt` / `dcast` | `tidyr::pivot_longer` / `pivot_wider` |
| `plyr::ddply` | `dplyr::group_by + summarise` |
| `magrittr::%>%` (in new scripts) | Native pipe `|>` |
| `summarise_at/_if/_all` | `summarise(across(...))` |
| `aggregate()` (base R, for analysis) | `dplyr::summarise` |
| `cbind.data.frame()` chains | `dplyr::bind_rows` / `bind_cols` |
| `ifelse()` for vectorized type-changing logic | `dplyr::case_when` or `dplyr::if_else` (type-stable) |
