# Forest plots — Cox multivariable and subgroup

Manuscript-grade forest plots for oncology: precision-scaled two-panel `ggplot2` + `patchwork` for Cox MV coefficients; subgroup forests across many small Cox models; explicit handling of the "HR = 1.00 (1.00–1.00)" non-estimable trap; `forestmodel` as a QC sanity check, not the manuscript figure.

> **When to use this doc.** Trigger phrases: "forest plot", "subgroup forest", "HR figure", "make a figure showing the adjusted HRs", reviewer asks for a figure that combines the Cox MV table with a visual, or any time you need a publication-ready forest with non-default styling. For the survival mechanics behind the Cox model, see `references/02-survival-analysis.md`. For general export conventions (PDF/EPS/TIFF/PNG, DPI, `dev.off()` rules), see `references/05-manuscript-figures.md`. For the companion `.docx` regression table, see `references/06-manuscript-tables.md`.

---

## 1. Decision tree — which forest are you building?

Answer this first. The three layouts are NOT interchangeable.

```
        ┌─────────────────────────────────────────────────────────┐
        │ Are you plotting...                                     │
        └─────────────────────────────────────────────────────────┘
                              │
        ┌─────────────────────┼──────────────────────┐
        ▼                     ▼                      ▼
  (a) ALL coefficients   (b) ONE treatment HR    (c) Pool HRs
      of ONE Cox MV          per subgroup            across published
      model                  across MANY small       studies
                             Cox fits
        │                     │                      │
        ▼                     ▼                      ▼
   §3 oncology layout    §5 subgroup forest    Out of scope —
   (precision-scaled,    (purrr::map_dfr        see `05-manuscript-
    two-panel ggplot)     builder + sibling      figures.md` →
                          plot fn)               Meta-analysis forest
```

If you are unsure which one you need, ask the user. A common confusion: a "subgroup analysis" in an oncology paper can mean either (a) — adjusted HRs of every covariate including a subgroup variable, in one MV model — or (b) — the treatment HR re-estimated within each subgroup level. They look superficially similar but require different code paths and different statistical reasoning (see §5.3 on interaction p-values).

---

## 2. Package stack

| Package | Role here | Key calls |
|---|---|---|
| **survival** | Fit the Cox model(s). | `coxph()`, `Surv()`, `cox.zph()` |
| **broom.helpers** | Tidy bridge from `coxph` to a manuscript-ready data frame with `reference_row` / `header_row` flags. **Without this you cannot do a custom layout.** | `tidy_plus_plus(model, exponentiate = TRUE, conf.int = TRUE, add_reference_rows = TRUE, add_header_rows = TRUE)` |
| **ggplot2** | The two ggplot panels (text + forest). | `ggplot()`, `geom_text()`, `geom_point()`, `geom_errorbar(orientation = "y")`, `geom_vline()`, `scale_x_log10()` |
| **patchwork** | Combine the text and forest panels with controlled widths. | `plot_layout(widths = c(1.4, 1))`, `+` |
| **scales** | Rescale `1 / SE` to a fixed visual square-size range; adaptive log breaks. | `rescale()`, `log_breaks()` |
| **gtsummary** | Companion `.docx` regression table; `modify_table_body()` for masking non-estimable rows. | `tbl_regression()`, `modify_table_body()` |
| **flextable** | Export the companion table to `.docx`. | `as_flex_table()`, `save_as_docx()` |
| **forestmodel** | QC sanity-check forest. **Not the manuscript figure.** See §6. | `forest_model()` |
| **forcats** | Set the reference level BEFORE fitting; collapse rare levels. | `fct_relevel()`, `fct_collapse()` |
| **purrr** | Loop over subgroups in §5. | `map_dfr()` |
| **here** | Resolve output paths relative to project root. | `here()` |

Libraries to load at the top of any script that uses the helpers in this doc:

```r
library(survival)
library(broom.helpers)
library(ggplot2)
library(patchwork)
library(scales)
library(gtsummary)
library(flextable)
library(forestmodel)
library(forcats)
library(purrr)
library(dplyr)
library(here)
```

---

## 3. The oncology layout — precision-scaled two-panel forest

This is the manuscript figure. It is the answer to "make me a forest plot from my multivariable Cox model" for any oncology paper.

### 3.1 Why precision-scaled squares (not uniform)

Each square's area is proportional to `1 / SE(coef)`, rescaled to a fixed visual range (default `c(2.5, 7)`) via `scales::rescale()`. This is the convention readers expect from Lancet / NEJM / JCO forest plots: the eye is drawn to precise estimates, not to whichever row happens to have the widest CI whisker. `forestmodel::forest_model()` uses uniform sizes — fine for QC, wrong for the figure.

### 3.2 Anatomy of the figure

```
┌──────────────────────────────────────┬──────────────────────────────┐
│ LEFT PANEL (text, no axes)           │ RIGHT PANEL (forest, axes)   │
│                                      │                              │
│ Histology                            │                              │
│   Atypical (G2)    —                 │              (no point)      │
│   Malignant (G3)   2.31 (1.42–3.74), │            ▣────────────     │
│                    p=0.001           │                              │
│                                      │                              │
│ Treatment                            │                              │
│   Surgery          —                 │              (no point)      │
│   Surgery + RT     0.62 (0.41–0.93), │        ──────▣──             │
│                    p=0.020           │                              │
│   CT alone         Not estimable     │              (no point)      │
│                                      │  ─────────┊─────────         │
│                                      │  0.1   0.5 1  2    10        │
│                                      │  HR (95% CI, log scale)      │
└──────────────────────────────────────┴──────────────────────────────┘
   widths = c(1.4, 1)        patchwork::plot_layout()
```

Key visual rules:

- Left panel: **bold** variable header rows (e.g., "Histology"), indented level labels, monospaced "HR (95% CI, p)" text column. Reference rows show `—`; degenerate rows show `Not estimable`.
- Right panel: dashed vertical line at HR = 1, log-scale x-axis with adaptive breaks (drops breaks outside the data range), horizontal error bars via `geom_errorbar(orientation = "y", width = 0)` (no whisker caps — they add visual noise without information).
- Composition: `patchwork::plot_layout(widths = c(1.4, 1))` — the text panel gets ~58% of the width.

### 3.3 Reference rows and degenerate rows in the visual

- **Reference level** (the baseline category for each factor): rendered as `—` in the text column, no point on the right panel. This is the row the other HRs of that variable are relative to.
- **Non-estimable row** (sparse stratum, n=1 cell, perfect separation): rendered as `Not estimable` in the text column, no point and no whisker on the right panel. **Do not** show `HR = 1.00 (1.00–1.00)` for these — that fake null is the bug §4 is about.

### 3.4 Code — `plot_forest_oncology()` (paste-and-adapt)

Drop this function into your analysis script (after the library block in §2). It takes a fitted `coxph` model and writes `<out_dir>/<file_base>.pdf` and `<out_dir>/<file_base>.png` (300 dpi). Returns the patchwork object invisibly.

```r
#' Precision-scaled two-panel oncology forest plot for a Cox MV model.
#'
#' @param model           A fitted `coxph` object.
#' @param file_base       File stem; output is "<out_dir>/<file_base>.pdf" + ".png".
#' @param endpoint_label  Title prefix, e.g. "Overall survival" or "LRFS".
#' @param ci_label        CI label string. Default "95% CI".
#' @param box_fill        Square fill color (hex). Default "#1f3b6f" (publication navy).
#' @param ci_extreme_hi   Upper CI bound above which a row is treated as non-estimable.
#' @param ci_extreme_lo   Lower CI bound below which a row is treated as non-estimable.
#' @param x_breaks        Candidate x-axis breaks; the function keeps only those in range.
#' @param w               Output width (inches). Default 11.
#' @param h               Output height (inches); NULL = auto from row count.
#' @param out_dir         Directory for PDF + PNG. Default here::here("output").
plot_forest_oncology <- function(model,
                                 file_base,
                                 endpoint_label  = "HR",
                                 ci_label        = "95% CI",
                                 box_fill        = "#1f3b6f",
                                 ci_extreme_hi   = 1e3,
                                 ci_extreme_lo   = 1e-4,
                                 x_breaks        = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 25),
                                 w               = 11,
                                 h               = NULL,
                                 out_dir         = here::here("output")) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  td <- broom.helpers::tidy_plus_plus(
    model,
    exponentiate       = TRUE,
    conf.int           = TRUE,
    add_reference_rows = TRUE,
    add_header_rows    = TRUE
  )

  td <- td |>
    dplyr::mutate(
      is_header = !is.na(header_row) & header_row,
      is_ref    = !is.na(reference_row) & reference_row,
      ## Continuous variables get no explicit header row from tidy_plus_plus
      ## -- their single row is BOTH the header and the estimate. Treat as
      ## a standalone entry so `age` renders at the group-header column
      ## (bold, flush-left) rather than indented under the previous factor.
      is_standalone = !is_header & !is_ref &
                      !is.na(var_type) & var_type == "continuous",
      group_lab = dplyr::case_when(is_header      ~ var_label,
                                   is_standalone  ~ var_label,
                                   TRUE           ~ ""),
      level_lab = dplyr::case_when(is_header      ~ "",
                                   is_standalone  ~ "",
                                   TRUE           ~ label),
      plottable = !is_header & !is_ref &
                  is.finite(estimate) & is.finite(conf.low) & is.finite(conf.high) &
                  conf.high < ci_extreme_hi & conf.low > ci_extreme_lo &
                  !is.na(std.error) & std.error > 0 &
                  !is.na(p.value),
      plot_est = ifelse(plottable, estimate, NA_real_),
      plot_lo  = ifelse(plottable, conf.low, NA_real_),
      plot_hi  = ifelse(plottable, conf.high, NA_real_),
      precision = ifelse(plottable & is.finite(std.error) & std.error > 0,
                         1 / std.error, NA_real_),
      hr_text = dplyr::case_when(
        is_header                     ~ "",
        is_ref                        ~ "-",  ## ASCII hyphen renders in any font; swap for "—" if using cairo_pdf
        !plottable                    ~ "Not estimable",
        TRUE ~ sprintf(
          "%.2f (%.2f-%.2f), p=%s",
          estimate, conf.low, conf.high,
          format.pval(p.value, digits = 2, eps = 0.001)
        )
      )
    )

  if (any(!is.na(td$precision))) {
    td$box_size <- scales::rescale(
      td$precision, to = c(2.5, 7),
      from = range(td$precision, na.rm = TRUE)
    )
  } else {
    td$box_size <- 3
  }
  td$box_size[td$is_header | td$is_ref | !td$plottable] <- NA

  td$row_id <- seq_len(nrow(td))
  td$y      <- max(td$row_id) - td$row_id + 1

  if (is.null(h)) h <- max(4.5, 0.32 * nrow(td) + 1.5)

  x_finite <- c(td$plot_lo, td$plot_hi)
  x_finite <- x_finite[is.finite(x_finite)]
  if (length(x_finite) == 0) x_finite <- c(0.5, 2)
  x_breaks_used <- x_breaks[x_breaks >= min(x_finite, 1) * 0.5 &
                            x_breaks <= max(x_finite, 1) * 2]
  if (length(x_breaks_used) < 3) x_breaks_used <- x_breaks
  x_lim <- range(c(x_breaks_used, x_finite, 1), na.rm = TRUE)

  p_text <- ggplot2::ggplot(td, ggplot2::aes(y = y)) +
    ggplot2::geom_text(ggplot2::aes(x = 0,    label = group_lab),
                       hjust = 0, fontface = "bold", size = 3.6) +
    ggplot2::geom_text(ggplot2::aes(x = 0.18, label = level_lab),
                       hjust = 0, size = 3.3) +
    ggplot2::geom_text(ggplot2::aes(x = 1.05, label = hr_text),
                       hjust = 0, size = 3.3, family = "mono") +
    ggplot2::scale_x_continuous(limits = c(0, 2.4), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, max(td$y) + 0.6), expand = c(0, 0)) +
    ggplot2::labs(title = paste0(endpoint_label, " ; HR (", ci_label, ", p value)")) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title  = ggplot2::element_text(face = "bold", hjust = 0, size = 11),
                   plot.margin = ggplot2::margin(8, 2, 8, 8))

  p_forest <- ggplot2::ggplot(td, ggplot2::aes(y = y)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = plot_lo, xmax = plot_hi),
                           orientation = "y", width = 0, na.rm = TRUE,
                           linewidth = 0.5, color = "black") +
    ggplot2::geom_point(ggplot2::aes(x = plot_est, size = box_size),
                        shape = 22, fill = box_fill, color = "black",
                        na.rm = TRUE, stroke = 0.4) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_log10(breaks = x_breaks_used,
                           labels = x_breaks_used,
                           limits = x_lim,
                           expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(0.4, max(td$y) + 0.6), expand = c(0, 0)) +
    ggplot2::labs(x = paste0(endpoint_label, ", ", ci_label, " (log scale)"), y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.y        = ggplot2::element_blank(),
                   axis.ticks.y       = ggplot2::element_blank(),
                   panel.border       = ggplot2::element_blank(),
                   axis.line.x        = ggplot2::element_line(color = "black"),
                   plot.margin        = ggplot2::margin(8, 8, 8, 2))

  combo <- p_text + p_forest + patchwork::plot_layout(widths = c(1.4, 1))

  grDevices::pdf(file.path(out_dir, paste0(file_base, ".pdf")),
                 width = w, height = h)
  print(combo); grDevices::dev.off()
  grDevices::png(file.path(out_dir, paste0(file_base, ".png")),
                 width = w, height = h, units = "in", res = 300)
  print(combo); grDevices::dev.off()

  invisible(combo)
}
```

### 3.5 Calling it

```r
mv_os <- coxph(Surv(os_months, os_event) ~ Histology + Treatment + Education,
               data = df_analysis)
plot_forest_oncology(mv_os, file_base = "fig_forest_mv_os",
                     endpoint_label = "Overall survival")
```

Two files appear in `here::here("output")`: `fig_forest_mv_os.pdf` and `fig_forest_mv_os.png`.

### 3.6 Output convention

- **PDF**: base `pdf()` device, vector, the manuscript figure. (If you need embedded Unicode fonts, swap in `grDevices::cairo_pdf()` — requires Cairo / XQuartz on macOS.)
- **PNG**: 300 dpi, the slide / supplement preview.
- Default `w = 11` inches; `h` auto-scales as `max(4.5, 0.32 * n_rows + 1.5)` so a 12-row forest is ~5.3 in and a 30-row forest is ~11.1 in. Override both if your journal has hard size limits.
- `out_dir` defaults to `here::here("output")`. Create the directory if needed (the function does this for you).

---

## 4. The "HR = 1.00 (1.00–1.00)" trap

> ❌ The single most dangerous failure mode of an oncology forest plot. The bug is **silent** — your figure will render cleanly and a degenerate row will look like a real null result. Reviewers will not catch it.

### 4.1 Three degenerate shapes — what they look like and why

`coxph` does not raise an error when a covariate level cannot be estimated. Instead it returns one of three shapes:

| Shape | What `coxph` returns | What `exp()` of it looks like | Why it happens in oncology data |
|---|---|---|---|
| **(a) Zero placeholder** | `coef = 0`, `SE = 0`, `p = NaN` | **HR = 1.00 (1.00–1.00)** — looks like a real null | Sparse stratum (e.g., n=1 patient in "CT alone"); single-cell pattern with no events |
| **(b) Plain NA** | `coef = NA`, `SE = NA` | NA — `broom.helpers` usually filters, but double-check | Perfect collinearity; categorical with one empty level after `na.omit` |
| **(c) Extreme finite** | `coef ≈ ±20`, `SE ≈ thousands` | HR like `4.85e8 (0, Inf)` | Quasi-separation: an n=1 event cell or a covariate that perfectly predicts the outcome |

Shape (a) is the killer. Shape (c) crashes `forestmodel` outright (§6). Shape (b) is mostly handled, but trust nothing.

The detector below catches shapes (a) and (b) directly; shape (c) is caught downstream by the extreme-CI filter in §4.3.

### 4.2 Detect — `nonest_terms()` helper

```r
#' Identify Cox terms with degenerate coefficients.
#'
#' Returns the row names of `summary(model)$coefficients` rows where
#' SE is NA or zero, p-value is NA, or coefficient is NA. These rows
#' MUST be masked in the forest data and the .docx table.
nonest_terms <- function(model) {
  s <- summary(model)$coefficients
  if (is.null(s) || nrow(s) == 0) return(character(0))
  bad <- is.na(s[, "se(coef)"]) | s[, "se(coef)"] == 0 |
         is.na(s[, "Pr(>|z|)"]) | is.na(s[, "coef"])
  rownames(s)[bad]
}
```

### 4.3 Filter in the forest data

The filter that `plot_forest_oncology()` (§3.4) applies inside its `plottable = ...` mutate is the canonical version. If you build your own forest data frame from `broom.helpers::tidy_plus_plus()`, apply BOTH the non-estimable filter AND the extreme-CI filter:

```r
plottable <- !is_header & !is_ref &
             is.finite(estimate) & is.finite(conf.low) & is.finite(conf.high) &
             conf.high < ci_extreme_hi & conf.low  > ci_extreme_lo &
             !is.na(std.error) & std.error > 0 &
             !is.na(p.value)
```

> ✅ All five clauses are required. Drop `!is.na(std.error) & std.error > 0` and the shape (a) HR=1.00 sneaks through. Drop the `conf.high < ci_extreme_hi` clause and shape (c) plots a square at HR = 4.85e8 that compresses every real HR into a vertical sliver.

### 4.4 Mark in the companion `.docx` table — `mark_nonest_tbl()`

The forest plot and the `gtsummary` regression table must agree. If the forest renders "Not estimable" but the `.docx` table still reads `HR = 1.00 (1.00, 1.00), p = 1.000`, the reader sees a contradiction and trusts the table. Mask the degenerate rows in the table too:

```r
#' Blank HR / CI / p in a gtsummary table for non-estimable Cox terms.
#'
#' Requires gtsummary >= 1.7 (column names: estimate, conf.low, conf.high, p.value).
mark_nonest_tbl <- function(tbl, model) {
  bad <- nonest_terms(model)
  if (length(bad) == 0) return(tbl)
  message("[mark_nonest_tbl] Non-estimable terms blanked: ",
          paste(bad, collapse = "; "))
  tbl |>
    gtsummary::modify_table_body(
      ~ dplyr::mutate(
          .x,
          estimate  = ifelse(!is.na(term) & term %in% bad, NA_real_, estimate),
          conf.low  = ifelse(!is.na(term) & term %in% bad, NA_real_, conf.low),
          conf.high = ifelse(!is.na(term) & term %in% bad, NA_real_, conf.high),
          p.value   = ifelse(!is.na(term) & term %in% bad, NA_real_, p.value)
        )
    )
}
```

Use it like:

```r
tbl_mv_os <- tbl_regression(mv_os, exponentiate = TRUE) |>
  mark_nonest_tbl(mv_os)
tbl_mv_os |> as_flex_table() |>
  flextable::save_as_docx(path = here::here("output", "table_mv_os.docx"))
```

### 4.5 Fix the data — `fct_collapse()` worked example

Masking is a display fix, not a statistical fix. The underlying problem is that some covariate levels are too sparse to estimate. If clinically meaningful, **collapse rare levels before fitting**. The trade-off: interpretability goes down (you lose granularity) but every coefficient becomes estimable. Collapses should be pre-specified in the analysis plan; ad-hoc collapsing after seeing the results is a form of p-hacking.

Example — a registry `Treatment` factor with 9 levels, three of which have n ≤ 2 patients:

```r
table(df_analysis$Treatment)
#> Surgery                       Surgery + RT                 Surgery + CT
#>     58                                 41                            2
#> Surgery + RT + CT             RT alone                     RT + CT
#>     12                                  8                            5
#> CT alone                      Best supportive care         Unknown
#>      1                                  2                            1

df_analysis <- df_analysis |>
  dplyr::mutate(
    Treatment_collapsed = forcats::fct_collapse(
      Treatment,
      "Surgery only"             = c("Surgery"),
      "Surgery + adjuvant"       = c("Surgery + RT", "Surgery + CT",
                                     "Surgery + RT + CT"),
      "Non-surgical / palliative" = c("RT alone", "RT + CT", "CT alone",
                                      "Best supportive care", "Unknown")
    ) |>
    forcats::fct_relevel("Surgery only")  # set the clinical reference
  )
```

When in doubt about the right collapse, **ask the user** — the right grouping is driven by clinical intent (curative vs palliative, surgery-first vs not), not by sample size alone.

---

## 5. Subgroup forest — one Cox per subgroup

### 5.1 When to use this layout (and when not)

**Use it when**: you want to show whether the treatment effect is consistent across pre-specified subgroups (age bands, histology, ECOG, stage). Each row is the treatment HR within one subgroup level, estimated from a Cox model fit on that subgroup only.

**Don't use it when**: you have run 12 post-hoc subgroup analyses and want to advertise the one with a "significant" HR. That is a fishing expedition. Use a single overall model with an interaction term and report ONE p-interaction per subgroup variable (§5.3) instead of 12 separate confidence intervals.

See SKILL.md "Exploratory vs confirmatory" for the broader principle.

### 5.2 Building the tidy data frame (`purrr::map_dfr`)

For each level of one subgroup variable, fit `coxph(Surv(time, event) ~ treatment)` on the subset and extract the treatment HR. Skip levels with `< 10` events (render as "Not estimable" rather than silently drop — the reader needs to see WHY a stratum is missing).

```r
#' Build a tidy data frame of treatment HRs across subgroup levels.
#'
#' @param data        Analysis data frame.
#' @param time_var    Name of the time-to-event column (string).
#' @param event_var   Name of the 0/1 event column (string).
#' @param treat_var   Name of the treatment column (string); the reference
#'                    level must already be set via forcats::fct_relevel().
#' @param subgroup_var Name of the subgroup column (string).
#' @param min_events  Minimum events per stratum to attempt a fit. Default 10.
build_subgroup_tidy <- function(data, time_var, event_var, treat_var,
                                subgroup_var, min_events = 10) {
  levels_use <- levels(droplevels(data[[subgroup_var]]))

  purrr::map_dfr(levels_use, function(lev) {
    d <- data[!is.na(data[[subgroup_var]]) & data[[subgroup_var]] == lev, ]
    n      <- nrow(d)
    events <- sum(d[[event_var]] == 1, na.rm = TRUE)

    if (events < min_events ||
        length(unique(d[[treat_var]])) < 2) {
      return(tibble::tibble(
        subgroup_label = lev, n = n, events = events,
        estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        std.error = NA_real_, p.value = NA_real_,
        note = sprintf("Not estimable (n=%d, events=%d)", n, events)
      ))
    }

    f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, treat_var))
    fit <- tryCatch(survival::coxph(f, data = d), error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble::tibble(
        subgroup_label = lev, n = n, events = events,
        estimate = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        std.error = NA_real_, p.value = NA_real_,
        note = "coxph error"
      ))
    }

    s <- summary(fit)$coefficients
    ci <- summary(fit)$conf.int
    tibble::tibble(
      subgroup_label = lev,
      n = n, events = events,
      estimate  = ci[1, "exp(coef)"],
      conf.low  = ci[1, "lower .95"],
      conf.high = ci[1, "upper .95"],
      std.error = s[1, "se(coef)"],
      p.value   = s[1, "Pr(>|z|)"],
      note      = NA_character_
    )
  })
}
```

Call site:

```r
sub_age <- build_subgroup_tidy(df_analysis,
                               time_var = "os_months", event_var = "os_event",
                               treat_var = "Treatment_collapsed",
                               subgroup_var = "age_band")
```

### 5.3 Interaction p — the right way to advertise heterogeneity

Do NOT compare the row-level p-values across subgroups and claim "the treatment works in subgroup X". Per-stratum p-values reflect the within-subgroup precision, not heterogeneity of the effect. The right test is a likelihood-ratio test (LRT) of `treatment * subgroup_var` against `treatment + subgroup_var` on the WHOLE dataset:

```r
interaction_p <- function(data, time_var, event_var, treat_var, subgroup_var) {
  f_add <- as.formula(sprintf("Surv(%s, %s) ~ %s + %s",
                              time_var, event_var, treat_var, subgroup_var))
  f_int <- as.formula(sprintf("Surv(%s, %s) ~ %s * %s",
                              time_var, event_var, treat_var, subgroup_var))
  m_add <- survival::coxph(f_add, data = data)
  m_int <- survival::coxph(f_int, data = data)
  anova(m_add, m_int, test = "LRT")[2, "Pr(>|Chi|)"]
}
```

Call it once per subgroup variable:

```r
p_int_age <- interaction_p(df_analysis, "os_months", "os_event",
                           "Treatment_collapsed", "age_band")
```

Display ONE `p_interaction` per subgroup variable next to the forest header (e.g., as a subtitle or in the bottom-right corner), not on each row.

### 5.4 Code — `plot_forest_subgroup()` (designed for this reference)

> ⚠ This function is **designed from first principles for this reference**; it has not yet been battle-tested in user projects the way `plot_forest_oncology()` has. Treat it as a starting template. Mirrors the visual style of §3.4 (precision-scaled squares, two-panel patchwork, log-scale x-axis, dashed HR=1 line) but takes a tidy data frame instead of a model.

Required input columns (the output of `build_subgroup_tidy()` already conforms):

| Column | Type | Required | Notes |
|---|---|---|---|
| `subgroup_label` | character | yes | One row per subgroup level. |
| `n`             | integer   | yes | Patients in the stratum. |
| `events`        | integer   | yes | Events in the stratum. |
| `estimate`      | double    | yes | Treatment HR; NA for non-estimable. |
| `conf.low`      | double    | yes | NA for non-estimable. |
| `conf.high`     | double    | yes | NA for non-estimable. |
| `std.error`     | double    | yes | Used for precision-scaled squares. |
| `p.value`       | double    | yes | Per-stratum; NOT used for the heterogeneity claim. |
| `p_interaction` | double    | no  | A single scalar; pass via the `p_interaction` argument instead of as a column. |

```r
#' Subgroup forest plot — one row per subgroup level.
plot_forest_subgroup <- function(tidy_df,
                                 file_base,
                                 endpoint_label  = "HR",
                                 ci_label        = "95% CI",
                                 p_interaction   = NA_real_,
                                 box_fill        = "#1f3b6f",
                                 ci_extreme_hi   = 1e3,
                                 ci_extreme_lo   = 1e-4,
                                 x_breaks        = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
                                 w               = 11,
                                 h               = NULL,
                                 out_dir         = here::here("output")) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  td <- tidy_df |>
    dplyr::mutate(
      plottable = is.finite(estimate) & is.finite(conf.low) & is.finite(conf.high) &
                  conf.high < ci_extreme_hi & conf.low > ci_extreme_lo &
                  !is.na(std.error) & std.error > 0 & !is.na(p.value),
      plot_est = ifelse(plottable, estimate, NA_real_),
      plot_lo  = ifelse(plottable, conf.low, NA_real_),
      plot_hi  = ifelse(plottable, conf.high, NA_real_),
      n_ev_lab = sprintf("%d / %d", events, n),
      hr_text = dplyr::case_when(
        !plottable ~ "Not estimable",
        TRUE ~ sprintf("%.2f (%.2f-%.2f), p=%s",
                       estimate, conf.low, conf.high,
                       format.pval(p.value, digits = 2, eps = 0.001))
      ),
      precision = ifelse(plottable, 1 / std.error, NA_real_)
    )

  if (any(!is.na(td$precision))) {
    td$box_size <- scales::rescale(td$precision, to = c(2.5, 7),
                                   from = range(td$precision, na.rm = TRUE))
  } else {
    td$box_size <- 3
  }
  td$box_size[!td$plottable] <- NA

  td$row_id <- seq_len(nrow(td))
  td$y      <- max(td$row_id) - td$row_id + 1

  if (is.null(h)) h <- max(4.5, 0.32 * nrow(td) + 1.5)

  x_finite <- c(td$plot_lo, td$plot_hi)
  x_finite <- x_finite[is.finite(x_finite)]
  if (length(x_finite) == 0) x_finite <- c(0.5, 2)
  x_breaks_used <- x_breaks[x_breaks >= min(x_finite, 1) * 0.5 &
                            x_breaks <= max(x_finite, 1) * 2]
  if (length(x_breaks_used) < 3) x_breaks_used <- x_breaks
  x_lim <- range(c(x_breaks_used, x_finite, 1), na.rm = TRUE)

  subtitle <- if (!is.na(p_interaction)) {
    sprintf("p (interaction) = %s",
            format.pval(p_interaction, digits = 2, eps = 0.001))
  } else {
    ""
  }

  p_text <- ggplot2::ggplot(td, ggplot2::aes(y = y)) +
    ggplot2::geom_text(ggplot2::aes(x = 0,    label = subgroup_label),
                       hjust = 0, fontface = "bold", size = 3.6) +
    ggplot2::geom_text(ggplot2::aes(x = 0.85, label = n_ev_lab),
                       hjust = 0, size = 3.3, family = "mono") +
    ggplot2::geom_text(ggplot2::aes(x = 1.30, label = hr_text),
                       hjust = 0, size = 3.3, family = "mono") +
    ggplot2::scale_x_continuous(limits = c(0, 2.7), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, max(td$y) + 0.6), expand = c(0, 0)) +
    ggplot2::labs(title    = paste0(endpoint_label, " by subgroup"),
                  subtitle = subtitle) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", hjust = 0, size = 11),
                   plot.subtitle = ggplot2::element_text(hjust = 0, size = 10),
                   plot.margin   = ggplot2::margin(8, 2, 8, 8))

  p_forest <- ggplot2::ggplot(td, ggplot2::aes(y = y)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = plot_lo, xmax = plot_hi),
                           orientation = "y", width = 0, na.rm = TRUE,
                           linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = plot_est, size = box_size),
                        shape = 22, fill = box_fill, color = "black",
                        na.rm = TRUE, stroke = 0.4) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_x_log10(breaks = x_breaks_used, labels = x_breaks_used,
                           limits = x_lim,
                           expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(0.4, max(td$y) + 0.6), expand = c(0, 0)) +
    ggplot2::labs(x = paste0(endpoint_label, " HR, ", ci_label, " (log scale)"), y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   axis.text.y        = ggplot2::element_blank(),
                   axis.ticks.y       = ggplot2::element_blank(),
                   panel.border       = ggplot2::element_blank(),
                   axis.line.x        = ggplot2::element_line(color = "black"),
                   plot.margin        = ggplot2::margin(8, 8, 8, 2))

  combo <- p_text + p_forest + patchwork::plot_layout(widths = c(1.4, 1))

  grDevices::pdf(file.path(out_dir, paste0(file_base, ".pdf")),
                 width = w, height = h)
  print(combo); grDevices::dev.off()
  grDevices::png(file.path(out_dir, paste0(file_base, ".png")),
                 width = w, height = h, units = "in", res = 300)
  print(combo); grDevices::dev.off()

  invisible(combo)
}
```

### 5.5 Calling it

```r
sub_age <- build_subgroup_tidy(df_analysis, "os_months", "os_event",
                               "Treatment_collapsed", "age_band")
p_int   <- interaction_p(df_analysis, "os_months", "os_event",
                         "Treatment_collapsed", "age_band")
plot_forest_subgroup(sub_age, file_base = "fig_forest_subgroup_age_os",
                     endpoint_label = "Overall survival",
                     p_interaction  = p_int)
```

---

## 6. `forestmodel` — when it helps, when it breaks

### 6.1 Use as a sanity check

Render `forestmodel::forest_model(cox_mv)` next to the custom plot from §3.4. Divergence between the two — a row that exists in one and not the other, an HR that differs, a variable that one renders and the other drops — almost always means there is a bug in your manual data wrangling (a missing filter, a typo in a column name, a `mutate` that overwrote `estimate` with `log(estimate)`). The custom plot is the manuscript figure; `forest_model()` is the second pair of eyes.

### 6.2 The base-R `axis()` crash

`forestmodel::forest_model()` renders log-scale ticks via base R's `axis(..., log = TRUE)`. If ANY covariate in the model has a coefficient that spans more than ~8 orders of magnitude — which happens whenever shape (c) from §4.1 sneaks through — the tick generator fails with:

```
Error in log - axis(), 'at' creation, _LARGE_ range: invalid {xy}axp or par
```

One degenerate row kills the entire plot. There is no warning until you call `print()` on the result.

### 6.3 Workaround — `save_forestmodel()` wrapper

Auto-exclude the degenerate covariates from the `forest_model()` call via the `covariates =` argument, using the SAME `nonest_terms()` + extreme-CI filter you already trust from §4. Log which variables got dropped so the reader of the run log understands why the QC plot is missing rows the manuscript figure shows as "Not estimable".

```r
#' Render forestmodel::forest_model with auto-exclusion of degenerate covariates.
#'
#' Drops any variable that contains a non-estimable term (per `nonest_terms()`)
#' or a level whose CI bound is outside [ci_extreme_lo, ci_extreme_hi]. Logs
#' the exclusions. Saves PDF + PNG.
save_forestmodel <- function(model, file_base, w = 11, h = NULL,
                             ci_extreme_hi = 1e3, ci_extreme_lo = 1e-4,
                             out_dir = here::here("output")) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  td <- tryCatch(
    broom.helpers::tidy_plus_plus(model, exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL
  )
  if (is.null(td)) {
    message("[save_forestmodel] tidy_plus_plus failed for ", file_base)
    return(invisible(NULL))
  }

  bad_terms_v <- nonest_terms(model)
  is_extreme <- (!is.na(td$conf.high) & td$conf.high > ci_extreme_hi) |
                (!is.na(td$conf.low)  & td$conf.low  < ci_extreme_lo)
  is_bad     <- (!is.na(td$term) & td$term %in% bad_terms_v) | is_extreme
  exclude_vars <- unique(td$variable[is_bad])
  exclude_vars <- exclude_vars[!is.na(exclude_vars)]

  all_vars <- unique(td$variable[!is.na(td$variable)])
  use_vars <- setdiff(all_vars, exclude_vars)

  if (length(exclude_vars) > 0) {
    message("[save_forestmodel] ", file_base,
            " - dropping variables with degenerate/extreme levels from the ",
            "forestmodel rendering (kept in the custom plot): ",
            paste(exclude_vars, collapse = ", "))
  }
  if (length(use_vars) == 0) {
    message("[save_forestmodel] ", file_base,
            " - no variables left after exclusion; skipping forestmodel plot.")
    return(invisible(NULL))
  }

  if (is.null(h)) {
    n_levels <- sum(td$variable %in% use_vars)
    h <- max(4.5, 0.45 * n_levels + 3)
  }

  out <- tryCatch({
    forestmodel::forest_model(model, covariates = use_vars)
  }, error = function(e) {
    message("[save_forestmodel] forestmodel::forest_model failed for ",
            file_base, " even after exclusion: ", conditionMessage(e))
    NULL
  })
  if (is.null(out)) return(invisible(NULL))

  ggplot2::ggsave(file.path(out_dir, paste0(file_base, ".pdf")),
                  out, width = w, height = h, device = "pdf")
  ggplot2::ggsave(file.path(out_dir, paste0(file_base, ".png")),
                  out, width = w, height = h, dpi = 300, device = "png")
  invisible(out)
}
```

### 6.4 Uniform-square caveat

`forest_model()` uses uniform square sizes — fine for QC, wrong for the manuscript figure (see §3.1). When the manuscript reviewer asks "why are the squares all the same size?", they have looked at the QC plot, not the manuscript figure. Make sure §3.4 is the one that goes into the paper and §6.3 is the one that goes into the supplement (if anywhere).

---

## 7. Practical guardrails

- **Set the reference level explicitly with `forcats::fct_relevel()` BEFORE `coxph()`.** The reference category drives the interpretation of every HR in the forest. When the right reference is not obvious (e.g., for a histology variable with both "Atypical (G2)" and "Malignant (G3)"), **ask the user** which one they want as the comparator. Never rely on alphabetical default.

- **Save both PDF (`pdf()`) and PNG (300 dpi).** PDF is the manuscript figure; PNG is the slide / supplement preview. Both `plot_forest_oncology()` and `plot_forest_subgroup()` do this for you. Swap `pdf()` → `grDevices::cairo_pdf()` if you need embedded Unicode fonts and have Cairo / XQuartz.

- **Print `summary(model)` to the run log** after every Cox fit so the table, the forest, and the raw Wald output can all be cross-checked. If three sources disagree, the bug is upstream of all of them.

- **Wrap `cox.zph()` in `tryCatch`** — singular models will throw, and you do not want one PH diagnostic to abort an entire analysis script:

  ```r
  zph_os <- tryCatch(cox.zph(mv_os), error = function(e) {
    message("cox.zph failed (likely sparse stratum / singular information matrix).")
    NULL
  })
  if (!is.null(zph_os)) print(zph_os)
  ```

- **For MV Cox with rare strata, run BOTH a primary sparse-stratum-aware MV AND a sensitivity collapsed-Treatment MV (§4.5).** Show the primary in the main figure and the sensitivity in the supplement; reviewers will ask. If the two analyses qualitatively agree, the conclusion is robust; if they diverge, that itself is the finding.

---

## 8. Common failure modes — quick reference

| Failure mode | Symptom | Fix |
|---|---|---|
| Fake null HR | A row in the forest shows `HR = 1.00 (1.00–1.00), p = 1.000` for a covariate level you know is sparse. | §4 — apply `nonest_terms()` + the §4.3 filter; render the row as "Not estimable". |
| `forestmodel` axis crash | `Error in log - axis(), 'at' creation, _LARGE_ range: invalid {xy}axp or par` when `print()`-ing `forest_model(cox_mv)`. | §6.3 — wrap the call via `save_forestmodel()`, which auto-excludes degenerate covariates. |
| Subgroup-fishing claim | "The treatment works in patients over 65 (HR 0.42, p = 0.03)" pulled from a 12-subgroup forest. | §5.3 — replace per-row p-values with ONE `p_interaction` per subgroup variable from an interaction LRT. |
| Forest and `.docx` table disagree | Forest shows "Not estimable"; companion table shows `HR = 1.00 (1.00, 1.00), p = 1.000`. | §4.4 — pipe the gtsummary table through `mark_nonest_tbl()`. |
| Squares all the same size | The manuscript figure has uniform squares. | §3.1 — you are looking at the `forestmodel` QC plot, not the precision-scaled `plot_forest_oncology()` output. |
| Wrong reference category | An HR seems flipped (e.g., "Surgery + RT" shows HR > 1 when you expected < 1). | §7 — check `levels(df$Treatment)[1]`; set explicitly with `forcats::fct_relevel()` and refit. |

---

## See also

- `references/02-survival-analysis.md` — Cox model fitting, PH diagnostics, when to switch to RMST.
- `references/05-manuscript-figures.md` — general export conventions (PDF/EPS/TIFF/PNG, `dev.off()`); meta-analysis forest plots for pooled HRs across studies.
- `references/06-manuscript-tables.md` — exporting the `gtsummary` companion table to `.docx`.
- `references/09-package-quickref.md` — one-line reminders for `broom.helpers`, `patchwork`, `scales`, `forestmodel`.
