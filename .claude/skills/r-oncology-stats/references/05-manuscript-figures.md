# Manuscript-quality figures (KM curves, forest plots, export)

This file covers: choosing `ggsurvplot` vs `ggsurvfit`, theme and label hygiene, forest plots, file formats, dimensions, and the `dev.off()` rule that trips up new R users.

---

## Kaplan-Meier curves

### `ggsurvplot` (survminer) vs `ggsurvfit`

| Need | Pick |
|---|---|
| Quick manuscript-ready figure with risk table + p-value + censor ticks | `survminer::ggsurvplot` |
| Tidier `ggplot2` composability (layer geoms, facet, theme) | `ggsurvfit::ggsurvfit` |
| Cumulative incidence / competing risks visual | `ggsurvfit` (better) or `cmprsk::cuminc` |
| Multi-panel comparison across endpoints | `ggsurvfit` + `patchwork` |

Both produce publication-quality output. `survminer` has more knobs out of the box; `ggsurvfit` is cleaner if you want to layer things.

### Manuscript template (`ggsurvplot`)

```r
library(survminer)
library(ggplot2)

km <- survminer::ggsurvplot(
  os_fit, data = df_analysis,

  # core
  risk.table        = TRUE,
  risk.table.height = 0.22,
  conf.int          = FALSE,            # CIs only when single-arm or specifically wanted
  censor            = TRUE,
  censor.shape      = "|",
  censor.size       = 2,

  # statistics overlay
  pval              = TRUE,
  pval.method       = TRUE,
  pval.coord        = c(3, 0.05),
  pval.size         = 4,

  # axes & labels
  xlab              = "Time since randomization (months)",
  ylab              = "Overall survival",     # ⚠ specific endpoint
  xlim              = c(0, 60),
  break.x.by        = 12,
  ylim              = c(0, 1),
  break.y.by        = 0.2,

  # legend
  legend.title      = "",
  legend.labs       = c("Standard of care", "Experimental"),
  legend            = c(0.8, 0.85),

  # theme & palette
  ggtheme           = ggplot2::theme_classic(base_size = 12),
  palette           = c("#1F77B4", "#D62728"),

  # risk table styling
  risk.table.col    = "strata",
  risk.table.y.text = FALSE,
  tables.theme      = survminer::theme_cleantable()
)
print(km)
```

### Label hygiene (⚠ IRON RULE)

- ✅ "Overall survival" — known endpoint.
- ✅ "Progression-free survival" — known endpoint.
- ✅ "Disease-free survival", "Recurrence-free survival", "Event-free survival" — known endpoints.
- ❌ "Cumulative survival" — meaningless generic; never use when the endpoint is known.
- ❌ "Survival probability" — only acceptable in a generic methods illustration, not a results figure.
- ❌ "Time (mo)" — write the full word "months" for the axis; abbreviations look amateurish in print.

The y-axis can be on probability scale (0–1) or percent (0–100). Most journals prefer probability with tick labels "0", "0.2", ..., "1.0" or percent "0%", "20%", ..., "100%". Pick one and stick to it across all figures in the paper.

### Risk table tips

- Show the risk table whenever you show a KM curve. It's standard for clinical journals and lets the reader gauge precision over time.
- Don't include the risk table on a single-arm curve unless follow-up is heterogeneous.
- If risk-table labels collide with axis labels, set `risk.table.y.text = FALSE` and label by color in the figure legend.

---

## Forest plots

Three flavors:

- **Cox multivariable coefficients (manuscript-grade)** — precision-scaled two-panel `ggplot2` + `patchwork` layout with explicit handling of non-estimable rows (the `HR = 1.00 (1.00–1.00)` trap). See `references/11-forest-plots.md` §3.
- **Subgroup forest (one Cox per subgroup level)** — `purrr::map_dfr` builder + sibling `plot_forest_subgroup()`. See `references/11-forest-plots.md` §5.
- **Meta-analysis forest (pooling published HRs)** — `meta::metagen()` + `forest()`. See below.

> `forestmodel::forest_model()` is documented in `references/11-forest-plots.md` §6 as a QC sanity check, not the manuscript figure.

### Meta-analysis forest (`meta` or `metafor`)

```r
library(meta)
m <- metagen(TE = log_hr, seTE = se_log_hr, studlab = study,
             data = pooled, sm = "HR")
forest(m)
```

---

## Export — formats and when

### Vector vs raster

| Format | Use for | Notes |
|---|---|---|
| **PDF** | Manuscript figures with text, lines, transparency | Default for everything. Scales infinitely, small file size. |
| **EPS** | Some old journals require EPS | Older PostScript-based vector; produced by `cairo_ps` device. |
| **SVG** | Web display, slides that need re-editing | Native to browsers; readable XML. |
| **TIFF** | Journal-required raster at 300+ DPI | Larger files; use `compression = "lzw"`. |
| **PNG** | Web display, draft sharing, presentation slides | Lossless raster; lighter than TIFF. |
| **JPEG** | Almost never | Lossy; avoid for any figure with text or sharp lines. |

**Default:** PDF for vector figures; TIFF 600 DPI when the journal requires raster.

### Export idioms

```r
# ggplot output (KM, forest, anything ggplot-based)
ggplot2::ggsave(
  filename = here::here("output", "fig1_km_os.pdf"),
  plot     = km$plot,                    # km$plot from ggsurvplot
  width    = 6.5, height = 5,
  units    = "in",
  device   = "pdf"                       # portable default; see cairo_pdf note below
)

# For survminer object with risk table, save the composite:
pdf(here::here("output", "fig1_km_os.pdf"), width = 6.5, height = 5)
print(km)
dev.off()                                # ⚠ REQUIRED — see below
```

### The `dev.off()` rule

When using **base graphics devices** (`pdf()`, `tiff()`, `png()`, `cairo_pdf()`), you MUST close the device with `dev.off()` after writing. Otherwise the file is incomplete or unreadable.

```r
pdf("fig.pdf", width = 6, height = 4)
plot(my_fit)
dev.off()                                # ⚠ MUST be present
```

When using **`ggsave()`**, you do NOT need `dev.off()` — it handles the device lifecycle internally.

```r
ggsave("fig.pdf", my_ggplot, width = 6, height = 4)   # no dev.off() needed
```

Common bug: opening `pdf()` inside a function, the function errors, and the device is left open — subsequent plots write to the corrupt PDF. Wrap critical exports in `tryCatch` + `on.exit(dev.off())` for production scripts.

### `cairo_pdf` vs `"pdf"` portability

`grDevices::cairo_pdf()` embeds Unicode fonts and produces cleaner text (en-dashes, Greek letters, non-ASCII glyphs), but it depends on the Cairo graphics stack:

- On **Linux servers** (Docker, CI, most HPCs): usually present — safe to use.
- On **macOS**: requires **XQuartz** to be installed. When XQuartz is missing, `ggsave(..., device = cairo_pdf)` fails **silently** (a warning is emitted but no file is written, and the return value looks fine to callers).
- On **Windows**: usually present when R was installed with the standard binary.

Default to `device = "pdf"` in reference/template code so scripts are portable. When a manuscript figure needs Unicode glyphs, switch to `device = cairo_pdf` on a machine you have verified has Cairo — and **check `file.exists(...)` after each `ggsave`** so a missing device doesn't slip past.

```r
# Portable default:
ggplot2::ggsave("fig.pdf", plot = p, width = 6.5, height = 5, device = "pdf")

# Cairo variant — verify the device works on your machine first:
if (capabilities("cairo")) {
  ggplot2::ggsave("fig.pdf", plot = p, width = 6.5, height = 5, device = cairo_pdf)
} else {
  ggplot2::ggsave("fig.pdf", plot = p, width = 6.5, height = 5, device = "pdf")
}
```

---

## Dimensions for journals

Most clinical oncology journals (JCO, JAMA Onc, Lancet Onc, NEJM) accept:

- **Single column**: ~3.5 in / 89 mm wide
- **1.5 column**: ~5 in / 127 mm wide
- **Double column** (full width): ~7.0 in / 178 mm wide

For KM curves with a risk table, **5–6 in wide × 4–5 in tall** is a good default. Forest plots are usually full-width (7 in) and proportionally tall.

Check the journal's "Information for authors" for exact specs and resolution requirements.

---

## Color choices

- **Two-arm**: a colorblind-safe pair like `#1F77B4` (blue) and `#D62728` (red), or `viridis::cividis(2)`. Avoid red/green pairs.
- **Three+ arms**: `RColorBrewer::brewer.pal(n, "Dark2")` or `viridis::viridis(n)`.
- **Print-friendly**: test in grayscale (`scales::show_col(gray.colors(n))`) — if arms become indistinguishable, also vary linetype.
- Match colors across all figures in the manuscript — same arm should always have the same color.

---

## Font sizes

- Axis title: 11–12 pt
- Axis ticks: 10–11 pt
- Legend: 10–11 pt
- Risk table: 9–10 pt
- p-value annotation: 10–11 pt

`base_size = 12` in `theme_classic()` is a reasonable global default; override individual elements as needed.

---

## Anti-patterns

- ❌ `theme_gray()` (the ggplot default with gray background) for manuscript figures — use `theme_classic()`.
- ❌ Jet/rainbow palettes — perceptually non-uniform.
- ❌ 3D effects, drop shadows, gradients on bars — out of place in clinical journals.
- ❌ Tick marks every month on a 60-month axis — set `break.x.by = 12`.
- ❌ Omitting the risk table — required by most reviewers.
- ❌ "Survival" with no units on the x-axis ("Time" with no unit).
- ❌ Saving as PNG when the journal asks for vector or 300+ DPI raster.
