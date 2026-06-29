# Data inspection — the first 5 minutes with a new user file

When the user hands you a `.xlsx` or `.csv` (or a data frame), do this BEFORE any analysis. Skipping these steps is the most common source of clinically wrong results.

---

## The inspection checklist

Run all of these, in order, and share the output (or a summary) with the user:

1. **Load the file** — use the right reader, don't auto-detect via the wrong package.
2. **Print variable names** — clean them with `janitor::clean_names()` to a new object; keep the original loaded for reference.
3. **Print structure** — variable types (numeric? factor? character?).
4. **Print missingness** — per-variable count and percent.
5. **Print a 6-row preview** — to confirm the shape matches the user's mental model.
6. **For survival data specifically** — verify event coding, time range, and reference level.

---

## 1. Loading

### Excel

```r
library(readxl)
library(here)

# Inspect sheets first
readxl::excel_sheets(here::here("data", "trial_data.xlsx"))

# Then read the target sheet
df_raw <- readxl::read_excel(
  path  = here::here("data", "trial_data.xlsx"),
  sheet = "patients",
  na    = c("", "NA", "N/A", "NULL", ".", "?")     # common sentinels
)
```

Common pitfalls:

- `read_excel` reads the first sheet by default — explicit `sheet =` prevents accidentally analyzing the wrong tab.
- Excel dates come in as numeric (days since 1899-12-30 on Windows, 1904-01-01 on old Mac). Convert with `janitor::excel_numeric_to_date()` or `lubridate::as_date(value, origin = "1899-12-30")`.
- Merged cells, hidden rows, and color-coded metadata are all lost — ask the user to confirm none of those carry information.

### CSV

```r
library(readr)

df_raw <- readr::read_csv(
  file  = here::here("data", "trial_data.csv"),
  na    = c("", "NA", "N/A", "NULL", ".", "?"),
  guess_max = 10000          # avoid type-guessing errors when early rows are NA
)

# If the file is European-formatted (semicolon-separated, comma decimal):
df_raw <- readr::read_csv2(here::here("data", "trial_data_eu.csv"))
```

Pitfalls:

- `read_csv` (American) vs `read_csv2` (European). Brazilian/European files often use `;` and `,` as decimal — wrong reader → garbage.
- Date columns parsed as character → convert explicitly with `lubridate::dmy()`, `mdy()`, `ymd()` per the source format.

---

## 2. Variable names

```r
library(janitor)

df_clean <- df_raw |> janitor::clean_names()

# Compare
tibble::tibble(original = names(df_raw), cleaned = names(df_clean))
```

`clean_names()` lowercases, snake_cases, removes special characters, and de-duplicates. Always create a NEW object (`df_clean`) — keep `df_raw` for reference.

---

## 3. Structure

```r
dplyr::glimpse(df_clean)
# OR
skimr::skim(df_clean)
```

`glimpse` is concise; `skimr` gives missingness and distribution summaries in one call. Use `skim` for a thorough first look, `glimpse` for follow-up checks.

For survival data, watch for:

- Time-to-event variable typed as character (often happens with Excel exports of integer-looking columns) → must convert to numeric.
- Date columns typed as character (or as numeric Excel serial) → must convert to Date.
- Event variable typed as character with "yes"/"no" or "alive"/"dead" → must recode to 0/1.

---

## 4. Missingness

```r
df_clean |>
  dplyr::summarise(dplyr::across(everything(), ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  dplyr::mutate(pct_missing = round(100 * n_missing / nrow(df_clean), 1)) |>
  dplyr::arrange(dplyr::desc(n_missing))
```

Or use `naniar`:

```r
naniar::miss_var_summary(df_clean)
naniar::vis_miss(df_clean)             # visual missingness map
```

For any variable in the planned analysis: ask the user whether missingness is informative (MAR vs MCAR vs MNAR) and what to do — complete-case, single-imputation, multiple imputation (`mice`), or treat as "Unknown" category.

---

## 5. Preview

```r
df_clean |> dplyr::slice_sample(n = 6)     # random 6 rows — better than head() for catching patterns
```

Confirm with the user: "Does this look like what you expect for a row?" Surprisingly often there's a one-row offset, a hidden header row, or an extra column.

---

## 6. Survival-specific sanity checks

### Event coding

```r
df_clean |> dplyr::count(os_event)
```

Expected: two levels (typically 0 and 1). If you see 1/2, "Alive"/"Dead", TRUE/FALSE, recode:

```r
df_clean <- df_clean |>
  dplyr::mutate(
    os_event = dplyr::case_when(
      os_event %in% c(1, "Dead", "Yes", TRUE)  ~ 1L,
      os_event %in% c(0, "Alive", "No", FALSE) ~ 0L,
      TRUE                                     ~ NA_integer_
    )
  )
```

Verify post-recode counts match the raw table.

### Time range

```r
summary(df_clean$os_months)
hist(df_clean$os_months, breaks = 40, main = "OS time distribution")
```

Sanity checks:

- Negative values → typo or wrong reference date; flag.
- Zeros → patients with event at randomization; usually means time was rounded, ask the user.
- Implausibly large values → time was in days but you assumed months (or vice versa). Confirm the time unit BEFORE plotting.

### Reference level

```r
levels(factor(df_clean$arm))
```

If the reference is wrong:

```r
df_clean$arm <- forcats::fct_relevel(df_clean$arm, "Control")
```

For multi-level categorical covariates, set the clinically meaningful reference (ECOG 0, Stage I, wild-type, etc.).

### Event rate per arm

```r
df_clean |>
  dplyr::group_by(arm) |>
  dplyr::summarise(
    n         = dplyr::n(),
    events    = sum(os_event, na.rm = TRUE),
    censored  = sum(os_event == 0, na.rm = TRUE),
    median_t  = median(os_months, na.rm = TRUE)
  )
```

This catches:

- Very low event counts → HR estimates will be unstable.
- All events in one arm → cross-arm comparison degenerate.
- Median time wildly different from what the user expected → wrong column or wrong unit.

---

## What to share with the user

After running the checklist, summarize in ~5 lines:

> Read `trial_data.xlsx` (sheet "patients"): 342 rows, 18 columns. Variables renamed to snake_case. Missingness > 5% in `ldh` (12.6%) and `bmi` (8.2%); everything else < 2%. `os_event` is coded 0/1 (286 censored, 56 events) — confirming "1 = death" by your description. `os_months` ranges 0.5–58.3, median 24.7 — looks like months, confirm. `arm` is a 2-level factor ("Control", "Experimental"); reference set to "Control".

Then ask any necessary follow-ups (event coding for ambiguous variables, time unit if unclear, missingness handling).

---

## Anti-patterns

- ❌ Running `survfit()` directly on a freshly-read dataset without inspection.
- ❌ Overwriting `df_raw` in place (`df_raw <- clean_names(df_raw)`).
- ❌ Coercing event with `as.numeric(event)` when the variable is character "Alive"/"Dead" — produces NAs silently.
- ❌ Assuming time is months without checking.
- ❌ Letting `read_excel` auto-pick sheet 1 when the user has multiple tabs.
- ❌ Skipping the missingness summary — surprises the user later.
