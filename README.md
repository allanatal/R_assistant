# R_assistant

A Claude Code skill for **clinical oncology statistical analysis in R**. Built for medical oncologists and clinical researchers who work in R on survival data, manuscript figures and tables, and reconstructed individual patient data (IPD) from published Kaplan–Meier curves.

The skill teaches Claude Code to:

- Ask the right clarifying questions (time unit, event coding, endpoint name, reference category) *before* generating any survival-analysis code.
- Produce copy-paste-ready R code using tidyverse-style, manuscript-quality idioms.
- Refuse silently-wrong analyses — flag PH violations, propose RMST when curves cross, use reverse Kaplan–Meier for median follow-up, label KM y-axes with the actual endpoint (not "Cumulative survival").
- Inspect Excel/CSV files before analysis (variable types, missingness, factor levels, event-coding sanity).
- Review and refactor existing R scripts conservatively.

## What's inside

```
.claude/skills/r-oncology-stats/
├── SKILL.md                              # the main skill — iron rules, routing, protocol
├── references/                           # methodology deep-dives Claude loads on demand
│   ├── 01-clarifying-questions.md
│   ├── 02-survival-analysis.md
│   ├── 03-ipd-reconstruction.md
│   ├── 04-baseline-tables.md
│   ├── 05-manuscript-figures.md
│   ├── 06-manuscript-tables.md
│   ├── 07-data-inspection.md
│   ├── 08-code-review-checklist.md
│   └── 09-package-quickref.md
├── scripts/                              # paste-and-adapt R templates
│   ├── km_curve.R
│   ├── cox_univariable.R
│   ├── cox_multivariable.R
│   ├── rmst.R
│   ├── reverse_km_followup.R
│   ├── ipd_from_km.R
│   ├── baseline_table1.R
│   └── forest_plot.R
└── templates/
    └── analysis_skeleton.R               # reproducible-project skeleton
```

## Preferred R packages

Survival: `survival`, `survminer`, `ggsurvfit`, `prodlim`, `survRM2`, `IPDfromKM`
Tables: `gtsummary`, `flextable`, `gt`, `officer`, `broom`, `broom.helpers`, `forestmodel`
Wrangling: `dplyr`, `tidyr`, `forcats`, `stringr`, `janitor`
I/O & paths: `readr`, `readxl`, `here`
Plots: `ggplot2`
Meta: `meta`, `metafor`

See `.claude/skills/r-oncology-stats/references/09-package-quickref.md` for one-line usage notes.

## Install

### Option 1 — Plugin marketplace (recommended)

In Claude Code:

```
/plugin marketplace add allanatal/R_assistant
/plugin install r-assistant@r-assistant
```

### Option 2 — User-level skill (symlink)

```bash
git clone https://github.com/allanatal/R_assistant.git ~/R_assistant
ln -s ~/R_assistant/.claude/skills/r-oncology-stats ~/.claude/skills/r-oncology-stats
```

Then restart Claude Code. The skill should appear in the available-skills list as `r-oncology-stats`.

### Option 3 — Project-level skill

Drop the `.claude/` folder into a project root and the skill becomes available in any Claude Code session opened from that project.

## How to use

Once installed, start a Claude Code session and ask anything oncology-survival-related in R. The skill triggers on terms like *Kaplan-Meier*, *Cox model*, *RMST*, *forest plot*, *baseline table*, *IPD reconstruction*, or in Portuguese: *sobrevida*, *razão de risco*, *curva de sobrevida*.

Example prompts:

> "I have an Excel file with OS data for two arms. Make me a KM curve with risk table and log-rank p-value, manuscript-ready."

> "Fit a multivariable Cox model adjusting for age, ECOG, and stage. Run the PH check and tell me if I should switch to RMST."

> "Reconstruct IPD from this published KM curve — I have the digitized coordinates and the numbers-at-risk table."

> "Review this R script and tell me what's wrong with the survival analysis."

Claude will start by confirming the time unit, event coding, endpoint, and reference category — then generate the R code.

## Status

Version `0.1.0` — focused on **frequentist** survival analysis, baseline comparisons, IPD reconstruction, and manuscript outputs.

Planned for later versions:

- Bayesian survival models (`brms` + Weibull/lognormal).
- Competing-risks analysis (`tidycmprsk`, Fine-Gray).
- Prediction-model workflows (`rms`, calibration).

## License

MIT — see [LICENSE](LICENSE).

## Acknowledgments

Skill structure inspired by [`ab604/claude-code-r-skills`](https://github.com/ab604/claude-code-r-skills) (general-purpose R skill bundle). This repo extends the pattern with clinical-oncology-specific content.
