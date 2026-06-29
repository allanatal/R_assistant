# Clarifying-questions protocol

Before any new survival analysis, ask these as **one grouped message**. Do not drip-feed. Mark plausible defaults in `[brackets]` so the user can answer fast.

If the user already provided answers in the prompt or in a dataset you've inspected, do not re-ask — preface with a one-line "I'm assuming X, Y, Z — say if anything's off" and proceed.

---

## The grouped checklist

Adapt the wording. Always cover these six categories:

```
Before I write the analysis I need to confirm a few things — answer
whichever aren't already obvious from the data:

  1. Time-to-event variable:  [name in your dataset, e.g. `os_months`]
     Time unit:               [days / months / years]

  2. Event/status variable:   [name in your dataset, e.g. `os_event`]
     Event coding:            [1 = event, 0 = censored — confirm]
     If coded differently (e.g. 0/1, 1/2, "Alive"/"Dead"), I'll recode.

  3. Endpoint name:           [Overall survival / Progression-free /
                               Disease-free / Recurrence-free /
                               Event-free / Time-to-failure / other]

  4. Group / exposure variable: [e.g. treatment_arm]
     Reference level:           [e.g. "Control" or "Standard"]
     If not stated, I'll use the first alphabetical level — please
     confirm or override with forcats::fct_relevel().

  5. Censoring definition: How are patients censored?
     [end of follow-up / lost to follow-up / death from other cause /
     administrative cutoff date / competing event treated as censoring]
     If competing risks matter (e.g. non-cancer death in elderly cohort),
     say so — Cox-on-cause-specific-hazard vs Fine-Gray is a different model.

  6. What outputs do you want?
     [KM curve with risk table / HR with 95% CI / median survival with CI /
     landmark survival at 1/2/5 years / RMST / median follow-up / all of these]
     Target: [exploration / manuscript figure / abstract / slide]
```

For follow-up requests on the same dataset, you don't need the full checklist again — just confirm what's new.

---

## Decisions that often come up after the checklist

These usually emerge AFTER the user answers; ask them once relevant:

- **Reverse-KM follow-up.** Always recommend reporting median follow-up via reverse KM. Ask if they want it overall, by group, or both.
- **Proportional hazards check.** For any Cox model, plan to run `cox.zph()` and a Schoenfeld plot. Ask if they want it inline or saved as a supplementary figure.
- **RMST.** If you suspect PH violation OR the user mentions crossing curves OR the endpoint is plateauing (e.g., immunotherapy long-term tail), proactively suggest RMST alongside the Cox HR.
- **Multivariable adjustment.** If the user wants a multivariable Cox, ask which covariates — clinically chosen, not automated. Discourage stepwise selection unless they justify it.
- **Subgroups.** If subgroups are mentioned, ask whether they are pre-specified or exploratory, and whether they want interaction tests (recommended) or just forest plots.

---

## What to skip the checklist for

- Simple data inspection / cleaning ("what's in this file?") — go straight to `references/07-data-inspection.md`.
- Code review / refactor of existing analysis — go to `references/08-code-review-checklist.md`.
- Already-fit-model post-processing (e.g., "make this KM prettier", "export this table to Word") — proceed.
- IPD reconstruction from a published curve — different checklist lives in `references/03-ipd-reconstruction.md`.

---

## Anti-patterns

- ❌ Asking the questions one at a time across multiple messages.
- ❌ Assuming `status = 1` means event without checking — many EHR exports use 1=alive / 2=dead, or text strings.
- ❌ Picking the reference category alphabetically when a clinically meaningful one exists (always reference is "Control" / "Standard of care" / "Wild-type", NOT alphabetically first).
- ❌ Defaulting time unit to "months" without asking — a Cox model is unit-agnostic for the HR but the figure axis and median survival are not.
