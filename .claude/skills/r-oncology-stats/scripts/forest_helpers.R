## forest_helpers.R -- 7 helper functions for oncology forest plots.
##
## MIRRORED from references/11-forest-plots.md (paste-and-adapt reference).
## Sections in the reference:
##   plot_forest_oncology   ->  reference doc §3.4
##   nonest_terms           ->  reference doc §4.2
##   mark_nonest_tbl        ->  reference doc §4.4
##   build_subgroup_tidy    ->  reference doc §5.2
##   interaction_p          ->  reference doc §5.3
##   plot_forest_subgroup   ->  reference doc §5.4
##   save_forestmodel       ->  reference doc §6.3
##
## If you edit the doc, regenerate this file via /tmp/gen_forest_helpers.R.


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


interaction_p <- function(data, time_var, event_var, treat_var, subgroup_var) {
  f_add <- as.formula(sprintf("Surv(%s, %s) ~ %s + %s",
                              time_var, event_var, treat_var, subgroup_var))
  f_int <- as.formula(sprintf("Surv(%s, %s) ~ %s * %s",
                              time_var, event_var, treat_var, subgroup_var))
  m_add <- survival::coxph(f_add, data = data)
  m_int <- survival::coxph(f_int, data = data)
  anova(m_add, m_int, test = "LRT")[2, "Pr(>|Chi|)"]
}


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


