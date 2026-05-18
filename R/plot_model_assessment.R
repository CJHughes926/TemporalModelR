#' Visualize Model Assessment Metrics Across Time
#'
#' Generates diagnostic plots from the per-timestep assessment table produced
#' by \code{\link{generate_spatiotemporal_predictions}}, optionally overlaying
#' overall reference values from the model result object.
#'
#' @param predictions List returned by
#'   \code{\link{generate_spatiotemporal_predictions}}, or a named list with
#'   at least a \code{timestep_metrics} element. The \code{timestep_metrics}
#'   element may also be a path to a \code{Timestep_Assessment_Metrics.csv}
#'   file produced by \code{\link{generate_spatiotemporal_predictions}}.
#' @param time_column Character. Name of the primary time column in
#'   \code{timestep_metrics} to use as the x axis (e.g. \code{"year"}).
#'   When predictions span multiple time columns (e.g. \code{"year"} and
#'   \code{"season"}), provide all relevant column names as a character
#'   vector and control how secondary columns are handled via
#'   \code{secondary_time_mode}.
#' @param secondary_time_mode Character. How to handle secondary time columns
#'   when \code{time_column} has length > 1. One of:
#'   \itemize{
#'     \item \code{"combine"} (default): secondary time values are appended to
#'       the primary value to form a single ordered x-axis label
#'       (e.g. \code{1_Spring}, \code{1_Summer}, \code{2_Spring}, ...).
#'     \item \code{"facet"}: a separate plot is produced for each unique
#'       combination of secondary time values, with the primary time column
#'       as the x axis on every panel.
#'   }
#' @param model_result List or character. Optional. Output from a
#'   \code{build_temporal_*()} function or path to its \code{.rds} file. When
#'   supplied, overall sensitivity and specificity from
#'   \code{model_result$fold_test_metrics} are added as per-fold reference
#'   lines. Default is \code{NULL}.
#' @param cbp_threshold Numeric. Significance threshold for CBP. Default is
#'   \code{0.05}.
#' @param plot_palette Character. Name of an HCL or RColorBrewer palette used
#'   to color folds in diagnostic plots. Accepts any HCL palette name (see
#'   \code{\link[grDevices]{hcl.pals}}) or, if \pkg{RColorBrewer} is installed,
#'   any Brewer palette name. Default is \code{"Dark 2"}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing.
#'
#' @return Invisibly returns a named list containing:
#' \itemize{
#'   \item \code{pct_suitable}: Recorded plot of proportion of study area
#'     predicted suitable per time step.
#'   \item \code{sensitivity}: Recorded plot of per-timestep sensitivity.
#'   \item \code{specificity}: Recorded plot of per-timestep specificity
#'     (only present when pseudoabsence data were used).
#'   \item \code{cbp}: Recorded plot of cumulative binomial probability per
#'     time step on a log scale.
#'   \item \code{tp_fn}: Recorded plot of true positives and false negatives
#'     per time step.
#'   \item \code{tn_fp}: Recorded plot of true negatives and false positives
#'     per time step (only present when pseudoabsence data were used).
#'   \item \code{timestep_summary}: Data frame of per-time-step cross-fold
#'     mean and SD for each metric.
#'   \item \code{overall_summary}: Data frame from
#'     \code{predictions$overall_summary}, when present.
#' }
#'
#'@details
#' Plots per-fold and per-timestep diagnostic plots for data produced by
#' \code{\link{generate_spatiotemporal_predictions}}. These quick visuals can
#' be used by users to assess model performance and significance and decide
#' if the model's performance warrants further interpretation of the results
#' through postprocessing analyses.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_partition}},
#'   \code{\link{generate_absences}}
#'
#' Modeling: \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_rf}}, \code{\link{build_temporal_hv}},
#'
#' Postprocessing: \code{\link{generate_spatiotemporal_predictions}},
#'   \code{\link{summarize_raster_outputs}}
#'
#' @examples
#' \donttest{
#'   data(tmr_predictions, package = "TemporalModelR")
#'
#'   plot_model_assessment(
#'     predictions         = tmr_predictions,
#'     time_column         = c("year", "season"),
#'     secondary_time_mode = "combine",
#'     verbose             = FALSE
#'   )
#' }

#' @export
#' @importFrom graphics abline axis legend lines mtext par plot points rect
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom stats sd
#' @importFrom utils capture.output read.csv
plot_model_assessment <- function(predictions,
                                  time_column,
                                  secondary_time_mode = "combine",
                                  model_result  = NULL,
                                  cbp_threshold = 0.05,
                                  plot_palette  = "Dark 2",
                                  verbose       = TRUE) {

  if (missing(predictions)) {
    stop(paste0("ERROR: 'predictions' is required. Please provide the list ",
                "output from generate_spatiotemporal_predictions()."))
  }
  if (missing(time_column) || !is.character(time_column) || length(time_column) == 0) {
    stop(paste0("ERROR: 'time_column' must be a character string (or vector of strings) ",
                "naming the time column(s) in timestep_metrics."))
  }

  if (is.character(predictions) && length(predictions) == 1) {
    if (file.exists(predictions)) {
      if (verbose) message("Reading timestep metrics from CSV: ", basename(predictions))
      an <- tryCatch(
        utils::read.csv(predictions, stringsAsFactors = FALSE),
        error = function(e) stop(paste0("ERROR reading 'predictions' CSV: ", e$message))
      )
      ov_sum <- NULL
    } else {
      stop(paste0("ERROR: 'predictions' CSV file not found: ", predictions))
    }
  } else if (is.data.frame(predictions)) {
    an <- predictions; ov_sum <- NULL
  } else if (is.list(predictions)) {
    if (!"timestep_metrics" %in% names(predictions))
      stop("ERROR: 'predictions' must contain 'timestep_metrics'.")
    tm <- predictions$timestep_metrics
    if (is.data.frame(tm)) {
      an <- tm
    } else if (is.character(tm) && file.exists(tm)) {
      an <- tryCatch(utils::read.csv(tm, stringsAsFactors = FALSE),
                     error = function(e) stop(paste0("ERROR reading 'timestep_metrics': ", e$message)))
    } else {
      stop(paste0("ERROR: 'timestep_metrics' file not found: ", tm))
    }
    ov_sum <- if ("overall_summary" %in% names(predictions) &&
                  !is.null(predictions$overall_summary) &&
                  nrow(predictions$overall_summary) > 0)
      predictions$overall_summary else NULL
  } else stop("ERROR: 'predictions' must be the list from generate_spatiotemporal_predictions(), a data frame, or a path to a Timestep_Assessment_Metrics CSV.")

  if (nrow(an) == 0) stop("ERROR: timestep_metrics is empty.")

  an_aliases <- list(
    Pct_Suitable = c("Pct_Suitable", "G_Volume"),
    TP           = c("TP",  "TP_G"),
    FN           = c("FN",  "FN_G"),
    TN           = c("TN",  "TN_G"),
    FP           = c("FP",  "FP_G"),
    Sensitivity  = c("Sensitivity",  "Sensitivity_G"),
    Specificity  = c("Specificity",  "Specificity_G"),
    CBP          = c("CBP", "CBP_G")
  )
  for (canon in names(an_aliases)) {
    if (!canon %in% names(an)) {
      for (alias in an_aliases[[canon]]) {
        if (alias %in% names(an)) { an[[canon]] <- an[[alias]]; break }
      }
    }
  }

  missing_cols <- setdiff(time_column, names(an))
  if (length(missing_cols) > 0)
    stop(paste0("ERROR: time_column value(s) not found in metrics: ",
                paste(missing_cols, collapse = ", "),
                ". Available columns: ", paste(names(an), collapse = ", ")))

  secondary_time_mode <- match.arg(secondary_time_mode, c("combine", "facet"))

  plot_time_col  <- time_column[1]
  secondary_cols <- if (length(time_column) > 1) time_column[-1] else character(0)

  secondary_combos <- if (length(secondary_cols) > 0) {
    unique(an[, secondary_cols, drop = FALSE])
  } else {
    NULL
  }

  if (length(secondary_cols) > 0 && secondary_time_mode == "combine") {

    an$..time_label.. <- apply(
      an[, time_column, drop = FALSE], 1,
      function(r) paste(r, collapse = "_")
    )
    all_labels  <- unique(an$..time_label..)
    label_order <- order(
      as.numeric(an[match(all_labels, an$..time_label..), plot_time_col]),
      as.character(an[match(all_labels, an$..time_label..), secondary_cols[1]])
    )
    an$..x_val..  <- match(an$..time_label.., all_labels[label_order])
    x_labels      <- all_labels[label_order]
    plot_time_col <- "..x_val.."
    secondary_cols <- character(0)
    attr(an, "x_labels") <- x_labels
  }

  has_spec <- "Specificity" %in% names(an) && !all(is.na(an$Specificity))
  has_tnfp <- all(c("TN","FP") %in% names(an)) &&
    !all(is.na(an$TN)) && !all(is.na(an$FP))

  fold_test <- NULL
  if (!is.null(model_result)) {
    if (is.character(model_result) && file.exists(model_result))
      model_result <- tryCatch(readRDS(model_result),
                               error = function(e) { warning(paste("Could not read model_result:", e$message)); NULL })
    if (is.list(model_result) && "fold_test_metrics" %in% names(model_result) &&
        nrow(model_result$fold_test_metrics) > 0) {
      fold_test <- model_result$fold_test_metrics
      if (verbose) message(paste("Loaded fold_test_metrics for", nrow(fold_test), "fold(s)."))
    }
  }

  all_folds   <- sort(unique(an$Fold))
  fold_colors <- .resolve_palette(length(all_folds), plot_palette)
  names(fold_colors) <- as.character(all_folds)

  col_tp <- "#0072B2"
  col_fn <- "#D55E00"
  col_tn <- "#009E73"
  col_fp <- "#CC79A7"
  col_cbp_good <- col_tp
  col_cbp_bad  <- col_fn
  col_mean <- "black"

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  if (secondary_time_mode == "facet" && !is.null(secondary_combos) &&
      nrow(secondary_combos) > 0) {
    orig_secondary_cols <- time_column[-1]
    primary_col         <- time_column[1]

    facet_results <- list()

    if ("Pct_Suitable" %in% names(an)) {
      .plot_facet_metric(an, "Pct_Suitable", "Prop. Suitable", "G-Space Suitability",
                         secondary_combos, primary_col, orig_secondary_cols,
                         fold_colors, col_mean, y_floor = 0, y_ceil = 1)
      facet_results$pct_suitable <- grDevices::recordPlot()
    }
    if ("Sensitivity" %in% names(an)) {
      .plot_facet_metric(an, "Sensitivity", "Sensitivity", "Sensitivity",
                         secondary_combos, primary_col, orig_secondary_cols,
                         fold_colors, col_mean, y_floor = 0, y_ceil = 1)
      facet_results$sensitivity <- grDevices::recordPlot()
    }
    if (has_spec && "Specificity" %in% names(an)) {
      .plot_facet_metric(an, "Specificity", "Specificity", "Specificity",
                         secondary_combos, primary_col, orig_secondary_cols,
                         fold_colors, col_mean, y_floor = 0, y_ceil = 1)
      facet_results$specificity <- grDevices::recordPlot()
    }
    if ("CBP" %in% names(an)) {
      .plot_facet_metric(an, "CBP", "CBP", "CBP",
                         secondary_combos, primary_col, orig_secondary_cols,
                         fold_colors, col_mean)
      facet_results$cbp <- grDevices::recordPlot()
    }
    facet_results$timestep_summary <- data.frame()
    n_plots <- length(facet_results) - 1
    if (verbose) message(paste0("Facet mode: produced ", n_plots, " stacked plot(s) across ",
                                nrow(secondary_combos), " ",
                                paste(orig_secondary_cols, collapse = "/"),
                                " value(s)."))
    return(invisible(facet_results))
  }

  agg_cols <- c("Pct_Suitable","Sensitivity","CBP","TP","FN")
  if (has_spec) agg_cols <- c(agg_cols, "Specificity")
  if (has_tnfp) agg_cols <- c(agg_cols, "TN","FP")
  agg_cols <- agg_cols[agg_cols %in% names(an)]

  time_vals <- sort(unique(an[[plot_time_col]]))
  an_summary_rows <- lapply(time_vals, function(tv) {
    grp <- an[!is.na(an[[plot_time_col]]) & an[[plot_time_col]] == tv, , drop = FALSE]
    row <- data.frame(tv, stringsAsFactors = FALSE); names(row) <- plot_time_col
    for (col in agg_cols) {
      row[[paste0("Mean_", col)]] <- mean(grp[[col]], na.rm = TRUE)
      row[[paste0("SD_",   col)]] <- stats::sd(grp[[col]], na.rm = TRUE)
    }
    row
  })
  an_summary <- do.call(rbind.data.frame, an_summary_rows)
  rownames(an_summary) <- NULL
  an_summary$CBP_good <- if ("Mean_CBP" %in% names(an_summary))
    !is.na(an_summary$Mean_CBP) & an_summary$Mean_CBP < cbp_threshold
  else rep(FALSE, nrow(an_summary))

  t_all <- an[[plot_time_col]]
  t_sum <- an_summary[[plot_time_col]]
  x_lim <- range(t_all, na.rm = TRUE)

  x_labels_attr <- attr(an, "x_labels")

  t_uniq <- sort(unique(t_all[!is.na(t_all)]))
  x_tck  <- if (length(t_uniq) <= 15) {
    t_uniq
  } else {
    t_uniq[seq(1, length(t_uniq), by = ceiling(length(t_uniq) / 15))]
  }

  x_axis_label <- if (!is.null(x_labels_attr)) {
    paste(time_column, collapse = " / ")
  } else {
    plot_time_col
  }

  plot_list <- list()

  if ("Pct_Suitable" %in% names(an) && "Mean_Pct_Suitable" %in% names(an_summary)) {

    ref_vals <- if (!is.null(ov_sum) && "Mean_Pct_Suitable" %in% names(ov_sum))
      ov_sum$Mean_Pct_Suitable else NULL

    .new_plot(an$Pct_Suitable, ref_vals,
              title = "G-Space Predictions per Time Step",
              ylab  = "Proportion of Study Area Predicted Suitable",
              x_lim = x_lim, x_axis_label = x_axis_label,
              x_tck = x_tck, x_labels_attr = x_labels_attr,
              y_floor = 0, y_ceil = 1)

    .ref_lines(ov_sum, "Mean_Pct_Suitable", all_folds, fold_colors)
    .tseries(t_all, an$Pct_Suitable, an$Fold, t_sum, an_summary$Mean_Pct_Suitable,
             all_folds, fold_colors, col_mean)

    leg_lab <- "Cross-fold mean"
    leg_col <- col_mean
    leg_lty <- 1
    leg_lwd <- 3
    if (!is.null(ov_sum) && "Mean_Pct_Suitable" %in% names(ov_sum)) {
      leg_lab <- c(leg_lab, "Overall mean (per fold)")
      leg_col <- c(leg_col, "gray50")
      leg_lty <- c(leg_lty, 2)
      leg_lwd <- c(leg_lwd, 1.5)
    }
    .line_legend(leg_lab, leg_col, leg_lty, leg_lwd)
    plot_list$pct_suitable <- grDevices::recordPlot()
  }

  if ("Sensitivity" %in% names(an) && "Mean_Sensitivity" %in% names(an_summary)) {

    ref_vals <- c(
      if (!is.null(ov_sum)    && "Overall_Sensitivity" %in% names(ov_sum))    ov_sum$Overall_Sensitivity,
      if (!is.null(fold_test) && "Sensitivity"          %in% names(fold_test)) fold_test$Sensitivity
    )

    .new_plot(an$Sensitivity, ref_vals,
              title = "Sensitivity per Time Step", ylab = "Sensitivity",
              x_lim = x_lim, x_axis_label = x_axis_label,
              x_tck = x_tck, x_labels_attr = x_labels_attr,
              y_floor = 0, y_ceil = 1)

    .ref_lines(ov_sum,    "Overall_Sensitivity", all_folds, fold_colors, lty = 2, lwd = 1.2)
    .ref_lines(fold_test, "Sensitivity",          all_folds, fold_colors, lty = 3, lwd = 1.5)
    .tseries(t_all, an$Sensitivity, an$Fold, t_sum, an_summary$Mean_Sensitivity,
             all_folds, fold_colors, col_mean)

    leg_lab <- "Cross-fold mean"
    leg_col <- col_mean
    leg_lty <- 1
    leg_lwd <- 3
    if (!is.null(ov_sum) && "Overall_Sensitivity" %in% names(ov_sum)) {
      leg_lab <- c(leg_lab, "Overall sensitivity (per fold)")
      leg_col <- c(leg_col, "gray50")
      leg_lty <- c(leg_lty, 2)
      leg_lwd <- c(leg_lwd, 1.5)
    }
    if (!is.null(fold_test) && "Sensitivity" %in% names(fold_test)) {
      leg_lab <- c(leg_lab, "Model sensitivity incl. pseudoabsences (per fold)")
      leg_col <- c(leg_col, "gray50")
      leg_lty <- c(leg_lty, 3)
      leg_lwd <- c(leg_lwd, 1.5)
    }
    .line_legend(leg_lab, leg_col, leg_lty, leg_lwd)
    plot_list$sensitivity <- grDevices::recordPlot()
  }

  if (has_spec && "Mean_Specificity" %in% names(an_summary)) {

    ref_vals <- c(
      if (!is.null(ov_sum)    && "Overall_Specificity" %in% names(ov_sum))    ov_sum$Overall_Specificity,
      if (!is.null(fold_test) && "Specificity"          %in% names(fold_test)) fold_test$Specificity
    )

    .new_plot(an$Specificity, ref_vals,
              title = "Specificity per Time Step", ylab = "Specificity",
              x_lim = x_lim, x_axis_label = x_axis_label,
              x_tck = x_tck, x_labels_attr = x_labels_attr,
              y_floor = 0, y_ceil = 1)

    .ref_lines(ov_sum,    "Overall_Specificity", all_folds, fold_colors, lty = 2, lwd = 1.2)
    .ref_lines(fold_test, "Specificity",          all_folds, fold_colors, lty = 3, lwd = 1.5)
    .tseries(t_all, an$Specificity, an$Fold, t_sum, an_summary$Mean_Specificity,
             all_folds, fold_colors, col_mean)

    leg_lab <- "Cross-fold mean"
    leg_col <- col_mean
    leg_lty <- 1
    leg_lwd <- 3
    if (!is.null(ov_sum) && "Overall_Specificity" %in% names(ov_sum)) {
      leg_lab <- c(leg_lab, "Overall specificity (per fold)")
      leg_col <- c(leg_col, "gray50")
      leg_lty <- c(leg_lty, 2)
      leg_lwd <- c(leg_lwd, 1.5)
    }
    if (!is.null(fold_test) && "Specificity" %in% names(fold_test)) {
      leg_lab <- c(leg_lab, "Model specificity incl. pseudoabsences (per fold)")
      leg_col <- c(leg_col, "gray50")
      leg_lty <- c(leg_lty, 3)
      leg_lwd <- c(leg_lwd, 1.5)
    }
    .line_legend(leg_lab, leg_col, leg_lty, leg_lwd)
    plot_list$specificity <- grDevices::recordPlot()
  }

  if ("CBP" %in% names(an) && "Mean_CBP" %in% names(an_summary)) {

    all_cbp <- c(an$CBP,
                 if (!is.null(ov_sum) && "Overall_CBP" %in% names(ov_sum)) ov_sum$Overall_CBP)
    pos_cbp <- all_cbp[!is.na(all_cbp) & all_cbp > 0]

    if (length(pos_cbp) > 0) {
      lo   <- min(pos_cbp); hi <- max(pos_cbp)
      lim  <- c(min(lo, cbp_threshold) * 0.8, max(hi, cbp_threshold) * 1.5)
      lg   <- floor(log10(lim[1])):ceiling(log10(lim[2]))
      brks <- 10^lg

      bot_mar <- if (!is.null(x_labels_attr)) 6 else 4
      graphics::par(mar = c(bot_mar, 5, 3.5, 12), xpd = FALSE)
      graphics::plot(NULL, xlim = x_lim, ylim = lim, log = "y",
                     xlab = x_axis_label, ylab = "CBP (log scale)",
                     main = "Cumulative Binomial Probability per Time Step",
                     las = 1, xaxt = "n", yaxt = "n")
      .draw_xaxis(x_tck, x_labels_attr)

      log_brks <- brks[brks > 0]
      log_labs <- vapply(log_brks, function(b) {
        if (b >= 0.001) formatC(b, format = "f", digits = max(0, -floor(log10(b))))
        else            formatC(b, format = "e", digits = 0)
      }, character(1))
      graphics::axis(2, at = log_brks, labels = log_labs, las = 1, cex.axis = 0.85)

      graphics::abline(h = cbp_threshold, lty = 1, col = "black", lwd = 1.8)

      .ref_lines(ov_sum, "Overall_CBP", all_folds, fold_colors, lty = 4, lwd = 1.0)
      .tseries(t_all, an$CBP, an$Fold, t_sum, an_summary$Mean_CBP,
               all_folds, fold_colors, col_mean,
               col_cbp_good = col_cbp_good, col_cbp_bad = col_cbp_bad,
               cbp_good_sum = an_summary$CBP_good)

      leg_lab <- c(paste0("p < ", cbp_threshold), paste0("p >= ", cbp_threshold), "Threshold")
      leg_col <- c(col_cbp_good, col_cbp_bad, "black")
      leg_lty <- c(1, 1, 1)
      leg_lwd <- c(3, 3, 1.8)
      if (!is.null(ov_sum) && "Overall_CBP" %in% names(ov_sum)) {
        leg_lab <- c(leg_lab, "Overall CBP (per fold)")
        leg_col <- c(leg_col, "gray50")
        leg_lty <- c(leg_lty, 4)
        leg_lwd <- c(leg_lwd, 1.2)
      }
      .line_legend(leg_lab, leg_col, leg_lty, leg_lwd)
      plot_list$cbp <- grDevices::recordPlot()
    }
  }

  if (all(c("TP","FN","Mean_TP","Mean_FN") %in% c(names(an), names(an_summary)))) {

    tp_all <- an$TP; fn_all <- an$FN
    tp_sum <- an_summary$Mean_TP; fn_sum <- an_summary$Mean_FN
    valid  <- !is.na(tp_all) & !is.na(fn_all) & (tp_all > 0 | fn_all > 0)

    if (any(valid)) {
      y_max <- max(c(tp_all[valid], tp_sum), na.rm = TRUE) * 1.15
      y_min <- -max(c(fn_all[valid], fn_sum), na.rm = TRUE) * 1.15

      bot_mar <- if (!is.null(x_labels_attr)) 6 else 4
      graphics::par(mar = c(bot_mar, 5, 3.5, 12), xpd = FALSE)
      graphics::plot(NULL, xlim = x_lim, ylim = c(y_min, y_max),
                     xlab = x_axis_label, ylab = "Count  (TP above 0 | FN below 0)",
                     main = "True Positives and False Negatives per Time Step",
                     las = 1, xaxt = "n")
      .draw_xaxis(x_tck, x_labels_attr)
      graphics::abline(h = 0, col = "gray40", lwd = 0.8)

      t_uniq <- sort(unique(t_sum))
      bw <- if (length(t_uniq) > 1) diff(range(t_uniq)) / (length(t_uniq) * 2.2) else 0.4

      for (ti in seq_along(t_sum)) {
        tx <- t_sum[ti]
        if (!is.na(tp_sum[ti])) graphics::rect(tx - bw, 0, tx + bw,  tp_sum[ti],
                                               col = grDevices::adjustcolor(col_tp, 0.7), border = "gray30", lwd = 0.4)
        if (!is.na(fn_sum[ti])) graphics::rect(tx - bw, 0, tx + bw, -fn_sum[ti],
                                               col = grDevices::adjustcolor(col_fn, 0.7), border = "gray30", lwd = 0.4)
      }

      pt_tp <- grDevices::adjustcolor(col_tp, 0.7)
      pt_fn <- grDevices::adjustcolor(col_fn, 0.7)
      graphics::points(t_all[valid],  tp_all[valid], pch = 21, bg = pt_tp, cex = 0.9, col = "gray30")
      graphics::points(t_all[valid], -fn_all[valid], pch = 21, bg = pt_fn, cex = 0.9, col = "gray30")

      .bar_legend(
        bar_labs = c("Mean TP", "Mean FN"),
        bar_cols = c(pt_tp, pt_fn),
        pt_lab   = "Individual folds",
        pt_col   = grDevices::adjustcolor("gray30", 0.6)
      )

      plot_list$tp_fn <- grDevices::recordPlot()
    }
  }

  if (has_tnfp && all(c("Mean_TN","Mean_FP") %in% names(an_summary))) {

    tn_all <- an$TN; fp_all <- an$FP
    tn_sum <- an_summary$Mean_TN; fp_sum <- an_summary$Mean_FP
    valid  <- !is.na(tn_all) & !is.na(fp_all) & (tn_all > 0 | fp_all > 0)

    if (any(valid)) {
      y_max <- max(c(tn_all[valid], tn_sum), na.rm = TRUE) * 1.15
      y_min <- -max(c(fp_all[valid], fp_sum), na.rm = TRUE) * 1.15

      bot_mar <- if (!is.null(x_labels_attr)) 6 else 4
      graphics::par(mar = c(bot_mar, 5, 3.5, 12), xpd = FALSE)
      graphics::plot(NULL, xlim = x_lim, ylim = c(y_min, y_max),
                     xlab = x_axis_label, ylab = "Count  (TN above 0 | FP below 0)",
                     main = "True Negatives and False Positives per Time Step",
                     las = 1, xaxt = "n")
      .draw_xaxis(x_tck, x_labels_attr)
      graphics::abline(h = 0, col = "gray40", lwd = 0.8)

      t_uniq <- sort(unique(t_sum))
      bw <- if (length(t_uniq) > 1) diff(range(t_uniq)) / (length(t_uniq) * 2.2) else 0.4

      for (ti in seq_along(t_sum)) {
        tx <- t_sum[ti]
        if (!is.na(tn_sum[ti])) graphics::rect(tx - bw, 0, tx + bw,  tn_sum[ti],
                                               col = grDevices::adjustcolor(col_tn, 0.7), border = "gray30", lwd = 0.4)
        if (!is.na(fp_sum[ti])) graphics::rect(tx - bw, 0, tx + bw, -fp_sum[ti],
                                               col = grDevices::adjustcolor(col_fp, 0.7), border = "gray30", lwd = 0.4)
      }

      pt_tn <- grDevices::adjustcolor(col_tn, 0.7)
      pt_fp <- grDevices::adjustcolor(col_fp, 0.7)
      graphics::points(t_all[valid],  tn_all[valid], pch = 21, bg = pt_tn, cex = 0.9, col = "gray30")
      graphics::points(t_all[valid], -fp_all[valid], pch = 21, bg = pt_fp, cex = 0.9, col = "gray30")

      .bar_legend(
        bar_labs = c("Mean TN", "Mean FP"),
        bar_cols = c(pt_tn, pt_fp),
        pt_lab   = "Individual folds",
        pt_col   = grDevices::adjustcolor("gray30", 0.6)
      )

      plot_list$tn_fp <- grDevices::recordPlot()
    }
  }

  summary_rows <- list()
  if ("Pct_Suitable" %in% names(an)) {
    summary_rows[["Pct_Suitable"]] <- data.frame(
      Metric = "Pct_Suitable",
      Mean   = round(mean(an$Pct_Suitable, na.rm = TRUE), 4),
      SD     = round(stats::sd(an$Pct_Suitable, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }
  if ("Sensitivity" %in% names(an)) {
    summary_rows[["Sensitivity"]] <- data.frame(
      Metric = "Sensitivity",
      Mean   = round(mean(an$Sensitivity, na.rm = TRUE), 4),
      SD     = round(stats::sd(an$Sensitivity, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }
  if (has_spec) {
    summary_rows[["Specificity"]] <- data.frame(
      Metric = "Specificity",
      Mean   = round(mean(an$Specificity, na.rm = TRUE), 4),
      SD     = round(stats::sd(an$Specificity, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }
  if ("CBP" %in% names(an)) {
    pct_good <- round(mean(an$CBP < cbp_threshold, na.rm = TRUE) * 100, 1)
    summary_rows[["CBP"]] <- data.frame(
      Metric = paste0("CBP < ", cbp_threshold, " (%)"),
      Mean   = pct_good,
      SD     = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (length(summary_rows) > 0) {
    summary_df <- do.call(rbind, summary_rows)
    rownames(summary_df) <- NULL
    if (verbose) {
      message("\nTimestep assessment summary.")
      message(paste(utils::capture.output(
        print(summary_df, row.names = FALSE)), collapse = "\n"))
    }
  }

  graphics::par(opar)
  plot_list$timestep_summary <- an_summary
  if (!is.null(ov_sum)) plot_list$overall_summary <- ov_sum
  invisible(plot_list)
}
