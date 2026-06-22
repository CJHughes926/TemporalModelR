#' Build Temporal GLM Models Across Cross-Validation Folds
#'
#' Modeling function that constructs binomial generalized linear models (GLMs)
#' for each cross-validation fold using presence and pseudoabsence data. Each
#' model reserves one fold as testing data and uses the remaining folds as
#' training data. The user supplies the model formula directly, giving full
#' control over predictor terms, polynomials, and interactions. The link
#' function can be set to logit, probit, complementary log-log, or cauchit.
#' Supports automatic or manual probability thresholding for converting
#' continuous predictions to binary suitability classifications necessary for
#' downstream analyses. The returned object follows the same
#' structure as \code{\link{build_temporal_hv}}, \code{\link{build_temporal_gam}},
#' and \code{\link{build_temporal_rf}}, and is accepted directly by
#' \code{\link{generate_spatiotemporal_predictions}}.
#'
#' @usage
#' build_temporal_glm(partition_result, pseudoabsence_result, model_formula,
#'           link = "logit", threshold_method = "tss",
#'           output_dir = file.path(tempdir(), "GLM_Models"),
#'           create_plot = TRUE, plot_palette = "Dark 2", overwrite = FALSE,
#'           time_cols = NULL, verbose = TRUE)
#'
#' @param partition_result List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to an \code{.rds}
#'   file containing that output.
#' @param pseudoabsence_result List or character. Output from
#'   \code{\link{generate_absences}} or path to an \code{.rds} file
#'   containing that output.
#' @param model_formula Formula or character. The right-hand side of the model
#'   formula supplied as either a formula object or a character string. The
#'   response variable (\code{presence}) is always added automatically on the
#'   left-hand side, so only the right-hand side needs to be provided. Both of
#'   the following are accepted and equivalent:
#'   \itemize{
#'     \item \code{~ Var1 + Var2 + I(Var1^2)}
#'     \item \code{"~ Var1 + Var2 + I(Var1^2)"}
#'   }
#'   Standard R formula syntax applies: \code{+} for additive terms, \code{*}
#'   for main effects plus interaction, \code{:} for interaction only,
#'   \code{I()} for arithmetic transformations, \code{poly()} for orthogonal
#'   polynomials, \code{log()}, \code{sqrt()}, and any other base R function
#'   that can appear in a formula. All predictor names referenced in the formula
#'   must be present as columns in both the presence and pseudoabsence data.
#' @param link Character. The link function for the binomial GLM. One of
#'   \code{"logit"} (default), \code{"probit"}, \code{"cloglog"}, or
#'   \code{"cauchit"}. See \code{\link[stats]{binomial}} for details on each
#'   link function.
#' @param threshold_method Character or numeric. Method used to convert
#'   continuous predicted probabilities to binary suitability. Accepted values:
#'   \itemize{
#'     \item \code{"prevalence"}: Sets threshold equal to the prevalence
#'       (proportion of presences) in the training data for that fold.
#'     \item \code{"tss"}: Selects the threshold that maximizes the True Skill
#'       Statistic (sensitivity + specificity - 1) on the training data.
#'       Default.
#'     \item A numeric value between 0 and 1 (e.g. \code{0.4}): Uses that
#'       value as a fixed threshold for all folds directly.
#'   }
#' @param output_dir Character. Directory to write output files including saved
#'   model objects and plots. Default is \code{file.path(tempdir(), "GLM_Models")}.
#' @param create_plot Logical. If \code{TRUE}, generates per-fold response
#'   curve plots and a combined ROC curve summary. Default is \code{TRUE}.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing saved model
#'   files. If \code{FALSE}, loads existing files when available. Default is
#'   \code{FALSE}.
#' @param plot_palette Character. Name of an HCL or RColorBrewer palette used
#'   to color folds in diagnostic plots. Accepts any HCL palette name (see
#'   \code{\link[grDevices]{hcl.pals}}) or, if \pkg{RColorBrewer} is installed,
#'   any Brewer palette name. Default is \code{"Dark 2"}.
#' @param time_cols Character. Name of the column(s) containing year or time
#'   step values in the occurrence data. Must match \code{time_cols} used in
#'   \code{\link{spatiotemporal_partition}}. Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes per-fold training summaries and
#'   file-saved messages. The completion summary and metrics table are
#'   always printed regardless of this setting.
#'
#' @return A list with class \code{"TemporalGLM"} containing:
#' \itemize{
#'   \item \code{models}: Named list of fitted \code{glm} objects, one per
#'     fold.
#'   \item \code{thresholds}: Named numeric vector of probability thresholds
#'     used for binary classification, one per fold.
#'   \item \code{threshold_method}: Character string recording the thresholding
#'     method used.
#'   \item \code{model_formula}: The formula object as passed to the fitting
#'     function.
#'   \item \code{link}: Character string recording the link function used.
#'   \item \code{model_vars}: Character vector of predictor names extracted
#'     from the formula right-hand side.
#'   \item \code{fold_training_data}: Named list of training data frames used
#'     to fit each fold model, retained for downstream prediction.
#'   \item \code{fold_test_metrics}: Data frame of held-out test fold metrics
#'     per fold: \code{Threshold}, \code{AUC}, \code{TSS}, \code{Kappa},
#'     \code{Sensitivity}, and \code{Specificity}. Also written to
#'     \code{Fold_Test_Metrics.csv} in \code{output_dir}.
#'   \item \code{output_dir}: Path to the output directory.
#'   \item \code{model_type}: Character string \code{"glm"}, used by
#'     \code{\link{generate_spatiotemporal_predictions}}.
#'   \item \code{plots}: Named list of recorded plot objects when
#'     \code{create_plot = TRUE}. Plots can be replayed with
#'     \code{grDevices::replayPlot()}.
#' }
#'
#' @details
#'
#' The \code{model_formula} argument accepts any standard R formula right-hand
#' side. The response (\code{presence}) is prepended automatically. All R
#' formula operators are valid, including \code{I()}, \code{poly()},
#' \code{log()}, \code{sqrt()}, \code{:}, and \code{*}. Variable names must
#' match column names in the data exactly. Predictor names for response curve
#' plots are extracted via \code{all.vars()}, which correctly unwraps terms
#' such as \code{I(Var1^2)} to the base variable \code{Var1}.
#'
#' All models are fit as \code{stats::glm(..., family = binomial(link = link))}.
#' Predicted values are always probabilities on the 0-1 scale.
#'
#' The returned object is recognized by
#' \code{\link{generate_spatiotemporal_predictions}}, which uses the
#' \code{model_type} field to use the correct prediction and evaluation
#' logic.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_partition}},
#'   \code{\link{generate_absences}}
#'
#' Modeling: \code{\link{build_temporal_gam}}, \code{\link{build_temporal_rf}},
#'   \code{\link{build_temporal_hv}},
#'   \code{\link{generate_spatiotemporal_predictions}}
#'
#' @examples
#' data(tmr_partition, package = "TemporalModelR")
#'
#' data(tmr_absences,  package = "TemporalModelR")
#'
#' build_temporal_glm(
#'   partition_result     = tmr_partition,
#'   pseudoabsence_result = tmr_absences,
#'   model_formula        = ~ elevation + forest_cover + prseas,
#'   threshold_method     = "tss",
#'   output_dir           = tempdir(),
#'   create_plot          = FALSE,
#'   time_cols            = c("year", "season"),
#'   verbose              = FALSE
#' )

#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom stats glm binomial predict var complete.cases as.formula
#' @importFrom grDevices adjustcolor dev.cur dev.off png recordPlot
#' @importFrom graphics plot lines abline legend par mtext points rug
#' @importFrom tools file_ext
#' @importFrom utils capture.output write.csv
build_temporal_glm <- function(partition_result,
                               pseudoabsence_result,
                               model_formula,
                               link              = "logit",
                               threshold_method  = "tss",
                               output_dir        = file.path(tempdir(), "GLM_Models"),
                               create_plot       = TRUE,
                               plot_palette      = "Dark 2",
                               overwrite         = FALSE,
                               time_cols         = NULL,
                               verbose           = TRUE) {

  if (missing(partition_result)) {
    stop(paste0("ERROR: 'partition_result' is required. Provide output from ",
                "spatiotemporal_partition() or a path to an .rds file."))
  }

  if (missing(pseudoabsence_result)) {
    stop(paste0("ERROR: 'pseudoabsence_result' is required. Provide output from ",
                "generate_absences() or a path to an .rds file."))
  }

  if (missing(model_formula)) {
    stop(paste0(
      "ERROR: 'model_formula' is required. Provide the right-hand side of the model formula as a ",
      "formula object or character string, e.g. ~ Var1 + Var2 + I(Var1^2)."
    ))
  }

  if (!link %in% c("logit", "probit", "cloglog", "cauchit")) {
    stop("ERROR: 'link' must be one of 'logit', 'probit', 'cloglog', or 'cauchit'.")
  }

  thr_resolved     <- .resolve_threshold_method(threshold_method)
  threshold_method <- thr_resolved$method
  manual_threshold <- thr_resolved$manual_value

  model_formula <- .parse_model_formula(model_formula)

  family_obj   <- stats::binomial(link = link)
  family_label <- paste0("binomial(", link, ")")

  partition_result     <- .load_partition_result(partition_result, verbose = verbose)
  pseudoabsence_result <- .load_pseudoabsence_result(pseudoabsence_result, verbose = verbose)

  rhs_vars <- all.vars(model_formula)

  prep             <- .prepare_combined_modeling_data(partition_result, pseudoabsence_result,
                                                      model_vars = rhs_vars, time_cols = time_cols,
                                                      verbose = verbose)
  combined_df      <- prep$combined_df
  all_folds        <- prep$all_folds
  single_fold_mode <- prep$single_fold_mode

  rhs_str      <- gsub("\\s+", " ", paste(deparse(model_formula[[2]]), collapse = " "))
  full_formula <- stats::as.formula(paste("presence ~", rhs_str))

  if (verbose) message(paste0("Building ", family_label, " model for ", length(all_folds), " fold(s)."))
  if (verbose) message(paste("Formula:", paste(deparse(full_formula), collapse = " ")))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  combined_file <- file.path(output_dir, paste0("all_models_binomial_", link, ".rds"))

  model_list      <- list()
  thresholds <- numeric(length(all_folds))
  names(thresholds) <- paste0("fold", all_folds)
  train_data_list <- list()
  metric_rows     <- list()

  if (file.exists(combined_file) && !overwrite) {

    if (verbose) message(paste("Loading existing model file:", basename(combined_file)))
    saved <- tryCatch(
      readRDS(combined_file),
      error = function(e) stop(paste0("ERROR loading saved models. Try overwrite = TRUE: ", e$message))
    )
    model_list      <- saved$models
    thresholds      <- saved$thresholds
    train_data_list <- saved$fold_training_data
    metric_rows     <- saved$fold_test_metrics_raw

  } else {

    for (fold in all_folds) {

      fold_key   <- paste0("fold", fold)
      train_data <- if (single_fold_mode) combined_df else combined_df[combined_df$fold != fold, , drop = FALSE]
      test_data  <- combined_df[combined_df$fold == fold, , drop = FALSE]

      n_pr_train <- sum(train_data$presence == 1, na.rm = TRUE)
      n_pa_train <- sum(train_data$presence == 0, na.rm = TRUE)
      n_pr_test  <- sum(test_data$presence == 1, na.rm = TRUE)
      n_pa_test  <- sum(test_data$presence == 0, na.rm = TRUE)

      if (verbose) message(paste0(
        "Fold ", fold, ": training on ", nrow(train_data), " points (",
        n_pr_train, " presences, ", n_pa_train, " pseudoabsences), testing on ",
        nrow(test_data), " (", n_pr_test, " presences, ", n_pa_test, " pseudoabsences)."
      ))

      if (n_pr_train < 5) {
        stop(paste0("ERROR: Fold ", fold, " has only ", n_pr_train, " training presences. Minimum 5 required."))
      }

      if (n_pa_train < 1) {
        stop(paste0("ERROR: Fold ", fold, " has no training pseudoabsences. Check generate_absences() output."))
      }

      var_check <- vapply(train_data[, rhs_vars, drop = FALSE],
                          function(x) stats::var(x, na.rm = TRUE), numeric(1))
      zero_var  <- names(var_check[!is.na(var_check) & var_check == 0])
      if (length(zero_var) > 0) {
        warning(paste0("Zero variance in Fold ", fold, " for: ", paste(zero_var, collapse = ", "),
                       ". These predictors will have no effect in this fold."))
      }

      model_fit <- tryCatch(
        stats::glm(full_formula, data = train_data, family = family_obj),
        error = function(e) stop(paste0("ERROR: glm() failed for fold ", fold, ": ", e$message))
      )

      thr <- {
        train_pred_prob_thr <- stats::predict(model_fit, newdata = train_data,
                                              type = "response")
        .select_threshold_prob(
          pred_prob    = train_pred_prob_thr,
          observed     = train_data$presence,
          method       = threshold_method,
          manual_value = manual_threshold
        )
      }
      thresholds[fold_key] <- thr

      model_list[[fold_key]]      <- model_fit
      train_data_list[[fold_key]] <- train_data

      if (nrow(test_data) > 0) {
        test_pred_prob <- stats::predict(model_fit, newdata = test_data, type = "response")
        test_pred_bin  <- as.integer(test_pred_prob >= thr)
        m <- .compute_confusion_metrics(
          observed  = test_data$presence,
          predicted = test_pred_bin,
          pred_prob = test_pred_prob
        )
      } else {
        m <- list(auc = NA, tss = NA, sensitivity = NA, specificity = NA,
                  kappa = NA, tp = NA, fn = NA, tn = NA, fp = NA)
      }

      metric_rows[[fold_key]] <- .build_fold_metric_row(fold, thr, m)
    }

    tryCatch({
      saveRDS(
        list(models = model_list, thresholds = thresholds,
             fold_training_data = train_data_list, fold_test_metrics_raw = metric_rows),
        combined_file
      )
      if (verbose) message(paste("Saved models to:", basename(combined_file)))
    }, error = function(e) {
      warning(paste0("Could not save model file: ", e$message))
    })
  }

  fold_metrics_df           <- do.call(rbind, metric_rows)
  rownames(fold_metrics_df) <- NULL

  metrics_file <- file.path(output_dir, "Fold_Test_Metrics.csv")
  tryCatch(
    utils::write.csv(fold_metrics_df, metrics_file, row.names = FALSE),
    error = function(e) warning(paste0("Could not save fold test metrics: ", e$message))
  )

  plot_list <- list()

  if (create_plot) {
    if (verbose) message("Generating diagnostic plots...")

    fold_colors <- .resolve_palette(length(all_folds), plot_palette)
    names(fold_colors) <- paste0("fold", all_folds)

    for (fold in all_folds) {
      fold_key   <- paste0("fold", fold)
      model_fit  <- model_list[[fold_key]]
      train_data <- train_data_list[[fold_key]]
      thr <- thresholds[fold_key]

      plot_file <- file.path(output_dir, paste0("ResponseCurves_Fold", fold, ".png"))
      n_vars    <- length(rhs_vars)
      n_cols    <- min(3, n_vars)
      n_rows    <- ceiling(n_vars / n_cols)
      var_means <- colMeans(train_data[, rhs_vars, drop = FALSE], na.rm = TRUE)

      tryCatch({
        grDevices::png(plot_file, width = 350 * n_cols, height = 320 * n_rows, res = 120)
        .plot_response_curves_parametric(model_fit, train_data, rhs_vars,
                                         var_means, fold_colors[fold_key], thr,
                                         fold, n_rows, n_cols, "GLM")
        grDevices::dev.off()
        if (verbose) message(paste("Saved response curve plot:", basename(plot_file)))
      }, error = function(e) {
        warning(paste0("Could not save response curve plot for Fold ", fold, ": ", e$message))
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })

      tryCatch({
        .plot_response_curves_parametric(model_fit, train_data, rhs_vars,
                                         var_means, fold_colors[fold_key], thr,
                                         fold, n_rows, n_cols, "GLM")
        plot_list[[paste0("response_curves_fold", fold)]] <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not record response curve plot for Fold ", fold, ": ", e$message))
      })
    }

    roc_file <- file.path(output_dir, "ROC_Summary.png")

    tryCatch({
      grDevices::png(roc_file, width = 600, height = 580, res = 120)
      .plot_roc_summary_parametric(model_list, train_data_list,
                                   thresholds, fold_metrics_df,
                                   all_folds, fold_colors, "GLM")
      grDevices::dev.off()
      if (verbose) message(paste("Saved ROC summary plot:", basename(roc_file)))
    }, error = function(e) {
      warning(paste0("Could not save ROC summary plot: ", e$message))
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    })

    tryCatch({
      .plot_roc_summary_parametric(model_list, train_data_list,
                                   thresholds, fold_metrics_df,
                                   all_folds, fold_colors, "GLM")
      plot_list[["roc_summary"]] <- grDevices::recordPlot()
    }, error = function(e) {
      warning(paste0("Could not record ROC plot: ", e$message))
    })
  }

  if (verbose) message("\nGLM modeling complete.")
  if (verbose) message(paste0("Link: ", link, " | Threshold: ",
                              threshold_method,
                              " | Folds: ", length(all_folds)))
  if (nrow(fold_metrics_df) > 0 && verbose) {
    message(paste("Metrics saved to:", basename(metrics_file)))
    message(paste(utils::capture.output(
      print(fold_metrics_df, row.names = FALSE)), collapse = "\n"))
  }

  if (create_plot && nrow(fold_metrics_df) > 0) {
    plot_rate_cols <- c("Sensitivity", "Specificity", "TSS", "Kappa", "AUC")
    plot_rate_cols <- plot_rate_cols[plot_rate_cols %in% names(fold_metrics_df) &
                                       !vapply(fold_metrics_df[plot_rate_cols],
                                               function(x) all(is.na(x)), logical(1))]
    if (length(plot_rate_cols) > 0) {
      bar_file <- file.path(output_dir, "Fold_Test_Metrics.png")
      tryCatch({
        grDevices::png(bar_file, width = 300 * length(plot_rate_cols) + 200,
                       height = 700, res = 150)
        .plot_fold_metrics_bars(
          metrics_df  = fold_metrics_df,
          rate_cols   = plot_rate_cols,
          fold_colors = fold_colors,
          title       = "GLM Time Independent Assessment Metrics"
        )
        grDevices::dev.off()
        if (verbose) message(paste("Saved fold metrics plot:", basename(bar_file)))
      }, error = function(e) {
        warning(paste0("Could not save fold metrics plot: ", e$message))
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })
      tryCatch({
        .plot_fold_metrics_bars(
          metrics_df  = fold_metrics_df,
          rate_cols   = plot_rate_cols,
          fold_colors = fold_colors,
          title       = "GLM Time Independent Assessment Metrics"
        )
        plot_list[["fold_test_metrics"]] <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not record fold metrics plot: ", e$message))
      })
    }
  }

  result <- list(
    models             = model_list,
    thresholds         = thresholds,
    threshold_method   = threshold_method,
    model_formula      = full_formula,
    link               = link,
    model_vars         = rhs_vars,
    fold_training_data = train_data_list,
    fold_test_metrics  = fold_metrics_df,
    output_dir         = output_dir,
    model_type         = "glm",
    plots              = plot_list
  )

  class(result) <- c("TemporalGLM", "list")
  invisible(result)
}
