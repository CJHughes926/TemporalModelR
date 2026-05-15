#' Build Temporal Random Forest Models Across Cross-Validation Folds
#'
#' Modeling function that constructs Random Forest classification models for
#' each cross-validation fold using presence and pseudoabsence data. Each model
#' reserves one fold as testing data and uses the remaining folds as training
#' data. The user specifies predictors as a character vector.
#' Predicted probabilities of presence are extracted from the
#' out-of-bag or in-bag vote fractions and thresholded to produce binary
#' suitability classifications. Variable importance is recorded for each fold.
#' The returned object follows the same
#' structure as \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#' and \code{\link{build_temporal_hv}}, and is accepted directly by
#' \code{\link{generate_spatiotemporal_predictions}}.
#'
#' @param partition_result List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to an \code{.rds}
#'   file containing that output.
#' @param pseudoabsence_result List or character. Output from
#'   \code{\link{generate_absences}} or path to an \code{.rds} file
#'   containing that output.
#' @param model_vars Character vector. Names of predictor columns to include in
#'   the Random Forest. All variables must be present as columns in both the
#'   presence and pseudoabsence data.
#' @param rf_params Named list. Additional arguments passed to
#'   \code{\link[randomForest]{randomForest}}, such as \code{ntree} (number of
#'   trees, default 500), \code{mtry} (number of variables tried at each split,
#'   default \code{floor(sqrt(length(model_vars)))}), and \code{nodesize}
#'   (minimum node size, default 1 for classification). Default is an empty
#'   list, which uses \pkg{randomForest} defaults.
#' @param threshold_method Character or numeric. Method used to convert
#'   continuous predicted probabilities to binary suitability. Accepted values:
#'   \itemize{
#'     \item \code{"prevalence"}: Sets threshold equal to the prevalence
#'       (proportion of presences) in the training data for that fold.
#'     \item \code{"tss"}: Selects the threshold that maximises the True Skill
#'       Statistic (sensitivity + specificity - 1) on the training data.
#'       Default.
#'     \item A numeric value between 0 and 1 (e.g. \code{0.4}): Uses that
#'       value as a fixed threshold for all folds directly.
#'   }
#' @param output_dir Character. Directory to write output files including saved
#'   model objects and plots. Default is \code{file.path(tempdir(), "RF_Models")}.
#' @param create_plot Logical. If \code{TRUE}, generates a per-fold variable
#'   importance plot, partial dependence curves for each predictor, and a
#'   combined ROC curve summary. Default is \code{TRUE}.
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
#'   file-saved messages.
#'
#' @return A list with class \code{"TemporalRF"} containing:
#' \itemize{
#'   \item \code{models}: Named list of fitted \code{randomForest} objects,
#'     one per fold.
#'   \item \code{thresholds}: Named numeric vector of probability thresholds
#'     used for binary classification, one per fold.
#'   \item \code{threshold_method}: Character string recording the thresholding
#'     method used.
#'   \item \code{model_vars}: Character vector of predictor names used.
#'   \item \code{variable_importance}: Named list of importance data frames,
#'     one per fold, with mean decrease in accuracy for each predictor.
#'   \item \code{fold_training_data}: Named list of training data frames used
#'     to fit each fold model, retained for downstream prediction.
#'   \item \code{fold_test_metrics}: Data frame of held-out test fold metrics
#'     per fold: \code{Threshold}, \code{AUC}, \code{TSS}, \code{Kappa},
#'     \code{Sensitivity}, and \code{Specificity}. Also written to
#'     \code{Fold_Test_Metrics.csv} in \code{output_dir}.
#'   \item \code{output_dir}: Path to the output directory.
#'   \item \code{model_type}: Character string \code{"rf"}, used by
#'     \code{\link{generate_spatiotemporal_predictions}}.
#'   \item \code{plots}: Named list of recorded plot objects when
#'     \code{create_plot = TRUE}. Plots can be replayed with
#'     \code{grDevices::replayPlot()}.
#' }
#'
#' @details
#'
#' Random Forests are fit using \code{\link[randomForest]{randomForest}} from
#' the \pkg{randomForest} package. The response is treated as a factor
#' (\code{0}/\code{1}) so the model runs in classification mode, which produces
#' class vote fractions used as predicted probabilities. Importance is computed
#' with \code{importance = TRUE} and \code{type = 1} (mean decrease in
#' accuracy).
#'
#' Predicted probabilities are the vote fraction for class \code{1} from
#' \code{predict(..., type = "prob")[, "1"]}. These are used for threshold
#' selection and ROC curve construction.
#'
#' Diagnostic plots include: a variable importance bar chart (mean decrease in
#' accuracy across folds), partial dependence curves for each predictor showing
#' the marginal effect of each variable while averaging over all others (with
#' rug marks for presences and pseudoabsences), and a combined ROC curve panel.
#'
#' The returned object is recognised by
#' \code{\link{generate_spatiotemporal_predictions}}, which uses the
#' \code{model_type} field to use the correct prediction and evaluation
#' logic.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_partition}},
#'   \code{\link{generate_absences}}
#'
#' Modeling: \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_hv}},
#'   \code{\link{generate_spatiotemporal_predictions}}
#'
#' External: \code{\link[randomForest]{randomForest}}
#'
#' @examples
#' \donttest{
#'   data(tmr_partition, package = "TemporalModelR")
#'
#'   data(tmr_absences,  package = "TemporalModelR")
#'
#'   build_temporal_rf(
#'     partition_result     = tmr_partition,
#'     pseudoabsence_result = tmr_absences,
#'     model_vars           = c("elevation", "forest_cover", "prseas"),
#'     rf_params            = list(ntree = 100),
#'     threshold_method     = "tss",
#'     output_dir           = tempdir(),
#'     create_plot          = FALSE,
#'     time_cols            = c("year", "season"),
#'     verbose              = FALSE
#'   )
#' }

#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom stats predict complete.cases
#' @importFrom grDevices adjustcolor dev.cur dev.off png recordPlot
#' @importFrom graphics plot barplot abline legend mtext par rug
#' @importFrom tools file_ext
#' @importFrom utils capture.output write.csv
build_temporal_rf <- function(partition_result,
                              pseudoabsence_result,
                              model_vars,
                              rf_params         = list(),
                              threshold_method  = "tss",
                              output_dir        = file.path(tempdir(), "RF_Models"),
                              create_plot       = TRUE,
                              plot_palette      = "Dark 2",
                              overwrite         = FALSE,
                              time_cols         = NULL,
                              verbose           = TRUE) {


  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop(paste0("ERROR: The 'randomForest' package is required for build_temporal_rf(). ",
                "Install it with install.packages('randomForest')."))
  }

  if (missing(partition_result)) {
    stop(paste0("ERROR: 'partition_result' is required. Provide output from ",
                "spatiotemporal_partition() or a path to an .rds file."))
  }

  if (missing(pseudoabsence_result)) {
    stop(paste0("ERROR: 'pseudoabsence_result' is required. Provide output from ",
                "generate_absences() or a path to an .rds file."))
  }

  if (missing(model_vars) || !is.character(model_vars) || length(model_vars) == 0) {
    stop("ERROR: 'model_vars' must be a non-empty character vector of predictor names.")
  }

  thr_resolved     <- .resolve_threshold_method(threshold_method)
  threshold_method <- thr_resolved$method
  manual_threshold <- thr_resolved$manual_value

  partition_result     <- .load_partition_result(partition_result, verbose = verbose)
  pseudoabsence_result <- .load_pseudoabsence_result(pseudoabsence_result, verbose = verbose)

  prep             <- .prepare_combined_modeling_data(partition_result, pseudoabsence_result,
                                                      model_vars = model_vars, time_cols = time_cols,
                                                      verbose = verbose)
  combined_df      <- prep$combined_df
  all_folds        <- prep$all_folds
  single_fold_mode <- prep$single_fold_mode

  if (verbose) message(paste0("Building Random Forest for ", length(all_folds), " fold(s)."))
  if (verbose) message(paste("Predictors:", paste(model_vars, collapse = ", ")))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  combined_file <- file.path(output_dir, "all_rf_models.rds")

  model_list      <- list()
  thresholds <- numeric(length(all_folds))
  names(thresholds) <- paste0("fold", all_folds)
  train_data_list <- list()
  importance_list <- list()
  metric_rows     <- list()

  if (file.exists(combined_file) && !overwrite) {

    if (verbose) message(paste("Loading existing RF model file:", basename(combined_file)))
    saved <- tryCatch(
      readRDS(combined_file),
      error = function(e) stop(paste0("ERROR loading saved models. Try overwrite = TRUE: ", e$message))
    )
    model_list      <- saved$models
    thresholds      <- saved$thresholds
    train_data_list <- saved$fold_training_data
    importance_list <- saved$variable_importance
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

      train_rf         <- train_data[, model_vars, drop = FALSE]
      train_rf$presence <- factor(train_data$presence, levels = c(0, 1))

      rf_args <- c(
        list(
          formula    = presence ~ .,
          data       = train_rf,
          importance = TRUE
        ),
        rf_params
      )

      model_fit <- tryCatch(
        do.call(randomForest::randomForest, rf_args),
        error = function(e) stop(paste0("ERROR: randomForest() failed for fold ", fold, ": ", e$message))
      )

      train_pred_mat  <- stats::predict(model_fit, newdata = train_rf, type = "prob")
      train_pred_prob <- train_pred_mat[, "1"]

      thr <- {
        .select_threshold_prob(
          pred_prob    = train_pred_prob,
          observed     = train_data$presence,
          method       = threshold_method,
          manual_value = manual_threshold
        )
      }
      thresholds[fold_key] <- thr

      model_list[[fold_key]]      <- model_fit
      train_data_list[[fold_key]] <- train_data

      imp_raw  <- randomForest::importance(model_fit, type = 1)
      imp_gini <- randomForest::importance(model_fit, type = 2)
      importance_list[[fold_key]] <- data.frame(
        variable         = rownames(imp_raw),
        mean_decr_acc    = imp_raw[, 1],
        mean_decr_gini   = imp_gini[, 1],
        row.names        = NULL,
        stringsAsFactors = FALSE
      )

      if (nrow(test_data) > 0) {
        test_rf <- test_data
        test_rf$presence <- factor(test_rf$presence, levels = c(0, 1))
        test_pred_mat  <- stats::predict(model_fit, newdata = test_rf, type = "prob")
        test_pred_prob <- test_pred_mat[, "1"]
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
             fold_training_data = train_data_list,
             variable_importance = importance_list,
             fold_test_metrics_raw = metric_rows),
        combined_file
      )
      if (verbose) message(paste("Saved RF models to:", basename(combined_file)))
    }, error = function(e) {
      warning(paste0("Could not save RF model file: ", e$message))
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

    imp_file <- file.path(output_dir, "VariableImportance.png")

    tryCatch({
      n_vars_imp  <- length(model_vars)
      n_folds_imp <- length(all_folds)
      n_cols_imp  <- min(n_vars_imp, 3)
      n_rows_imp  <- ceiling(n_vars_imp / n_cols_imp)
      png_w <- 280 * n_cols_imp + 60
      png_h <- (220 + n_folds_imp * 25) * n_rows_imp + 80
      grDevices::png(imp_file, width = png_w, height = png_h, res = 120)
      .plot_rf_importance(model_vars, all_folds, importance_list, fold_colors)
      grDevices::dev.off()
      if (verbose) message(paste("Saved variable importance plot:", basename(imp_file)))
    }, error = function(e) {
      warning(paste0("Could not save variable importance plot: ", e$message))
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    })

    tryCatch({
      .plot_rf_importance(model_vars, all_folds, importance_list, fold_colors)
      plot_list[["variable_importance"]] <- grDevices::recordPlot()
    }, error = function(e) {
      warning(paste0("Could not record variable importance plot: ", e$message))
    })

    for (fold in all_folds) {
      fold_key   <- paste0("fold", fold)
      model_fit  <- model_list[[fold_key]]
      train_data <- train_data_list[[fold_key]]
      thr <- thresholds[fold_key]

      plot_file <- file.path(output_dir, paste0("MarginalCurves_Fold", fold, ".png"))
      n_vars    <- length(model_vars)
      n_cols    <- min(3, n_vars)
      n_rows    <- ceiling(n_vars / n_cols)

      var_means <- colMeans(train_data[, model_vars, drop = FALSE], na.rm = TRUE)

      tryCatch({
        grDevices::png(plot_file, width = 350 * n_cols, height = 320 * n_rows, res = 120)
        .plot_rf_marginal_curves(model_fit, train_data, model_vars,
                                 var_means, fold_colors[fold_key], thr,
                                 fold, n_rows, n_cols)
        grDevices::dev.off()
        if (verbose) message(paste("Saved marginal curve plot:", basename(plot_file)))
      }, error = function(e) {
        warning(paste0("Could not save marginal curve plot for Fold ", fold, ": ", e$message))
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })

      tryCatch({
        .plot_rf_marginal_curves(model_fit, train_data, model_vars,
                                 var_means, fold_colors[fold_key], thr,
                                 fold, n_rows, n_cols)
        plot_list[[paste0("marginal_curves_fold", fold)]] <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not record marginal curve plot for Fold ", fold, ": ", e$message))
      })
    }

  }

  if (verbose) message("\nRandom Forest modeling complete.")
  if (verbose) message(paste0("Threshold: ",
                              threshold_method,
                              " | Folds: ", length(all_folds)))
  if (nrow(fold_metrics_df) > 0 && verbose) {
    message(paste("Metrics saved to:", basename(metrics_file)))
    message(paste(utils::capture.output(
      print(fold_metrics_df, row.names = FALSE)), collapse = "\n"))
  }

  if (create_plot && nrow(fold_metrics_df) > 0) {
    plot_rate_cols <- c("Sensitivity", "Specificity", "TSS", "Kappa")
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
          title       = "Random Forest Time Independent Assessment Metrics"
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
          title       = "Random Forest Time Independent Assessment Metrics"
        )
        plot_list[["fold_test_metrics"]] <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not record fold metrics plot: ", e$message))
      })
    }
  }

  result <- list(
    models              = model_list,
    thresholds          = thresholds,
    threshold_method    = threshold_method,
    model_vars          = model_vars,
    variable_importance = importance_list,
    fold_training_data  = train_data_list,
    fold_test_metrics   = fold_metrics_df,
    output_dir          = output_dir,
    model_type          = "rf",
    plots               = plot_list
  )

  class(result) <- c("TemporalRF", "list")
  invisible(result)
}
