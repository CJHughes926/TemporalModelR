#' Generate Spatiotemporal Predictions from Temporal Models
#'
#' Generates temporally explicit habitat suitability predictions by projecting
#' fitted models onto environmental raster stacks matching each time period.
#' Accepts output from any of the four TemporalModelR modeling functions.
#' Produces per-time-step assessment metrics showing how predictions and
#' test point coverage vary over time.
#'
#'
#' @usage
#' generate_spatiotemporal_predictions(partition_result, model_result,
#'                                     pseudoabsence_result = NULL, raster_dir,
#'                                     variable_patterns, time_cols, time_steps,
#'                                     output_dir = file.path(tempdir(), "Predictions"),
#'                                     overwrite = FALSE,
#'                                     verbose = TRUE)
#'
#' @param partition_result List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to an \code{.rds}
#'   file containing that output.
#' @param model_result List or character. Output from any of
#'   \code{\link{build_temporal_hv}}, \code{\link{build_temporal_glm}},
#'   \code{\link{build_temporal_gam}}, or \code{\link{build_temporal_rf}}, or
#'   a path to an \code{.rds} file. Model type is detected automatically from
#'   the \code{model_type} field.
#' @param pseudoabsence_result List, character, or \code{NULL}. Optional.
#'   Output from \code{\link{generate_absences}} or path to an
#'   \code{.rds} file. When supplied for presence/absence models (GLM, GAM, RF),
#'   the held-out pseudoabsence test points for each fold are filtered to the
#'   current time step and used alongside presence test points to compute
#'   per-timestep TN, FP, Specificity, and TSS. These columns are added to the
#'   timestep metrics table when pseudoabsences are available. Ignored for
#'   hypervolume models. Default is \code{NULL}.
#' @param raster_dir Character. Directory containing environmental raster
#'   files (\code{.tif}), typically the output of
#'   \code{\link{raster_align}} or \code{\link{scale_rasters}}. File names
#'   must follow the patterns supplied in \code{variable_patterns}, with any
#'   time placeholder substituted for the corresponding value from
#'   \code{time_cols}.
#' @param variable_patterns Named character vector mapping clean variable names
#'   to raster filename patterns. For time-varying variables include the time
#'   placeholder in the pattern (e.g. \code{"forest_cover" = "forest_cover_YEAR"});
#'   for static variables omit it (e.g. \code{"elevation" = "elevation"}). Time
#'   placeholders must match entries in \code{time_cols}.
#' @param time_cols Character. Name of the column(s) containing year or time
#'   step values in the occurrence data. Must match \code{time_cols} used in
#'   \code{\link{spatiotemporal_partition}} and the time placeholders used in
#'   \code{variable_patterns}.
#' @param time_steps Vector, data frame, or matrix of time periods for which
#'   to generate predictions.
#' @param output_dir Character. Directory to write prediction rasters and the
#'   assessment metrics CSV. Default is \code{file.path(tempdir(), "Predictions")}.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing output files.
#'   If \code{FALSE} (default), existing files are skipped.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes per-time-step prediction progress.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{timestep_metrics}: data frame of per-fold per-time-step
#'     metrics. Columns always present: \code{Fold}, the time column(s),
#'     \code{Pct_Suitable}, \code{N_Pres}, \code{TP}, \code{FN},
#'     \code{Sensitivity}, \code{CBP}. When \code{pseudoabsence_result} is
#'     supplied for parametric models, additionally: \code{N_Abs}, \code{TN},
#'     \code{FP}, \code{Specificity}, \code{TSS}. Saved as
#'     \code{Timestep_Assessment_Metrics.csv}.
#'   \item \code{overall_summary}: data frame of pooled metrics per fold.
#'     Columns: \code{Fold}, \code{N_Timesteps}, \code{Mean_Pct_Suitable},
#'     \code{Total_TP}, \code{Total_FN}, \code{Overall_Sensitivity},
#'     \code{Overall_CBP}. When pseudoabsences available: additionally
#'     \code{Total_TN}, \code{Total_FP}, \code{Overall_Specificity},
#'     \code{Overall_TSS}.
#'   \item \code{prediction_files}: character vector of paths to saved
#'     prediction rasters.
#'   \item \code{model_type}: character string recording the model type used.
#' }
#'
#' @details
#' G-space predictions are produced as rasters for each time-step and fold of an
#' input model.
#'
#' Per-timestep metrics are computed per fold per time step by extracting
#' raster predictions at the held-out presence points that fall within that
#' time step. \code{Pct_Suitable} records the proportion of the study area
#' predicted suitable. \code{Sensitivity} shows how consistently test points are
#' captured. \code{CBP} (cumulative binomial probability) tests whether the
#' observed number of correctly predicted test points is better than expected
#' by random placement: under the null, each point has \code{Pct_Suitable}
#' probability of falling in suitable area, so
#' \code{dbinom(TP, N_Test_Pts, Pct_Suitable)} gives the probability of
#' observing exactly \code{TP} correct by chance. Small values indicate
#' predictions are better than random.
#'
#' @seealso
#' Preprocessing: \code{\link{scale_rasters}},
#'   \code{\link{spatiotemporal_partition}}
#'
#' Modeling: \code{\link{build_temporal_hv}},
#'   \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_rf}}
#'
#' Post-processing: \code{\link{summarize_raster_outputs}},
#'   \code{\link{plot_model_assessment}}
#'
#' @examples
#' data(tmr_partition, package = "TemporalModelR")
#'
#' data(tmr_glm,       package = "TemporalModelR")
#'
#' data(tmr_absences,  package = "TemporalModelR")
#'
#' scl_dir    <- system.file("extdata/rasters_scaled",
#'                           package = "TemporalModelR")
#'
#' time_steps <- expand.grid(
#'   year             = 1:15,
#'   season           = "Spring",
#'   stringsAsFactors = FALSE
#' )
#'
#' generate_spatiotemporal_predictions(
#'   partition_result     = tmr_partition,
#'   model_result         = tmr_glm,
#'   pseudoabsence_result = tmr_absences,
#'   raster_dir           = scl_dir,
#'   variable_patterns    = c(
#'     "elevation"    = "elevation",
#'     "forest_cover" = "forest_cover_YEAR",
#'     "prseas"       = "prseas_YEAR_SEASON"
#'   ),
#'   time_cols            = c("year", "season"),
#'   time_steps           = time_steps,
#'   output_dir           = tempdir(),
#'   overwrite            = TRUE,
#'   verbose              = FALSE
#' )

#' @export
#' @importFrom terra rast extract values vect ncell nlyr writeRaster
#' @importFrom sf st_drop_geometry
#' @importFrom stats complete.cases predict dbinom setNames
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils read.csv write.csv
generate_spatiotemporal_predictions <- function(partition_result,
                                                model_result,
                                                pseudoabsence_result = NULL,
                                                raster_dir,
                                                variable_patterns,
                                                time_cols,
                                                time_steps,
                                                output_dir = file.path(tempdir(), "Predictions"),
                                                overwrite  = FALSE,
                                                verbose    = TRUE) {

  if (missing(model_result)) {
    stop("ERROR: 'model_result' is required.")
  }

  partition_result  <- .load_partition_result(partition_result)
  occurrence_points <- partition_result$points_sf

  if (is.null(occurrence_points) || nrow(occurrence_points) == 0) {
    stop("ERROR: partition_result$points_sf is empty.")
  }
  if (!"fold" %in% names(occurrence_points)) {
    stop("ERROR: partition_result$points_sf is missing 'fold' column.")
  }

  occ_df    <- sf::st_drop_geometry(occurrence_points)
  all_folds <- sort(unique(occ_df$fold[!is.na(occ_df$fold)]))
  points_sf <- occurrence_points
  if (verbose) message(paste("Loaded", nrow(occ_df), "points across", length(all_folds), "folds."))

  pseudoabs_sf <- NULL
  if (!is.null(pseudoabsence_result)) {
    pseudoabsence_result <- .load_pseudoabsence_result(pseudoabsence_result)
    pseudoabs_sf <- pseudoabsence_result$pseudoabsences
    if (is.null(pseudoabs_sf) || nrow(pseudoabs_sf) == 0) {
      warning("pseudoabsence_result$pseudoabsences is empty  -  per-timestep TN/FP/Specificity/TSS will not be computed.")
      pseudoabs_sf <- NULL
    } else {
      if (verbose) message(paste("Loaded", nrow(pseudoabs_sf), "pseudoabsence points for per-timestep evaluation."))
    }
  }

  if (is.character(model_result)) {
    if (!file.exists(model_result)) stop(paste0("ERROR: 'model_result' file not found: ", model_result))
    if (tolower(tools::file_ext(model_result)) != "rds") stop("ERROR: 'model_result' must be .rds format.")
    if (verbose) message(paste("Reading model results from:", basename(model_result)))
    model_result <- tryCatch(readRDS(model_result),
                             error = function(e) stop(paste0("ERROR reading model_result: ", e$message)))
  }

  if (!is.list(model_result) || is.null(model_result$model_type)) {
    stop("ERROR: 'model_result' must be output from a build_temporal_*() function.")
  }

  model_type <- model_result$model_type
  if (!model_type %in% c("hypervolume", "glm", "gam", "rf")) {
    stop(paste0("ERROR: Unrecognized model_type '", model_type, "'."))
  }
  if (verbose) message(paste("Model type:", model_type))

  if (model_type == "hypervolume" && !requireNamespace("hypervolume", quietly = TRUE)) {
    stop("ERROR: The 'hypervolume' package is required. Install with: install.packages('hypervolume')")
  }

  model_list    <- model_result$models
  model_vars    <- model_result$model_vars
  is_parametric <- model_type != "hypervolume"
  thresholds    <- NULL

  if (is.null(model_list) || length(model_list) == 0) stop("ERROR: model_result$models is empty.")
  if (is.null(model_vars) || length(model_vars) == 0) stop("ERROR: model_result$model_vars is missing.")

  if (is_parametric) {
    thresholds <- model_result$thresholds
    if (is.null(thresholds) || all(is.na(thresholds))) {
      stop("ERROR: model_result$thresholds is missing or all NA. Re-run the model with a valid threshold_method.")
    }
  }

  if (verbose) message(paste("Loaded", length(model_list), "model(s)."))

  if (missing(variable_patterns)) {
    stop("ERROR: 'variable_patterns' is required.")
  }
  .validate_variable_patterns(variable_patterns)
  if (missing(time_cols) || is.null(time_cols)) stop("ERROR: 'time_cols' is required.")
  time_cols <- as.character(time_cols)
  if (any(!time_cols %in% names(occ_df))) {
    stop(paste0("ERROR: 'time_cols' '", paste(setdiff(time_cols, names(occ_df)), collapse = ", "),
                "' not found in occurrence data."))
  }
  if (missing(time_steps)) stop("ERROR: 'time_steps' is required.")
  if (is.data.frame(time_steps) || is.matrix(time_steps)) {
    time_steps_df <- as.data.frame(time_steps, stringsAsFactors = FALSE)
    missing_cols  <- setdiff(time_cols, names(time_steps_df))
    if (length(missing_cols) > 0) stop(paste0("ERROR: 'time_steps' missing columns: ", paste(missing_cols, collapse = ", ")))
  } else if (is.vector(time_steps) && !is.list(time_steps)) {
    if (length(time_cols) == 1) {
      time_steps_df <- stats::setNames(data.frame(time_steps, stringsAsFactors = FALSE), time_cols)
    } else {
      other_cols   <- time_cols[-1]
      other_combos <- unique(occ_df[, other_cols, drop = FALSE])
      cross <- expand.grid(
        c(list(time_steps), as.list(other_combos)),
        stringsAsFactors = FALSE
      )
      names(cross)  <- time_cols
      time_steps_df <- cross
    }
  } else {
    stop("ERROR: 'time_steps' must be a plain vector of time step values, or a data frame with one column per time_cols.")
  }
  if (missing(raster_dir) || !dir.exists(raster_dir)) {
    stop(paste0("ERROR: 'raster_dir' does not exist: ", raster_dir))
  }

  if (verbose) message(paste("Processing", nrow(time_steps_df), "time steps."))

  var_classes  <- .classify_variable_patterns(variable_patterns, time_cols)
  dynamic_vars <- var_classes$dynamic
  static_vars  <- var_classes$static

  if (verbose) message(paste("Dynamic variables:", if (length(dynamic_vars) > 0) paste(dynamic_vars, collapse = ", ") else "none"))
  if (verbose) message(paste("Static variables:",  if (length(static_vars)  > 0) paste(static_vars,  collapse = ", ") else "none"))

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fold_keys <- paste0("fold", all_folds)

  timestep_metrics_file <- file.path(output_dir, "Timestep_Assessment_Metrics.csv")
  timestep_metrics <- if (!overwrite && file.exists(timestep_metrics_file)) {
    tryCatch(utils::read.csv(timestep_metrics_file, stringsAsFactors = FALSE),
             error = function(e) data.frame())
  } else data.frame()

  prediction_files <- character(0)

  for (i in seq_len(nrow(time_steps_df))) {

    time_values <- time_steps_df[i, , drop = FALSE]
    time_label  <- paste(vapply(time_cols, function(tc) as.character(time_values[[tc]]),
                                character(1)), collapse = "_")
    output_path <- file.path(output_dir, paste0("Prediction_", time_label, ".tif"))

    if (file.exists(output_path) && !overwrite) {
      if (verbose) message(paste("Skipping", time_label, " - output already exists."))
      prediction_files <- c(prediction_files, output_path)
      next
    }

    if (verbose) message(paste0("--- Time step: ", time_label, " ---"))

    raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                          time_cols, time_values, raster_dir)
    if (is.null(raster_paths)) {
      warning(paste("Missing rasters for", time_label, " - skipping."))
      next
    }

    raster_stack <- tryCatch(terra::rast(raster_paths),
                             error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL })
    if (is.null(raster_stack)) next

    expected_names <- vapply(raster_paths, function(p) {
      tools::file_path_sans_ext(basename(p))
    }, character(1))
    if (terra::nlyr(raster_stack) == length(expected_names)) {
      names(raster_stack) <- expected_names
    }

    time_filter <- rep(TRUE, nrow(points_sf))
    for (tc in time_cols) {
      time_filter <- time_filter & !is.na(points_sf[[tc]]) & (points_sf[[tc]] == time_values[[tc]])
    }

    fold_projections <- vector("list", length(all_folds))
    names(fold_projections) <- fold_keys

    for (fi in seq_along(all_folds)) {

      fold     <- all_folds[fi]
      fold_key <- fold_keys[fi]

      if (verbose) message(paste("  Fold", fold, "..."))

      test_pts_time <- points_sf[!is.na(points_sf$fold) & points_sf$fold == fold & time_filter, ]

      pseudo_pts_time <- if (!is.null(pseudoabs_sf) && is_parametric) {
        pa_fold_mask <- !is.na(pseudoabs_sf$fold) & pseudoabs_sf$fold == fold
        pa_time_mask <- rep(TRUE, nrow(pseudoabs_sf))
        for (tc in time_cols) {
          if (tc %in% names(pseudoabs_sf)) {
            pa_time_mask <- pa_time_mask &
              !is.na(pseudoabs_sf[[tc]]) & (pseudoabs_sf[[tc]] == time_values[[tc]])
          }
        }
        pseudoabs_sf[pa_fold_mask & pa_time_mask, ]
      } else NULL

      binary_raster <- if (model_type == "hypervolume") {
        .predict_hypervolume_raster(model_list[[fold_key]], raster_stack, fold, time_label)
      } else {
        thr <- thresholds[fold_key]
        if (is.null(thr) || is.na(thr)) { warning(paste("No threshold for", fold_key)); NULL }
        else .predict_parametric_raster(model_list[[fold_key]], raster_stack, model_vars,
                                        thr, model_type, variable_patterns, time_cols, time_values)
      }

      fold_projections[[fold_key]] <- binary_raster
      if (is.null(binary_raster)) next

      timestep_row <- .compute_timestep_metrics(
        binary_raster   = binary_raster,
        test_pts_time   = test_pts_time,
        pseudo_pts_time = pseudo_pts_time,
        model           = if (is_parametric) model_list[[fold_key]] else NULL,
        model_vars      = model_vars,
        threshold       = if (is_parametric) thresholds[fold_key] else NULL,
        model_type      = model_type,
        time_values     = time_values,
        time_cols       = time_cols,
        fold            = fold
      )
      timestep_metrics <- rbind(timestep_metrics, timestep_row)
    }

    valid_proj <- fold_projections[!vapply(fold_projections, is.null, logical(1))]
    if (length(valid_proj) == 0) { warning(paste("No valid projections for", time_label)); next }

    combined_raster <- Reduce("+", valid_proj)

    tryCatch({
      terra::writeRaster(combined_raster, output_path, overwrite = TRUE)
      if (verbose) message(paste("  Saved:", basename(output_path)))
      prediction_files <- c(prediction_files, output_path)
    }, error = function(e) warning(paste("Could not save raster for", time_label, ":", e$message)))

    tryCatch(utils::write.csv(timestep_metrics, timestep_metrics_file, row.names = FALSE),
             error = function(e) warning(paste("Could not save timestep metrics:", e$message)))
  }

  if (verbose) message("--- Predictions complete ---")
  if (verbose) message(paste("Timestep metrics:", basename(timestep_metrics_file)))

  overall_summary <- data.frame()

  if (nrow(timestep_metrics) > 0) {
    has_pa_metrics <- all(c("TN", "FP") %in% names(timestep_metrics))

    for (fold in all_folds) {
      fm <- timestep_metrics[!is.na(timestep_metrics$Fold) & timestep_metrics$Fold == fold, ]
      if (nrow(fm) == 0) next

      sum_tp   <- sum(fm$TP,           na.rm = TRUE)
      sum_fn   <- sum(fm$FN,           na.rm = TRUE)
      tot_p    <- sum_tp + sum_fn
      mean_pct <- mean(fm$Pct_Suitable, na.rm = TRUE)
      sens     <- if (tot_p > 0) sum_tp / tot_p else NA_real_
      cbp_ov   <- if (tot_p > 0 && !is.na(mean_pct)) {
        stats::dbinom(sum_tp, size = tot_p, prob = mean_pct)
      } else NA_real_

      row <- data.frame(
        Fold               = fold,
        N_Timesteps        = nrow(fm),
        Mean_Pct_Suitable  = round(mean_pct, 4),
        Total_TP           = sum_tp,
        Total_FN           = sum_fn,
        Overall_Sensitivity = round(sens,    4),
        Overall_CBP        = cbp_ov,
        stringsAsFactors   = FALSE
      )

      if (has_pa_metrics) {
        sum_tn  <- sum(fm$TN, na.rm = TRUE)
        sum_fp  <- sum(fm$FP, na.rm = TRUE)
        tot_n   <- sum_tn + sum_fp
        spec    <- if (tot_n > 0) sum_tn / tot_n else NA_real_
        tss     <- if (!is.na(sens) && !is.na(spec)) sens + spec - 1 else NA_real_
        row$Total_TN            <- sum_tn
        row$Total_FP            <- sum_fp
        row$Overall_Specificity <- round(spec, 4)
        row$Overall_TSS         <- round(tss,  4)
      }

      overall_summary <- rbind(overall_summary, row)
    }

    if (nrow(overall_summary) > 0) {
      rate_cols <- c("Overall_Sensitivity", if (has_pa_metrics) c("Overall_Specificity", "Overall_TSS") else NULL)
      hdr_cols  <- c("Fold", "N_Steps", "MeanPct", "Tot_TP", "Tot_FN",
                     if (has_pa_metrics) c("Tot_TN", "Tot_FP") else NULL,
                     "Sens", if (has_pa_metrics) c("Spec", "TSS") else NULL, "CBP")
      if (verbose) message("  Overall pooled metrics per fold:")
      if (verbose) message(paste0("  ", paste(sprintf("%-8s", hdr_cols), collapse = " ")))

      for (i in seq_len(nrow(overall_summary))) {
        r <- overall_summary[i, ]
        vals <- c(r$Fold, r$N_Timesteps,
                  round(r$Mean_Pct_Suitable, 3),
                  r$Total_TP, r$Total_FN,
                  if (has_pa_metrics) c(r$Total_TN, r$Total_FP) else NULL,
                  round(r$Overall_Sensitivity, 3),
                  if (has_pa_metrics) c(round(r$Overall_Specificity, 3), round(r$Overall_TSS, 3)) else NULL,
                  formatC(r$Overall_CBP, format = "e", digits = 2))
        if (verbose) message(paste0("  ", paste(sprintf("%-8s", vals), collapse = " ")))
      }

      means <- c("",
                 round(mean(overall_summary$N_Timesteps, na.rm = TRUE), 0),
                 round(mean(overall_summary$Mean_Pct_Suitable, na.rm = TRUE), 3),
                 round(mean(overall_summary$Total_TP, na.rm = TRUE), 1),
                 round(mean(overall_summary$Total_FN, na.rm = TRUE), 1),
                 if (has_pa_metrics) c(
                   round(mean(overall_summary$Total_TN, na.rm = TRUE), 1),
                   round(mean(overall_summary$Total_FP, na.rm = TRUE), 1)
                 ) else NULL,
                 round(mean(overall_summary$Overall_Sensitivity, na.rm = TRUE), 3),
                 if (has_pa_metrics) c(
                   round(mean(overall_summary$Overall_Specificity, na.rm = TRUE), 3),
                   round(mean(overall_summary$Overall_TSS,         na.rm = TRUE), 3)
                 ) else NULL,
                 formatC(mean(overall_summary$Overall_CBP, na.rm = TRUE), format = "e", digits = 2))
      if (verbose) message(paste0("  Means    ", paste(sprintf("%-8s", means), collapse = " ")))
    }
  }

  invisible(list(
    timestep_metrics = timestep_metrics,
    overall_summary  = overall_summary,
    prediction_files = prediction_files,
    model_type       = model_type
  ))
}
