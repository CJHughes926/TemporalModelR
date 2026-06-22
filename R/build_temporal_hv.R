#' Build Temporal Hypervolume Models Across Cross-Validation Folds
#'
#' Modeling function that constructs hypervolume models for each
#' cross-validation fold using either Gaussian kernel density estimation or
#' one-class SVM. Each hypervolume reserves one fold as testing data and uses
#' the remaining folds as training data. The returned object follows the same
#' structure as \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#' and \code{\link{build_temporal_rf}}, and is accepted directly by
#' \code{\link{generate_spatiotemporal_predictions}}.
#'
#' @usage
#' build_temporal_hv(partition_result, model_vars, method,
#'                   hypervolume_params = list(),
#'                   output_dir = file.path(tempdir(), "Hypervolume_Models"),
#'                   create_plot = TRUE, overwrite = FALSE,
#'                   plot_palette = "Dark 2", verbose = TRUE)
#'
#' @param partition_result List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to an \code{.rds}
#'   file containing that output.
#' @param model_vars Character vector. Names of predictor columns to use in
#'   hypervolume construction. All variables must be present as columns in
#'   the occurrence data produced by \code{\link{temporally_explicit_extraction}}.
#' @param method Character. Hypervolume method. One of \code{"gaussian"} for
#'   Gaussian kernel density estimation or \code{"svm"} for one-class support
#'   vector machine. Required.
#' @param hypervolume_params Named list. Additional parameters passed to
#'   \code{\link[hypervolume]{hypervolume_gaussian}} or
#'   \code{\link[hypervolume]{hypervolume_svm}}. For Gaussian models, valid
#'   keys are \code{kde.bandwidth}, \code{quantile.requested},
#'   \code{quantile.requested.type}, \code{chunk.size}, \code{verbose}, and
#'   \code{samples.per.point}. For SVM models, valid keys are \code{svm.nu},
#'   \code{svm.gamma}, \code{chunk.size}, \code{verbose}, and
#'   \code{samples.per.point}. Default is an empty list, which uses the
#'   built-in defaults.
#' @param output_dir Character. Directory to write output files including saved
#'   hypervolume objects and plots. Default is \code{file.path(tempdir(), "Hypervolume_Models")}.
#' @param create_plot Logical. If \code{TRUE}, generates pairplot visualisations
#'   of each fold's hypervolume and a combined comparison plot. Default is
#'   \code{TRUE}.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing saved
#'   hypervolume files. If \code{FALSE}, loads existing files when available.
#'   Default is \code{FALSE}.
#' @param plot_palette Character. Name of an HCL or RColorBrewer palette used
#'   to color folds in diagnostic plots. Accepts any HCL palette name (see
#'   \code{\link[grDevices]{hcl.pals}}) or, if \pkg{RColorBrewer} is installed,
#'   any Brewer palette name. Default is \code{"Dark 2"}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes per-fold hypervolume construction
#'   progress, volume summaries, and overlap statistics.
#'
#' @return A list with class \code{"TemporalHypervolume"} containing:
#' \itemize{
#'   \item \code{models}: Named list of fitted \code{Hypervolume} objects, one
#'     per fold, named \code{fold1}, \code{fold2}, etc. This naming convention
#'     matches the other \code{build_temporal_*()} functions and is required by
#'     \code{\link{generate_spatiotemporal_predictions}}.
#'   \item \code{volumes}: Named numeric vector of hypervolume sizes (units
#'     depend on the number of dimensions and bandwidth).
#'   \item \code{overlaps}: Named list of pairwise percent volume overlaps
#'     between all fold combinations.
#'   \item \code{method}: Character string recording the method used.
#'   \item \code{model_vars}: Character vector of predictor names used.
#'   \item \code{fold_training_data}: Named list of training data frames used
#'     to fit each fold model, retained for consistency with the other model
#'     types and for downstream use.
#'   \item \code{fold_test_metrics}: Data frame of E-space inclusion metrics
#'     computed on the held-out test points for each fold. Columns: \code{fold},
#'     \code{n_test}, \code{volume}, \code{tp}, \code{fn}, \code{sensitivity},
#'     \code{Sensitivity}. Printed as a summary table at the end of
#'     model building.
#'   \item \code{output_dir}: Path to the output directory.
#'   \item \code{model_type}: Character string \code{"hypervolume"}, used by
#'     \code{\link{generate_spatiotemporal_predictions}}.
#'   \item \code{plots}: List of recorded plot objects (if
#'     \code{create_plot = TRUE}). Plots can be replayed with
#'     \code{grDevices::replayPlot()}.
#' }
#'
#' @details
#' For N folds, constructs N hypervolumes where each fold reserves one group
#' of points for testing and uses the remaining N-1 groups for training.
#' Pairwise overlap statistics quantify similarity between hypervolumes across
#' folds, with low overlap indicating that the environmental space sampled
#' differs substantially across folds.
#'
#' Hypervolume objects are saved as a combined \code{.rds} file in
#' \code{output_dir}. If the file already exists and \code{overwrite = FALSE},
#' the saved file is loaded rather than re-fitting.
#'
#' The returned object is accepted directly by
#' \code{\link{generate_spatiotemporal_predictions}}, which uses the
#' \code{model_type} field to dispatch hypervolume-specific projection logic
#' via \code{\link[hypervolume]{hypervolume_project}}.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_partition}}
#'
#' Modeling: \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_rf}},
#'   \code{\link{generate_spatiotemporal_predictions}}
#'
#' External: \code{\link[hypervolume]{hypervolume_gaussian}},
#'   \code{\link[hypervolume]{hypervolume_svm}}
#'
#' @examples
#' data(tmr_partition, package = "TemporalModelR")

#' build_temporal_hv(
#'   partition_result = tmr_partition,
#'   model_vars       = c("elevation", "forest_cover", "prseas"),
#'   method           = "svm",
#'   output_dir       = tempdir(),
#'   create_plot      = FALSE,
#'   verbose          = FALSE
#' )

#' @export
#' @importFrom sf st_drop_geometry
#' @importFrom stats complete.cases var
#' @importFrom grDevices dev.cur dev.off png recordPlot
#' @importFrom graphics plot
#' @importFrom utils combn modifyList write.csv
build_temporal_hv <- function(partition_result,
                              model_vars,
                              method,
                              hypervolume_params = list(),
                              output_dir         = file.path(tempdir(), "Hypervolume_Models"),
                              create_plot        = TRUE,
                              overwrite          = FALSE,
                              plot_palette       = "Dark 2",
                              verbose            = TRUE) {

  ### Package check

  if (!requireNamespace("hypervolume", quietly = TRUE)) {
    stop(paste0(
      "ERROR: The 'hypervolume' package is required to use build_temporal_hv(). ",
      "Install it with: install.packages('hypervolume')"
    ))
  }

  if (missing(partition_result)) {
    stop(paste0("ERROR: 'partition_result' is required. Provide output from ",
                "spatiotemporal_partition() or a path to an .rds file."))
  }

  if (missing(method) || !method %in% c("gaussian", "svm")) {
    stop("ERROR: 'method' is required and must be 'gaussian' or 'svm'.")
  }

  if (missing(model_vars) || !is.character(model_vars) || length(model_vars) == 0) {
    stop("ERROR: 'model_vars' must be a non-empty character vector of predictor names.")
  }

  partition_result <- .load_partition_result(partition_result)
  occurrence_points <- partition_result$points_sf

  if (is.null(occurrence_points) || nrow(occurrence_points) == 0) {
    stop("ERROR: 'partition_result$points_sf' is empty.")
  }
  if (!"fold" %in% names(occurrence_points)) {
    stop(paste0(
      "ERROR: Missing 'fold' column in 'partition_result$points_sf'. Available columns: ",
      paste(names(occurrence_points)[names(occurrence_points) != "geometry"], collapse = ", ")
    ))
  }

  occ_df   <- sf::st_drop_geometry(occurrence_points)
  all_folds <- sort(unique(occ_df$fold[!is.na(occ_df$fold)]))

  if (length(all_folds) == 0) {
    stop("ERROR: All fold values are NA. Re-run spatiotemporal_partition().")
  }

  single_fold_mode <- length(all_folds) == 1
  if (single_fold_mode) {
    if (verbose) message("Single-fold mode: all points used for both training and testing.")
  } else if (length(all_folds) < 2) {
    warning(paste0("Only ", length(all_folds), " unique fold detected."))
  }

  if (verbose) message(paste("Loaded", nrow(occ_df), "points across", length(all_folds), "folds."))

  available_vars <- names(occ_df)[names(occ_df) != "geometry"]
  missing_vars   <- model_vars[!model_vars %in% available_vars]
  if (length(missing_vars) > 0) {
    stop(paste0(
      "ERROR: The following 'model_vars' are not present in the occurrence data: ",
      paste(missing_vars, collapse = ", "),
      ". Available: ", paste(available_vars, collapse = ", ")
    ))
  }

  if (length(hypervolume_params) > 0) {
    valid_gaussian <- c("kde.bandwidth", "quantile.requested", "quantile.requested.type",
                        "chunk.size", "verbose", "samples.per.point")
    valid_svm      <- c("svm.nu", "svm.gamma", "chunk.size", "verbose", "samples.per.point")
    valid_params   <- if (method == "gaussian") valid_gaussian else valid_svm
    invalid        <- setdiff(names(hypervolume_params), valid_params)
    if (length(invalid) > 0) {
      warning(paste0(
        "Unrecognized hypervolume_params for method '", method, "': ",
        paste(invalid, collapse = ", "),
        ". Valid keys: ", paste(valid_params, collapse = ", ")
      ))
    }
    if (method == "gaussian" && "quantile.requested" %in% names(hypervolume_params)) {
      qr <- hypervolume_params$quantile.requested
      if (!is.numeric(qr) || qr <= 0 || qr > 1) {
        warning(paste0("quantile.requested should be between 0 and 1, got: ", qr))
      }
    }
    if (method == "svm" && "svm.nu" %in% names(hypervolume_params)) {
      nu <- hypervolume_params$svm.nu
      if (!is.numeric(nu) || nu <= 0 || nu > 1) {
        warning(paste0("svm.nu should be between 0 and 1, got: ", nu))
      }
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  n_original    <- nrow(occ_df)
  na_counts     <- vapply(model_vars, function(v) sum(is.na(occ_df[[v]])), integer(1))
  if (any(na_counts > 0)) {
    warning(paste0(
      "NA values in: ",
      paste(paste0(names(na_counts[na_counts > 0]), " (", na_counts[na_counts > 0], ")"), collapse = ", ")
    ))
  }
  occ_df    <- occ_df[stats::complete.cases(occ_df[, model_vars, drop = FALSE]), ]
  n_removed <- n_original - nrow(occ_df)
  if (n_removed > 0) {
    warning(paste0("Removed ", n_removed, " rows with NAs. Remaining: ", nrow(occ_df)))
  }
  if (nrow(occ_df) == 0) {
    stop("ERROR: No complete cases remain after removing NA values in model predictors.")
  }

  fold_ids <- sort(unique(occ_df$fold))
  if (verbose) message(paste0(
    "Building hypervolume models (", method, ") for ",
    length(fold_ids), " folds: ", paste(fold_ids, collapse = ", ")
  ))

  default_gaussian <- list(
    quantile.requested      = 0.95,
    quantile.requested.type = "probability",
    chunk.size              = 1000,
    verbose                 = FALSE
  )
  default_svm <- list(
    svm.nu     = 0.01,
    svm.gamma  = 0.5,
    chunk.size = 1000,
    verbose    = FALSE
  )

  combined_file <- file.path(output_dir, paste0("all_hypervolumes_", method, ".rds"))

  model_list      <- list()
  train_data_list <- list()

  if (file.exists(combined_file) && !overwrite) {

    if (verbose) message(paste("Loading existing hypervolume file:", basename(combined_file)))
    saved <- tryCatch(
      readRDS(combined_file),
      error = function(e) stop(paste0("ERROR loading saved hypervolumes. Try overwrite = TRUE: ", e$message))
    )

    if (is.list(saved) && all(c("models", "fold_test_metrics", "model_type") %in% names(saved))) {
      if (verbose) message(paste("Loaded complete result with", length(saved$models), "hypervolumes."))
      class(saved) <- c("TemporalHypervolume", "list")
      return(invisible(saved))
    }

    if (is.list(saved) && "models" %in% names(saved)) {
      model_list      <- saved$models
      train_data_list <- saved$fold_training_data
    } else {
      model_list <- saved
      for (fold in fold_ids) {
        fold_key <- paste0("fold", fold)
        train_data_list[[fold_key]] <- occ_df[occ_df$fold != fold, model_vars, drop = FALSE]
      }
    }

    for (fold in fold_ids) {
      fold_key <- paste0("fold", fold)
      if (!is.null(model_list[[fold_key]]) &&
          (is.null(model_list[[fold_key]]@Name) || model_list[[fold_key]]@Name == "untitled")) {
        model_list[[fold_key]]@Name <- paste("Fold", fold)
      }
    }
    if (verbose) message(paste("Loaded", length(model_list), "hypervolumes from file."))

  } else {

    hv_fn <- if (method == "gaussian") hypervolume::hypervolume_gaussian else hypervolume::hypervolume_svm
    defaults <- if (method == "gaussian") default_gaussian else default_svm

    for (fold in fold_ids) {

      fold_key   <- paste0("fold", fold)
      train_data <- if (single_fold_mode) occ_df[, model_vars, drop = FALSE] else occ_df[occ_df$fold != fold, model_vars, drop = FALSE]
      n_train    <- nrow(train_data)

      if (verbose) message(paste0("Fold ", fold, ": training on ", n_train, " points."))

      if (n_train < 5) {
        stop(paste0("ERROR: Fold ", fold, " has only ", n_train, " training points. Minimum 5 required."))
      }
      if (n_train < 10) {
        warning(paste0("Fold ", fold, " has only ", n_train, " training points. Results may be unreliable."))
      }

      var_check <- vapply(train_data[, model_vars, drop = FALSE], function(x) stats::var(x, na.rm = TRUE), numeric(1))
      zero_var  <- names(var_check[!is.na(var_check) & var_check == 0])
      if (length(zero_var) > 0) {
        warning(paste0("Zero variance in Fold ", fold, " for: ", paste(zero_var, collapse = ", "),
                       ". These predictors will have no effect in this fold."))
      }

      hv_params <- modifyList(defaults, hypervolume_params)

      hv <- tryCatch({
        h <- suppressMessages(do.call(hv_fn, c(list(data = train_data), hv_params)))
        h@Name <- paste("Fold", fold)
        h
      }, error = function(e) {
        stop(paste0("ERROR: Hypervolume construction failed for Fold ", fold, ": ", e$message))
      })

      vol <- tryCatch(
        hypervolume::get_volume(hv),
        error = function(e) NA_real_
      )
      if (verbose) message(paste0("  Volume: ", round(vol, 4)))

      model_list[[fold_key]]      <- hv
      train_data_list[[fold_key]] <- train_data
    }

    tryCatch({
      saveRDS(
        list(models = model_list, fold_training_data = train_data_list),
        combined_file
      )
      if (verbose) message(paste("Saved hypervolumes to:", basename(combined_file)))
    }, error = function(e) {
      warning(paste0("Could not save hypervolume file: ", e$message))
    })
  }

  volumes <- vapply(model_list, function(h) {
    tryCatch(hypervolume::get_volume(h), error = function(e) NA_real_)
  }, numeric(1))

  if (verbose) message("--- Fold Volumes ---")
  for (fold in fold_ids) {
    fold_key <- paste0("fold", fold)
    if (verbose) message(paste0("  Fold ", fold, ": ", round(volumes[fold_key], 4)))
  }

  overlap_stats <- list()

  if (length(model_list) > 1) {
    if (verbose) message("Calculating pairwise hypervolume overlaps...")
    fold_pairs <- utils::combn(names(model_list), 2, simplify = FALSE)

    for (pair in fold_pairs) {
      pair_label <- paste(
        "Fold", gsub("fold", "", pair[1]),
        "vs Fold", gsub("fold", "", pair[2])
      )
      pct <- tryCatch({
        hv_set  <- suppressMessages(hypervolume::hypervolume_set(
          model_list[[pair[1]]], model_list[[pair[2]]],
          check.memory = FALSE, verbose = FALSE
        ))
        stats_out <- suppressMessages(hypervolume::hypervolume_overlap_statistics(hv_set))
        if (is.list(stats_out) && "percent_volume_overlap" %in% names(stats_out)) {
          stats_out$percent_volume_overlap
        } else if (is.numeric(stats_out)) {
          stats_out[1]
        } else NA_real_
      }, error = function(e) {
        warning(paste0("Could not calculate overlap for ", pair_label, ": ", e$message))
        NA_real_
      })
      overlap_stats[[pair_label]] <- pct
      if (verbose) message(paste0("  ", pair_label, ": ", round(pct * 100, 2), "%"))
    }
  }

  plot_list   <- list()
  fold_colors <- .resolve_palette(length(fold_ids), plot_palette)
  names(fold_colors) <- paste0("fold", fold_ids)

  if (create_plot && length(model_vars) >= 2) {
    if (verbose) message("Generating hypervolume plots...")

    for (fold in fold_ids) {
      fold_key  <- paste0("fold", fold)
      plot_file <- file.path(output_dir, paste0("Hypervolume_Fold", fold, ".png"))

      tryCatch({
        grDevices::png(plot_file, width = 1200, height = 1000, res = 150)
        suppressMessages(graphics::plot(
          model_list[[fold_key]],
          pairplot = TRUE, show.3d = FALSE,
          main     = paste("Fold", fold, "Hypervolume"),
          colors   = fold_colors[fold_key]
        ))
        grDevices::dev.off()
        if (verbose) message(paste("Saved Fold", fold, "plot:", basename(plot_file)))
      }, error = function(e) {
        warning(paste0("Could not save plot for Fold ", fold, ": ", e$message))
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })

      tryCatch({
        suppressMessages(graphics::plot(
          model_list[[fold_key]],
          pairplot = TRUE, show.3d = FALSE,
          main     = paste("Fold", fold, "Hypervolume"),
          colors   = fold_colors[fold_key]
        ))
        plot_list[[paste0("fold", fold)]] <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not record plot for Fold ", fold, ": ", e$message))
      })
    }

    if (length(model_list) > 1) {
      comp_file <- file.path(output_dir, "Hypervolume_Comparison.png")
      hv_joined <- tryCatch(
        suppressMessages(do.call(hypervolume::hypervolume_join, model_list)),
        error = function(e) {
          warning(paste0("Could not join hypervolumes for comparison plot: ", e$message))
          NULL
        }
      )
      if (!is.null(hv_joined)) {
        tryCatch({
          grDevices::png(comp_file, width = 1400, height = 1200, res = 150)
          suppressMessages(graphics::plot(
            hv_joined,
            pairplot = TRUE, show.3d = FALSE,
            main     = paste("Hypervolume Comparison:", length(model_list), "Folds"),
            colors   = fold_colors
          ))
          grDevices::dev.off()
          if (verbose) message(paste("Saved comparison plot:", basename(comp_file)))
        }, error = function(e) {
          warning(paste0("Could not save comparison plot: ", e$message))
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
        })
        tryCatch({
          suppressMessages(graphics::plot(
            hv_joined,
            pairplot = TRUE, show.3d = FALSE,
            main     = paste("Hypervolume Comparison:", length(model_list), "Folds"),
            colors   = fold_colors
          ))
          plot_list[["comparison"]] <- grDevices::recordPlot()
        }, error = function(e) {
          warning(paste0("Could not record comparison plot: ", e$message))
        })
      }
    }
  }

  if (verbose) message("Computing E-space inclusion metrics on held-out test folds...")

  metric_rows <- list()

  for (fold in fold_ids) {
    fold_key  <- paste0("fold", fold)
    hv_model  <- model_list[[fold_key]]
    test_data <- occ_df[occ_df$fold == fold, model_vars, drop = FALSE]
    test_data <- test_data[stats::complete.cases(test_data), , drop = FALSE]

    tp <- NA_integer_; fn <- NA_integer_
    sensitivity <- NA_real_

    if (!is.null(hv_model) && nrow(test_data) > 0) {
      inclusion <- tryCatch(
        suppressMessages(suppressWarnings(
          hypervolume::hypervolume_inclusion_test(hv_model, test_data)
        )),
        error = function(e) {
          warning(paste0("Inclusion test failed for Fold ", fold, ": ", e$message))
          NULL
        }
      )
      if (!is.null(inclusion)) {
        tp  <- sum(inclusion, na.rm = TRUE)
        fn  <- sum(!inclusion, na.rm = TRUE)
        tot <- tp + fn
        sensitivity <- if (tot > 0) tp / tot else NA_real_
      }
    }

    vol <- tryCatch(hypervolume::get_volume(hv_model), error = function(e) NA_real_)

    metric_rows[[fold_key]] <- data.frame(
      Fold          = fold,
      N_Test        = nrow(test_data),
      E_Volume      = round(vol, 4),
      Testing_TP    = tp,
      Testing_FN    = fn,
      Sensitivity   = round(sensitivity, 4),
      stringsAsFactors = FALSE
    )
  }

  fold_metrics_df           <- do.call(rbind, metric_rows)
  rownames(fold_metrics_df) <- NULL

  metrics_file <- file.path(output_dir, "Fold_Test_Metrics.csv")
  tryCatch(
    utils::write.csv(fold_metrics_df, metrics_file, row.names = FALSE),
    error = function(e) warning(paste0("Could not save fold test metrics: ", e$message))
  )

  if (verbose) message("--- Hypervolume Modeling Complete ---")
  if (verbose) message(paste0("  Method: ", method))
  if (nrow(fold_metrics_df) > 0) {
    if (verbose) message("  Held-out E-space inclusion metrics:")
    if (verbose) message(paste0("  ", paste(sprintf("%-10s",
                                                    c("Fold", "N_Test", "E_Volume", "Testing_TP", "Testing_FN", "Sensitivity")),
                                            collapse = " ")))
    for (i in seq_len(nrow(fold_metrics_df))) {
      r <- fold_metrics_df[i, ]
      if (verbose) message(paste0("  ", paste(sprintf("%-10s", c(
        r$Fold,
        r$N_Test,
        round(r$E_Volume,    3),
        r$Testing_TP,
        r$Testing_FN,
        round(r$Sensitivity, 3)
      )), collapse = " ")))
    }
    if (verbose) message(paste0("  Means      ",
                                paste(sprintf("%-10s", c(
                                  "",
                                  round(mean(fold_metrics_df$N_Test,      na.rm = TRUE), 0),
                                  round(mean(fold_metrics_df$E_Volume,    na.rm = TRUE), 3),
                                  round(mean(fold_metrics_df$Testing_TP,  na.rm = TRUE), 1),
                                  round(mean(fold_metrics_df$Testing_FN,  na.rm = TRUE), 1),
                                  round(mean(fold_metrics_df$Sensitivity, na.rm = TRUE), 3)
                                )), collapse = " ")
    ))
    if (verbose) message(paste("  Saved to:", basename(metrics_file)))
  }

  if (create_plot && nrow(fold_metrics_df) > 0) {
    plot_rate_cols <- c("Sensitivity")
    bar_file <- file.path(output_dir, "Fold_Test_Metrics.png")
    tryCatch({
      grDevices::png(bar_file, width = 300 * length(plot_rate_cols) + 200,
                     height = 700, res = 150)
      .plot_fold_metrics_bars(
        metrics_df  = fold_metrics_df,
        rate_cols   = plot_rate_cols,
        fold_colors = fold_colors,
        title       = "Hypervolume Time Independent Assessment Metrics"
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
        title       = "Hypervolume Time Independent Assessment Metrics"
      )
      plot_list[["fold_test_metrics"]] <- grDevices::recordPlot()
    }, error = function(e) {
      warning(paste0("Could not record fold metrics plot: ", e$message))
    })
  }

  result <- list(
    models             = model_list,
    volumes            = volumes,
    overlaps           = overlap_stats,
    method             = method,
    model_vars         = model_vars,
    fold_training_data = train_data_list,
    fold_test_metrics  = fold_metrics_df,
    output_dir         = output_dir,
    model_type         = "hypervolume",
    plots              = plot_list
  )

  class(result) <- c("TemporalHypervolume", "list")

  tryCatch({
    saveRDS(result, combined_file)
    if (verbose) message(paste("Updated saved object with metrics and plots:", basename(combined_file)))
  }, error = function(e) {
    warning(paste0("Could not update saved hypervolume file: ", e$message))
  })

  invisible(result)
}
