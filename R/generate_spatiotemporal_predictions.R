#' Generate Spatiotemporal Predictions from Hypervolume Models
#'
#' Modeling function that generates temporally explicit habitat suitability
#' predictions by projecting hypervolume models onto environmental raster stacks
#' for each time period. Produces binary prediction rasters and comprehensive
#' evaluation metrics.
#'
#' @param partition_results List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to .rds file.
#' @param hypervolume_results List or character. Output from
#'   \code{\link{build_hypervolume_models}} or path to .rds file.
#' @param raster_dir Character. Directory containing scaled environmental
#'   rasters. Should contain output from \code{\link{scale_rasters}}.
#' @param variable_patterns Named character vector mapping variable names to
#'   raster filename patterns. Must match patterns used in
#'   \code{\link{temporally_explicit_extraction}}.
#' @param time_cols Character vector. Time column names for temporal matching.
#' @param output_dir Character. Directory to write prediction rasters and
#'   evaluation tables.
#' @param output_prefix Character. Prefix for output filenames. Default is
#'   "predictions".
#' @param threshold Numeric. Threshold for binary classification. Default is
#'   0.5.
#' @param overwrite Logical. If TRUE, overwrites existing prediction files. If
#'   FALSE, skips files that already exist. Default is FALSE.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item assessment_table: data frame of evaluation metrics for all
#'     fold-time period combinations
#'   \item prediction_files: character vector of paths to prediction rasters
#' }
#'
#' @details
#' For each cross-validation fold and time period, projects the hypervolume
#' model onto scaled environmental rasters using
#' \code{\link[hypervolume]{hypervolume_project}}. Generates binary predictions
#' and calculates assessment metrics via
#' \code{\link{model_assessment_metrics}}.
#'
#' Output prediction rasters serve as input for
#' \code{\link{summarize_raster_outputs}} and subsequent temporal pattern
#' analyses.
#'
#' @seealso
#' Preprocessing: \code{\link{scale_rasters}},
#' \code{\link{spatiotemporal_partition}}
#'
#' Modeling: \code{\link{build_hypervolume_models}},
#' \code{\link{model_assessment_metrics}}
#'
#' Postprocessing: \code{\link{summarize_raster_outputs}},
#' \code{\link{plot_model_assessment}}
#'
#' External: \code{\link[hypervolume]{hypervolume_project}}
#'
#' @examples
#' \dontrun{
#' predictions <- generate_spatiotemporal_predictions(
#'   partition_results = "partition_results.rds",
#'   hypervolume_results = "hypervolume_models.rds",
#'   raster_dir = "scaled_rasters/",
#'   variable_patterns = c("bio1" = "bio1_YEAR", "bio12" = "bio12_YEAR"),
#'   time_cols = "YEAR",
#'   output_dir = "predictions/"
#' )
#' }
#'
#' @export
#' @importFrom hypervolume hypervolume_project
#' @importFrom terra rast nlyr writeRaster mask
#' @importFrom raster stack
#' @importFrom sp SpatialPointsDataFrame CRS spTransform
#' @importFrom sf st_coordinates st_drop_geometry st_as_sf
#' @importFrom utils write.csv
#' @importFrom tools file_ext
generate_spatiotemporal_predictions <- function(partition_results,
                                                hypervolume_results,
                                                time_cols,
                                                time_steps,
                                                variable_patterns,
                                                raster_dir,
                                                output_dir = "./Predictions",
                                                overwrite = FALSE) {

  require(hypervolume)
  require(terra)
  require(dplyr)
  require(sf)

  print("Generating spatiotemporal predictions")

  ### VALIDATE AND LOAD partition_results INPUT

  if (missing(partition_results)) {
    stop("partition_results is required. Provide output from spatiotemporal_partition() as list or .rds file path")
  }

  if (is.character(partition_results)) {
    if (!file.exists(partition_results)) {
      stop(paste("partition_results file not found:", partition_results, "Working directory:", getwd()))
    }
    file_ext <- tolower(tools::file_ext(partition_results))
    if (file_ext != "rds") {
      stop(paste("partition_results must be .rds format, got:", file_ext))
    }
    print(paste("Reading partition results from:", basename(partition_results)))
    partition_data <- readRDS(partition_results)
  } else if (is.list(partition_results)) {
    print("Using provided partition results list...")
    partition_data <- partition_results
  } else {
    stop(paste("partition_results must be file path or list, got:", class(partition_results)[1]))
  }

  if (!"points_sf" %in% names(partition_data)) {
    stop(paste("partition_results missing 'points_sf' element. Available:", paste(names(partition_data), collapse = ", ")))
  }

  occurrence_points <- partition_data$points_sf

  if (is.null(occurrence_points) || nrow(occurrence_points) == 0) {
    stop("partition_results$points_sf is empty")
  }

  if (!"fold" %in% names(occurrence_points)) {
    stop(paste("Missing 'fold' column in partition_results$points_sf. Available columns:", paste(names(occurrence_points)[names(occurrence_points) != "geometry"], collapse = ", ")))
  }

  if (inherits(occurrence_points, "sf")) {
    folds <- st_drop_geometry(occurrence_points)$fold
  } else {
    folds <- occurrence_points$fold
  }

  if (all(is.na(folds))) {
    stop("All fold values are NA in partition_results")
  }

  unique_folds <- sort(unique(folds[!is.na(folds)]))
  print(paste("Loaded", nrow(occurrence_points), "points with", length(unique_folds), "folds"))

  ### VALIDATE AND LOAD hypervolume_results INPUT

  if (missing(hypervolume_results)) {
    stop("hypervolume_results is required. Provide output from build_hypervolume_models() as list or .rds file path")
  }

  if (is.character(hypervolume_results)) {
    if (!file.exists(hypervolume_results)) {
      stop(paste("hypervolume_results file not found:", hypervolume_results, "Working directory:", getwd()))
    }
    file_ext <- tolower(tools::file_ext(hypervolume_results))
    if (file_ext != "rds") {
      stop(paste("hypervolume_results must be .rds format, got:", file_ext))
    }
    print(paste("Reading hypervolume results from:", basename(hypervolume_results)))
    hv_data <- readRDS(hypervolume_results)
  } else if (is.list(hypervolume_results)) {
    print("Using provided hypervolume results list...")
    hv_data <- hypervolume_results
  } else {
    stop(paste("hypervolume_results must be file path or list, got:", class(hypervolume_results)[1]))
  }

  if ("hypervolumes" %in% names(hv_data)) {
    hypervolume_models <- hv_data$hypervolumes
    print(paste("Extracted hypervolumes from hypervolume results"))
  } else if (all(sapply(hv_data, function(x) inherits(x, "Hypervolume")))) {
    hypervolume_models <- hv_data
    print("Using hypervolume list directly")
  } else {
    stop("hypervolume_results must be output from build_hypervolume_models() or a list of Hypervolume objects")
  }

  if (!is.list(hypervolume_models) || !all(sapply(hypervolume_models, function(x) inherits(x, "Hypervolume")))) {
    stop("All elements in hypervolume_models must be Hypervolume objects")
  }

  n_folds <- length(hypervolume_models)
  if (n_folds != length(unique_folds)) {
    stop(paste("Number of hypervolume models (", n_folds, ") does not match number of unique folds (", length(unique_folds), ")", sep = ""))
  }
  print(paste("Loaded", n_folds, "hypervolume models"))

  ### INHERIT model_vars FROM variable_patterns
  if (missing(variable_patterns)) {
    stop("variable_patterns is required. Provide a named vector mapping model variables to raster file patterns.")
  }
  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop("variable_patterns must be a named vector.")
  }

  model_vars <- names(variable_patterns)
  clean_model_vars <- model_vars
  print(paste("Model variables (inherited from variable_patterns):", paste(model_vars, collapse = ", ")))

  ### VALIDATE time_cols
  if (missing(time_cols)) {
    stop("time_cols is required. Provide a character vector of time column names.")
  }
  if (!is.character(time_cols) || length(time_cols) == 0) {
    stop("time_cols must be a character vector with at least one time column name")
  }
  missing_time_cols <- setdiff(time_cols, names(occurrence_points))
  if (length(missing_time_cols) > 0) {
    stop(paste("time_cols not found in occurrence_points:", paste(missing_time_cols, collapse = ", ")))
  }
  print(paste("Time columns:", paste(time_cols, collapse = ", ")))

  ### VALIDATE time_steps
  if (missing(time_steps)) stop("time_steps is required.")
  if (is.vector(time_steps) && length(time_cols) == 1) {
    time_steps_df <- data.frame(x = time_steps)
    names(time_steps_df) <- time_cols
  } else if (is.data.frame(time_steps) || is.matrix(time_steps)) {
    time_steps_df <- as.data.frame(time_steps)
    missing_cols <- setdiff(time_cols, names(time_steps_df))
    if (length(missing_cols) > 0) stop(paste("time_steps missing columns:", paste(missing_cols, collapse = ", ")))
  } else stop("time_steps must be a vector (single column), data frame, or matrix")
  print(paste("Processing", nrow(time_steps_df), "time steps"))

  ### VALIDATE raster_dir
  if (missing(raster_dir)) stop("raster_dir is required. Provide directory containing scaled raster files")
  if (!dir.exists(raster_dir)) stop(paste("raster_dir does not exist:", raster_dir))
  print(paste("Reading rasters from:", raster_dir))

  ### CREATE OUTPUT DIRECTORY
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  }

  ### CONVERT OCCURRENCE POINTS TO SF-compatible
  if (inherits(occurrence_points, "sf")) {
    occ_df <- st_drop_geometry(occurrence_points)
    coords <- st_coordinates(occurrence_points)
    points_sf <- st_as_sf(data.frame(occ_df, coords), coords = c("X","Y"), crs = 4326)
  } else if (inherits(occurrence_points, "SpatialPointsDataFrame")) {
    points_sf <- st_as_sf(occurrence_points)
  } else {
    stop("occurrence_points must be an sf object or SpatialPointsDataFrame")
  }

  points_sf$fold <- folds

  ### IDENTIFY VARIABLE TYPES
  dynamic_vars <- c()
  static_vars <- c()
  for (var_name in names(variable_patterns)) {
    pattern <- variable_patterns[var_name]
    has_time_component <- any(sapply(time_cols, function(tc) grepl(tc, pattern, ignore.case = TRUE)))
    if (has_time_component) dynamic_vars <- c(dynamic_vars, var_name) else static_vars <- c(static_vars, var_name)
  }
  print(paste("Dynamic variables:", ifelse(length(dynamic_vars) > 0, paste(dynamic_vars, collapse = ", "), "none")))
  print(paste("Static variables:", ifelse(length(static_vars) > 0, paste(static_vars, collapse = ", "), "none")))

  metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")  # define first
  if (overwrite || !file.exists(metrics_file)) {
    all_model_data <- data.frame()
  } else {
    # Load existing metrics so we can append new results
    all_model_data <- tryCatch({
      read.csv(metrics_file)
    }, error = function(e) {
      data.frame()
    })
  }

  ### PROCESS EACH TIME STEP
  for (i in 1:nrow(time_steps_df)) {
    time_values <- time_steps_df[i, , drop = FALSE]
    time_label <- paste(sapply(time_cols, function(tc) as.character(time_values[[tc]])), collapse = "_")
    output_path <- file.path(output_dir, paste0("Prediction_", time_label, ".tif"))
    if (file.exists(output_path) && !overwrite) {
      print(paste("Skipping time step", time_label, "- output already exists"))
      next
    }
    print(paste("========================================"))
    print(paste("Processing time step:", time_label))
    print(paste("========================================"))

    ### LOAD AND STACK RASTERS
    raster_files <- c()
    for (var in dynamic_vars) {
      pattern <- variable_patterns[var]
      file_name <- pattern

      # Only substitute time columns that are actually present in this pattern
      for (tc in time_cols) {
        if (grepl(tc, file_name, ignore.case = TRUE)) {
          file_name <- gsub(tc, as.character(time_values[[tc]]), file_name, ignore.case = TRUE)
        }
      }

      raster_files <- c(raster_files, paste0(file_name, "_Scaled.tif"))
    }
    for (var in static_vars) raster_files <- c(raster_files, paste0(variable_patterns[var], "_Scaled.tif"))

    raster_paths <- file.path(raster_dir, raster_files)
    missing_files <- raster_paths[!file.exists(raster_paths)]
    if (length(missing_files) > 0) {
      warning(paste("Missing raster files for time step", time_label, ":", paste(basename(missing_files), collapse = ", ")))
      next
    }

    s2 <- rast(raster_paths)

    ### PROCESS EACH FOLD
    fold_projections <- list()
    for (j in 1:n_folds) {
      current_fold <- unique_folds[j]
      print(paste("  Processing Fold", current_fold, "..."))
      hypervolume_model <- hypervolume_models[[j]]

      dims <- hypervolume_model@Dimensionality
      n_rasters <- nlyr(s2)
      if (dims != n_rasters) {
        stop(paste0("Dimensionality mismatch: hypervolume models have ",
                    paste(dims, collapse = ", "),
                    " dimensions, but raster stack has ", n_rasters, " layers."))
      }



      test_points_all <- points_sf[points_sf$fold == current_fold, ]
      time_filter <- rep(TRUE, nrow(points_sf))
      for (tc in time_cols) time_filter <- time_filter & (points_sf[[tc]] == time_values[[tc]])
      test_points_year <- points_sf[points_sf$fold == current_fold & time_filter, ]
      hv_map <- tryCatch({
        hypervolume_project(hypervolume_model, rasters = s2, type = 'inclusion', fast.or.accurate = 'fast', verbose = FALSE)
      }, error = function(e) {
        warning(paste("Error projecting hypervolume for Fold", current_fold, "at time step", time_label, ":", e$message))
        NULL
      })
      fold_projections[[j]] <- hv_map
      if (!is.null(hv_map)) {
        metrics <- model_assessment_metrics(hypervolume_model = hypervolume_model,
                                            projected_raster = hv_map,
                                            test_points_current_year = test_points_year,
                                            test_points_all_years = test_points_all,
                                            variable_patterns = variable_patterns)
        time_data <- lapply(time_cols, function(tc) time_values[[tc]])
        names(time_data) <- time_cols
        model_data <- data.frame(time_data, Fold = current_fold,
                                 G_volume = round(metrics$G_volume, 3),
                                 E_volume = round(metrics$E_volume, 3),
                                 CBP_test_G = metrics$CBP_test_G,
                                 TP_test_G = metrics$TP_test_G,
                                 FN_test_G = metrics$FN_test_G,
                                 sensitivity_test_G = round(metrics$sensitivity_test_G, 3),
                                 omission_test_G = round(metrics$omission_test_G, 3),
                                 CBP_test_E = metrics$CBP_test_E,
                                 TP_test_E = metrics$TP_test_E,
                                 FN_test_E = metrics$FN_test_E,
                                 sensitivity_test_E = round(metrics$sensitivity_test_E, 3),
                                 omission_test_E = round(metrics$omission_test_E, 3),
                                 stringsAsFactors = FALSE)
        all_model_data <- rbind(all_model_data, model_data)
      }
    }

    ### COMBINE FOLDS AND SAVE
    valid_projections <- fold_projections[!sapply(fold_projections, is.null)]
    if (length(valid_projections) == 0) {
      warning(paste("No valid projections for time step", time_label, "- skipping"))
      next
    }
    result_raster <- valid_projections[[1]]
    if (length(valid_projections) > 1) for (k in 2:length(valid_projections)) result_raster <- result_raster + valid_projections[[k]]
    tryCatch({
      writeRaster(result_raster, output_path, overwrite = TRUE)
      print(paste("  Saved combined prediction to:", basename(output_path)))
    }, error = function(e) warning(paste("Error saving combined prediction for time step", time_label, ":", e$message)))

    # Save metrics after each time step
    metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")
    tryCatch({
      write.csv(all_model_data, metrics_file, row.names = FALSE)
      print(paste("  Metrics saved after time step", time_label))
    }, error = function(e) warning(paste("Error saving metrics after time step", time_label, ":", e$message)))

    print(paste("  Completed time step:", time_label))
  }

  print("========================================")
  print("All predictions complete!")
  print(paste("Metrics saved to:", file.path(output_dir, "Model_Assessment_Metrics.csv")))
  print("========================================")

  return(all_model_data)
}
