generate_spatiotemporal_predictions <- function(hypervolume_models,
                                                occurrence_points,
                                                folds,
                                                time_col = "Year",
                                                time_steps,
                                                variable_patterns,
                                                raster_dir,
                                                mask = NULL,
                                                output_dir = "./Predictions",
                                                model_vars,
                                                combine_folds = TRUE,
                                                overwrite = FALSE) {

  require(hypervolume)
  require(raster)
  require(terra)
  require(sp)
  require(dplyr)
  require(sf)

  # Validate inputs
  if (!is.list(hypervolume_models)) {
    stop("hypervolume_models must be a list of hypervolume objects")
  }

  n_folds <- length(hypervolume_models)
  unique_folds <- sort(unique(folds))

  if (n_folds != length(unique_folds)) {
    stop(paste("Number of hypervolume models (", n_folds,
               ") does not match number of unique folds (",
               length(unique_folds), ")", sep = ""))
  }

  print(paste("Generating spatiotemporal predictions with", n_folds, "folds"))
  print(paste("Processing", length(time_steps), "time steps"))

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  }

  # Convert occurrence points to SpatialPointsDataFrame if needed
  if (inherits(occurrence_points, "sf")) {
    coords <- st_coordinates(occurrence_points)
    occ_df <- st_drop_geometry(occurrence_points)
    points_sp <- SpatialPointsDataFrame(
      coords, occ_df,
      proj4string = CRS("+proj=longlat +datum=WGS84")
    )
  } else if (inherits(occurrence_points, "SpatialPointsDataFrame")) {
    points_sp <- occurrence_points
  } else {
    stop("occurrence_points must be an sf object or SpatialPointsDataFrame")
  }

  # Add fold assignments
  points_sp$fold <- folds

  # Load mask if provided
  if (!is.null(mask)) {
    if (is.character(mask)) {
      mask_raster <- raster(mask)
    } else {
      mask_raster <- mask
    }
    # Transform points to mask CRS
    points_sp <- spTransform(points_sp, crs(mask_raster))
  } else {
    mask_raster <- NULL
  }

  # Identify static vs dynamic variables
  dynamic_vars <- names(variable_patterns)[grepl("YEAR", variable_patterns)]
  static_vars <- names(variable_patterns)[!grepl("YEAR", variable_patterns)]

  print(paste("Dynamic variables:", paste(dynamic_vars, collapse = ", ")))
  print(paste("Static variables:", paste(static_vars, collapse = ", ")))

  # Clean model variable names
  clean_model_vars <- gsub("^X", "", model_vars)

  # Initialize results
  all_model_data <- data.frame()

  # Process each time step
  for (time_step in time_steps) {

    output_path <- file.path(output_dir, paste0("Prediction_", time_step, ".tif"))

    if (file.exists(output_path) && !overwrite) {
      print(paste("Skipping time step", time_step, "- output already exists"))
      next
    }

    print(paste("\n========================================"))
    print(paste("Processing time step:", time_step))
    print(paste("========================================"))

    # =============================
    # LOAD AND STACK RASTERS
    # =============================
    raster_files <- c()

    # Dynamic variables - substitute YEAR with time_step
    for (var in dynamic_vars) {
      pattern <- variable_patterns[var]
      file_name <- gsub("YEAR", time_step, pattern)
      file_name <- paste0(file_name, "_Scaled.tif")
      raster_files <- c(raster_files, file_name)
    }

    # Static variables - no year substitution
    for (var in static_vars) {
      pattern <- variable_patterns[var]
      file_name <- paste0(pattern, "_Scaled.tif")
      raster_files <- c(raster_files, file_name)
    }

    print(paste("  Loading", length(raster_files), "raster files..."))

    # Stack rasters
    raster_paths <- file.path(raster_dir, raster_files)

    # Check if all files exist
    missing_files <- raster_paths[!file.exists(raster_paths)]
    if (length(missing_files) > 0) {
      warning(paste("Missing raster files for time step", time_step, ":",
                    paste(basename(missing_files), collapse = ", ")))
      next
    }

    s <- stack(raster_paths)

    # Apply mask if provided
    if (!is.null(mask_raster)) {
      s_masked <- terra::mask(s, mask = sum(s))
      s2 <- stack(rast(s_masked))
    } else {
      s2 <- s
    }

    # Sort layers to match model variable order
    s2 <- s2[[sort(names(s2))]]

    print(paste("  Raster stack created with", nlayers(s2), "layers"))

    # =============================
    # PROCESS EACH FOLD
    # =============================
    fold_projections <- list()

    for (i in 1:n_folds) {
      current_fold <- unique_folds[i]
      print(paste("  Processing Fold", current_fold, "..."))

      # Get hypervolume model for this fold
      hypervolume_model <- hypervolume_models[[i]]

      # Subset points
      # Testing: all points in current fold
      test_points_all <- subset(points_sp, points_sp$fold == current_fold)
      test_points_year <- subset(points_sp,
                                 points_sp$fold == current_fold &
                                   points_sp[[time_col]] == time_step)

      # Project hypervolume to raster
      hv_map <- hypervolume_project(
        hypervolume_model,
        rasters = s2,
        type = 'inclusion',
        fast.or.accurate = 'fast',
        verbose = FALSE
      )

      # Store projection for later combination
      fold_projections[[i]] <- hv_map

      # Calculate metrics
      metrics <- model_assessment_metrics(
        hypervolume_model = hypervolume_model,
        projected_raster = hv_map,
        test_points_current_year = test_points_year,
        test_points_all_years = test_points_all,
        model_vars = model_vars
      )

      # Store results
      model_data <- data.frame(
        Time_Step = time_step,
        Fold = current_fold,
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
        stringsAsFactors = FALSE
      )

      all_model_data <- rbind(all_model_data, model_data)
    }

    # =============================
    # COMBINE FOLDS AND SAVE
    # =============================
    if (combine_folds) {
      print(paste("  Combining all", n_folds, "folds..."))

      # Sum all fold projections
      result_raster <- fold_projections[[1]]
      if (n_folds > 1) {
        for (i in 2:n_folds) {
          result_raster <- result_raster + fold_projections[[i]]
        }
      }

      # Save combined projection
      writeRaster(result_raster, output_path, overwrite = TRUE)
      print(paste("  Saved combined prediction to:", output_path))

    } else {
      # Save separately for each fold
      for (i in 1:n_folds) {
        current_fold <- unique_folds[i]
        output_path_fold <- file.path(output_dir,
                                      paste0("Prediction_Fold", current_fold,
                                             "_", time_step, ".tif"))
        writeRaster(fold_projections[[i]], output_path_fold, overwrite = TRUE)
        print(paste("  Saved Fold", current_fold, "prediction to:", output_path_fold))
      }
    }

    # Save metrics after each time step
    metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")
    write.csv(all_model_data, metrics_file, row.names = FALSE)

    print(paste("  Completed time step:", time_step))
  }

  # Final save of metrics
  metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")
  write.csv(all_model_data, metrics_file, row.names = FALSE)

  print("\n========================================")
  print("All predictions complete!")
  print(paste("Metrics saved to:", metrics_file))
  print("========================================")

  return(all_model_data)
}
