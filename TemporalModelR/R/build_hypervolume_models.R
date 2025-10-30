build_hypervolume_models <- function(partition_results,
                                     model_vars,
                                     method = "gaussian",
                                     output_dir = "./Hypervolume_Models",
                                     hypervolume_params = list(),
                                     create_plot = TRUE,
                                     overwrite = FALSE) {

  require(hypervolume)
  require(dplyr)
  require(sf)

  print(paste("Building hypervolume models using method:", method))

  if (missing(partition_results)) {
    stop(paste(
      "Error: partition_results is required.",
      "\nThis should be the output from spatiotemporal_partition().",
      "\nProvide either:",
      "\n  - The list object directly from spatiotemporal_partition()",
      "\n  - A file path to the saved .rds file"
    ))
  }

  # Handle partition_results input
  if (is.character(partition_results)) {
    if (!file.exists(partition_results)) {
      stop(paste("Error: File does not exist:", partition_results))
    }

    file_ext <- tolower(tools::file_ext(partition_results))
    if (file_ext == "rds") {
      print(paste("Reading partition results from:", basename(partition_results)))
      partition_data <- readRDS(partition_results)
    } else {
      stop(paste("Error: partition_results file must be .rds format. Got:", file_ext))
    }

  } else if (is.list(partition_results)) {
    print("Using provided partition results list...")
    partition_data <- partition_results
  } else {
    stop("Error: partition_results must be either a file path to an .rds file or a list object from spatiotemporal_partition()")
  }

  # Validate partition_data structure
  if (!"points_sf" %in% names(partition_data)) {
    stop(paste(
      "Error: partition_results must contain 'points_sf' element.",
      "\nThis should be the output from spatiotemporal_partition().",
      "\nAvailable elements:", paste(names(partition_data), collapse = ", ")
    ))
  }

  # Extract occurrence points
  occurrence_points <- partition_data$points_sf

  # Check if fold column exists in points_sf
  if (!"fold" %in% names(occurrence_points)) {
    stop("Error: partition_results$points_sf must contain a 'fold' column")
  }

  # Extract fold column from points_sf
  if (inherits(occurrence_points, "sf")) {
    folds <- st_drop_geometry(occurrence_points)$fold
  } else {
    folds <- occurrence_points$fold
  }

  print(paste("Loaded", nrow(occurrence_points), "points from partition results"))

  ### VALIDATE REMAINING INPUTS

  if (length(folds) != nrow(occurrence_points)) {
    stop(paste(
      "Error: Length of folds (", length(folds),
      ") must match number of occurrence points (", nrow(occurrence_points), ")",
      sep = ""
    ))
  }

  if (!all(model_vars %in% colnames(occurrence_points))) {
    missing_vars <- model_vars[!model_vars %in% colnames(occurrence_points)]
    available_vars <- colnames(occurrence_points)
    available_vars <- available_vars[available_vars != "geometry"]
    stop(paste(
      "Error: Missing variables:", paste(missing_vars, collapse = ", "),
      "\nAvailable variables:", paste(available_vars, collapse = ", ")
    ))
  }

  if (!method %in% c("gaussian", "svm")) {
    stop("Error: method must be either 'gaussian' or 'svm'")
  }

  ### CREATE OUTPUT DIRECTORY

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  }

  ### PREPARE DATA

  occ_df <- st_drop_geometry(occurrence_points)
  occ_df$fold <- folds

  # Check for missing values in model variables
  n_original <- nrow(occ_df)
  occ_df <- occ_df[complete.cases(occ_df[, model_vars]), ]
  n_removed <- n_original - nrow(occ_df)

  if (n_removed > 0) {
    pct_removed <- round(n_removed / n_original * 100, 2)
    warning(paste0("Removed ", n_removed, " rows (", pct_removed,
                   "%) with missing values in model variables"))
  }

  if (nrow(occ_df) == 0) {
    stop("Error: No complete cases remaining after removing missing values")
  }

  fold_ids <- sort(unique(occ_df$fold))
  print(paste("Detected", length(fold_ids), "folds:", paste(fold_ids, collapse = ", ")))

  ### DEFAULT HYPERVOLUME PARAMETERS

  default_params_gaussian <- list(
    quantile.requested = 0.95,
    quantile.requested.type = "probability",
    chunk.size = 1000,
    verbose = TRUE
  )
  default_params_svm <- list(
    svm.nu = 0.01,
    svm.gamma = 0.5,
    chunk.size = 1000,
    verbose = TRUE
  )

  ### BUILD HYPERVOLUMES FOR EACH FOLD

  hv_list <- list()
  for (fold in fold_ids) {
    hv_file <- file.path(output_dir, paste0("hypervolume_model_fold", fold, "_", method, ".rds"))

    if (file.exists(hv_file) && !overwrite) {
      print(paste("Loading existing hypervolume for Fold", fold, "..."))
      hv <- readRDS(hv_file)
    } else {
      print(paste("\nBuilding hypervolume for Fold", fold))

      train_data <- occ_df[occ_df$fold != fold, model_vars]
      print(paste("Training with", nrow(train_data), "points"))

      if (nrow(train_data) < 10) {
        warning(paste("Fold", fold, "has only", nrow(train_data),
                      "training points. Results may be unreliable."))
      }

      hv_params <- modifyList(
        if (method == "gaussian") default_params_gaussian else default_params_svm,
        hypervolume_params
      )

      fn <- if (method == "gaussian") hypervolume_gaussian else hypervolume_svm
      hv <- do.call(fn, c(list(data = train_data), hv_params))

      saveRDS(hv, hv_file)
      print(paste("Saved Fold", fold, "hypervolume to:", basename(hv_file)))
    }

    hv_list[[paste0("fold", fold)]] <- hv
    print(paste("Fold", fold, "volume:", round(get_volume(hv), 4)))
  }

  ### CREATE PLOTS

  if (create_plot && length(model_vars) >= 2) {
    print("\nCreating hypervolume plots")

    for (fold in fold_ids) {
      plot_file <- file.path(output_dir, paste0("hypervolume_plot_fold", fold, "_", method, ".png"))
      png(plot_file, width = 10, height = 10, units = "in", res = 300)
      plot(hv_list[[paste0("fold", fold)]], pairplot = TRUE, show.3d = FALSE,
           main = paste("Fold", fold, "Hypervolume"))
      dev.off()
      print(paste("Saved:", basename(plot_file)))
    }

    # Combined comparison plot (if >1 fold)
    if (length(hv_list) > 1) {
      hv_joined <- do.call(hypervolume_join, hv_list)
      comparison_file <- file.path(output_dir, paste0("hypervolume_comparison_", method, ".png"))
      png(comparison_file, width = 10, height = 10, units = "in", res = 300)
      plot(hv_joined, pairplot = TRUE, show.3d = FALSE,
           main = paste("Hypervolume Comparison across", length(hv_list), "Folds"))
      dev.off()
      print(paste("Comparison plot saved to:", basename(comparison_file)))
    }
  } else if (create_plot) {
    print("Skipping plots (need ≥2 variables)")
  }

  ### SUMMARY STATISTICS

  print("\nSUMMARY STATISTICS")

  volumes <- sapply(hv_list, get_volume)
  print(volumes)

  ### PAIRWISE OVERLAPS

  overlap_stats <- list()
  if (length(hv_list) > 1) {
    print("\nCalculating pairwise overlaps...")
    fold_pairs <- combn(names(hv_list), 2, simplify = FALSE)
    for (pair in fold_pairs) {
      hv_set <- hypervolume_set(hv_list[[pair[1]]], hv_list[[pair[2]]],
                                check.memory = FALSE, verbose = FALSE)
      stats <- hypervolume_overlap_statistics(hv_set)
      overlap_stats[[paste(pair, collapse = "_vs_")]] <- stats
    }
  }

  print("\nProcessing complete!")

  ### RETURN

  return(list(
    hypervolumes = hv_list,
    volumes = volumes,
    overlaps = overlap_stats,
    method = method,
    model_vars = model_vars,
    output_dir = output_dir
  ))
}
