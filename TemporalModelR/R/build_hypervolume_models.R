build_hypervolume_models <- function(occurrence_points,
                                     folds,
                                     model_vars,
                                     method = "gaussian",
                                     output_dir = "./Hypervolume_Models",
                                     hypervolume_params = list(),
                                     create_plot = TRUE,
                                     overwrite = FALSE) {

  # --- Load required libraries ---
  require(hypervolume)
  require(dplyr)
  require(sf)

  print(paste("Building hypervolume models using method:", method))

  # --- Validate inputs ---
  if (!inherits(occurrence_points, "sf")) stop("occurrence_points must be an sf object")
  if (length(folds) != nrow(occurrence_points)) stop("Length of folds must match number of occurrence points")
  if (!all(model_vars %in% colnames(occurrence_points))) {
    missing_vars <- model_vars[!model_vars %in% colnames(occurrence_points)]
    stop(paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  }
  if (!method %in% c("gaussian", "svm")) stop("method must be either 'gaussian' or 'svm'")

  # --- Create output directory ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    print(paste("Created output directory:", output_dir))
  }

  # --- Prepare data ---
  occ_df <- st_drop_geometry(occurrence_points)
  occ_df$fold <- folds
  occ_df <- occ_df[complete.cases(occ_df[, model_vars]), ]

  fold_ids <- sort(unique(occ_df$fold))
  print(paste("Detected", length(fold_ids), "folds:", paste(fold_ids, collapse = ", ")))

  # --- Default hypervolume parameters ---
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

  # --- Build hypervolumes for each fold ---
  hv_list <- list()
  for (fold in fold_ids) {
    hv_file <- file.path(output_dir, paste0("hypervolume_model_fold", fold, "_", method, ".rds"))

    if (file.exists(hv_file) && !overwrite) {
      print(paste("Loading existing hypervolume for Fold", fold, "..."))
      hv <- readRDS(hv_file)
    } else {
      print("\n========================================")
      print(paste("Building hypervolume for Fold", fold))
      print("========================================")

      train_data <- occ_df[occ_df$fold != fold, model_vars]
      print(paste("Training with", nrow(train_data), "points"))

      hv_params <- modifyList(
        if (method == "gaussian") default_params_gaussian else default_params_svm,
        hypervolume_params
      )

      fn <- if (method == "gaussian") hypervolume_gaussian else hypervolume_svm
      hv <- do.call(fn, c(list(data = train_data), hv_params))

      saveRDS(hv, hv_file)
      print(paste("Saved Fold", fold, "hypervolume to:", hv_file))
    }

    hv_list[[paste0("fold", fold)]] <- hv
    print(paste("Fold", fold, "volume:", round(get_volume(hv), 4)))
  }

  # --- Create plots ---
  if (create_plot && length(model_vars) >= 2) {
    print("\n========================================")
    print("Creating hypervolume plots")
    print("========================================")

    for (fold in fold_ids) {
      plot_file <- file.path(output_dir, paste0("hypervolume_plot_fold", fold, "_", method, ".png"))
      png(plot_file, width = 10, height = 10, units = "in", res = 300)
      plot(hv_list[[paste0("fold", fold)]], pairplot = TRUE, show.3d = FALSE,
           main = paste("Fold", fold, "Hypervolume"))
      dev.off()
      print(paste("Saved:", plot_file))
    }

    # Combined comparison plot (if >1 fold)
    if (length(hv_list) > 1) {
      hv_joined <- do.call(hypervolume_join, hv_list)
      comparison_file <- file.path(output_dir, paste0("hypervolume_comparison_", method, ".png"))
      png(comparison_file, width = 10, height = 10, units = "in", res = 300)
      plot(hv_joined, pairplot = TRUE, show.3d = FALSE,
           main = paste("Hypervolume Comparison across", length(hv_list), "Folds"))
      dev.off()
      print(paste("Comparison plot saved to:", comparison_file))
    }
  } else {
    print("Skipping plots (need ≥2 variables)")
  }

  # --- Summary Statistics ---
  print("\n========================================")
  print("SUMMARY STATISTICS")
  print("========================================")

  volumes <- sapply(hv_list, get_volume)
  print(volumes)

  # --- Pairwise Overlaps ---
  overlap_stats <- list()
  if (length(hv_list) > 1) {
    fold_pairs <- combn(names(hv_list), 2, simplify = FALSE)
    for (pair in fold_pairs) {
      hv_set <- hypervolume_set(hv_list[[pair[1]]], hv_list[[pair[2]]], check.memory = FALSE, verbose = FALSE)
      stats <- hypervolume_overlap_statistics(hv_set)
      overlap_stats[[paste(pair, collapse = "_vs_")]] <- stats
    }
  }

  # --- Return ---
  return(list(
    hypervolumes = hv_list,
    volumes = volumes,
    overlaps = overlap_stats,
    method = method,
    model_vars = model_vars,
    output_dir = output_dir
  ))
}
