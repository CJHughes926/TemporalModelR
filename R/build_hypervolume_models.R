#' Build hypervolume models across folds
#'
#' Fits per-fold hypervolumes (Gaussian KDE or one-class SVM), saves/loads a combined RDS, plots (optional), and returns volumes and overlaps.
#'
#' @param partition_results List or .rds path from spatiotemporal_partition().
#' @param model_vars Character vector of predictor column names.
#' @param method "gaussian" or "svm".
#' @param output_dir Output directory.
#' @param hypervolume_params Named list of args passed to the selected hypervolume function.
#' @param create_plot Logical.
#' @param overwrite Logical.
#'
#' @return A list with elements: \code{hypervolumes}, \code{volumes}, \code{overlaps}, \code{method}, \code{model_vars}, \code{output_dir}, \code{combined_file}.
#'
#' @export
#' @importFrom hypervolume hypervolume_gaussian hypervolume_svm hypervolume_join
#' @importFrom hypervolume hypervolume_set hypervolume_overlap_statistics get_volume
#' @importFrom sf st_drop_geometry
#' @importFrom utils readRDS saveRDS combn
#' @importFrom tools file_ext
#' @importFrom grDevices rainbow
#' @importFrom graphics plot
build_hypervolume_models <- function(partition_results,
                                     model_vars,
                                     method,
                                     output_dir = "./Hypervolume_Models",
                                     hypervolume_params = list(),
                                     create_plot = TRUE,
                                     overwrite = FALSE) {

  require(hypervolume)
  require(dplyr)
  require(sf)

  ### VALIDATE method INPUT

  if (missing(method)) {
    stop("method is required. Expecting 'svm' or 'gaussian'")
  }

  if (!method %in% c("gaussian", "svm")) {
    stop(paste("method must be 'gaussian' or 'svm', got:", method))
  }

  print(paste("Building hypervolume models using method:", method))

  ### VALIDATE partition_results INPUT

  if (missing(partition_results)) {
    stop("partition_results is required. Provide output from spatiotemporal_partition() as list or .rds file path")
  }

  if (is.character(partition_results)) {
    if (!file.exists(partition_results)) {
      stop(paste("File not found:", partition_results, "Working directory:", getwd()))
    }

    file_ext <- tolower(tools::file_ext(partition_results))
    if (file_ext != "rds") {
      stop(paste("partition_results must be .rds format, got:", file_ext))
    }

    print(paste("Reading partition results from:", basename(partition_results)))
    tryCatch({
      partition_data <- readRDS(partition_results)
    }, error = function(e) {
      stop(paste("Error reading .rds file:", e$message))
    })

  } else if (is.list(partition_results)) {
    print("Using provided partition results list...")
    partition_data <- partition_results
  } else {
    stop(paste("partition_results must be file path or list, got:", class(partition_results)[1]))
  }

  ### VALIDATE partition_data STRUCTURE

  if (!"points_sf" %in% names(partition_data)) {
    stop(paste("partition_results missing 'points_sf' element. Available:",
               paste(names(partition_data), collapse = ", ")))
  }

  occurrence_points <- partition_data$points_sf

  if (is.null(occurrence_points) || nrow(occurrence_points) == 0) {
    stop("partition_results$points_sf is empty")
  }

  if (!"fold" %in% names(occurrence_points)) {
    stop(paste("Missing 'fold' column. Available columns:",
               paste(names(occurrence_points)[names(occurrence_points) != "geometry"], collapse = ", ")))
  }

  if (inherits(occurrence_points, "sf")) {
    folds <- st_drop_geometry(occurrence_points)$fold
  } else {
    folds <- occurrence_points$fold
  }

  if (all(is.na(folds))) {
    stop("All fold values are NA. Re-run spatiotemporal_partition()")
  }

  if (length(unique(folds[!is.na(folds)])) < 2) {
    warning(paste("Only", length(unique(folds[!is.na(folds)])), "unique fold detected"))
  }

  print(paste("Loaded", nrow(occurrence_points), "points from partition results"))

  ### VALIDATE model_vars INPUT

  if (missing(model_vars)) {
    available_vars <- colnames(occurrence_points)
    available_vars <- available_vars[available_vars != "geometry"]
    stop(paste("model_vars required. Available:", paste(available_vars, collapse = ", ")))
  }

  if (!is.character(model_vars) || length(model_vars) == 0) {
    stop("model_vars must be character vector with at least one variable")
  }

  if (!all(model_vars %in% colnames(occurrence_points))) {
    missing_vars <- model_vars[!model_vars %in% colnames(occurrence_points)]
    available_vars <- colnames(occurrence_points)
    available_vars <- available_vars[available_vars != "geometry"]
    stop(paste("Missing variables:", paste(missing_vars, collapse = ", "),
               "Available:", paste(available_vars, collapse = ", ")))
  }

  if (length(folds) != nrow(occurrence_points)) {
    stop(paste("Fold length mismatch:", length(folds), "vs", nrow(occurrence_points)))
  }

  ### VALIDATE hypervolume_params

  if (length(hypervolume_params) > 0) {
    if (method == "gaussian") {
      valid_params <- c("kde.bandwidth", "quantile.requested", "quantile.requested.type",
                        "chunk.size", "verbose", "samples.per.point")
      invalid_params <- setdiff(names(hypervolume_params), valid_params)
      if (length(invalid_params) > 0) {
        warning(paste("Invalid parameters for gaussian method:", paste(invalid_params, collapse = ", "),
                      "Valid parameters:", paste(valid_params, collapse = ", ")))
      }

      if ("quantile.requested" %in% names(hypervolume_params)) {
        qr <- hypervolume_params$quantile.requested
        if (!is.numeric(qr) || qr <= 0 || qr > 1) {
          warning(paste("quantile.requested should be between 0 and 1, got:", qr))
        }
      }

      if ("quantile.requested.type" %in% names(hypervolume_params)) {
        qrt <- hypervolume_params$quantile.requested.type
        if (!qrt %in% c("probability", "volume")) {
          warning(paste("quantile.requested.type should be 'probability' or 'volume', got:", qrt))
        }
      }

    } else if (method == "svm") {
      valid_params <- c("svm.nu", "svm.gamma", "chunk.size", "verbose", "samples.per.point")
      invalid_params <- setdiff(names(hypervolume_params), valid_params)
      if (length(invalid_params) > 0) {
        warning(paste("Invalid parameters for svm method:", paste(invalid_params, collapse = ", "),
                      "Valid parameters:", paste(valid_params, collapse = ", ")))
      }

      if ("svm.nu" %in% names(hypervolume_params)) {
        nu <- hypervolume_params$svm.nu
        if (!is.numeric(nu) || nu <= 0 || nu > 1) {
          warning(paste("svm.nu should be between 0 and 1, got:", nu))
        }
      }

      if ("svm.gamma" %in% names(hypervolume_params)) {
        gamma <- hypervolume_params$svm.gamma
        if (!is.numeric(gamma) || gamma <= 0) {
          warning(paste("svm.gamma should be positive, got:", gamma))
        }
      }
    }

    if ("samples.per.point" %in% names(hypervolume_params)) {
      spp <- hypervolume_params$samples.per.point
      if (!is.numeric(spp) || spp < 1) {
        warning(paste("samples.per.point should be positive integer, got:", spp))
      }
    }
  }

  ### CREATE OUTPUT DIRECTORY

  tryCatch({
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      print(paste("Created output directory:", output_dir))
    }
  }, error = function(e) {
    stop(paste("Could not create output directory:", output_dir, "-", e$message))
  })

  ### PREPARE DATA

  occ_df <- st_drop_geometry(occurrence_points)
  occ_df$fold <- folds

  n_original <- nrow(occ_df)
  na_summary <- sapply(model_vars, function(v) sum(is.na(occ_df[[v]])))

  if (any(na_summary > 0)) {
    na_vars <- names(na_summary[na_summary > 0])
    warning(paste("NAs in:", paste(paste(na_vars, na_summary[na_summary > 0], sep = "="), collapse = ", ")))
  }

  occ_df <- occ_df[complete.cases(occ_df[, model_vars]), ]
  n_removed <- n_original - nrow(occ_df)

  if (n_removed > 0) {
    pct_removed <- round(n_removed / n_original * 100, 2)
    warning(paste0("Removed ", n_removed, " rows (", pct_removed, "%) with NAs. Remaining: ", nrow(occ_df)))
  }

  if (nrow(occ_df) == 0) {
    stop("No complete cases after removing NAs. Check environmental data extraction")
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

  ### CHECK FOR EXISTING COMBINED FILE

  combined_file <- file.path(output_dir, paste0("all_hypervolumes_", method, ".rds"))

  if (file.exists(combined_file) && !overwrite) {
    print(paste("Loading existing combined hypervolume file:", basename(combined_file)))
    tryCatch({
      hv_list <- readRDS(combined_file)
      print(paste("Loaded", length(hv_list), "hypervolumes from file"))

      for (i in seq_along(hv_list)) {
        fold_name <- names(hv_list)[i]
        fold_id <- gsub("fold", "", fold_name)
        if (is.null(hv_list[[i]]@Name) || hv_list[[i]]@Name == "untitled") {
          hv_list[[i]]@Name <- paste("Fold", fold_id)
        }
      }

    }, error = function(e) {
      stop(paste("Error loading combined hypervolume file. Try overwrite = TRUE:", e$message))
    })

  } else {

    ### BUILD HYPERVOLUMES FOR EACH FOLD

    hv_list <- list()
    for (fold in fold_ids) {

      print(paste("Building hypervolume for Fold", fold))

      train_data <- occ_df[occ_df$fold != fold, model_vars]
      print(paste("Training with", nrow(train_data), "points"))

      if (nrow(train_data) < 5) {
        stop(paste("Fold", fold, "has only", nrow(train_data), "training points. Minimum 5 required"))
      }

      if (nrow(train_data) < 10) {
        warning(paste("Fold", fold, "has only", nrow(train_data), "training points. Results may be unreliable"))
      }

      var_check <- sapply(train_data, var)
      if (any(var_check == 0)) {
        zero_var <- names(var_check[var_check == 0])
        stop(paste("Zero variance in Fold", fold, "for:", paste(zero_var, collapse = ", ")))
      }

      hv_params <- modifyList(
        if (method == "gaussian") default_params_gaussian else default_params_svm,
        hypervolume_params
      )

      fn <- if (method == "gaussian") hypervolume_gaussian else hypervolume_svm

      tryCatch({
        hv <- suppressMessages(do.call(fn, c(list(data = train_data), hv_params)))
        hv@Name <- paste("Fold", fold)
      }, error = function(e) {
        stop(paste("Hypervolume construction failed for Fold", fold, ":", e$message))
      })

      hv_list[[paste0("fold", fold)]] <- hv

      tryCatch({
        vol <- get_volume(hv)
        print(paste("Fold", fold, "volume:", round(vol, 4)))
      }, error = function(e) {
        warning(paste("Could not calculate volume for Fold", fold))
      })
    }

    ### SAVE COMBINED HYPERVOLUME FILE

    tryCatch({
      saveRDS(hv_list, combined_file)
      print(paste("Saved combined hypervolumes to:", basename(combined_file)))
    }, error = function(e) {
      warning(paste("Could not save combined hypervolume file:", e$message))
    })
  }

  ### CREATE PLOTS

  if (create_plot && length(model_vars) >= 2) {
    print("Creating hypervolume plots")

    for (fold in fold_ids) {
      tryCatch({
        suppressMessages(
          plot(hv_list[[paste0("fold", fold)]], pairplot = TRUE, show.3d = FALSE,
               main = paste("Fold", fold, "Hypervolume"))
        )
      }, error = function(e) {
        warning(paste("Could not create plot for Fold", fold))
      })
    }

    if (length(hv_list) > 1) {
      tryCatch({
        hv_joined <- suppressMessages(do.call(hypervolume_join, hv_list))

        suppressMessages(
          plot(hv_joined, pairplot = TRUE, show.3d = FALSE,
               main = paste("Hypervolume Comparison:", length(hv_list), "Folds"),
               colors = rainbow(length(hv_list)))
        )
      }, error = function(e) {
        warning("Could not create comparison plot")
      })
    }
  } else if (create_plot && length(model_vars) < 2) {
    print("Skipping plots: need at least 2 variables")
  }

  ### SUMMARY STATISTICS

  print("SUMMARY STATISTICS")

  volumes <- sapply(hv_list, function(h) {
    tryCatch(get_volume(h), error = function(e) NA)
  })
  print(volumes)

  ### PAIRWISE OVERLAPS

  overlap_stats <- list()
  if (length(hv_list) > 1) {
    print("Calculating pairwise overlaps...")
    fold_pairs <- combn(names(hv_list), 2, simplify = FALSE)

    for (pair in fold_pairs) {
      pair_label <- paste("Fold", gsub("fold", "", pair[1]), "vs Fold", gsub("fold", "", pair[2]))

      tryCatch({
        hv_set <- suppressMessages(
          hypervolume_set(hv_list[[pair[1]]], hv_list[[pair[2]]],
                          check.memory = FALSE, verbose = FALSE)
        )
        stats <- suppressMessages(hypervolume_overlap_statistics(hv_set))
        overlap_stats[[pair_label]] <- stats
      }, error = function(e) {
        warning(paste("Could not calculate overlap for", pair_label))
      })
    }
  }

  print("Processing complete!")

  return(list(
    hypervolumes = hv_list,
    volumes = volumes,
    overlaps = overlap_stats,
    method = method,
    model_vars = model_vars,
    output_dir = output_dir,
    combined_file = combined_file
  ))
}
