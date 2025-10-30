spatiotemporal_partition <- function(
    reference_shapefile_path,
    points_file_path,
    time_col = NULL,
    xcol = NULL,
    ycol = NULL,
    points_crs = NULL,
    total_folds = 4,
    n_temporal = 2,
    n_spatial = 8,
    blocking_priority = NULL,
    max_imbalance = 0.2,
    max_retries = 10,
    generate_plots = TRUE,
    output_file = NULL
) {

  require(sf)
  require(dplyr)
  require(ggplot2)
  require(viridis)
  require(gridExtra)
  require(deldir)

  ### INPUT VALIDATION

  if (is.null(reference_shapefile_path)) {
    stop("ERROR: 'reference_shapefile_path' is required.")
  }

  if (is.null(points_file_path)) {
    stop("ERROR: 'points_file_path' is required.")
  }

  if (is.null(blocking_priority)) {
    stop("ERROR: 'blocking_priority' is required. Must be one of: 'balanced', 'spatial', or 'temporal'.")
  }

  # Handle missing time_col
  if (is.null(time_col)) {
    warning("WARNING: 'time_col' not provided. Setting n_temporal to 1.")
    n_temporal <- 1
    temporal_partitioning <- FALSE
  } else {
    temporal_partitioning <- TRUE
  }

  if (!blocking_priority %in% c("balanced", "spatial", "temporal")) {
    stop("ERROR: 'blocking_priority' must be one of: 'balanced', 'spatial', or 'temporal'")
  }

  if (!temporal_partitioning && blocking_priority == "temporal") {
    stop("ERROR: Cannot use blocking_priority='temporal' without specifying a time_col.")
  }

  if (!temporal_partitioning && blocking_priority == "balanced") {
    warning("WARNING: blocking_priority='balanced' requires time_col. Switching to 'spatial'.")
    blocking_priority <- "spatial"
  }

  if (is.null(n_spatial)) {
    n_spatial <- total_folds * n_temporal
  }

  if (n_spatial <= 0) n_spatial <- 1
  if (n_temporal <= 0) n_temporal <- 1

  ### LOAD REFERENCE SHAPEFILE

  if (inherits(reference_shapefile_path, "sf")) {
    reference_shapefile <- reference_shapefile_path
  } else {
    if (!file.exists(reference_shapefile_path)) {
      stop(paste0("ERROR: Reference shapefile not found at: ", reference_shapefile_path))
    }
    reference_shapefile <- st_read(reference_shapefile_path, quiet = TRUE)
  }

  ### LOAD AND CONVERT POINTS DATA

  # Handle points_file_path input
  if (is.character(points_file_path)) {
    if (!file.exists(points_file_path)) {
      stop(paste("ERROR: File does not exist:", points_file_path))
    }

    file_ext <- tolower(tools::file_ext(points_file_path))

    if (file_ext == "csv") {

      # Require xcol, ycol, and points_crs for CSV files
      if (is.null(xcol)) {
        stop("ERROR: 'xcol' is required when reading CSV files.")
      }
      if (is.null(ycol)) {
        stop("ERROR: 'ycol' is required when reading CSV files.")
      }
      if (is.null(points_crs)) {
        stop("ERROR: 'points_crs' is required when reading CSV files.")
      }

      print(paste("Reading CSV file:", basename(points_file_path)))
      pts <- read.csv(points_file_path, stringsAsFactors = FALSE)

      if (!xcol %in% names(pts)) {
        stop(paste0("ERROR: Column '", xcol, "' not found in CSV. Available columns: ",
                    paste(names(pts), collapse = ", ")))
      }
      if (!ycol %in% names(pts)) {
        stop(paste0("ERROR: Column '", ycol, "' not found in CSV. Available columns: ",
                    paste(names(pts), collapse = ", ")))
      }

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_file_path)))
      pts_sf_raw <- st_read(points_file_path, quiet = TRUE)
      pts <- st_drop_geometry(pts_sf_raw)

      # Extract coordinates
      coords <- st_coordinates(pts_sf_raw)
      if (is.null(xcol)) xcol <- "X"
      if (is.null(ycol)) ycol <- "Y"
      pts[[xcol]] <- coords[, 1]
      pts[[ycol]] <- coords[, 2]

      # Get CRS if not provided
      if (is.null(points_crs)) {
        points_crs <- st_crs(pts_sf_raw)
      }

    } else {
      stop(paste("ERROR: Unsupported file format:", file_ext,
                 "\nSupported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (inherits(points_file_path, "sf")) {

    # If it's already an sf object
    print("Using provided sf object...")
    pts <- st_drop_geometry(points_file_path)

    # Extract coordinates
    coords <- st_coordinates(points_file_path)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]

    # Get CRS if not provided
    if (is.null(points_crs)) {
      points_crs <- st_crs(points_file_path)
    }

  } else if (is.data.frame(points_file_path)) {

    # Require xcol, ycol, and points_crs for data frames
    if (is.null(xcol)) {
      stop("ERROR: 'xcol' is required when providing a data frame.")
    }
    if (is.null(ycol)) {
      stop("ERROR: 'ycol' is required when providing a data frame.")
    }
    if (is.null(points_crs)) {
      stop("ERROR: 'points_crs' is required when providing a data frame.")
    }

    print("Using provided data frame...")
    pts <- points_file_path

    if (!xcol %in% names(pts)) {
      stop(paste0("ERROR: Column '", xcol, "' not found in data frame. Available columns: ",
                  paste(names(pts), collapse = ", ")))
    }
    if (!ycol %in% names(pts)) {
      stop(paste0("ERROR: Column '", ycol, "' not found in data frame. Available columns: ",
                  paste(names(pts), collapse = ", ")))
    }

  } else {
    stop("ERROR: points_file_path must be an sf object, data frame with x/y columns, or file path")
  }

  ### DATA CLEANING

  n_original <- nrow(pts)
  pts_complete <- pts[complete.cases(pts), ]
  n_removed <- n_original - nrow(pts_complete)

  if (n_removed > 0) {
    pct_removed <- round(n_removed/n_original * 100, 2)
    warning(paste0("Removed ", n_removed, " incomplete rows (", pct_removed, "%)"))
  }

  if (nrow(pts_complete) == 0) {
    stop("ERROR: No complete rows remaining.")
  }

  ### VALIDATE REQUIRED COLUMNS

  if (!xcol %in% names(pts_complete)) {
    stop(paste0("ERROR: Column '", xcol, "' not found. Available: ", paste(names(pts_complete), collapse = ", ")))
  }
  if (!ycol %in% names(pts_complete)) {
    stop(paste0("ERROR: Column '", ycol, "' not found. Available: ", paste(names(pts_complete), collapse = ", ")))
  }
  if (temporal_partitioning && !time_col %in% names(pts_complete)) {
    stop(paste0("ERROR: Column '", time_col, "' not found. Available: ", paste(names(pts_complete), collapse = ", ")))
  }

  ### CREATE SF OBJECT

  pts_sf <- st_as_sf(pts_complete, coords = c(xcol, ycol), crs = points_crs)
  pts_sf <- st_transform(pts_sf, crs = st_crs(reference_shapefile))

  total_points <- nrow(pts_sf)

  if (total_points < total_folds) {
    stop(paste0("ERROR: Not enough points (", total_points, ") for ", total_folds, " folds."))
  }

  ### RETRY LOOP

  best_imbalance <- Inf
  best_results <- NULL

  for (attempt in 1:max_retries) {

    if (attempt > 1) {
      cat(paste0("\n--- Attempt ", attempt, "/", max_retries, " ---\n"))
    }

    pts_sf_attempt <- pts_sf
    pts_sf_attempt$spatial_block <- NULL
    pts_sf_attempt$temporal_block <- NULL
    pts_sf_attempt$fold <- NULL
    pts_sf_attempt$block_type <- NULL

    ### FOLD STRUCTURE

    if (blocking_priority == "spatial") {
      n_spatially_exclusive_folds <- total_folds
      n_temporally_exclusive_folds <- 0
      min_exclusive_blocks <- max(1, floor(n_spatial / 2))
      n_shared_blocks <- n_spatial - min_exclusive_blocks

      if (attempt == 1) {
        cat("BLOCKING PRIORITY: SPATIAL\n")
        cat(paste0("  ", total_folds, " spatially exclusive folds\n"))
        if (temporal_partitioning) {
          cat(paste0("  ", min_exclusive_blocks, " dedicated + ", n_shared_blocks, " shared blocks\n\n"))
        } else {
          cat(paste0("  ", n_spatial, " spatial blocks\n\n"))
        }
      }

    } else if (blocking_priority == "temporal") {
      n_temporally_exclusive_folds <- total_folds
      n_spatially_exclusive_folds <- 0
      min_exclusive_blocks <- 0
      n_shared_blocks <- n_spatial

      if (attempt == 1) {
        cat("BLOCKING PRIORITY: TEMPORAL\n")
        cat(paste0("  ", total_folds, " temporally exclusive folds\n"))
        cat(paste0("  ", n_spatial, " shared spatial blocks\n\n"))
      }

    } else {
      n_spatially_exclusive_folds <- max(0, total_folds - n_temporal)
      n_temporally_exclusive_folds <- min(n_temporal, total_folds)
      min_exclusive_blocks <- n_spatially_exclusive_folds
      n_shared_blocks <- n_spatial - min_exclusive_blocks

      if (attempt == 1) {
        cat("BLOCKING PRIORITY: BALANCED\n")
        cat(paste0("  ", n_spatially_exclusive_folds, " spatially exclusive folds\n"))
        cat(paste0("  ", n_temporally_exclusive_folds, " temporally exclusive folds\n\n"))
      }
    }

    ### SPATIAL PARTITIONING

    if (n_spatial == 1) {
      pts_sf_attempt$spatial_block <- 1
      voronoi_sf <- st_sf(
        spatial_block = 1,
        geometry = st_geometry(st_union(reference_shapefile))
      )
    } else {
      coords <- st_coordinates(pts_sf_attempt)
      kmeans_result <- kmeans(coords, centers = n_spatial, nstart = 50, iter.max = 100)
      pts_sf_attempt$spatial_block <- kmeans_result$cluster

      spatial_centers <- data.frame(
        spatial_block = 1:n_spatial,
        center_x = kmeans_result$centers[, 1],
        center_y = kmeans_result$centers[, 2]
      )

      study_bbox <- st_bbox(reference_shapefile)

      voronoi_deldir <- deldir(
        spatial_centers$center_x,
        spatial_centers$center_y,
        rw = c(study_bbox["xmin"], study_bbox["xmax"],
               study_bbox["ymin"], study_bbox["ymax"])
      )
      voronoi_tiles <- tile.list(voronoi_deldir)

      voronoi_polygons <- lapply(seq_along(voronoi_tiles), function(i) {
        tile <- voronoi_tiles[[i]]
        coords <- cbind(c(tile$x, tile$x[1]), c(tile$y, tile$y[1]))
        st_polygon(list(coords))
      })

      voronoi_sf <- st_sf(
        spatial_block = 1:n_spatial,
        geometry = st_sfc(voronoi_polygons, crs = st_crs(pts_sf_attempt))
      )

      voronoi_sf <- suppressWarnings(st_intersection(voronoi_sf, st_union(reference_shapefile)))
    }

    ### TEMPORAL PARTITIONING

    if (temporal_partitioning) {
      temporal_values <- pts_sf_attempt[[time_col]]

      if (n_temporal == 1) {
        pts_sf_attempt$temporal_block <- 1
      } else if (n_temporal == 2) {
        temporal_threshold <- mean(temporal_values, na.rm = TRUE)
        pts_sf_attempt$temporal_block <- ifelse(temporal_values <= temporal_threshold, 1, 2)
      } else {
        temporal_breaks <- quantile(temporal_values,
                                    probs = seq(0, 1, length.out = n_temporal + 1),
                                    na.rm = TRUE)
        if (any(duplicated(temporal_breaks))) {
          temporal_breaks <- unique(temporal_breaks)
          warning("Duplicate temporal breaks detected")
        }
        pts_sf_attempt$temporal_block <- cut(temporal_values,
                                             breaks = temporal_breaks,
                                             labels = FALSE,
                                             include.lowest = TRUE)
      }
    } else {
      pts_sf_attempt$temporal_block <- 1
    }

    ### FOLD ASSIGNMENT

    pts_sf_attempt$fold <- NA
    pts_sf_attempt$block_type <- NA

    if (blocking_priority == "spatial") {

      all_spatial_blocks <- 1:n_spatial

      if (temporal_partitioning && n_shared_blocks > 0) {
        dedicated_spatial_blocks <- sample(all_spatial_blocks, min_exclusive_blocks, replace = FALSE)
        shared_spatial_blocks <- setdiff(all_spatial_blocks, dedicated_spatial_blocks)
      } else {
        dedicated_spatial_blocks <- all_spatial_blocks
        shared_spatial_blocks <- c()
      }

      dedicated_block_to_fold <- data.frame(
        spatial_block = dedicated_spatial_blocks,
        fold = rep(1:total_folds, length.out = length(dedicated_spatial_blocks))
      )
      dedicated_block_to_fold <- dedicated_block_to_fold[sample(nrow(dedicated_block_to_fold)), ]

      for (i in 1:nrow(dedicated_block_to_fold)) {
        s <- dedicated_block_to_fold$spatial_block[i]
        f <- dedicated_block_to_fold$fold[i]

        block_idx <- which(pts_sf_attempt$spatial_block == s)
        pts_sf_attempt$fold[block_idx] <- f
        pts_sf_attempt$block_type[block_idx] <- "dedicated_spatial"
      }

      if (temporal_partitioning && length(shared_spatial_blocks) > 0) {
        for (s in shared_spatial_blocks) {
          for (t in 1:n_temporal) {
            temp_idx <- which(pts_sf_attempt$spatial_block == s &
                                pts_sf_attempt$temporal_block == t &
                                is.na(pts_sf_attempt$fold))

            if (length(temp_idx) > 0) {
              chosen_fold <- ((t - 1) %% total_folds) + 1
              pts_sf_attempt$fold[temp_idx] <- chosen_fold
              pts_sf_attempt$block_type[temp_idx] <- "shared_temporal_split"
            }
          }
        }
      }

    } else if (blocking_priority == "temporal") {

      temporal_fold_mapping <- data.frame(
        fold = 1:total_folds,
        temporal_block = rep(1:n_temporal, length.out = total_folds)
      )

      for (t in 1:n_temporal) {
        folds_for_temporal <- temporal_fold_mapping$fold[temporal_fold_mapping$temporal_block == t]

        temp_points <- which(pts_sf_attempt$temporal_block == t & is.na(pts_sf_attempt$fold))

        if (length(temp_points) > 0) {
          fold_assignments <- rep(folds_for_temporal, length.out = length(temp_points))
          pts_sf_attempt$fold[temp_points] <- fold_assignments
          pts_sf_attempt$block_type[temp_points] <- "temporal_exclusive"
        }
      }

    } else {

      all_spatial_blocks <- 1:n_spatial

      if (min_exclusive_blocks > 0) {
        dedicated_spatial_blocks <- sample(all_spatial_blocks, min_exclusive_blocks, replace = FALSE)
        shared_spatial_blocks <- setdiff(all_spatial_blocks, dedicated_spatial_blocks)

        dedicated_block_to_fold <- data.frame(
          spatial_block = dedicated_spatial_blocks,
          fold = 1:n_spatially_exclusive_folds
        )
        dedicated_block_to_fold <- dedicated_block_to_fold[sample(nrow(dedicated_block_to_fold)), ]

        for (i in 1:nrow(dedicated_block_to_fold)) {
          s <- dedicated_block_to_fold$spatial_block[i]
          f <- dedicated_block_to_fold$fold[i]
          block_idx <- which(pts_sf_attempt$spatial_block == s)
          pts_sf_attempt$fold[block_idx] <- f
          pts_sf_attempt$block_type[block_idx] <- "dedicated_spatial"
        }
      } else {
        shared_spatial_blocks <- all_spatial_blocks
      }

      temporally_exclusive_folds <- (total_folds - n_temporally_exclusive_folds + 1):total_folds
      temporal_fold_mapping <- data.frame(
        fold = temporally_exclusive_folds,
        temporal_block = 1:n_temporally_exclusive_folds
      )

      for (s in shared_spatial_blocks) {
        for (t in 1:n_temporal) {
          temp_idx <- which(pts_sf_attempt$spatial_block == s &
                              pts_sf_attempt$temporal_block == t &
                              is.na(pts_sf_attempt$fold))

          if (length(temp_idx) > 0) {
            eligible_folds <- 1:total_folds

            for (tex_fold in temporally_exclusive_folds) {
              designated_temporal <- temporal_fold_mapping$temporal_block[temporal_fold_mapping$fold == tex_fold]
              if (designated_temporal != t) {
                eligible_folds <- setdiff(eligible_folds, tex_fold)
              }
            }

            current_counts <- table(factor(pts_sf_attempt$fold[!is.na(pts_sf_attempt$fold)], levels = 1:total_folds))
            eligible_counts <- current_counts[eligible_folds]
            chosen_fold <- eligible_folds[which.min(eligible_counts)]

            pts_sf_attempt$fold[temp_idx] <- chosen_fold
            pts_sf_attempt$block_type[temp_idx] <- "shared_balanced"
          }
        }
      }
    }

    unassigned_idx <- which(is.na(pts_sf_attempt$fold))
    if (length(unassigned_idx) > 0) {
      for (idx in unassigned_idx) {
        current_counts <- table(factor(pts_sf_attempt$fold[!is.na(pts_sf_attempt$fold)], levels = 1:total_folds))
        chosen_fold <- which.min(current_counts)[1]
        pts_sf_attempt$fold[idx] <- chosen_fold
        pts_sf_attempt$block_type[idx] <- "remainder"
      }
    }

    final_fold_counts <- table(factor(pts_sf_attempt$fold, levels = 1:total_folds))

    ### BALANCE ASSESSMENT

    mean_per_fold <- mean(final_fold_counts)
    max_deviation <- max(abs(final_fold_counts - mean_per_fold))
    imbalance <- max_deviation / mean_per_fold

    if (imbalance < best_imbalance) {
      best_imbalance <- imbalance
      best_results <- list(
        pts_sf = pts_sf_attempt,
        voronoi_sf = voronoi_sf,
        final_fold_counts = final_fold_counts,
        mean_per_fold = mean_per_fold,
        imbalance = imbalance,
        n_spatially_exclusive_folds = n_spatially_exclusive_folds,
        n_temporally_exclusive_folds = n_temporally_exclusive_folds,
        min_exclusive_blocks = min_exclusive_blocks,
        n_shared_blocks = n_shared_blocks,
        best_attempt = attempt
      )
    }

    cat(paste0("  Imbalance: ", round(imbalance * 100, 2), "%\n"))
  }

  pts_sf <- best_results$pts_sf
  voronoi_sf <- best_results$voronoi_sf
  final_fold_counts <- best_results$final_fold_counts
  mean_per_fold <- best_results$mean_per_fold
  imbalance <- best_results$imbalance
  n_spatially_exclusive_folds <- best_results$n_spatially_exclusive_folds
  n_temporally_exclusive_folds <- best_results$n_temporally_exclusive_folds
  min_exclusive_blocks <- best_results$min_exclusive_blocks
  n_shared_blocks <- best_results$n_shared_blocks
  best_attempt <- best_results$best_attempt

  cat(paste0("\n✓ Best balance: attempt ", best_attempt, "/", max_retries,
             " (", round(imbalance * 100, 2), "%)\n"))

  if (imbalance > max_imbalance) {
    warning(paste0("Imbalance ", round(imbalance * 100, 2),
                   "% exceeds threshold ", round(max_imbalance * 100, 2), "%"))
  }

  ### REPORTING

  cat("\n=== FOLD STRUCTURE ===\n")
  cat(paste0(total_folds, " folds | ", n_spatial, " spatial blocks"))
  if (temporal_partitioning) {
    cat(paste0(" | ", n_temporal, " temporal blocks\n"))
  } else {
    cat(" | No temporal blocks\n")
  }
  cat(paste0("Priority: ", toupper(blocking_priority), "\n\n"))

  cat("=== FOLD SIZES ===\n")
  for (f in 1:total_folds) {
    points <- as.numeric(final_fold_counts[f])
    percent <- round(points / total_points * 100, 2)
    cat(paste0("Fold ", f, ": ", points, " (", percent, "%)\n"))
  }
  cat("\n")

  plot_list <- list()

  if (generate_plots) {

    if (!exists("study_bbox")) {
      study_bbox <- st_bbox(reference_shapefile)
    }
    bbox <- study_bbox

    if (temporal_partitioning) {
      subtitle_text <- paste0(n_spatial, " spatial | ", n_temporal, " temporal | ", toupper(blocking_priority))
    } else {
      subtitle_text <- paste0(n_spatial, " spatial | ", toupper(blocking_priority))
    }

    if (n_spatial > 1) {
      plot_list$main <- ggplot() +
        geom_sf(data = reference_shapefile, fill = "gray98", color = "gray40", linewidth = 0.5) +
        geom_sf(data = voronoi_sf, fill = NA, color = "gray60", linewidth = 0.8, linetype = "dashed") +
        geom_sf(data = pts_sf, aes(color = factor(fold)), size = 1.5, alpha = 0.7) +
        coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
                 ylim = c(bbox["ymin"], bbox["ymax"])) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray40"),
          legend.position = "right",
          panel.grid = element_line(color = "gray90", linewidth = 0.3),
          panel.background = element_rect(fill = "aliceblue", color = NA)
        ) +
        labs(title = paste0("Spatiotemporal Partitioning: ", total_folds, " Folds"),
             subtitle = subtitle_text)
    } else {
      plot_list$main <- ggplot() +
        geom_sf(data = reference_shapefile, fill = "gray98", color = "gray40", linewidth = 0.5) +
        geom_sf(data = pts_sf, aes(color = factor(fold)), size = 1.5, alpha = 0.7) +
        coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
                 ylim = c(bbox["ymin"], bbox["ymax"])) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(color = "gray40"),
          legend.position = "right",
          panel.grid = element_line(color = "gray90", linewidth = 0.3),
          panel.background = element_rect(fill = "aliceblue", color = NA)
        ) +
        labs(title = paste0("Temporal Partitioning: ", total_folds, " Folds"),
             subtitle = subtitle_text)
    }

    if (temporal_partitioning && n_temporal > 1) {
      plot_list$temporal <- ggplot(pts_sf, aes(x = get(time_col), fill = factor(temporal_block))) +
        geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
        scale_fill_viridis_d(option = "mako", name = "Temporal\nBlock") +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right") +
        labs(title = "Temporal Distribution", x = time_col, y = "Count")
    }

    fold_balance_df <- data.frame(
      fold = factor(names(final_fold_counts), levels = sort(as.numeric(names(final_fold_counts)))),
      count = as.numeric(final_fold_counts)
    )

    plot_list$balance <- ggplot(fold_balance_df, aes(x = fold, y = count)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_hline(yintercept = mean_per_fold, linetype = "dashed", color = "red", linewidth = 0.8) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold")) +
      labs(title = "Fold Balance", x = "Fold", y = "Points")

    if (n_spatial > 1 && temporal_partitioning && n_temporal > 1) {
      plot_list$combined <- ggplot() +
        geom_sf(data = reference_shapefile, fill = "gray98", color = "gray40", linewidth = 0.5) +
        geom_sf(data = voronoi_sf, fill = NA, color = "black", linewidth = 1) +
        geom_sf(data = pts_sf, aes(color = factor(fold), shape = factor(temporal_block)),
                size = 2, alpha = 0.8) +
        scale_shape_manual(values = c(16, 17, 15, 18)[1:n_temporal],
                           name = "Temporal\nBlock") +
        coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
                 ylim = c(bbox["ymin"], bbox["ymax"])) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid = element_line(color = "gray90", linewidth = 0.3),
          panel.background = element_rect(fill = "aliceblue", color = NA)
        ) +
        labs(title = "Complete Partitioning Structure")
    } else {
      plot_list$combined <- plot_list$main
    }

    print(plot_list$combined)
  }

  summary_stats <- data.frame(
    parameter = c("spatial_blocks", "temporal_blocks", "total_folds",
                  "blocking_priority", "spatially_exclusive_folds",
                  "temporally_exclusive_folds", "dedicated_spatial_blocks",
                  "shared_spatial_blocks", "total_points", "points_removed",
                  "pct_rows_removed", "max_imbalance_pct", "total_attempts",
                  "best_attempt", "temporal_partitioning_enabled"),
    value = c(n_spatial, n_temporal, total_folds, blocking_priority,
              n_spatially_exclusive_folds, n_temporally_exclusive_folds,
              min_exclusive_blocks, n_shared_blocks, nrow(pts_sf), n_removed,
              ifelse(n_removed > 0, round(n_removed/n_original * 100, 2), 0),
              round(imbalance * 100, 2), max_retries, best_attempt,
              as.character(temporal_partitioning))
  )

  folds_cols <- c("fold", "spatial_block", "temporal_block", "block_type")
  if (temporal_partitioning) {
    folds_cols <- c(folds_cols, time_col)
  }

  results <- list(
    folds = pts_sf %>% st_drop_geometry() %>% dplyr::select(all_of(folds_cols)),
    points_sf = pts_sf,
    voronoi = voronoi_sf,
    summary = summary_stats,
    plots = plot_list
  )

  if (!is.null(output_file)) {
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    saveRDS(results, output_file)
    cat(paste0("Results saved to: ", output_file, "\n\n"))
  }

  return(results)
}
