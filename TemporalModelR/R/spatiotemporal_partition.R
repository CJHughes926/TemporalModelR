spatiotemporal_partition <- function(
    reference_shapefile_path,
    points_file_path,
    time_col = NULL,
    total_folds = 4,
    n_temporal = 2,
    n_spatial = 8,
    blocking_priority = "balanced",
    max_imbalance = 0.2,
    generate_plots = TRUE,
    output_file = NULL
) {

  require(sf)
  require(dplyr)
  require(ggplot2)
  require(viridis)
  require(gridExtra)
  require(deldir)

  if (!blocking_priority %in% c("balanced", "spatial", "temporal")) {
    stop("blocking_priority must be one of: 'balanced', 'spatial', or 'temporal'")
  }

  if (is.null(n_spatial)) {
    n_spatial <- total_folds * n_temporal
  }

  if (n_spatial <= 0) n_spatial <- 1
  if (n_temporal <= 0) n_temporal <- 1

  # ============================================================================
  # CALCULATE FOLD STRUCTURE BY PRIORITY
  # ============================================================================

  if (blocking_priority == "spatial") {
    # SPATIAL: All folds spatially exclusive, none temporally exclusive
    n_spatially_exclusive_folds <- total_folds
    n_temporally_exclusive_folds <- 0
    min_exclusive_blocks <- max(1, floor(n_spatial / 2))
    n_shared_blocks <- n_spatial - min_exclusive_blocks

    cat("BLOCKING PRIORITY: SPATIAL\n")
    cat(paste0("  All ", total_folds, " folds are spatially exclusive\n"))
    cat(paste0("  ", min_exclusive_blocks, " dedicated blocks + ",
               n_shared_blocks, " shared blocks (split by temporal)\n\n"))

  } else if (blocking_priority == "temporal") {
    # TEMPORAL: All folds temporally exclusive, none spatially exclusive
    n_temporally_exclusive_folds <- total_folds
    n_spatially_exclusive_folds <- 0
    min_exclusive_blocks <- 0
    n_shared_blocks <- n_spatial

    cat("BLOCKING PRIORITY: TEMPORAL\n")
    cat(paste0("  All ", total_folds, " folds are temporally exclusive\n"))
    cat(paste0("  All ", n_spatial, " spatial blocks are shared across temporal folds\n\n"))

  } else {
    # BALANCED: Mix of both types
    n_spatially_exclusive_folds <- max(0, total_folds - n_temporal)
    n_temporally_exclusive_folds <- min(n_temporal, total_folds)
    min_exclusive_blocks <- n_spatially_exclusive_folds
    n_shared_blocks <- n_spatial - min_exclusive_blocks

    cat("BLOCKING PRIORITY: BALANCED\n")
    cat(paste0("  ", n_spatially_exclusive_folds, " spatially exclusive folds\n"))
    cat(paste0("  ", n_temporally_exclusive_folds, " temporally exclusive folds\n\n"))
  }

  # ============================================================================
  # LOAD DATA
  # ============================================================================

  reference_shapefile <- st_read(reference_shapefile_path, quiet = TRUE)
  pts <- read.csv(points_file_path)
  pts_sf <- st_as_sf(pts, coords = c("LONGDD", "LATDD"), crs = 4326)

  if (!time_col %in% names(pts_sf)) {
    stop(paste("Time column", time_col, "not found in data!"))
  }

  pts_sf <- st_transform(pts_sf, crs = st_crs(reference_shapefile))

  # ============================================================================
  # SPATIAL PARTITIONING
  # ============================================================================

  if (n_spatial == 1) {
    pts_sf$spatial_block <- 1
    voronoi_sf <- st_sf(
      spatial_block = 1,
      geometry = st_geometry(st_union(reference_shapefile))
    )
  } else {
    coords <- st_coordinates(pts_sf)
    kmeans_result <- kmeans(coords, centers = n_spatial, nstart = 50, iter.max = 100)
    pts_sf$spatial_block <- kmeans_result$cluster

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
      geometry = st_sfc(voronoi_polygons, crs = st_crs(pts_sf))
    )

    voronoi_sf <- suppressWarnings(st_intersection(voronoi_sf, st_union(reference_shapefile)))
  }

  # ============================================================================
  # TEMPORAL PARTITIONING
  # ============================================================================

  temporal_values <- pts_sf[[time_col]]

  if (n_temporal == 1) {
    pts_sf$temporal_block <- 1
  } else if (n_temporal == 2) {
    temporal_threshold <- mean(temporal_values, na.rm = TRUE)
    pts_sf$temporal_block <- ifelse(temporal_values <= temporal_threshold, 1, 2)
  } else {
    temporal_breaks <- quantile(temporal_values,
                                probs = seq(0, 1, length.out = n_temporal + 1),
                                na.rm = TRUE)
    if (any(duplicated(temporal_breaks))) {
      temporal_breaks <- unique(temporal_breaks)
      warning("Duplicate temporal breaks detected - using unique values")
    }
    pts_sf$temporal_block <- cut(temporal_values,
                                 breaks = temporal_breaks,
                                 labels = FALSE,
                                 include.lowest = TRUE)
  }

  # ============================================================================
  # FOLD ASSIGNMENT
  # ============================================================================

  pts_sf$fold <- NA
  pts_sf$block_type <- NA

  if (blocking_priority == "spatial") {

    # Select dedicated vs shared blocks
    all_spatial_blocks <- 1:n_spatial
    dedicated_spatial_blocks <- sample(all_spatial_blocks, min_exclusive_blocks, replace = FALSE)
    shared_spatial_blocks <- setdiff(all_spatial_blocks, dedicated_spatial_blocks)

    # Map dedicated blocks to folds evenly
    dedicated_block_to_fold <- data.frame(
      spatial_block = dedicated_spatial_blocks,
      fold = rep(1:total_folds, length.out = length(dedicated_spatial_blocks))
    )
    dedicated_block_to_fold <- dedicated_block_to_fold[sample(nrow(dedicated_block_to_fold)), ]

    # Assign dedicated blocks (100% to one fold)
    for (i in 1:nrow(dedicated_block_to_fold)) {
      s <- dedicated_block_to_fold$spatial_block[i]
      f <- dedicated_block_to_fold$fold[i]

      block_idx <- which(pts_sf$spatial_block == s)
      pts_sf$fold[block_idx] <- f
      pts_sf$block_type[block_idx] <- "dedicated_spatial"
    }

    # Assign shared blocks split by temporal period
    for (s in shared_spatial_blocks) {
      for (t in 1:n_temporal) {
        temp_idx <- which(pts_sf$spatial_block == s &
                            pts_sf$temporal_block == t &
                            is.na(pts_sf$fold))

        if (length(temp_idx) > 0) {
          chosen_fold <- ((t - 1) %% total_folds) + 1
          pts_sf$fold[temp_idx] <- chosen_fold
          pts_sf$block_type[temp_idx] <- "shared_temporal_split"
        }
      }
    }

  } else if (blocking_priority == "temporal") {

    # Assign each temporal block to a fold
    temporal_fold_mapping <- data.frame(
      fold = 1:total_folds,
      temporal_block = rep(1:n_temporal, length.out = total_folds)
    )

    # All spatial blocks are shared - split by temporal
    for (t in 1:n_temporal) {
      folds_for_temporal <- temporal_fold_mapping$fold[temporal_fold_mapping$temporal_block == t]

      temp_points <- which(pts_sf$temporal_block == t & is.na(pts_sf$fold))

      if (length(temp_points) > 0) {
        # Distribute points from this temporal block across assigned folds
        fold_assignments <- rep(folds_for_temporal, length.out = length(temp_points))
        pts_sf$fold[temp_points] <- fold_assignments
        pts_sf$block_type[temp_points] <- "temporal_exclusive"
      }
    }

  } else {
    # BALANCED: Mix of dedicated spatial and temporally exclusive folds

    all_spatial_blocks <- 1:n_spatial

    # Select dedicated spatial blocks
    if (min_exclusive_blocks > 0) {
      dedicated_spatial_blocks <- sample(all_spatial_blocks, min_exclusive_blocks, replace = FALSE)
      shared_spatial_blocks <- setdiff(all_spatial_blocks, dedicated_spatial_blocks)

      # Map dedicated blocks to spatially exclusive folds
      dedicated_block_to_fold <- data.frame(
        spatial_block = dedicated_spatial_blocks,
        fold = 1:n_spatially_exclusive_folds
      )
      dedicated_block_to_fold <- dedicated_block_to_fold[sample(nrow(dedicated_block_to_fold)), ]

      # Assign dedicated blocks
      for (i in 1:nrow(dedicated_block_to_fold)) {
        s <- dedicated_block_to_fold$spatial_block[i]
        f <- dedicated_block_to_fold$fold[i]
        block_idx <- which(pts_sf$spatial_block == s)
        pts_sf$fold[block_idx] <- f
        pts_sf$block_type[block_idx] <- "dedicated_spatial"
      }
    } else {
      shared_spatial_blocks <- all_spatial_blocks
    }

    # Temporal fold mapping
    temporally_exclusive_folds <- (total_folds - n_temporally_exclusive_folds + 1):total_folds
    temporal_fold_mapping <- data.frame(
      fold = temporally_exclusive_folds,
      temporal_block = 1:n_temporally_exclusive_folds
    )

    # Assign shared blocks with balancing
    for (s in shared_spatial_blocks) {
      for (t in 1:n_temporal) {
        temp_idx <- which(pts_sf$spatial_block == s &
                            pts_sf$temporal_block == t &
                            is.na(pts_sf$fold))

        if (length(temp_idx) > 0) {
          eligible_folds <- 1:total_folds

          # Temporally exclusive folds only accept their designated temporal period
          for (tex_fold in temporally_exclusive_folds) {
            designated_temporal <- temporal_fold_mapping$temporal_block[temporal_fold_mapping$fold == tex_fold]
            if (designated_temporal != t) {
              eligible_folds <- setdiff(eligible_folds, tex_fold)
            }
          }

          # Choose fold with fewest points
          current_counts <- table(factor(pts_sf$fold[!is.na(pts_sf$fold)], levels = 1:total_folds))
          eligible_counts <- current_counts[eligible_folds]
          chosen_fold <- eligible_folds[which.min(eligible_counts)]

          pts_sf$fold[temp_idx] <- chosen_fold
          pts_sf$block_type[temp_idx] <- "shared_balanced"
        }
      }
    }
  }

  # Assign any remaining unassigned points
  unassigned_idx <- which(is.na(pts_sf$fold))
  if (length(unassigned_idx) > 0) {
    for (idx in unassigned_idx) {
      current_counts <- table(factor(pts_sf$fold[!is.na(pts_sf$fold)], levels = 1:total_folds))
      chosen_fold <- which.min(current_counts)[1]
      pts_sf$fold[idx] <- chosen_fold
      pts_sf$block_type[idx] <- "remainder"
    }
  }

  final_fold_counts <- table(factor(pts_sf$fold, levels = 1:total_folds))

  # ============================================================================
  # BALANCE ASSESSMENT
  # ============================================================================

  mean_per_fold <- mean(final_fold_counts)
  max_deviation <- max(abs(final_fold_counts - mean_per_fold))
  imbalance <- max_deviation / mean_per_fold

  if (imbalance > max_imbalance) {
    warning(paste0("IMBALANCE WARNING: ", round(imbalance * 100, 2),
                   "% exceeds threshold of ", round(max_imbalance * 100, 2), "%,
                   try increasing the number of spatial blocks..."))
  }

  # ============================================================================
  # REPORTING
  # ============================================================================

  cat("\n=== FOLD STRUCTURE ===\n")
  cat(paste0("Total folds: ", total_folds, " | ",
             "Spatial blocks: ", n_spatial, " | ",
             "Temporal blocks: ", n_temporal, "\n"))
  cat(paste0("Priority: ", toupper(blocking_priority), "\n\n"))

  if (blocking_priority == "spatial") {
    cat(paste0("All ", total_folds, " folds are SPATIALLY EXCLUSIVE\n"))
    cat(paste0("  → ", min_exclusive_blocks, " dedicated spatial blocks\n"))
    cat(paste0("  → ", n_shared_blocks, " shared blocks (split by temporal period)\n\n"))

  } else if (blocking_priority == "temporal") {
    cat(paste0("All ", total_folds, " folds are TEMPORALLY EXCLUSIVE\n"))
    cat(paste0("  → All ", n_spatial, " spatial blocks shared\n\n"))

  } else {
    if (n_spatially_exclusive_folds > 0) {
      cat(paste0("SPATIALLY EXCLUSIVE FOLDS: ", paste(1:n_spatially_exclusive_folds, collapse = ", "), "\n"))
      cat(paste0("  → ", min_exclusive_blocks, " dedicated spatial blocks\n\n"))
    }

    if (n_temporally_exclusive_folds > 0) {
      temporally_exclusive_folds <- (total_folds - n_temporally_exclusive_folds + 1):total_folds
      cat(paste0("TEMPORALLY EXCLUSIVE FOLDS: ", paste(temporally_exclusive_folds, collapse = ", "), "\n"))
      cat(paste0("  → Use ", n_shared_blocks, " shared spatial blocks\n\n"))
    }
  }

  cat("=== FOLD SIZES ===\n")
  total_points <- nrow(pts_sf)
  for (f in 1:total_folds) {
    points <- as.numeric(final_fold_counts[f])
    percent <- round(points / total_points * 100, 2)
    cat(paste0("Fold ", f, ": ", points, " points (", percent, "%)\n"))
  }
  cat("\n")

  plot_list <- list()

  if (generate_plots) {

    if (!exists("study_bbox")) {
      study_bbox <- st_bbox(reference_shapefile)
    }
    bbox <- study_bbox

    subtitle_text <- paste0(n_spatial, " spatial blocks | ",
                            n_temporal, " temporal blocks | ",
                            toupper(blocking_priority), " priority")

    # Main fold map
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
        labs(
          title = paste0("Spatiotemporal Partitioning: ", total_folds, " Folds"),
          subtitle = subtitle_text
        )
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
        labs(
          title = paste0("Temporal Partitioning: ", total_folds, " Folds"),
          subtitle = subtitle_text
        )
    }

    # Temporal distribution
    if (n_temporal > 1) {
      plot_list$temporal <- ggplot(pts_sf, aes(x = get(time_col), fill = factor(temporal_block))) +
        geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
        scale_fill_viridis_d(option = "mako", name = "Temporal\nBlock") +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right") +
        labs(title = "Temporal Distribution", x = time_col, y = "Count")
    }

    # Fold balance
    fold_balance_df <- data.frame(
      fold = factor(names(final_fold_counts), levels = sort(as.numeric(names(final_fold_counts)))),
      count = as.numeric(final_fold_counts)
    )

    plot_list$balance <- ggplot(fold_balance_df, aes(x = fold, y = count)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_hline(yintercept = mean_per_fold, linetype = "dashed", color = "red", linewidth = 0.8) +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold")) +
      labs(title = "Fold Balance", x = "Fold", y = "Number of Points")

    # Combined map
    if (n_spatial > 1 && n_temporal > 1) {
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
                  "shared_spatial_blocks", "total_points", "max_imbalance_pct"),
    value = c(n_spatial, n_temporal, total_folds, blocking_priority,
              n_spatially_exclusive_folds, n_temporally_exclusive_folds,
              min_exclusive_blocks, n_shared_blocks, nrow(pts_sf),
              round(imbalance * 100, 2))
  )

  results <- list(
    folds = pts_sf %>% st_drop_geometry() %>% select(fold, spatial_block, temporal_block, block_type, all_of(time_col)),
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

    # Save as RDS
    saveRDS(results, output_file)
    cat(paste0("Results saved to: ", output_file, "\n\n"))
  }

  return(results)
}
