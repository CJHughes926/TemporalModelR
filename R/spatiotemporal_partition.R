#' Spatiotemporal Cross-Validation Partitioning
#'
#' Preprocessing function that partitions species occurrence data into
#' spatially and temporally balanced folds for cross-validation. Generates
#' Voronoi tessellations for spatial blocks and distributes points into folds
#' based on user-specified priorities.
#'
#' @param reference_shapefile_path Character or sf object. Path to study area
#'   polygon or sf polygon object.
#' @param points_file_path Character, sf object, or data frame. Path to
#'   occurrence data (.csv/.shp/.geojson/.gpkg) or spatial object.
#' @param time_col Character. Column name for temporal blocking. Default is
#'   NULL for spatial-only partitioning.
#' @param xcol Character. Coordinate column name when reading CSV or data frame.
#'   Optional.
#' @param ycol Character. Coordinate column name when reading CSV or data frame.
#'   Optional.
#' @param points_crs Character or CRS object. CRS for input points if not
#'   embedded. Optional.
#' @param total_folds Integer. Total number of folds for cross-validation.
#'   Default is 4.
#' @param n_temporal Integer. Number of temporal blocks. Default is 2.
#' @param n_spatial Integer. Number of spatial blocks. Default is 8.
#' @param blocking_priority Character. One of "balanced", "spatial", or
#'   "temporal". Determines fold assignment strategy. Currently only "balanced"
#'   is fully implemented.
#' @param max_imbalance Numeric. Maximum allowed fold size imbalance as a
#'   proportion. Default is 0.05.
#' @param generate_plots Logical. If TRUE, draws summary plots of fold
#'   distribution. If FALSE, skips plot generation. Default is TRUE.
#' @param output_file Character. Optional path to save results as .rds file.
#'   Default is NULL.
#'
#' @return A list containing:
#' \itemize{
#'   \item folds: data frame with fold assignments
#'   \item points_sf: sf object of occurrence points
#'   \item voronoi: sf object of Voronoi polygons
#'   \item summary: data frame of fold statistics
#'   \item plots: list of ggplot objects
#' }
#'
#' @details
#' The function creates spatial blocks using k-means clustering and optional
#' temporal blocks from time_col. Points are distributed into folds based on
#' blocking_priority to balance spatial and temporal coverage.
#'
#' Partitioned data are used as input for \code{\link{build_hypervolume_models}}
#' to ensure spatially and temporally independent cross-validation.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_rarefication}},
#' \code{\link{temporally_explicit_extraction}}
#'
#' Modeling: \code{\link{build_hypervolume_models}}
#'
#' @examples
#' \dontrun{
#' partition_results <- spatiotemporal_partition(
#'   reference_shapefile_path = "study_area.shp",
#'   points_file_path = "occurrences.csv",
#'   time_col = "Year",
#'   xcol = "longitude",
#'   ycol = "latitude",
#'   points_crs = "EPSG:4326",
#'   total_folds = 4,
#'   n_temporal = 2,
#'   n_spatial = 8,
#'   total_folds = 5,
#'   blocking_priority = "balanced"
#' )
#' }
#'
#' @export
#' @importFrom sf st_read st_as_sf st_transform st_coordinates st_crs st_bbox
#'   st_union st_intersection st_drop_geometry st_sfc
#' @importFrom dplyr select all_of
#' @importFrom ggplot2 ggplot geom_sf aes coord_sf theme_minimal theme
#'   element_text labs geom_histogram geom_col geom_hline scale_fill_viridis_d
#' @importFrom deldir deldir tile.list
#' @importFrom stats kmeans quantile
#' @importFrom tools file_ext
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
    max_imbalance = 0.05,
    generate_plots = TRUE,
    output_file = NULL
) {
  require(sf)
  require(dplyr)
  require(ggplot2)
  require(viridis)
  require(gridExtra)
  require(deldir)
  require(RColorBrewer)

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

  if (blocking_priority == "spatial") {
    stop("ERROR: 'blocking_priority' = 'spatial' is currently under construction, try with balanced")
  }

  if (blocking_priority == "temporal") {
    stop("ERROR: 'blocking_priority' = 'temporal' is currently under construction, try with balanced")
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

  if (blocking_priority == "balanced") {
    if (total_folds < (n_temporal + 1)) {
      stop(paste0("ERROR: 'blocking_priority' = 'balanced' requires at least (n_temporal + 1) folds. ",
                  "You provided total_folds = ", total_folds, " and n_temporal = ", n_temporal, ". ",
                  "Try increasing 'total_folds' or use blocking_priority = 'spatial' or 'temporal'."))
    }
    if (n_spatial < 2) {
      stop(paste0("ERROR: 'blocking_priority' = 'balanced' requires at least 2 spatial blocks. ",
                  "You provided n_spatial = ", n_spatial, ". ",
                  "Try increasing 'n_spatial' or choose blocking_priority = 'spatial' or 'temporal'."))
    }
  }

  ### LOAD REFERENCE SHAPEFILE
  if (inherits(reference_shapefile_path, "sf")) {
    print("Using provided sf object for reference shapefile...")
    reference_shapefile <- reference_shapefile_path
  } else if (inherits(reference_shapefile_path, "sfc")) {
    print("Converting sfc object to sf for reference shapefile...")
    reference_shapefile <- st_sf(geometry = reference_shapefile_path)
  } else if (inherits(reference_shapefile_path, "Spatial")) {
    print("Converting Spatial object to sf for reference shapefile...")
    reference_shapefile <- st_as_sf(reference_shapefile_path)
  } else if (is.character(reference_shapefile_path)) {
    if (!file.exists(reference_shapefile_path)) {
      stop(paste0("ERROR: Reference shapefile not found at: ", reference_shapefile_path))
    }
    print(paste("Reading reference shapefile:", basename(reference_shapefile_path)))
    reference_shapefile <- st_read(reference_shapefile_path, quiet = TRUE)
  } else {
    stop("ERROR: reference_shapefile_path must be an sf object, sfc object, Spatial object, or file path")
  }

  ### LOAD AND CONVERT POINTS DATA
  if (is.character(points_file_path)) {
    if (!file.exists(points_file_path)) {
      stop(paste("ERROR: File does not exist:", points_file_path))
    }

    file_ext <- tolower(tools::file_ext(points_file_path))

    if (file_ext == "csv") {
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
        stop(paste0("ERROR: Column '", xcol, "' not found in CSV. Available columns: ", paste(names(pts), collapse = ", ")))
      }
      if (!ycol %in% names(pts)) {
        stop(paste0("ERROR: Column '", ycol, "' not found in CSV. Available columns: ", paste(names(pts), collapse = ", ")))
      }
    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_file_path)))
      pts_sf_raw <- st_read(points_file_path, quiet = TRUE)
      pts <- st_drop_geometry(pts_sf_raw)

      coords <- st_coordinates(pts_sf_raw)
      if (is.null(xcol)) xcol <- "X"
      if (is.null(ycol)) ycol <- "Y"
      pts[[xcol]] <- coords[, 1]
      pts[[ycol]] <- coords[, 2]

      if (is.null(points_crs)) {
        points_crs <- st_crs(pts_sf_raw)
      }
    } else {
      stop(paste("ERROR: Unsupported file format:", file_ext, "Supported formats: .csv, .shp, .geojson, .gpkg"))
    }
  } else if (inherits(points_file_path, "sf")) {
    print("Using provided sf object...")
    pts <- st_drop_geometry(points_file_path)

    coords <- st_coordinates(points_file_path)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]

    if (is.null(points_crs)) {
      points_crs <- st_crs(points_file_path)
    }
  } else if (inherits(points_file_path, "sfc")) {
    print("Converting sfc object to sf...")
    points_sf_temp <- st_sf(geometry = points_file_path)
    pts <- data.frame(row.names = 1:length(points_file_path))

    coords <- st_coordinates(points_sf_temp)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]

    if (is.null(points_crs)) {
      points_crs <- st_crs(points_sf_temp)
    }
  } else if (inherits(points_file_path, "Spatial")) {
    print("Converting Spatial object to sf...")
    points_sf_temp <- st_as_sf(points_file_path)
    pts <- st_drop_geometry(points_sf_temp)

    coords <- st_coordinates(points_sf_temp)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]

    if (is.null(points_crs)) {
      points_crs <- st_crs(points_sf_temp)
    }
  } else if (is.data.frame(points_file_path)) {
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
      stop(paste0("ERROR: Column '", xcol, "' not found in data frame. Available columns: ", paste(names(pts), collapse = ", ")))
    }
    if (!ycol %in% names(pts)) {
      stop(paste0("ERROR: Column '", ycol, "' not found in data frame. Available columns: ", paste(names(pts), collapse = ", ")))
    }
  } else {
    stop("ERROR: points_file_path must be an sf object, sfc object, Spatial object, data frame with x/y columns, or file path")
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

  ### INTERNAL RETRY LOOP
  max_internal_retries <- 10
  best_imbalance <- Inf
  best_results <- NULL

  for (attempt in 1:max_internal_retries) {

    pts_sf_attempt <- pts_sf
    pts_sf_attempt$spatial_block <- NULL
    pts_sf_attempt$temporal_block <- NULL
    pts_sf_attempt$fold <- NULL
    pts_sf_attempt$block_type <- NULL

    ### FOLD STRUCTURE
    n_spatially_exclusive_folds <- max(0, total_folds - n_temporal)
    n_temporally_exclusive_folds <- min(n_temporal, total_folds)
    n_spatial_exclusive_blocks <- n_spatially_exclusive_folds
    n_shared_blocks <- n_spatial - n_spatial_exclusive_blocks

    if (attempt == 1) {
      print("BLOCKING PRIORITY: BALANCED")
      print(paste0(" ", n_spatially_exclusive_folds, " spatially exclusive folds"))
      print(paste0(" ", n_temporally_exclusive_folds, " temporally exclusive folds"))
    }

    ### SPATIAL PARTITIONING
    if (n_spatial == 1) {
      pts_sf_attempt$spatial_block <- 1
      voronoi_sf <- st_sf(
        spatial_block = 1,
        geometry = st_geometry(st_union(reference_shapefile))
      )
      adjacency_matrix <- matrix(TRUE, 1, 1)
      spatial_centers <- data.frame(spatial_block = 1, center_x = mean(st_coordinates(pts_sf_attempt)[,1]),
                                    center_y = mean(st_coordinates(pts_sf_attempt)[,2]))
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
        rw = c(study_bbox["xmin"], study_bbox["xmax"], study_bbox["ymin"], study_bbox["ymax"])
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

      n_blocks <- nrow(voronoi_sf)
      adjacency_matrix <- matrix(FALSE, n_blocks, n_blocks)

      for (i in 1:(n_blocks - 1)) {
        for (j in (i + 1):n_blocks) {
          intersection <- suppressWarnings(st_intersection(voronoi_sf$geometry[i], voronoi_sf$geometry[j]))
          if (length(intersection) > 0 && !st_is_empty(intersection)) {
            geom_type <- st_geometry_type(intersection)
            if (geom_type %in% c("LINESTRING", "MULTILINESTRING", "GEOMETRYCOLLECTION")) {
              adjacency_matrix[i, j] <- TRUE
              adjacency_matrix[j, i] <- TRUE
            }
          }
        }
      }
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
        temporal_breaks <- quantile(temporal_values, probs = seq(0, 1, length.out = n_temporal + 1), na.rm = TRUE)
        if (any(duplicated(temporal_breaks))) {
          temporal_breaks <- unique(temporal_breaks)
          warning("Duplicate temporal breaks detected")
        }
        pts_sf_attempt$temporal_block <- cut(temporal_values, breaks = temporal_breaks, labels = FALSE, include.lowest = TRUE)
      }
    } else {
      pts_sf_attempt$temporal_block <- 1
    }

    ### FOLD ASSIGNMENT LOGIC
    pts_sf_attempt$fold <- NA
    pts_sf_attempt$block_type <- NA

    center_coords <- cbind(spatial_centers$center_x, spatial_centers$center_y)
    distance_matrix <- as.matrix(dist(center_coords))

    spatial_exclusive_blocks <- numeric(0)
    shared_blocks <- 1:n_spatial

    if (n_spatial_exclusive_blocks > 0) {
      all_blocks <- 1:n_spatial
      block_coords <- center_coords

      selected_blocks <- numeric(n_spatial_exclusive_blocks)

      first_block_candidates <- all_blocks
      best_min_dist <- -Inf
      best_first_block <- NULL

      for (candidate_first in first_block_candidates) {
        temp_selected <- candidate_first

        for (i in 2:n_spatial_exclusive_blocks) {
          remaining <- setdiff(all_blocks, temp_selected)
          min_distances <- apply(distance_matrix[remaining, temp_selected, drop = FALSE], 1, min)
          farthest <- remaining[which.max(min_distances)]
          temp_selected <- c(temp_selected, farthest)
        }

        min_pairwise_dist <- min(distance_matrix[temp_selected, temp_selected][upper.tri(distance_matrix[temp_selected, temp_selected])])

        if (min_pairwise_dist > best_min_dist) {
          best_min_dist <- min_pairwise_dist
          best_first_block <- temp_selected
        }
      }

      selected_blocks <- best_first_block

      spatial_exclusive_blocks <- selected_blocks
      shared_blocks <- setdiff(1:n_spatial, spatial_exclusive_blocks)
    }

    spatial_fold_assignment <- data.frame(
      fold = 1:n_spatially_exclusive_folds,
      spatial_block = spatial_exclusive_blocks
    )

    fold_core_centroids <- list()
    for (i in 1:nrow(spatial_fold_assignment)) {
      block_id <- spatial_fold_assignment$spatial_block[i]
      fold_id <- spatial_fold_assignment$fold[i]

      block_points <- which(pts_sf_attempt$spatial_block == block_id)
      pts_sf_attempt$fold[block_points] <- fold_id
      pts_sf_attempt$block_type[block_points] <- "spatial_exclusive"

      fold_core_centroids[[fold_id]] <- c(
        spatial_centers$center_x[spatial_centers$spatial_block == block_id],
        spatial_centers$center_y[spatial_centers$spatial_block == block_id]
      )
    }

    block_ownership <- rep(NA, n_spatial)
    block_ownership[spatial_exclusive_blocks] <- spatial_fold_assignment$fold

    fold_priority_blocks <- list()
    for (fold_id in 1:n_spatially_exclusive_folds) {
      core_block <- spatial_fold_assignment$spatial_block[spatial_fold_assignment$fold == fold_id]
      adjacent_to_core <- which(adjacency_matrix[core_block, ])
      adjacent_shared <- intersect(adjacent_to_core, shared_blocks)
      fold_priority_blocks[[fold_id]] <- adjacent_shared
    }

    if (temporal_partitioning && n_temporally_exclusive_folds > 0) {
      temporal_folds <- (n_spatially_exclusive_folds + 1):total_folds

      for (s in shared_blocks) {
        for (t in 1:n_temporal) {
          block_points_indices <- which(pts_sf_attempt$spatial_block == s &
                                          pts_sf_attempt$temporal_block == t)

          if (length(block_points_indices) > 0) {
            temporal_fold_id <- temporal_folds[t]
            pts_sf_attempt$fold[block_points_indices] <- temporal_fold_id
            pts_sf_attempt$block_type[block_points_indices] <- "temporal_exclusive"
          }
        }
      }
    }

    initial_fold_counts <- table(factor(pts_sf_attempt$fold, levels = 1:total_folds))
    mean_per_fold_target <- mean(initial_fold_counts)

    spatial_fold_deficits <- sapply(1:n_spatially_exclusive_folds, function(f) {
      max(0, mean_per_fold_target - initial_fold_counts[f])
    })

    total_deficit <- sum(spatial_fold_deficits)

    if (total_deficit > 0 && n_temporally_exclusive_folds > 0) {
      temporal_folds <- (n_spatially_exclusive_folds + 1):total_folds

      temporal_fold_limits <- rep(0, length(temporal_folds))
      for (i in seq_along(temporal_folds)) {
        temporal_fold_id <- temporal_folds[i]
        current_count <- initial_fold_counts[temporal_fold_id]
        min_acceptable <- mean_per_fold_target * (1 - max_imbalance)
        max_transferable <- floor(current_count - min_acceptable)
        temporal_fold_limits[i] <- max(0, max_transferable)
      }

      point_coords <- st_coordinates(pts_sf_attempt)

      fold_temporal_counts <- matrix(0, nrow = n_spatially_exclusive_folds, ncol = n_temporal)
      fold_targets <- ceiling(spatial_fold_deficits)
      temporal_remaining <- temporal_fold_limits

      max_round_robin_iterations <- sum(fold_targets) * 2
      round_robin_iteration <- 0

      while (sum(fold_targets) > 0 && sum(temporal_remaining) > 0 && round_robin_iteration < max_round_robin_iterations) {
        round_robin_iteration <- round_robin_iteration + 1

        for (spatial_fold_id in 1:n_spatially_exclusive_folds) {
          if (fold_targets[spatial_fold_id] <= 0) next

          for (temporal_idx in seq_along(temporal_folds)) {
            if (temporal_remaining[temporal_idx] <= 0) next

            temporal_fold_id <- temporal_folds[temporal_idx]
            temporal_block <- temporal_idx
            fold_centroid <- fold_core_centroids[[spatial_fold_id]]

            priority_blocks <- fold_priority_blocks[[spatial_fold_id]]

            point_moved <- FALSE

            for (priority_block in priority_blocks) {
              if (!is.na(block_ownership[priority_block]) && block_ownership[priority_block] != spatial_fold_id) {
                next
              }

              candidate_indices <- which(
                pts_sf_attempt$fold == temporal_fold_id &
                  pts_sf_attempt$spatial_block == priority_block &
                  pts_sf_attempt$temporal_block == temporal_block
              )

              if (length(candidate_indices) > 0) {
                candidate_coords <- point_coords[candidate_indices, , drop = FALSE]
                distances_to_centroid <- sqrt((candidate_coords[,1] - fold_centroid[1])^2 +
                                                (candidate_coords[,2] - fold_centroid[2])^2)

                point_to_move <- candidate_indices[which.min(distances_to_centroid)]

                pts_sf_attempt$fold[point_to_move] <- spatial_fold_id
                pts_sf_attempt$block_type[point_to_move] <- "rebalanced"

                if (is.na(block_ownership[priority_block])) {
                  block_ownership[priority_block] <- spatial_fold_id
                }

                fold_temporal_counts[spatial_fold_id, temporal_block] <- fold_temporal_counts[spatial_fold_id, temporal_block] + 1
                fold_targets[spatial_fold_id] <- fold_targets[spatial_fold_id] - 1
                temporal_remaining[temporal_idx] <- temporal_remaining[temporal_idx] - 1

                point_moved <- TRUE
                break
              }
            }

            if (!point_moved) {
              available_shared <- shared_blocks[is.na(block_ownership[shared_blocks])]

              if (length(available_shared) > 0) {
                core_block <- spatial_fold_assignment$spatial_block[spatial_fold_assignment$fold == spatial_fold_id]
                block_distances <- distance_matrix[core_block, available_shared]
                closest_block <- available_shared[which.min(block_distances)]

                candidate_indices <- which(
                  pts_sf_attempt$fold == temporal_fold_id &
                    pts_sf_attempt$spatial_block == closest_block &
                    pts_sf_attempt$temporal_block == temporal_block
                )

                if (length(candidate_indices) > 0) {
                  candidate_coords <- point_coords[candidate_indices, , drop = FALSE]
                  distances_to_centroid <- sqrt((candidate_coords[,1] - fold_centroid[1])^2 +
                                                  (candidate_coords[,2] - fold_centroid[2])^2)

                  point_to_move <- candidate_indices[which.min(distances_to_centroid)]

                  pts_sf_attempt$fold[point_to_move] <- spatial_fold_id
                  pts_sf_attempt$block_type[point_to_move] <- "rebalanced"
                  block_ownership[closest_block] <- spatial_fold_id
                  fold_priority_blocks[[spatial_fold_id]] <- c(fold_priority_blocks[[spatial_fold_id]], closest_block)

                  fold_temporal_counts[spatial_fold_id, temporal_block] <- fold_temporal_counts[spatial_fold_id, temporal_block] + 1
                  fold_targets[spatial_fold_id] <- fold_targets[spatial_fold_id] - 1
                  temporal_remaining[temporal_idx] <- temporal_remaining[temporal_idx] - 1

                  point_moved <- TRUE
                }
              }
            }

            if (fold_targets[spatial_fold_id] <= 0) break
          }
        }
      }
    }

    max_balance_iterations <- 100000
    balance_iteration <- 0

    point_coords <- st_coordinates(pts_sf_attempt)

    while (balance_iteration < max_balance_iterations) {
      balance_iteration <- balance_iteration + 1

      current_fold_counts <- table(factor(pts_sf_attempt$fold, levels = 1:total_folds))
      mean_per_fold <- mean(current_fold_counts)
      max_deviation <- max(abs(current_fold_counts - mean_per_fold))
      current_imbalance <- max_deviation / mean_per_fold

      if (current_imbalance <= max_imbalance) {
        break
      }

      if (balance_iteration %% 5000 == 0) {
        print(paste0("  Balance iteration ", balance_iteration, ": ", round(current_imbalance * 100, 2), "%"))
      }

      largest_fold <- as.numeric(names(which.max(current_fold_counts)))
      smallest_fold <- as.numeric(names(which.min(current_fold_counts)))

      if (current_fold_counts[largest_fold] <= current_fold_counts[smallest_fold]) {
        break
      }

      largest_is_temporal <- largest_fold > n_spatially_exclusive_folds
      smallest_is_spatial <- smallest_fold <= n_spatially_exclusive_folds

      if (!largest_is_temporal || !smallest_is_spatial) {
        break
      }

      priority_blocks <- fold_priority_blocks[[smallest_fold]]
      fold_centroid <- fold_core_centroids[[smallest_fold]]

      moved <- FALSE
      for (priority_block in priority_blocks) {
        if (!is.na(block_ownership[priority_block]) && block_ownership[priority_block] != smallest_fold) {
          next
        }

        candidate_indices <- which(
          pts_sf_attempt$fold == largest_fold &
            pts_sf_attempt$spatial_block == priority_block
        )

        if (length(candidate_indices) > 0) {
          candidate_coords <- point_coords[candidate_indices, , drop = FALSE]
          distances_to_centroid <- sqrt((candidate_coords[,1] - fold_centroid[1])^2 +
                                          (candidate_coords[,2] - fold_centroid[2])^2)

          point_to_move <- candidate_indices[which.min(distances_to_centroid)]

          pts_sf_attempt$fold[point_to_move] <- smallest_fold
          pts_sf_attempt$block_type[point_to_move] <- "rebalanced"

          if (is.na(block_ownership[priority_block])) {
            block_ownership[priority_block] <- smallest_fold
          }

          moved <- TRUE
          break
        }
      }

      if (!moved) {
        available_shared <- shared_blocks[is.na(block_ownership[shared_blocks])]

        if (length(available_shared) > 0) {
          core_block <- spatial_fold_assignment$spatial_block[spatial_fold_assignment$fold == smallest_fold]
          block_distances <- distance_matrix[core_block, available_shared]
          closest_block <- available_shared[which.min(block_distances)]

          candidate_indices <- which(
            pts_sf_attempt$fold == largest_fold &
              pts_sf_attempt$spatial_block == closest_block
          )

          if (length(candidate_indices) > 0) {
            candidate_coords <- point_coords[candidate_indices, , drop = FALSE]
            distances_to_centroid <- sqrt((candidate_coords[,1] - fold_centroid[1])^2 +
                                            (candidate_coords[,2] - fold_centroid[2])^2)

            point_to_move <- candidate_indices[which.min(distances_to_centroid)]

            pts_sf_attempt$fold[point_to_move] <- smallest_fold
            pts_sf_attempt$block_type[point_to_move] <- "rebalanced"
            block_ownership[closest_block] <- smallest_fold
            fold_priority_blocks[[smallest_fold]] <- c(fold_priority_blocks[[smallest_fold]], closest_block)
            moved <- TRUE
          }
        }
      }

      if (!moved) {
        break
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
        n_spatial_exclusive_blocks = n_spatial_exclusive_blocks,
        n_shared_blocks = n_shared_blocks,
        spatial_fold_assignment = spatial_fold_assignment,
        best_attempt = attempt
      )
    }

    if (imbalance <= max_imbalance) {
      break
    }
  }

  if (best_imbalance > max_imbalance) {
    stop(paste0("Failed to achieve balance within ", max_internal_retries, " attempts. ",
                "Final imbalance: ", round(best_imbalance * 100, 2), "%. ",
                "Try increasing n_spatial or increasing max_imbalance threshold."))
  }

  pts_sf <- best_results$pts_sf
  voronoi_sf <- best_results$voronoi_sf
  final_fold_counts <- best_results$final_fold_counts
  mean_per_fold <- best_results$mean_per_fold
  imbalance <- best_results$imbalance
  n_spatially_exclusive_folds <- best_results$n_spatially_exclusive_folds
  n_temporally_exclusive_folds <- best_results$n_temporally_exclusive_folds
  n_spatial_exclusive_blocks <- best_results$n_spatial_exclusive_blocks
  n_shared_blocks <- best_results$n_shared_blocks

  ### REPORTING
  print("=== FOLD STRUCTURE ===")
  print(paste0(total_folds, " folds | ", n_spatial, " spatial blocks"))
  if (temporal_partitioning) {
    print(paste0(" | ", n_temporal, " temporal blocks"))
  } else {
    print(" | No temporal blocks")
  }
  print(paste0("Priority: ", toupper(blocking_priority)))
  print("=== FOLD SIZES ===")
  for (f in 1:total_folds) {
    points <- as.numeric(final_fold_counts[f])
    percent <- round(points / total_points * 100, 2)
    print(paste0("Fold ", f, ": ", points, " (", percent, "%)"))
  }
  print("")

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

    if (temporal_partitioning && n_temporal > 1) {
      fold_colors <- viridis::viridis(n = total_folds, option = "turbo")
      plot_list$temporal <- ggplot(pts_sf, aes(x = get(time_col), fill = factor(fold))) +
        geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
        scale_fill_manual(values = fold_colors, name = "Fold") +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold"), legend.position = "right") +
        labs(title = "Temporal Distribution by Fold", x = time_col, y = "Count")
    }

    fold_balance_df <- data.frame(
      fold = factor(names(final_fold_counts), levels = sort(as.numeric(names(final_fold_counts)))),
      count = as.numeric(final_fold_counts)
    )

    fold_colors <- viridis::viridis(n = total_folds, option = "turbo")
    plot_list$balance <- ggplot(fold_balance_df, aes(x = fold, y = count, fill = fold)) +
      geom_col(alpha = 0.8) +
      geom_hline(yintercept = mean_per_fold, linetype = "dashed", color = "red", linewidth = 0.8) +
      scale_fill_manual(values = fold_colors, name = "Fold") +
      theme_minimal(base_size = 10) +
      theme(plot.title = element_text(face = "bold")) +
      labs(title = "Fold Balance", x = "Fold", y = "Points")

    if (n_spatial > 1 && temporal_partitioning && n_temporal > 1) {
      fold_colors <- viridis::viridis(n = total_folds, option = "turbo")
      plot_list$combined <- ggplot() +
        geom_sf(data = reference_shapefile, fill = "gray98", color = "gray40", linewidth = 0.5) +
        geom_sf(data = voronoi_sf, fill = NA, color = "black", linewidth = 1) +
        geom_sf(data = pts_sf, aes(color = factor(fold), shape = factor(temporal_block)), size = 2, alpha = 0.8) +
        scale_color_manual(values = fold_colors, name = "Fold") +
        scale_shape_manual(values = c(16, 17, 15, 18)[1:n_temporal], name = "Temporal Block") +
        coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid = element_line(color = "gray90", linewidth = 0.3),
          panel.background = element_rect(fill = "aliceblue", color = NA)
        ) +
        labs(title = "Partitioning Structure by Fold")
    } else {
      plot_list$combined <- plot_list$balance
    }

    print(plot_list$combined)
  }

  summary_stats <- data.frame(
    parameter = c("spatial_blocks", "temporal_blocks", "total_folds", "blocking_priority",
                  "spatially_exclusive_folds", "temporally_exclusive_folds",
                  "dedicated_spatial_blocks", "shared_spatial_blocks",
                  "total_points", "points_removed", "pct_rows_removed",
                  "final_imbalance_pct",
                  "temporal_partitioning_enabled"),
    value = c(n_spatial, n_temporal, total_folds, blocking_priority,
              n_spatially_exclusive_folds, n_temporally_exclusive_folds,
              n_spatial_exclusive_blocks, n_shared_blocks,
              nrow(pts_sf), n_removed,
              ifelse(n_removed > 0, round(n_removed/n_original * 100, 2), 0),
              round(imbalance * 100, 2),
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
    print(paste0("Results saved to: ", output_file))
  }

  return(results)
}
