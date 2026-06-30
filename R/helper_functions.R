#' @keywords internal
#' @noRd
.natural_sort <- function(paths) {
  keys <- lapply(basename(paths), function(f) {
    parts <- strsplit(f, "[^0-9]+", perl = TRUE)[[1]]
    parts <- suppressWarnings(as.numeric(parts))
    parts[!is.na(parts)]
  })
  paths[order(vapply(keys,
                     function(k) paste(sprintf("%020.0f", k), collapse = "_"),
                     character(1)))]
}

#' @keywords internal
#' @noRd
#' @importFrom sf sf_use_s2 st_as_sf st_as_sfc st_bbox st_buffer st_cast
#' @importFrom sf st_collection_extract st_contains st_convex_hull st_coordinates
#' @importFrom sf st_crs st_distance st_drop_geometry st_geometry st_geometry_type
#' @importFrom sf st_intersection st_intersects st_is_empty st_join st_polygon
#' @importFrom sf st_read st_sample st_sf st_sfc st_union st_voronoi
#' @importFrom deldir deldir tile.list
#' @importFrom stats as.formula binomial chisq.test coef dist glm kmeans setNames
#' @importFrom stats logLik pchisq pnorm
#' @importFrom grDevices adjustcolor colorRampPalette hcl.colors hcl.pals recordPlot
#' @importFrom graphics abline barplot legend mtext par plot points
#' @importFrom tools file_ext toTitleCase
#' @importFrom utils read.csv
.test_cp_likelihood <- function(data, cp, alpha = 0.05, use_neighbor = TRUE) {
  n <- nrow(data)
  seg1_data <- data[1:cp, ]
  seg2_data <- data[(cp + 1):n, ]

  if (nrow(seg1_data) < 3 || nrow(seg2_data) < 3) {
    return(list(p_value = NA, significant = FALSE, test_statistic = NA))
  }

  tryCatch({
    if (use_neighbor) {
      model_full <- stats::glm(y ~ lag1 + neighbor, data = data, family = stats::binomial())
      model_seg1 <- stats::glm(y ~ lag1 + neighbor, data = seg1_data, family = stats::binomial())
      model_seg2 <- stats::glm(y ~ lag1 + neighbor, data = seg2_data, family = stats::binomial())
    } else {
      model_full <- stats::glm(y ~ lag1, data = data, family = stats::binomial())
      model_seg1 <- stats::glm(y ~ lag1, data = seg1_data, family = stats::binomial())
      model_seg2 <- stats::glm(y ~ lag1, data = seg2_data, family = stats::binomial())
    }

    loglik_full <- stats::logLik(model_full)[1]
    loglik_seg <- stats::logLik(model_seg1)[1] + stats::logLik(model_seg2)[1]
    lr_stat <- 2 * (loglik_seg - loglik_full)
    df <- length(stats::coef(model_seg1)) + length(stats::coef(model_seg2)) - length(stats::coef(model_full))
    p_value <- 1 - stats::pchisq(lr_stat, df = max(df, 1))

    list(p_value = p_value, significant = p_value < alpha, test_statistic = lr_stat, df = df)
  }, error = function(e) {
    list(p_value = NA, significant = FALSE, test_statistic = NA)
  })
}

#' @keywords internal
#' @noRd
.test_cp_permutation <- function(data, cp, n_perm = 1000, alpha = 0.05) {
  n <- nrow(data)
  seg1_mean <- mean(data$y[1:cp])
  seg2_mean <- mean(data$y[(cp + 1):n])
  obs_stat <- abs(seg1_mean - seg2_mean)

  perm_stats <- replicate(n_perm, {
    perm_y <- sample(data$y)
    abs(mean(perm_y[1:cp]) - mean(perm_y[(cp + 1):n]))
  })

  p_value <- mean(perm_stats >= obs_stat)
  list(p_value = p_value, significant = p_value < alpha,
       test_statistic = obs_stat, effect_size = obs_stat)
}

#' @keywords internal
#' @noRd
#' @importFrom stats pnorm
.test_cp_proportion <- function(data, cp, alpha = 0.05) {
  n <- nrow(data)
  seg1_y <- data$y[1:cp]
  seg2_y <- data$y[(cp + 1):n]

  n1 <- length(seg1_y)
  n2 <- length(seg2_y)
  p1 <- mean(seg1_y)
  p2 <- mean(seg2_y)
  p_pooled <- (sum(seg1_y) + sum(seg2_y)) / (n1 + n2)

  if (p_pooled == 0 || p_pooled == 1 || p1 == p2) {
    return(list(p_value = 1, significant = FALSE, test_statistic = 0, effect_size = abs(p1 - p2)))
  }

  se <- sqrt(p_pooled * (1 - p_pooled) * (1 / n1 + 1 / n2))
  if (se == 0) {
    return(list(p_value = 1, significant = FALSE, test_statistic = 0, effect_size = abs(p1 - p2)))
  }

  z_stat <- (p1 - p2) / se
  p_value <- 2 * (1 - stats::pnorm(abs(z_stat)))

  list(p_value = p_value, significant = p_value < alpha,
       test_statistic = z_stat, effect_size = abs(p1 - p2))
}

#' @keywords internal
#' @noRd
#' @importFrom stats chisq.test
.test_cp_chisquare <- function(data, cp, alpha = 0.05) {
  n <- nrow(data)
  seg1_y <- data$y[1:cp]
  seg2_y <- data$y[(cp + 1):n]

  cont_table <- matrix(c(
    sum(seg1_y == 1), sum(seg1_y == 0),
    sum(seg2_y == 1), sum(seg2_y == 0)
  ), nrow = 2, byrow = TRUE)

  test_result <- stats::chisq.test(cont_table, correct = FALSE)

  list(p_value = test_result$p.value, significant = test_result$p.value < alpha,
       test_statistic = test_result$statistic, effect_size = abs(mean(seg1_y) - mean(seg2_y)))
}

#' @keywords internal
#' @noRd
.assess_changepoint_significance <- function(data, cp_set, alpha = 0.05, n_perm = 1000, use_neighbor = TRUE) {
  if (length(cp_set) == 0) return(data.frame())

  cp_set <- sort(cp_set)
  significant_cps <- integer(0)
  results_list <- list()

  for (i in seq_along(cp_set)) {
    cp <- cp_set[i]

    baseline_start <- if (length(significant_cps) == 0) 1 else max(significant_cps) + 1
    current_end <- if (i < length(cp_set)) cp_set[i + 1] else nrow(data)

    baseline_segment <- baseline_start:cp
    current_segment <- (cp + 1):current_end

    if (length(baseline_segment) < 3 || length(current_segment) < 3) {
      results_list[[i]] <- data.frame(
        ChangePoint = cp,
        LR_PValue = NA, LR_Significant = FALSE,
        Perm_PValue = NA, Perm_Significant = FALSE, Perm_EffectSize = NA,
        Prop_PValue = NA, Prop_Significant = FALSE,
        Chi_PValue = NA, Chi_Significant = FALSE,
        Seg1_Proportion = NA, Seg2_Proportion = NA,
        Overall_Significant = FALSE
      )
      next
    }

    baseline_prop <- mean(data$y[baseline_segment])
    current_prop <- mean(data$y[current_segment])

    temp_combined <- data[c(baseline_segment, current_segment), ]
    cp_temp <- length(baseline_segment)

    lr_test <- suppressWarnings(.test_cp_likelihood(temp_combined, cp_temp, alpha, use_neighbor))
    lr_p_value <- ifelse(is.null(lr_test$p.value), lr_test$p_value, lr_test$p.value)
    lr_significant <- ifelse(is.na(lr_p_value), FALSE, lr_p_value < alpha)

    perm_test <- suppressWarnings(.test_cp_permutation(temp_combined, cp_temp, n_perm, alpha))
    prop_test <- suppressWarnings(.test_cp_proportion(temp_combined, cp_temp, alpha))

    chi_test <- suppressWarnings(tryCatch({
      .test_cp_chisquare(temp_combined, cp_temp, alpha)
    }, error = function(e) {
      list(p_value = NA, significant = FALSE, test_statistic = NA)
    }))
    chi_significant <- ifelse(is.na(chi_test$p_value), FALSE, chi_test$p_value < alpha)

    sig_votes <- sum(c(lr_significant, perm_test$significant, prop_test$significant, chi_significant), na.rm = TRUE)
    overall_significant <- sig_votes >= 2

    if (overall_significant) significant_cps <- c(significant_cps, cp)

    results_list[[i]] <- data.frame(
      ChangePoint = cp,
      LR_PValue = lr_p_value, LR_Significant = lr_significant,
      Perm_PValue = perm_test$p_value, Perm_Significant = perm_test$significant,
      Perm_EffectSize = perm_test$effect_size,
      Prop_PValue = prop_test$p_value, Prop_Significant = prop_test$significant,
      Chi_PValue = chi_test$p_value, Chi_Significant = chi_significant,
      Seg1_Proportion = baseline_prop, Seg2_Proportion = current_prop,
      Overall_Significant = overall_significant
    )
  }

  do.call(rbind, results_list)
}

#' @keywords internal
#' @noRd
.classify_pattern <- function(sig_results) {
  sig <- sig_results[sig_results$Overall_Significant == TRUE, ]
  if (nrow(sig) == 0) return("No Pattern")

  tolerance <- 1e-10
  prop_diffs <- sig$Seg2_Proportion - sig$Seg1_Proportion
  has_inc <- any(prop_diffs > tolerance)
  has_dec <- any(prop_diffs < -tolerance)

  if (has_inc & has_dec) return("Fluctuating/Intermittent")
  if (has_inc) return("Increasing")
  if (has_dec) return("Decreasing")
  return("Failed Classification")
}

#' @keywords internal
#' @noRd
.classify_pixel_with_times <- function(pixel_vals, n_middle, time_steps,
                                       fastcpd_params = list(), alpha = 0.05,
                                       n_perm = 1000, use_neighbor = TRUE) {

  y <- pixel_vals[1:n_middle]
  lag <- pixel_vals[(n_middle + 1):(2 * n_middle)]

  if (use_neighbor) {
    neighbor <- pixel_vals[(2 * n_middle + 1):(3 * n_middle)]
    mean_val <- pixel_vals[3 * n_middle + 1]
  } else {
    neighbor <- NULL
    mean_val <- pixel_vals[2 * n_middle + 1]
  }

  if (is.na(mean_val)) return(c(NA, NA, NA))
  if (mean_val < 0.01) return(c(1, NA, NA))
  if (mean_val > 0.99) return(c(2, NA, NA))
  if (any(is.na(c(y, lag)))) return(c(NA, NA, NA))
  if (use_neighbor && any(is.na(neighbor))) return(c(NA, NA, NA))

  data_matrix <- if (use_neighbor) cbind(y, lag, neighbor) else cbind(y, lag)

  cp_result <- tryCatch({
    suppressWarnings(do.call(
      fastcpd::fastcpd.binomial,
      c(list(data = data_matrix, r.progress = FALSE), fastcpd_params)
    ))
  }, error = function(e) NULL)

  if (is.null(cp_result)) return(c(7, NA, NA))

  cp_set <- cp_result@cp_set
  if (length(cp_set) == 0) return(c(3, NA, NA))

  data_df <- if (use_neighbor) {
    data.frame(y = y, lag1 = lag, neighbor = neighbor)
  } else {
    data.frame(y = y, lag1 = lag)
  }

  sig_results <- tryCatch({
    suppressWarnings(.assess_changepoint_significance(data_df, cp_set, alpha, n_perm, use_neighbor))
  }, error = function(e) NULL)

  if (is.null(sig_results) || nrow(sig_results) == 0) return(c(7, NA, NA))

  pattern <- .classify_pattern(sig_results)
  pattern_codes <- c("Always Absent" = 1, "Always Present" = 2, "No Pattern" = 3,
                     "Increasing" = 4, "Decreasing" = 5, "Fluctuating/Intermittent" = 6,
                     "Failed Classification" = 7)
  classification_code <- pattern_codes[pattern]

  time_decrease <- NA
  time_increase <- NA
  sig_cps <- sig_results[sig_results$Overall_Significant == TRUE, ]

  if (nrow(sig_cps) > 0) {
    if (pattern == "Decreasing") {
      dec_cps <- sig_cps[sig_cps$Seg2_Proportion < sig_cps$Seg1_Proportion, ]
      if (nrow(dec_cps) > 0) time_decrease <- time_steps[min(dec_cps$ChangePoint) + 1]
    }
    if (pattern == "Increasing") {
      inc_cps <- sig_cps[sig_cps$Seg2_Proportion > sig_cps$Seg1_Proportion, ]
      if (nrow(inc_cps) > 0) time_increase <- time_steps[min(inc_cps$ChangePoint) + 1]
    }
  }

  c(classification_code, time_decrease, time_increase)
}

#' @keywords internal
#' @noRd
.load_points_data <- function(points_file_path, xcol, ycol, points_crs,
                              verbose = TRUE) {

  if (is.character(points_file_path)) {
    if (!file.exists(points_file_path)) stop(paste("ERROR: File does not exist:", points_file_path))
    file_ext <- tolower(tools::file_ext(points_file_path))

    if (file_ext == "csv") {
      if (is.null(xcol)) stop("ERROR: 'xcol' is required when reading CSV files.")
      if (is.null(ycol)) stop("ERROR: 'ycol' is required when reading CSV files.")
      if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when reading CSV files.")
      if (verbose) message(paste("Reading CSV file:", basename(points_file_path)))
      pts <- utils::read.csv(points_file_path, stringsAsFactors = FALSE)
      if (!xcol %in% names(pts)) stop(paste0("ERROR: Column '", xcol, "' not found in CSV."))
      if (!ycol %in% names(pts)) stop(paste0("ERROR: Column '", ycol, "' not found in CSV."))

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      if (verbose) message(paste("Reading spatial file:", basename(points_file_path)))
      pts_sf_raw <- sf::st_read(points_file_path, quiet = TRUE)
      pts <- sf::st_drop_geometry(pts_sf_raw)
      coords <- sf::st_coordinates(pts_sf_raw)
      if (is.null(xcol)) xcol <- "X"
      if (is.null(ycol)) ycol <- "Y"
      pts[[xcol]] <- coords[, 1]
      pts[[ycol]] <- coords[, 2]
      if (is.null(points_crs)) points_crs <- sf::st_crs(pts_sf_raw)
    } else {
      stop(paste("ERROR: Unsupported file format:", file_ext))
    }

  } else if (inherits(points_file_path, "sf")) {
    if (verbose) message("Using provided sf object...")
    pts <- sf::st_drop_geometry(points_file_path)
    coords <- sf::st_coordinates(points_file_path)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]
    if (is.null(points_crs)) points_crs <- sf::st_crs(points_file_path)

  } else if (inherits(points_file_path, "sfc")) {
    if (verbose) message("Converting sfc to sf...")
    pts_sf_tmp <- sf::st_sf(geometry = points_file_path)
    pts <- data.frame(row.names = seq_along(points_file_path))
    coords <- sf::st_coordinates(pts_sf_tmp)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]
    if (is.null(points_crs)) points_crs <- sf::st_crs(pts_sf_tmp)

  } else if (inherits(points_file_path, "Spatial")) {
    if (verbose) message("Converting Spatial object to sf...")
    pts_sf_tmp <- sf::st_as_sf(points_file_path)
    pts <- sf::st_drop_geometry(pts_sf_tmp)
    coords <- sf::st_coordinates(pts_sf_tmp)
    if (is.null(xcol)) xcol <- "X"
    if (is.null(ycol)) ycol <- "Y"
    pts[[xcol]] <- coords[, 1]
    pts[[ycol]] <- coords[, 2]
    if (is.null(points_crs)) points_crs <- sf::st_crs(pts_sf_tmp)

  } else if (is.data.frame(points_file_path)) {
    if (is.null(xcol)) stop("ERROR: 'xcol' is required when providing a data frame.")
    if (is.null(ycol)) stop("ERROR: 'ycol' is required when providing a data frame.")
    if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when providing a data frame.")
    if (verbose) message("Using provided data frame...")
    pts <- points_file_path
    if (!xcol %in% names(pts)) stop(paste0("ERROR: Column '", xcol, "' not found."))
    if (!ycol %in% names(pts)) stop(paste0("ERROR: Column '", ycol, "' not found."))
  } else {
    stop("ERROR: points_file_path must be an sf object, sfc, Spatial, data frame, or file path.")
  }

  attr(pts, "xcol") <- xcol
  attr(pts, "ycol") <- ycol
  attr(pts, "crs") <- points_crs
  pts
}

#' @keywords internal
#' @noRd
#' @importFrom deldir deldir tile.list
.build_spatial_blocks <- function(pts_sf, reference_shapefile, n_spatial) {

  coords <- sf::st_coordinates(pts_sf)

  if (n_spatial == 1) {
    centers <- data.frame(spatial_block = 1,
                          center_x = mean(coords[, 1]),
                          center_y = mean(coords[, 2]))
    voronoi_sf <- sf::st_sf(fold_group = 1,
                            geometry = sf::st_geometry(sf::st_union(reference_shapefile)))
    return(list(clusters = rep(1, nrow(pts_sf)), voronoi_sf = voronoi_sf,
                adjacency = matrix(TRUE, 1, 1), centers = centers, dist_matrix = matrix(0, 1, 1)))
  }

  km <- stats::kmeans(coords, centers = n_spatial, nstart = 50, iter.max = 100)
  centers <- data.frame(spatial_block = seq_len(n_spatial),
                        center_x = km$centers[, 1],
                        center_y = km$centers[, 2])

  bbox <- sf::st_bbox(reference_shapefile)
  vor <- deldir::deldir(centers$center_x, centers$center_y,
                        rw = c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]))
  tiles <- deldir::tile.list(vor)
  polys <- lapply(seq_along(tiles), function(i) {
    tile <- tiles[[i]]
    sf::st_polygon(list(cbind(c(tile$x, tile$x[1]), c(tile$y, tile$y[1]))))
  })

  voronoi_sf <- sf::st_sf(fold_group = seq_len(n_spatial),
                          geometry = sf::st_sfc(polys, crs = sf::st_crs(pts_sf)))
  voronoi_sf <- suppressWarnings(sf::st_intersection(voronoi_sf, sf::st_union(reference_shapefile)))

  n_blocks <- nrow(voronoi_sf)
  adj <- matrix(FALSE, n_blocks, n_blocks)
  for (i in seq_len(n_blocks - 1)) {
    for (j in (i + 1):n_blocks) {
      isect <- suppressWarnings(sf::st_intersection(voronoi_sf$geometry[i], voronoi_sf$geometry[j]))
      if (length(isect) > 0 && !sf::st_is_empty(isect)) {
        gtype <- sf::st_geometry_type(isect)
        if (gtype %in% c("LINESTRING", "MULTILINESTRING", "GEOMETRYCOLLECTION")) {
          adj[i, j] <- TRUE
          adj[j, i] <- TRUE
        }
      }
    }
  }

  dm <- as.matrix(stats::dist(cbind(centers$center_x, centers$center_y)))
  list(clusters    = km$cluster,
       voronoi_sf  = voronoi_sf,
       adjacency   = adj,
       centers     = centers,
       dist_matrix = dm)
}

#' @keywords internal
#' @noRd
.select_dispersed_blocks <- function(distance_matrix, n_select, all_blocks) {
  if (n_select == 1) return(all_blocks[1])

  best_min_dist <- -Inf
  best_selection <- NULL
  for (first in all_blocks) {
    selected <- integer(n_select)
    selected[1] <- first
    for (k in 2:n_select) {
      remaining <- setdiff(all_blocks, selected[seq_len(k - 1)])
      dists <- apply(distance_matrix[remaining, selected[seq_len(k - 1)], drop = FALSE], 1, min)
      selected[k] <- remaining[which.max(dists)]
    }
    min_pw <- min(distance_matrix[selected, selected][upper.tri(distance_matrix[selected, selected])])
    if (min_pw > best_min_dist) {
      best_min_dist  <- min_pw
      best_selection <- selected
    }
  }
  best_selection
}

#' @keywords internal
#' @noRd
.assign_balanced_folds <- function(pts_work, n_folds, n_spatial, spatial_centers,
                                   distance_matrix, adjacency_matrix, max_imbalance, total_folds) {

  pts_work$fold       <- NA_integer_
  pts_work$block_type <- NA_character_
  all_blocks          <- seq_len(n_spatial)

  core_blocks <- .select_dispersed_blocks(distance_matrix, n_folds, all_blocks)
  assigned_blocks <- core_blocks
  fold_block_pairs <- vector("list", n_folds)

  for (fold_id in seq_len(n_folds)) {
    core <- core_blocks[fold_id]
    adjacent <- which(adjacency_matrix[core, ])
    available_adj <- setdiff(adjacent, assigned_blocks)

    second <- if (length(available_adj) > 0) {
      available_adj[1]
    } else {
      available_any <- setdiff(all_blocks, assigned_blocks)
      if (length(available_any) > 0) available_any[which.min(distance_matrix[core, available_any])] else core
    }
    fold_block_pairs[[fold_id]] <- c(core, second)
    assigned_blocks             <- c(assigned_blocks, second)

    for (block_id in fold_block_pairs[[fold_id]]) {
      for (temp_id in seq_len(2)) {
        idx                      <- which(pts_work$spatial_block == block_id & pts_work$temporal_block == temp_id)
        pts_work$fold[idx]       <- fold_id
        pts_work$block_type[idx] <- "balanced_core"
      }
    }
  }

  unassigned <- setdiff(all_blocks, assigned_blocks)
  for (block_id in unassigned) {
    bc <- c(spatial_centers$center_x[block_id], spatial_centers$center_y[block_id])
    fold_dists <- vapply(seq_len(n_folds), function(f) {
      fb <- fold_block_pairs[[f]]
      fc <- c(mean(spatial_centers$center_x[fb]), mean(spatial_centers$center_y[fb]))
      sqrt(sum((bc - fc)^2))
    }, numeric(1))
    sorted <- order(fold_dists)

    idx_t1 <- which(pts_work$spatial_block == block_id & pts_work$temporal_block == 1)
    idx_t2 <- which(pts_work$spatial_block == block_id & pts_work$temporal_block == 2)
    if (length(idx_t1) > 0) {
      pts_work$fold[idx_t1]       <- sorted[1]
      pts_work$block_type[idx_t1] <- "balanced_shared"
    }
    if (length(idx_t2) > 0) {
      pts_work$fold[idx_t2]       <- sorted[min(2, length(sorted))]
      pts_work$block_type[idx_t2] <- "balanced_shared"
    }
  }

  .rebalance_folds(pts_work, total_folds, max_imbalance, max_iterations = 25, moveable_type = "balanced_shared")
}

#' @keywords internal
#' @noRd
.kd_split_st <- function(idx, n_folds, target_sizes, coords) {
  if (n_folds == 1) return(list(idx))
  if (length(idx) <= 1) return(list(idx))
  sub_coords  <- coords[idx, , drop = FALSE]
  x_range     <- diff(range(sub_coords[, 1]))
  y_range     <- diff(range(sub_coords[, 2]))
  axis        <- if (x_range >= y_range) 1 else 2
  order_along <- order(sub_coords[, axis])
  idx_sorted  <- idx[order_along]
  left_n      <- floor(n_folds / 2)
  right_n     <- n_folds - left_n
  split_at    <- sum(target_sizes[seq_len(left_n)])
  split_at    <- max(1, min(split_at, length(idx_sorted) - 1))
  left_idx    <- idx_sorted[seq_len(split_at)]
  right_idx   <- idx_sorted[(split_at + 1):length(idx_sorted)]
  c(.kd_split_st(left_idx,  left_n,  target_sizes[seq_len(left_n)], coords),
    .kd_split_st(right_idx, right_n, target_sizes[(left_n + 1):n_folds], coords))
}

#' @keywords internal
#' @noRd
.assign_spatiotemporal_folds <- function(pts_work, n_spatial_folds, n_temporal_folds,
                                         total_folds, n_spatial, n_temporal,
                                         spatial_centers, distance_matrix,
                                         adjacency_matrix, temporal_partitioning,
                                         max_imbalance, use_balanced) {

  pts_work$fold       <- NA_integer_
  pts_work$block_type <- NA_character_
  coords              <- sf::st_coordinates(pts_work)
  n_pts  <- nrow(coords)

  if (n_spatial_folds > 0 && n_temporal_folds == 0) {

    target_size  <- floor(n_pts / n_spatial_folds)
    n_larger     <- n_pts - target_size * n_spatial_folds
    target_sizes <- rep(target_size, n_spatial_folds)
    if (n_larger > 0) target_sizes[seq_len(n_larger)] <- target_size + 1

    groups   <- .kd_split_st(seq_len(n_pts), n_spatial_folds, target_sizes, coords)
    fold_vec <- integer(n_pts)
    for (f in seq_along(groups)) fold_vec[groups[[f]]] <- f

    group_centroids <- do.call(rbind, lapply(seq_along(groups), function(f) {
      colMeans(coords[groups[[f]], , drop = FALSE])
    }))
    dist_to_centroids <- matrix(NA_real_, nrow = n_pts, ncol = n_spatial_folds)
    for (f in seq_len(n_spatial_folds)) {
      dx <- coords[, 1] - group_centroids[f, 1]
      dy <- coords[, 2] - group_centroids[f, 2]
      dist_to_centroids[, f] <- sqrt(dx^2 + dy^2)
    }
    fold_vec <- apply(dist_to_centroids, 1, which.min)

    for (iter in seq_len(n_pts * 2)) {
      counts <- tabulate(fold_vec, nbins = n_spatial_folds)
      over   <- which(counts > target_sizes)
      under  <- which(counts < target_sizes)
      if (length(over) == 0 || length(under) == 0) break
      of <- over[1]; uf <- under[1]
      candidates <- which(fold_vec == of)
      best <- candidates[which.min(dist_to_centroids[candidates, uf])]
      fold_vec[best] <- uf
    }

    pts_work$fold       <- fold_vec
    pts_work$block_type <- "spatial_exclusive"
    return(pts_work)
  }

  temporal_folds <- as.integer((n_spatial_folds + 1):total_folds)
  n_groups       <- total_folds

  fold_size  <- floor(n_pts / n_groups)
  remainder  <- n_pts - fold_size * n_groups
  group_target_sizes <- rep(fold_size, n_groups)
  if (remainder > 0) group_target_sizes[seq_len(remainder)] <- fold_size + 1

  groups    <- .kd_split_st(seq_len(n_pts), n_groups, group_target_sizes, coords)
  group_vec <- integer(n_pts)
  for (g in seq_along(groups)) group_vec[groups[[g]]] <- g

  group_centroids <- do.call(rbind, lapply(seq_len(n_groups), function(g) {
    idx <- which(group_vec == g)
    if (length(idx) == 0) return(c(NA_real_, NA_real_))
    colMeans(coords[idx, , drop = FALSE])
  }))

  dist_to_groups <- matrix(NA_real_, nrow = n_pts, ncol = n_groups)
  for (g in seq_len(n_groups)) {
    if (any(is.na(group_centroids[g, ]))) next
    dx <- coords[, 1] - group_centroids[g, 1]
    dy <- coords[, 2] - group_centroids[g, 2]
    dist_to_groups[, g] <- sqrt(dx^2 + dy^2)
  }

  if (!any(is.na(dist_to_groups))) {
    group_vec <- apply(dist_to_groups, 1, which.min)
  }

  for (iter in seq_len(n_pts * 2)) {
    counts <- tabulate(group_vec, nbins = n_groups)
    over   <- which(counts > group_target_sizes)
    under  <- which(counts < group_target_sizes)
    if (length(over) == 0 || length(under) == 0) break
    of <- over[1]; uf <- under[1]
    candidates <- which(group_vec == of)
    if (length(candidates) == 0) break
    best <- candidates[which.min(dist_to_groups[candidates, uf])]
    group_vec[best] <- uf
  }
  centroid_dists <- as.matrix(stats::dist(group_centroids))

  if (temporal_partitioning && "temporal_block" %in% names(pts_work)) {

    group_tb_score <- vapply(seq_len(n_groups), function(g) {
      idx <- which(group_vec == g)
      tb  <- pts_work$temporal_block[idx]
      tb  <- tb[!is.na(tb)]
      if (length(tb) == 0) return(Inf)
      props <- as.numeric(table(factor(tb, levels = seq_len(n_temporal)))) / length(tb)
      sum((props - 1 / n_temporal) ^ 2)
    }, numeric(1))

    if (n_temporal_folds == 1) {
      pool_groups <- which.min(group_tb_score)
    } else {
      seed        <- which.min(group_tb_score)
      pool_groups <- seed
      remaining   <- setdiff(seq_len(n_groups), seed)
      for (step in seq_len(n_temporal_folds - 1)) {
        if (length(remaining) == 0) break
        if (length(pool_groups) == 1) {
          dists_to_pool <- centroid_dists[remaining, pool_groups]
        } else {
          dists_to_pool <- apply(centroid_dists[remaining, pool_groups, drop = FALSE], 1, min)
        }
        next_group  <- remaining[which.min(dists_to_pool)]
        pool_groups <- c(pool_groups, next_group)
        remaining   <- setdiff(remaining, next_group)
      }
    }

  } else {
    all_centroid <- colMeans(group_centroids, na.rm = TRUE)
    dists_to_center <- sqrt((group_centroids[, 1] - all_centroid[1])^2 +
                              (group_centroids[, 2] - all_centroid[2])^2)
    seed        <- which.min(dists_to_center)
    pool_groups <- seed
    remaining   <- setdiff(seq_len(n_groups), seed)
    for (step in seq_len(n_temporal_folds - 1)) {
      if (length(remaining) == 0) break
      if (length(pool_groups) == 1) {
        dists_to_pool <- centroid_dists[remaining, pool_groups]
      } else {
        dists_to_pool <- apply(centroid_dists[remaining, pool_groups, drop = FALSE], 1, min)
      }
      next_group  <- remaining[which.min(dists_to_pool)]
      pool_groups <- c(pool_groups, next_group)
      remaining   <- setdiff(remaining, next_group)
    }
  }

  spatial_groups <- setdiff(seq_len(n_groups), pool_groups)
  for (f in seq_len(n_spatial_folds)) {
    g   <- spatial_groups[f]
    idx <- which(group_vec == g)
    pts_work$fold[idx]       <- f
    pts_work$block_type[idx] <- "spatial_exclusive"
  }

  temporal_pool_idx <- which(group_vec %in% pool_groups)

  if (length(temporal_pool_idx) > 0) {

    if (temporal_partitioning &&
        "temporal_block" %in% names(pts_work) &&
        any(!is.na(pts_work$temporal_block[temporal_pool_idx]))) {

      t_blocks  <- pts_work$temporal_block[temporal_pool_idx]
      unique_tb <- sort(unique(t_blocks[!is.na(t_blocks)]))
      n_unique  <- length(unique_tb)

      tb_fold_map <- stats::setNames(
        temporal_folds[((seq_len(n_unique) - 1) %% n_temporal_folds) + 1],
        as.character(unique_tb)
      )

      for (k in seq_along(temporal_pool_idx)) {
        gi <- temporal_pool_idx[k]
        tb <- pts_work$temporal_block[gi]
        if (!is.na(tb)) {
          pts_work$fold[gi]       <- tb_fold_map[as.character(tb)]
          pts_work$block_type[gi] <- "temporal_exclusive"
        } else {
          t_counts <- tabulate(
            pts_work$fold[pts_work$fold %in% temporal_folds & !is.na(pts_work$fold)],
            nbins = max(temporal_folds)
          )[temporal_folds]
          pts_work$fold[gi]       <- temporal_folds[which.min(t_counts)]
          pts_work$block_type[gi] <- "temporal_exclusive"
        }
      }

    } else {
      pool_coords <- coords[temporal_pool_idx, , drop = FALSE]
      x_rng       <- diff(range(pool_coords[, 1]))
      y_rng       <- diff(range(pool_coords[, 2]))
      sort_axis   <- if (x_rng >= y_rng) 1 else 2
      pool_order  <- order(pool_coords[, sort_axis])
      n_pool      <- length(temporal_pool_idx)

      chunk_size      <- floor(n_pool / n_temporal_folds)
      n_larger_chunks <- n_pool - chunk_size * n_temporal_folds
      chunk_sizes     <- rep(chunk_size, n_temporal_folds)
      if (n_larger_chunks > 0)
        chunk_sizes[seq_len(n_larger_chunks)] <- chunk_size + 1

      pos <- 1
      for (ti in seq_len(n_temporal_folds)) {
        end     <- pos + chunk_sizes[ti] - 1
        sel_idx <- temporal_pool_idx[pool_order[pos:end]]
        pts_work$fold[sel_idx]       <- temporal_folds[ti]
        pts_work$block_type[sel_idx] <- "temporal_exclusive"
        pos <- end + 1
      }
    }
  }

  still_na <- which(is.na(pts_work$fold))
  if (length(still_na) > 0) {
    counts_all <- tabulate(pts_work$fold[!is.na(pts_work$fold)], nbins = total_folds)
    for (idx in still_na) {
      pts_work$fold[idx]       <- as.integer(which.min(counts_all))
      pts_work$block_type[idx] <- "remainder"
      counts_all[pts_work$fold[idx]] <- counts_all[pts_work$fold[idx]] + 1
    }
  }

  pts_work
}

#' @keywords internal
#' @noRd
.rebalance_folds <- function(pts_work, total_folds, max_imbalance,
                             max_iterations = 25, moveable_type = NULL) {

  point_coords <- sf::st_coordinates(pts_work)

  for (iter in seq_len(max_iterations)) {
    counts <- table(factor(pts_work$fold, levels = seq_len(total_folds)))
    mean_pf <- mean(counts)
    imb <- max(abs(counts - mean_pf)) / mean_pf
    if (imb <= max_imbalance) break

    largest <- as.integer(names(which.max(counts)))
    smallest <- as.integer(names(which.min(counts)))
    if (counts[largest] <= counts[smallest]) break

    candidates <- if (!is.null(moveable_type)) {
      which(pts_work$fold == largest & pts_work$block_type == moveable_type)
    } else {
      which(pts_work$fold == largest)
    }
    if (length(candidates) == 0) break

    smallest_pts <- which(pts_work$fold == smallest)
    if (length(smallest_pts) == 0) break
    centroid <- c(mean(point_coords[smallest_pts, 1]), mean(point_coords[smallest_pts, 2]))
    cc <- point_coords[candidates, , drop = FALSE]
    dists <- sqrt((cc[, 1] - centroid[1])^2 + (cc[, 2] - centroid[2])^2)
    pts_work$fold[candidates[which.min(dists)]] <- smallest
  }
  pts_work
}

#' @keywords internal
#' @noRd
.resolve_palette <- function(n, palette) {
  if (palette %in% grDevices::hcl.pals()) {
    return(grDevices::hcl.colors(n, palette = palette))
  }

  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    brewer_pals <- rownames(RColorBrewer::brewer.pal.info)
    if (palette %in% brewer_pals) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      pal_colors <- RColorBrewer::brewer.pal(min(n, max_n), palette)
      if (n > max_n) {
        pal_colors <- grDevices::colorRampPalette(pal_colors)(n)
      }
      return(pal_colors)
    }
  }

  stop(paste0(
    "ERROR: Palette '", palette, "' not recognized. ",
    "Use grDevices::hcl.pals() to see available HCL palettes, or install ",
    "RColorBrewer for additional palette options."
  ))
}

#' @keywords internal
#' @noRd
.build_fold_polygons <- function(pts_sf, reference_shapefile, total_folds, partition_mode) {

  if (is.null(reference_shapefile) || partition_mode == "random") return(NULL)

  pts_coords_map <- sf::st_coordinates(pts_sf)
  fold_vec_pts   <- pts_sf$fold
  ref_sf_union   <- sf::st_union(reference_shapefile)
  bbox           <- sf::st_bbox(reference_shapefile)

  is_temporal_mode <- partition_mode %in% c("spatiotemporal", "balanced", "temporal_only")

  if (is_temporal_mode && "block_type" %in% names(pts_sf)) {
    pts_df_bt <- sf::st_drop_geometry(pts_sf)
    spatial_folds <- sort(unique(pts_df_bt$fold[
      pts_df_bt$block_type %in% c("spatial_exclusive", "rebalanced")
    ]))
    temporal_folds <- sort(unique(pts_df_bt$fold[
      pts_df_bt$block_type == "temporal_exclusive"
    ]))
    fold_vec_render <- fold_vec_pts
    fold_vec_render[fold_vec_pts %in% temporal_folds] <- total_folds + 1
    render_folds <- as.integer(c(spatial_folds, total_folds + 1))
  } else {
    spatial_folds   <- as.integer(seq_len(total_folds))
    temporal_folds  <- integer(0)
    fold_vec_render <- fold_vec_pts
    render_folds    <- as.integer(seq_len(total_folds))
  }

  max_render_idx <- max(render_folds)
  fold_mcps      <- vector("list", max_render_idx)
  valid_render   <- integer(0)

  for (f in render_folds) {
    fold_pts <- pts_coords_map[fold_vec_render == f, , drop = FALSE]
    if (nrow(fold_pts) == 0) next
    fold_sfc <- sf::st_as_sf(as.data.frame(fold_pts), coords = c("X", "Y"),
                             crs = sf::st_crs(pts_sf))
    fold_mcps[[f]] <- if (nrow(fold_pts) >= 3) {
      sf::st_convex_hull(sf::st_union(fold_sfc))
    } else {
      sf::st_buffer(sf::st_union(fold_sfc), dist = 1e-6)
    }
    valid_render <- c(valid_render, f)
  }
  valid_render <- as.integer(valid_render)

  fold_real_pts <- lapply(seq_len(max_render_idx), function(f) {
    if (!f %in% valid_render) return(NULL)
    fold_pts <- pts_coords_map[fold_vec_render == f, , drop = FALSE]
    if (nrow(fold_pts) == 0) return(NULL)
    sf::st_as_sf(as.data.frame(fold_pts), coords = c("X", "Y"),
                 crs = sf::st_crs(pts_sf))
  })

  set.seed(42)
  sample_pts <- sf::st_sample(ref_sf_union, size = 1000, type = "random")
  sample_pts <- sf::st_as_sf(sample_pts)
  sample_pts <- sample_pts[!sf::st_is_empty(sample_pts), ]
  sample_xy  <- sf::st_coordinates(sample_pts)
  n_sample   <- nrow(sample_xy)

  mcp_geoms      <- do.call(c, lapply(valid_render, function(f) fold_mcps[[f]]))
  contain_sparse <- suppressMessages(suppressWarnings(
    sf::st_intersects(sample_pts, mcp_geoms, sparse = TRUE)
  ))

  n_hits      <- lengths(contain_sparse)
  idx_one     <- which(n_hits == 1)
  idx_multi   <- which(n_hits > 1)
  idx_outside <- which(n_hits == 0)
  sample_fold <- integer(n_sample)

  if (length(idx_one) > 0) {
    sample_fold[idx_one] <- valid_render[vapply(contain_sparse[idx_one],
                                                `[[`, integer(1), 1)]
  }

  if (length(idx_multi) > 0) {
    sample_fold[idx_multi] <- vapply(idx_multi, function(i) {
      hit_folds <- valid_render[contain_sparse[[i]]]
      best_fold <- hit_folds[1]
      best_dist <- Inf
      pt_i      <- sample_pts[i, ]
      for (hf in hit_folds) {
        if (is.null(fold_real_pts[[hf]])) next
        d <- min(as.numeric(sf::st_distance(pt_i, fold_real_pts[[hf]])))
        if (d < best_dist) {
          best_dist <- d
          best_fold <- hf
        }
      }
      as.integer(best_fold)
    }, integer(1))
  }

  if (length(idx_outside) > 0) {
    dist_mat <- suppressMessages(suppressWarnings(
      sf::st_distance(sample_pts[idx_outside, ], mcp_geoms)
    ))
    nearest_col              <- apply(dist_mat, 1, which.min)
    sample_fold[idx_outside] <- valid_render[nearest_col]
  }

  sample_assigned <- sf::st_as_sf(
    data.frame(fold = sample_fold, sample_xy),
    coords = c("X", "Y"), crs = sf::st_crs(pts_sf)
  )

  s2_state <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  vor_raw <- tryCatch({
    suppressMessages(suppressWarnings(
      sf::st_voronoi(sf::st_union(sample_assigned),
                     envelope = sf::st_as_sfc(bbox))
    ))
  }, error = function(e) NULL)
  suppressMessages(sf::sf_use_s2(s2_state))

  if (is.null(vor_raw)) return(NULL)

  vor_polys <- sf::st_collection_extract(vor_raw, "POLYGON")
  vor_sf    <- sf::st_sf(geometry = vor_polys, crs = sf::st_crs(pts_sf))

  joined <- suppressMessages(suppressWarnings(
    sf::st_join(vor_sf, sample_assigned, join = sf::st_contains)
  ))

  real_pts_render <- sf::st_as_sf(
    data.frame(fold = fold_vec_render, pts_coords_map),
    coords = c("X", "Y"), crs = sf::st_crs(pts_sf)
  )

  sentinel     <- total_folds + 1
  orphan_geoms <- NULL
  keep_list    <- list()

  unique_render_folds <- unique(joined$fold[!is.na(joined$fold)])

  for (f in unique_render_folds) {
    sub <- joined[!is.na(joined$fold) & joined$fold == f, ]
    if (nrow(sub) == 0) next

    merged  <- suppressWarnings(sf::st_union(sub))
    clipped <- suppressWarnings(sf::st_intersection(merged, ref_sf_union))
    clipped <- clipped[!sf::st_is_empty(clipped)]
    if (length(clipped) == 0) next

    parts_sf <- suppressWarnings(
      sf::st_cast(sf::st_sf(fold = f, geometry = clipped), "POLYGON")
    )
    parts_sf <- parts_sf[!sf::st_is_empty(parts_sf), ]
    if (nrow(parts_sf) == 0) next

    if (f %in% spatial_folds) {
      fold_real <- real_pts_render[real_pts_render$fold == f, ]

      if (nrow(fold_real) == 0) {
        orphan_geoms <- c(orphan_geoms, sf::st_geometry(parts_sf))
        next
      }

      bbox_diag <- sqrt((bbox["xmax"] - bbox["xmin"])^2 +
                          (bbox["ymax"] - bbox["ymin"])^2)
      proximity_threshold <- bbox_diag * 0.001

      part_min_dists <- vapply(seq_len(nrow(parts_sf)), function(p) {
        suppressMessages(suppressWarnings(
          min(as.numeric(sf::st_distance(fold_real, parts_sf[p, ])))
        ))
      }, numeric(1))

      keep_mask <- part_min_dists <= proximity_threshold

      if (any(keep_mask)) {
        kept_geom <- suppressWarnings(sf::st_union(parts_sf[keep_mask, ]))
        keep_list[[length(keep_list) + 1]] <- sf::st_sf(
          fold     = f,
          geometry = sf::st_sfc(kept_geom, crs = sf::st_crs(pts_sf))
        )
      }
      if (any(!keep_mask)) {
        orphan_geoms <- c(orphan_geoms, sf::st_geometry(parts_sf[!keep_mask, ]))
      }

    } else {
      keep_list[[length(keep_list) + 1]] <- sf::st_sf(
        fold     = f,
        geometry = sf::st_sfc(suppressWarnings(sf::st_union(parts_sf)),
                              crs = sf::st_crs(pts_sf))
      )
    }
  }

  dissolved <- if (length(keep_list) > 0) do.call(rbind, keep_list) else NULL

  if (!is.null(orphan_geoms) && length(orphan_geoms) > 0 &&
      length(temporal_folds) > 0 && !is.null(dissolved)) {
    orphan_sfc <- sf::st_sfc(orphan_geoms, crs = sf::st_crs(pts_sf))
    existing_t <- dissolved[!is.na(dissolved$fold) & dissolved$fold == sentinel, ]
    if (nrow(existing_t) > 0) {
      merged_t <- suppressWarnings(
        sf::st_union(c(sf::st_geometry(existing_t), orphan_sfc))
      )
      dissolved <- dissolved[is.na(dissolved$fold) | dissolved$fold != sentinel, ]
    } else {
      merged_t <- suppressWarnings(sf::st_union(orphan_sfc))
    }
    dissolved <- rbind(
      dissolved,
      sf::st_sf(fold     = sentinel,
                geometry = sf::st_sfc(merged_t, crs = sf::st_crs(pts_sf)))
    )
  }

  dissolved
}

#' @keywords internal
#' @noRd
.plot_partitions_base <- function(pts_sf, reference_shapefile, final_fold_counts,
                                  mean_per_fold, total_folds, partition_mode,
                                  time_cols, temporal_partitioning, n_temporal,
                                  voronoi_sf, plot_palette = "Dark 2") {

  fold_colors <- .resolve_palette(total_folds, plot_palette)
  plot_list <- list()

  graphics::barplot(
    as.numeric(final_fold_counts), names.arg = names(final_fold_counts),
    col = fold_colors, main = "Fold Balance", xlab = "Fold", ylab = "Points", border = "gray30"
  )
  graphics::abline(h = mean_per_fold, lty = 2, col = "red", lwd = 2)
  plot_list$balance <- grDevices::recordPlot()

  if (!is.null(time_cols) && time_cols %in% names(pts_sf)) {
    time_vals <- sf::st_drop_geometry(pts_sf)[[time_cols]]
    fold_vals <- pts_sf$fold
    time_range <- range(time_vals, na.rm = TRUE)
    breaks_seq <- seq(time_range[1], time_range[2], length.out = 31)
    bin_labels <- round(breaks_seq[-length(breaks_seq)])

    bin_ids <- findInterval(time_vals, breaks_seq, rightmost.closed = TRUE)
    bin_ids[bin_ids == 0] <- 1
    n_bins <- length(breaks_seq) - 1

    count_matrix <- matrix(0, nrow = total_folds, ncol = n_bins)
    for (f in seq_len(total_folds)) {
      fold_bins <- bin_ids[fold_vals == f]
      tab <- table(factor(fold_bins, levels = seq_len(n_bins)))
      count_matrix[f, ] <- as.numeric(tab)
    }
    rownames(count_matrix) <- paste("Fold", seq_len(total_folds))

    label_interval <- max(1, round(n_bins / 15))
    display_labels <- rep("", n_bins)
    display_labels[seq(1, n_bins, by = label_interval)] <- bin_labels[seq(1, n_bins, by = label_interval)]

    graphics::barplot(
      count_matrix, beside = FALSE, col = fold_colors,
      names.arg = display_labels, las = 2, border = NA,
      main = "Temporal Distribution by Fold", xlab = time_cols, ylab = "Count"
    )
    graphics::legend("topright", legend = paste("Fold", seq_len(total_folds)),
                     fill = fold_colors, border = NA, cex = 0.8, bty = "n")
    plot_list$temporal <- grDevices::recordPlot()
  }

  if (!is.null(reference_shapefile)) {
    map_title <- paste(gsub("_", " ", tools::toTitleCase(partition_mode)), "Partitioning")
    graphics::plot(sf::st_geometry(reference_shapefile), col = "gray98", border = "gray40",
                   main = map_title)

    if (partition_mode != "random") {

      dissolved <- .build_fold_polygons(pts_sf, reference_shapefile, total_folds, partition_mode)

      is_temporal_mode <- partition_mode %in% c("spatiotemporal", "balanced", "temporal_only")

      if (is_temporal_mode && "block_type" %in% names(pts_sf)) {
        pts_df_bt      <- sf::st_drop_geometry(pts_sf)
        spatial_folds  <- sort(unique(pts_df_bt$fold[
          pts_df_bt$block_type %in% c("spatial_exclusive", "rebalanced")
        ]))
        temporal_folds <- sort(unique(pts_df_bt$fold[
          pts_df_bt$block_type == "temporal_exclusive"
        ]))
      } else {
        spatial_folds  <- as.integer(seq_len(total_folds))
        temporal_folds <- integer(0)
      }
      sentinel <- total_folds + 1

      if (!is.null(dissolved)) {
        for (f in spatial_folds) {
          sub <- dissolved[!is.na(dissolved$fold) & dissolved$fold == f, ]
          if (nrow(sub) == 0) next
          graphics::plot(sf::st_geometry(sub), add = TRUE,
                         col = grDevices::adjustcolor(fold_colors[f], alpha.f = 0.25),
                         border = fold_colors[f], lwd = 0.8)
        }

        if (length(temporal_folds) > 0) {
          sub_t <- dissolved[!is.na(dissolved$fold) & dissolved$fold == sentinel, ]
          if (nrow(sub_t) > 0) {
            graphics::plot(sf::st_geometry(sub_t), add = TRUE,
                           col = grDevices::adjustcolor("gray60", alpha.f = 0.20),
                           border = "gray50", lwd = 0.8, lty = 2)
          }
        }
      }

      graphics::plot(sf::st_geometry(reference_shapefile), add = TRUE,
                     col = NA, border = "gray40")
    }

    pts_coords <- sf::st_coordinates(pts_sf)

    use_temporal_shapes <- n_temporal > 1 && !any(is.na(pts_sf$temporal_block))
    pch_options <- c(16, 17, 15, 18)
    point_pch <- if (use_temporal_shapes) {
      pch_options[pmin(pts_sf$temporal_block, length(pch_options))]
    } else {
      16
    }

    graphics::points(pts_coords[, 1], pts_coords[, 2], col = fold_colors[pts_sf$fold],
                     pch = point_pch, cex = 0.8)

    legend_labels <- paste("Fold", seq_len(total_folds))
    graphics::legend("topright", legend = legend_labels,
                     col = fold_colors, pch = 16, cex = 0.8, bty = "n")
    plot_list$combined <- grDevices::recordPlot()
  }

  list(plots = plot_list, fold_polygons = if (exists("dissolved")) dissolved else NULL)
}

#' @keywords internal
#' @noRd
.build_partition_result <- function(pts_sf, voronoi_sf,
                                    final_fold_counts, mean_per_fold, imbalance,
                                    n_spatial, n_temporal, total_folds, partition_mode,
                                    temporal_partitioning, use_balanced,
                                    use_random, n_spatial_folds, n_temporal_folds,
                                    n_balanced_folds, n_random_folds,
                                    n_removed, n_original, plot_list, output_file,
                                    reference_shapefile = NULL,
                                    verbose = TRUE) {

  total_points <- nrow(pts_sf)

  if (verbose) message(paste0("--- Fold Structure (",
                              paste(gsub("_", " ", tools::toTitleCase(partition_mode))), ") ---"))
  if (verbose) message(paste0("  ", total_folds, " fold",
                              if (total_folds != 1) "s" else ""))
  for (f in seq_len(total_folds)) {
    n_pts <- as.numeric(final_fold_counts[f])
    if (verbose) message(paste0("  Fold ", f, ": ", n_pts, " (", round(n_pts / total_points * 100, 2), "%)"))
  }

  summary_stats <- data.frame(
    parameter = c("total_folds", "n_spatial_folds", "n_temporal_folds",
                  "n_balanced_folds", "n_random_folds", "partition_mode",
                  "total_points", "points_removed", "pct_rows_removed",
                  "final_imbalance_pct", "temporal_partitioning_enabled"),
    value = c(total_folds,
              n_spatial_folds, n_temporal_folds, n_balanced_folds,
              n_random_folds, partition_mode,
              total_points, n_removed,
              ifelse(n_removed > 0, round(n_removed / n_original * 100, 2), 0),
              round(imbalance * 100, 2),
              as.character(temporal_partitioning)),
    stringsAsFactors = FALSE
  )

  internal_cols <- c("spatial_block", "temporal_block", "block_type")
  pts_sf_public <- pts_sf[, !names(pts_sf) %in% internal_cols, drop = FALSE]

  fold_polygons <- .build_fold_polygons(pts_sf, reference_shapefile, total_folds, partition_mode)

  results <- list(
    folds         = sf::st_drop_geometry(pts_sf)[, "fold", drop = FALSE],
    points_sf     = pts_sf_public,
    voronoi_folds = if (!is.null(fold_polygons)) fold_polygons else voronoi_sf,
    summary       = summary_stats,
    plots         = plot_list
  )

  if (!is.null(output_file)) {
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(results, output_file)
    if (verbose) message(paste0("Results saved to: ", output_file))
  }

  invisible(results)
}

#' @keywords internal
#' @noRd
.plot_trend_timeseries <- function(timestep_summary, all_units, unit_colors,
                                   ts_x_rng, ts_y_max, x_ticks) {
  right_mar <- max(2, ceiling(max(nchar(all_units)) / 2.5) + 1)
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(4, 7.5, 3.5, right_mar),
                mgp = c(1.8, 0.7, 0))
  graphics::plot(NULL,
                 xlim = ts_x_rng,
                 ylim = c(0, ts_y_max * 1.1),
                 xlab = "Time Step", ylab = "",
                 main = "Suitable Habitat Pixels Over Time",
                 las = 1, xaxt = "n")
  graphics::mtext("Suitable Pixels", side = 2, line = 4.5, cex = 1)
  graphics::axis(1, at = x_ticks)
  for (u in all_units) {
    ud <- timestep_summary[timestep_summary$Spatial_Unit == u, ]
    ud <- ud[order(ud$Time_Step), ]
    graphics::lines(ud$Time_Step, ud$Pixels_Suitable,
                    col = unit_colors[u], lwd = 1.8)
    graphics::points(ud$Time_Step, ud$Pixels_Suitable,
                     col = unit_colors[u], pch = 19, cex = 0.7)
  }
  usr <- graphics::par("usr")
  x_legend <- usr[2] + (usr[2] - usr[1]) * 0.02
  y_legend <- usr[4]
  graphics::legend(x_legend, y_legend,
                   legend = all_units,
                   col    = unit_colors[all_units],
                   lwd    = 2, pch = 19,
                   bty    = "n", cex = 0.78,
                   xpd    = TRUE)
}

#' @keywords internal
#' @noRd
.plot_trend_change_per_step <- function(gain_v, loss_v, yrs_chg, ac_ylim) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(5, 8, 3.5, 2),
                mgp = c(1.8, 0.7, 0))
  graphics::barplot(gain_v,
                    names.arg = yrs_chg, las = 2,
                    ylim      = ac_ylim,
                    col       = "#267300",
                    border    = NA,
                    ylab      = "",
                    main      = "Annual Habitat Gains and Losses")
  graphics::mtext("Pixels", side = 2, line = 5, cex = 1)
  graphics::barplot(-loss_v,
                    add    = TRUE,
                    col    = "#E31A1C",
                    border = NA,
                    axes   = FALSE,
                    names.arg = rep("", length(yrs_chg)))
  graphics::abline(h = 0, col = "black", lwd = 0.8)
  graphics::legend("topright",
                   legend = c("Gain", "Loss"),
                   fill   = c("#267300", "#E31A1C"),
                   bty    = "n", cex = 0.85)
}

#' @keywords internal
#' @noRd
.plot_trend_total_unit <- function(g_tot, l_tot, shared_units, time_steps, tcu_ymax) {
  left_mar <- max(4, ceiling(max(nchar(shared_units)) / 2.5))
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(4, left_mar, 3.5, 1),
                mgp = c(1.8, 0.7, 0))

  bar_mat <- matrix(c(g_tot, l_tot), nrow = 2, byrow = TRUE,
                    dimnames = list(c("Gain", "Loss"), shared_units))

  graphics::barplot(
    bar_mat,
    beside    = TRUE,
    horiz     = TRUE,
    col       = c("#267300", "#E31A1C"),
    border    = NA,
    xlim      = c(0, tcu_ymax * 1.25),
    las       = 1,
    xlab      = "Total Pixels",
    main      = paste0("Total Gains and Losses by Unit (",
                       min(time_steps), "-", max(time_steps), ")")
  )
  graphics::legend("bottomright",
                   legend = c("Gain", "Loss"),
                   fill   = c("#267300", "#E31A1C"),
                   bty    = "n", cex = 0.9)
}

#' @keywords internal
#' @noRd
.plot_trend_facet_units <- function(change_by_timestep, all_units, n_rows, n_cols) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mfrow = c(n_rows, n_cols),
                mar   = c(3.5, 7, 2.5, 0.5),
                oma   = c(0, 0, 3, 0),
                mgp   = c(2.2, 0.6, 0))
  for (u in all_units) {
    ud    <- change_by_timestep[change_by_timestep$Spatial_Unit == u, ]
    ud    <- ud[order(ud$Time_Step), ]
    g_v   <- as.numeric(ud$Increase_Pixels)
    l_v   <- as.numeric(ud$Decrease_Pixels)
    y_max_raw <- suppressWarnings(max(c(g_v, l_v), na.rm = TRUE))
    y_max <- if (is.finite(y_max_raw) && y_max_raw > 0) y_max_raw * 1.2 else 1
    y_lim <- c(-y_max, y_max)

    graphics::barplot(g_v,
                      names.arg = ud$Time_Step,
                      ylim      = y_lim,
                      col       = "#267300",
                      border    = NA,
                      las       = 2,
                      cex.names = 0.55,
                      main      = u,
                      ylab      = "")
    graphics::mtext("Pixels", side = 2, line = 4.5, cex = 0.7)
    graphics::barplot(-l_v,
                      add       = TRUE,
                      col       = "#E31A1C",
                      border    = NA,
                      axes      = FALSE,
                      names.arg = rep("", nrow(ud)))
    graphics::abline(h = 0, col = "gray40", lwd = 0.6)
  }
  graphics::mtext("Gains and Losses by Unit",
                  outer = TRUE, side = 3, line = 1, cex = 1.1, font = 2)
}

#' @keywords internal
#' @noRd
.validate_variable_patterns <- function(variable_patterns) {
  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop(paste0(
      "ERROR: 'variable_patterns' must be a named character vector mapping variable names ",
      "to filename patterns (e.g. c(\"Var1\" = \"Var1_YEAR\", \"StaticVar\" = \"StaticVar\"))."
    ))
  }
  if (any(names(variable_patterns) == "")) {
    stop("ERROR: All elements in 'variable_patterns' must be named.")
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.classify_variable_patterns <- function(variable_patterns, time_cols) {
  if (is.null(time_cols)) time_cols <- character(0)
  time_cols <- as.character(time_cols)

  dynamic_vars <- character(0)
  for (v in names(variable_patterns)) {
    pat <- variable_patterns[v]
    is_dynamic <- FALSE
    for (tc in time_cols) {
      if (grepl(tc, pat, ignore.case = TRUE)) {
        is_dynamic <- TRUE
        break
      }
    }
    if (is_dynamic) dynamic_vars <- c(dynamic_vars, v)
  }
  static_vars <- setdiff(names(variable_patterns), dynamic_vars)

  if (length(dynamic_vars) > 0 && length(time_cols) > 0) {
    placeholders <- character(0)
    for (v in dynamic_vars) {
      parts <- strsplit(variable_patterns[v], "_")[[1]]
      placeholders <- c(placeholders, parts[toupper(parts) %in% toupper(time_cols)])
    }
    placeholders <- unique(placeholders)
    missing_in_patterns <- toupper(time_cols)[!toupper(time_cols) %in% toupper(placeholders)]
    if (length(missing_in_patterns) > 0) {
      warning("time_cols includes columns not found in variable_patterns: ",
              paste(missing_in_patterns, collapse = ", "))
    }
  }

  list(dynamic = dynamic_vars, static = static_vars)
}

#' @keywords internal
#' @noRd
.load_raster_input <- function(x, arg_name = "raster") {
  if (is.character(x)) {
    if (length(x) == 1 && dir.exists(x)) {
      tif_files <- list.files(x, pattern = "\\.tif$", full.names = TRUE)
      if (length(tif_files) == 0) {
        stop(paste0("ERROR: No .tif files found in directory '", x, "' for '", arg_name, "'."))
      }
      tif_files <- .natural_sort(tif_files)
      out <- terra::rast(tif_files)
      names(out) <- tools::file_path_sans_ext(basename(tif_files))
      return(out)
    }
    if (!file.exists(x)) {
      stop(paste0("ERROR: '", arg_name, "' file does not exist: '", x, "'."))
    }
    return(terra::rast(x))
  }
  if (inherits(x, "RasterLayer")) {
    warning(paste0("Converting '", arg_name, "' from 'RasterLayer' to 'SpatRaster' for terra compatibility."))
    return(terra::rast(x))
  }
  if (inherits(x, "SpatRaster")) {
    return(x)
  }
  stop(paste0("ERROR: '", arg_name, "' must be a file path, 'RasterLayer', or 'SpatRaster'. ",
              "Provided object is of class: ", class(x)[1]))
}

#' @keywords internal
#' @noRd
.load_shapefile_input <- function(x, arg_name = "shapefile") {
  if (inherits(x, "sf")) return(x)
  if (inherits(x, "sfc")) return(sf::st_sf(geometry = x))
  if (inherits(x, "Spatial")) return(sf::st_as_sf(x))
  if (is.character(x)) {
    if (!file.exists(x)) {
      stop(paste0("ERROR: '", arg_name, "' file not found: '", x, "'."))
    }
    return(sf::st_read(x, quiet = TRUE))
  }
  stop(paste0("ERROR: '", arg_name, "' must be an sf object, sfc, Spatial object, or file path. ",
              "Provided object is of class: ", class(x)[1]))
}

#' @keywords internal
#' @noRd
.resolve_threshold_method <- function(threshold_method) {
  if (is.numeric(threshold_method)) {
    if (length(threshold_method) != 1 || threshold_method <= 0 || threshold_method >= 1) {
      stop("ERROR: When supplying a numeric threshold, it must be a single value between 0 and 1.")
    }
    return(list(method = "manual", manual_value = threshold_method))
  }
  if (!threshold_method %in% c("prevalence", "tss", "manual")) {
    stop(paste0("ERROR: 'threshold_method' must be 'prevalence', 'tss', ",
                "or a numeric value between 0 and 1."))
  }
  list(method = threshold_method, manual_value = NULL)
}

#' @keywords internal
#' @noRd
.normalize_time_steps <- function(time_steps, fn_name) {
  secondary_filters <- list()
  if (is.data.frame(time_steps) || is.matrix(time_steps)) {
    time_steps_df <- as.data.frame(time_steps, stringsAsFactors = FALSE)
    if (ncol(time_steps_df) == 1) {
      time_steps <- time_steps_df[[1]]
    } else {
      primary_col    <- names(time_steps_df)[1]
      secondary_cols <- names(time_steps_df)[-1]
      for (sc in secondary_cols) {
        u_vals <- unique(time_steps_df[[sc]])
        if (length(u_vals) > 1) {
          stop(paste0(
            "ERROR: Secondary time column '", sc, "' has ", length(u_vals),
            " unique values (", paste(u_vals, collapse = ", "), "). ",
            fn_name, "() requires each secondary column to have ",
            "exactly one unique value so that raster layers are unambiguous. ",
            "Filter time_steps to a single value of '", sc, "' before calling this function. ",
            "For example: time_steps <- expand.grid(time_step = 1:15, season = \"Spring\", ",
            "stringsAsFactors = FALSE)"
          ))
        }
        secondary_filters[[sc]] <- u_vals
      }
      time_steps <- time_steps_df[[primary_col]]
    }
  }
  if (!is.atomic(time_steps) || is.list(time_steps)) {
    stop("ERROR: 'time_steps' must be a plain atomic vector of time step values.")
  }
  if (anyDuplicated(time_steps) > 0) {
    stop(paste0(
      "ERROR: 'time_steps' contains duplicate values. ",
      "Each raster layer must correspond to a unique time step. ",
      "If your predictions have repeated time steps across a secondary time column (e.g. season), ",
      "filter time_steps to a single value of that column before calling this function. ",
      "For example: time_steps <- expand.grid(time_step = 1:15, season = \"Spring\", ",
      "stringsAsFactors = FALSE)"
    ))
  }
  list(time_steps = time_steps, secondary_filters = secondary_filters)
}

#' @keywords internal
#' @noRd
.parse_model_formula <- function(f) {
  if (inherits(f, "formula")) {
    if (length(f) == 3) {
      warning(paste0(
        "model_formula appears to include a response variable on the left-hand side. ",
        "Only the right-hand side is used; 'presence' is always added automatically."
      ))
      return(stats::as.formula(
        paste("~", gsub("\\s+", " ", paste(deparse(f[[3]]), collapse = " ")))
      ))
    }
    return(f)
  }
  if (is.character(f) && length(f) == 1) {
    f <- trimws(f)
    if (!startsWith(f, "~")) {
      stop(paste0(
        "ERROR: model_formula as a character string must start with '~', ",
        "e.g. '~ Var1 + Var2'. Got: '", f, "'."
      ))
    }
    return(tryCatch(
      stats::as.formula(f),
      error = function(e) stop(paste0("ERROR: Could not parse model_formula '", f, "': ", e$message))
    ))
  }
  stop("ERROR: model_formula must be a formula object or a character string starting with '~'.")
}

#' @keywords internal
#' @noRd
.select_threshold_prob <- function(pred_prob, observed, method, manual_value) {
  if (method == "manual") return(manual_value)
  if (method == "prevalence") return(mean(observed, na.rm = TRUE))
  candidates <- seq(0.01, 0.99, by = 0.01)
  tss_vals <- vapply(candidates, function(thr) {
    pred_bin <- as.integer(pred_prob >= thr)
    tp   <- sum(pred_bin == 1 & observed == 1, na.rm = TRUE)
    fn   <- sum(pred_bin == 0 & observed == 1, na.rm = TRUE)
    tn   <- sum(pred_bin == 0 & observed == 0, na.rm = TRUE)
    fp   <- sum(pred_bin == 1 & observed == 0, na.rm = TRUE)
    sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
    sens + spec - 1
  }, numeric(1))
  candidates[which.max(tss_vals)]
}

#' @keywords internal
#' @noRd
.compute_confusion_metrics <- function(observed, predicted, pred_prob) {
  tp <- sum(predicted == 1 & observed == 1, na.rm = TRUE)
  fn <- sum(predicted == 0 & observed == 1, na.rm = TRUE)
  tn <- sum(predicted == 0 & observed == 0, na.rm = TRUE)
  fp <- sum(predicted == 1 & observed == 0, na.rm = TRUE)
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  tss         <- if (!is.na(sensitivity) && !is.na(specificity)) sensitivity + specificity - 1 else NA
  n     <- tp + tn + fp + fn
  po    <- (tp + tn) / n
  pe    <- ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / (n * n)
  kappa <- if ((1 - pe) > 0) (po - pe) / (1 - pe) else NA
  auc <- tryCatch({
    ord   <- order(pred_prob, decreasing = TRUE)
    o_obs <- observed[ord]
    n_pos <- sum(observed == 1, na.rm = TRUE)
    n_neg <- sum(observed == 0, na.rm = TRUE)
    if (n_pos == 0 || n_neg == 0) return(NA)
    tpr_v <- cumsum(o_obs == 1) / n_pos
    fpr_v <- cumsum(o_obs == 0) / n_neg
    sum(diff(fpr_v) * (tpr_v[-1] + tpr_v[-length(tpr_v)]) / 2)
  }, error = function(e) NA)
  list(tp = tp, fn = fn, tn = tn, fp = fp,
       sensitivity = sensitivity, specificity = specificity,
       tss = tss, kappa = kappa, auc = auc)
}

#' @keywords internal
#' @noRd
.compute_roc_curve <- function(observed, predicted) {
  thresholds <- sort(unique(predicted), decreasing = TRUE)
  tpr   <- numeric(length(thresholds))
  fpr   <- numeric(length(thresholds))
  n_pos <- sum(observed == 1, na.rm = TRUE)
  n_neg <- sum(observed == 0, na.rm = TRUE)
  for (i in seq_along(thresholds)) {
    pred_bin <- as.integer(predicted >= thresholds[i])
    tp <- sum(pred_bin == 1 & observed == 1, na.rm = TRUE)
    fp <- sum(pred_bin == 1 & observed == 0, na.rm = TRUE)
    tpr[i] <- if (n_pos > 0) tp / n_pos else 0
    fpr[i] <- if (n_neg > 0) fp / n_neg else 0
  }
  list(tpr = tpr, fpr = fpr, thresholds = thresholds)
}

#' @keywords internal
#' @noRd
.plot_response_curves_parametric <- function(model_fit, train_data, rhs_vars,
                                             var_means, fold_color, threshold,
                                             fold, n_rows, n_cols, model_label) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2.5, 1),
                oma = c(0, 0, 3, 0))

  pr_data <- train_data[train_data$presence == 1, , drop = FALSE]
  pa_data <- train_data[train_data$presence == 0, , drop = FALSE]
  for (v in rhs_vars) {
    v_range   <- seq(min(train_data[[v]], na.rm = TRUE),
                     max(train_data[[v]], na.rm = TRUE), length.out = 100)
    pred_grid <- as.data.frame(matrix(
      rep(var_means, each = 100), nrow = 100,
      dimnames = list(NULL, rhs_vars)
    ))
    pred_grid[[v]] <- v_range
    pred_prob <- stats::predict(model_fit, newdata = pred_grid, type = "response")
    graphics::plot(v_range, pred_prob, type = "l", lwd = 2,
                   col = fold_color, ylim = c(0, 1),
                   xlab = v, ylab = "P(presence)", main = v, las = 1)
    graphics::rug(pr_data[[v]], side = 3,
                  col = grDevices::adjustcolor("steelblue", alpha.f = 0.5),
                  ticksize = 0.03, lwd = 0.8)
    graphics::rug(pa_data[[v]], side = 1,
                  col = grDevices::adjustcolor("tomato", alpha.f = 0.5),
                  ticksize = 0.03, lwd = 0.8)
    if (!is.null(threshold)) {
      graphics::abline(h = threshold, lty = 2, col = "gray40", lwd = 1)
      graphics::legend("topright", legend = paste0("thr = ", round(threshold, 2)),
                       lty = 2, col = "gray40", bty = "n", cex = 0.8)
    }
  }
  graphics::mtext(
    paste0(model_label, " Fold ", fold, " Response Curves"),
    side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2
  )
}

#' @keywords internal
#' @noRd
.plot_roc_summary_parametric <- function(model_list, train_data_list,
                                         thresholds, fold_metrics_df,
                                         all_folds, fold_colors, model_label) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(4.5, 4.5, 3, 1.5))

  roc_curves <- lapply(all_folds, function(fold) {
    fold_key  <- paste0("fold", fold)
    pred_prob <- stats::predict(model_list[[fold_key]],
                                newdata = train_data_list[[fold_key]],
                                type = "response")
    .compute_roc_curve(observed = train_data_list[[fold_key]]$presence,
                       predicted = pred_prob)
  })
  graphics::plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "1 - Specificity (FPR)", ylab = "Sensitivity (TPR)",
                 main = paste0(model_label, " ROC Curves - All Folds"), las = 1)
  graphics::abline(a = 0, b = 1, lty = 2, col = "gray60")
  for (i in seq_along(all_folds)) {
    rc       <- roc_curves[[i]]
    fold_key <- paste0("fold", all_folds[i])
    graphics::lines(rc$fpr, rc$tpr, col = fold_colors[fold_key], lwd = 2)
    thr         <- thresholds[fold_key]
    closest_idx <- which.min(abs(rc$thresholds - thr))
    graphics::points(rc$fpr[closest_idx], rc$tpr[closest_idx],
                     pch = 19, cex = 1.3, col = fold_colors[fold_key])
  }
  graphics::legend(
    "bottomright",
    legend = paste0("Fold ", all_folds,
                    " (AUC = ", round(fold_metrics_df$AUC, 2), ")"),
    col = fold_colors, lwd = 2, pch = 19,
    bty = "n", cex = 0.85
  )
}

#' @keywords internal
#' @noRd
.plot_rf_importance <- function(model_vars, all_folds, importance_list, fold_colors) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  n_vars_loc  <- length(model_vars)
  n_cols_loc  <- min(n_vars_loc, 3)
  n_rows_loc  <- ceiling(n_vars_loc / n_cols_loc)
  bot_mar     <- max(nchar(paste0("Fold ", all_folds))) * 0.5 + 1.5

  graphics::par(mfrow = c(n_rows_loc, n_cols_loc),
                mar   = c(bot_mar, 4, 3, 1),
                oma   = c(0, 0, 3, 0))

  for (v in model_vars) {
    imp_vals <- vapply(all_folds, function(fold) {
      fk  <- paste0("fold", fold)
      df  <- importance_list[[fk]]
      if (is.null(df)) return(NA_real_)
      idx <- match(v, df$variable)
      if (is.na(idx)) NA_real_ else df$mean_decr_acc[idx]
    }, numeric(1))

    bar_cols  <- fold_colors[paste0("fold", all_folds)]
    fold_labs <- paste0("Fold ", all_folds)
    y_rng     <- c(0, max(imp_vals, na.rm = TRUE) * 1.15)

    graphics::barplot(
      imp_vals,
      names.arg = fold_labs,
      col       = bar_cols,
      border    = NA,
      horiz     = FALSE,
      las       = 2,
      ylim      = y_rng,
      ylab      = "Mean Decrease in Accuracy",
      main      = v,
      cex.names = 0.8,
      cex.main  = 0.95
    )
  }
  graphics::mtext(
    "RF Variable Importance by Fold",
    side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2
  )
}

#' @keywords internal
#' @noRd
.plot_rf_marginal_curves <- function(model_fit, train_data, model_vars,
                                     var_means, fold_color, threshold,
                                     fold, n_rows, n_cols) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 2.5, 1),
                oma = c(0, 0, 3, 0))

  pr_data <- train_data[train_data$presence == 1, , drop = FALSE]
  pa_data <- train_data[train_data$presence == 0, , drop = FALSE]
  for (v in model_vars) {
    v_range   <- seq(min(train_data[[v]], na.rm = TRUE),
                     max(train_data[[v]], na.rm = TRUE), length.out = 100)
    pred_grid <- as.data.frame(matrix(
      rep(var_means, each = 100), nrow = 100,
      dimnames = list(NULL, model_vars)
    ))
    pred_grid[[v]]     <- v_range
    pred_grid$presence <- factor(0, levels = c(0, 1))

    pred_prob <- tryCatch(
      stats::predict(model_fit, newdata = pred_grid, type = "prob")[, "1"],
      error = function(e) rep(NA_real_, 100)
    )

    graphics::plot(v_range, pred_prob, type = "l", lwd = 2,
                   col = fold_color, ylim = c(0, 1),
                   xlab = v, ylab = "P(presence)", main = v, las = 1)
    graphics::rug(pr_data[[v]], side = 3,
                  col = grDevices::adjustcolor("steelblue", alpha.f = 0.5),
                  ticksize = 0.03, lwd = 0.8)
    graphics::rug(pa_data[[v]], side = 1,
                  col = grDevices::adjustcolor("tomato", alpha.f = 0.5),
                  ticksize = 0.03, lwd = 0.8)
    if (!is.null(threshold)) {
      graphics::abline(h = threshold, lty = 2, col = "gray40", lwd = 1)
      graphics::legend("topright", legend = paste0("thr = ", round(threshold, 2)),
                       lty = 2, col = "gray40", bty = "n", cex = 0.8)
    }
  }
  graphics::mtext(
    paste0("RF Fold ", fold, " Marginal Prediction Curves"),
    side = 3, outer = TRUE, line = 1, cex = 1.1, font = 2
  )
}

#' @keywords internal
#' @noRd
#' @importFrom tools file_ext
.load_pseudoabsence_result <- function(pseudoabsence_result, verbose = TRUE) {
  if (is.character(pseudoabsence_result)) {
    if (!file.exists(pseudoabsence_result)) {
      stop(paste0("ERROR: File not found: ", pseudoabsence_result))
    }
    if (tolower(tools::file_ext(pseudoabsence_result)) != "rds") {
      stop("ERROR: pseudoabsence_result file must be .rds format.")
    }
    if (verbose) message(paste("Reading pseudoabsence results from:", basename(pseudoabsence_result)))
    pseudoabsence_result <- tryCatch(
      readRDS(pseudoabsence_result),
      error = function(e) stop(paste0("ERROR reading pseudoabsence_result .rds: ", e$message))
    )
  }
  if (!is.list(pseudoabsence_result) || !"pseudoabsences" %in% names(pseudoabsence_result)) {
    stop(paste0("ERROR: pseudoabsence_result must contain a 'pseudoabsences' element. ",
                "Re-run generate_absences()."))
  }
  pseudoabsence_result
}

#' @keywords internal
#' @noRd
#' @importFrom tools file_ext
.load_partition_result <- function(partition_result, verbose = TRUE) {
  if (is.character(partition_result)) {
    if (!file.exists(partition_result)) {
      stop(paste0("ERROR: File not found: ", partition_result))
    }
    if (tolower(tools::file_ext(partition_result)) != "rds") {
      stop("ERROR: partition_result file must be .rds format.")
    }
    if (verbose) message(paste("Reading partition results from:", basename(partition_result)))
    partition_result <- tryCatch(
      readRDS(partition_result),
      error = function(e) stop(paste0("ERROR reading partition_result .rds: ", e$message))
    )
  }
  if (!is.list(partition_result) || !"points_sf" %in% names(partition_result)) {
    stop(paste0("ERROR: partition_result must contain a 'points_sf' element. ",
                "Re-run spatiotemporal_partition()."))
  }
  partition_result
}

#' @keywords internal
#' @noRd
.prepare_combined_modeling_data <- function(partition_result, pseudoabsence_result,
                                            model_vars, time_cols = NULL,
                                            verbose = TRUE) {

  presence_sf <- partition_result$points_sf
  if (is.null(presence_sf) || nrow(presence_sf) == 0) {
    stop("ERROR: 'partition_result$points_sf' is empty.")
  }
  if (!"fold" %in% names(presence_sf)) {
    stop("ERROR: 'partition_result$points_sf' is missing 'fold' column. Re-run spatiotemporal_partition().")
  }
  presence_df <- sf::st_drop_geometry(presence_sf)
  if (!"presence" %in% names(presence_df)) {
    if (verbose) message("No 'presence' column found in partition_result$points_sf. Adding presence = 1 for all rows.")
    presence_df$presence <- 1
  }

  pseudoabs_sf <- pseudoabsence_result$pseudoabsences
  if (is.null(pseudoabs_sf) || nrow(pseudoabs_sf) == 0) {
    stop("ERROR: 'pseudoabsence_result$pseudoabsences' is empty.")
  }
  pseudoabs_df <- sf::st_drop_geometry(pseudoabs_sf)
  if (!"presence" %in% names(pseudoabs_df)) {
    pseudoabs_df$presence <- 0
  }
  if (!"fold" %in% names(pseudoabs_df)) {
    stop("ERROR: 'pseudoabsence_result$pseudoabsences' is missing 'fold' column. Re-run generate_absences().")
  }

  missing_pr <- model_vars[!model_vars %in% names(presence_df)]
  missing_pa <- model_vars[!model_vars %in% names(pseudoabs_df)]
  if (length(missing_pr) > 0) {
    stop(paste0(
      "ERROR: The following predictors are missing from the presence data: ",
      paste(missing_pr, collapse = ", "),
      ". Run temporally_explicit_extraction() before partitioning, or check variable names."
    ))
  }
  if (length(missing_pa) > 0) {
    stop(paste0(
      "ERROR: The following predictors are missing from the pseudoabsence data: ",
      paste(missing_pa, collapse = ", "),
      ". Ensure generate_absences() was run with a raster_dir."
    ))
  }

  keep_cols    <- unique(c("fold", "presence", model_vars, time_cols))
  keep_cols_pr <- keep_cols[keep_cols %in% names(presence_df)]
  keep_cols_pa <- keep_cols[keep_cols %in% names(pseudoabs_df)]

  combined_df <- rbind(
    presence_df[, keep_cols_pr, drop = FALSE],
    pseudoabs_df[, keep_cols_pa, drop = FALSE]
  )

  combined_df          <- combined_df[!is.na(combined_df$fold), ]
  combined_df$presence <- as.integer(combined_df$presence)

  complete_rows <- stats::complete.cases(combined_df[, model_vars, drop = FALSE])
  if (sum(!complete_rows) > 0) {
    if (verbose) message(paste0("Removing ", sum(!complete_rows), " rows with NA values in model predictors."))
    combined_df <- combined_df[complete_rows, ]
  }
  if (nrow(combined_df) == 0) {
    stop("ERROR: No complete cases remain after removing NA values in model predictors.")
  }

  all_folds <- sort(unique(combined_df$fold))
  if (length(all_folds) == 0) {
    stop("ERROR: No valid fold assignments found. Re-run spatiotemporal_partition().")
  }
  single_fold_mode <- length(all_folds) == 1
  if (single_fold_mode && verbose) {
    message("Single-fold mode: all points used for both training and testing.")
  }

  list(combined_df = combined_df, all_folds = all_folds,
       single_fold_mode = single_fold_mode)
}

#' @keywords internal
#' @noRd
.build_fold_metric_row <- function(fold, thr, metrics) {
  data.frame(
    Fold          = fold,
    Threshold     = round(thr, 4),
    Testing_TP    = metrics$tp,
    Testing_FN    = metrics$fn,
    Testing_TN    = metrics$tn,
    Testing_FP    = metrics$fp,
    Sensitivity   = round(metrics$sensitivity, 4),
    Specificity   = round(metrics$specificity, 4),
    TSS           = round(metrics$tss,         4),
    Kappa         = round(metrics$kappa,       4),
    AUC           = round(metrics$auc,         4),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.plot_fold_metrics_bars <- function(metrics_df, rate_cols, fold_colors, title) {

  n_metrics <- length(rate_cols)

  ref_h <- vapply(rate_cols, function(m) {
    if (m %in% c("TSS", "Kappa")) 0 else NA_real_
  }, numeric(1))

  fold_ids <- metrics_df$Fold

  bar_cols <- vapply(fold_ids, function(f) {
    key <- paste0("fold", f)
    if (key %in% names(fold_colors)) fold_colors[[key]]
    else if (as.character(f) %in% names(fold_colors)) fold_colors[[as.character(f)]]
    else "gray70"
  }, character(1))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  graphics::par(
    mfrow  = c(1, n_metrics),
    mar    = c(4, 4, 4, 1),
    oma    = c(0, 0, 3, 0)
  )

  for (mi in seq_along(rate_cols)) {
    col_name <- rate_cols[mi]
    vals     <- metrics_df[[col_name]]
    ref      <- ref_h[mi]

    y_lim <- c(
      max(min(c(0, if (!is.na(ref)) ref, vals), na.rm = TRUE) - 0.05, -0.1),
      min(max(c(1, if (!is.na(ref)) ref, vals), na.rm = TRUE) + 0.05,  1.1)
    )

    graphics::barplot(
      vals,
      names.arg = paste0("Fold ", fold_ids),
      col       = bar_cols,
      border    = "gray30",
      ylim      = y_lim,
      las       = 1,
      main      = col_name,
      ylab      = if (mi == 1) "Value" else "",
      cex.names = 0.85,
      cex.main  = 1.1,
      cex.axis  = 0.9
    )

    if (!is.na(ref)) {
      graphics::abline(h = ref, col = "gray50", lty = 3, lwd = 0.9)
    }

    mean_val <- mean(vals, na.rm = TRUE)
    graphics::abline(h = mean_val, col = "black", lty = 2, lwd = 1.8)
    graphics::mtext(paste0("mean = ", round(mean_val, 3)),
                    side = 3, line = 0, adj = 1, cex = 0.72, col = "gray30")
  }

  graphics::mtext(title, side = 3, outer = TRUE, line = 1, cex = 1.05, font = 2)
}

#' @keywords internal
#' @noRd
.compute_timestep_metrics <- function(binary_raster, test_pts_time,
                                      pseudo_pts_time = NULL,
                                      model = NULL, model_vars = NULL,
                                      threshold = NULL, model_type = NULL,
                                      time_values, time_cols, fold) {

  total_cells    <- terra::ncell(binary_raster) -
    sum(is.na(terra::values(binary_raster)), na.rm = TRUE)
  suitable_cells <- sum(terra::values(binary_raster) == 1, na.rm = TRUE)
  pct_suitable   <- if (total_cells > 0) suitable_cells / total_cells else NA_real_

  n_pres <- nrow(test_pts_time)
  n_abs  <- if (!is.null(pseudo_pts_time)) nrow(pseudo_pts_time) else 0
  has_pa <- !is.null(model) && !is.null(threshold) && !is.null(model_vars) &&
    !is.null(pseudo_pts_time) && n_abs > 0

  tp <- NA_integer_; fn <- NA_integer_
  tn <- NA_integer_; fp <- NA_integer_
  sensitivity <- NA_real_; specificity <- NA_real_
  tss <- NA_real_;  cbp <- NA_real_

  if (n_pres > 0 || has_pa) {

    if (has_pa) {

      pr_df <- sf::st_drop_geometry(test_pts_time)
      pa_df <- sf::st_drop_geometry(pseudo_pts_time)

      pr_has_vars <- all(model_vars %in% names(pr_df))
      pa_has_vars <- all(model_vars %in% names(pa_df))

      if (pr_has_vars || pa_has_vars) {
        rows <- list()
        if (pr_has_vars && n_pres > 0) {
          r <- pr_df[, model_vars, drop = FALSE]; r$obs <- 1; rows[["pr"]] <- r
        }
        if (pa_has_vars && n_abs > 0) {
          r <- pa_df[, model_vars, drop = FALSE]; r$obs <- 0; rows[["pa"]] <- r
        }
        combined      <- do.call(rbind, rows)
        complete_rows <- stats::complete.cases(combined[, model_vars, drop = FALSE])
        combined      <- combined[complete_rows, , drop = FALSE]

        if (nrow(combined) > 0) {
          pred_sub <- combined[, model_vars, drop = FALSE]
          raw <- tryCatch({
            if (model_type == "rf") {
              pred_sub$presence <- factor(0, levels = c(0, 1))
              stats::predict(model, newdata = pred_sub, type = "prob")[, "1"]
            } else {
              stats::predict(model, newdata = pred_sub, type = "response")
            }
          }, error = function(e) {
            warning(paste("Timestep prediction failed:", e$message))
            NULL
          })

          if (!is.null(raw)) {
            pred_bin <- as.integer(raw >= threshold)
            obs      <- combined$obs
            tp  <- sum(pred_bin == 1 & obs == 1, na.rm = TRUE)
            fn  <- sum(pred_bin == 0 & obs == 1, na.rm = TRUE)
            tn  <- sum(pred_bin == 0 & obs == 0, na.rm = TRUE)
            fp  <- sum(pred_bin == 1 & obs == 0, na.rm = TRUE)
            n_p <- tp + fn; n_n <- tn + fp
            sensitivity <- if (n_p > 0) tp / n_p else NA_real_
            specificity <- if (n_n > 0) tn / n_n else NA_real_
            tss         <- if (!is.na(sensitivity) && !is.na(specificity)) {
              sensitivity + specificity - 1
            } else NA_real_
            cbp <- if (n_p > 0 && total_cells > 0) {
              stats::dbinom(tp, size = n_p, prob = suitable_cells / total_cells)
            } else NA_real_
          }
        }
      }

    } else if (n_pres > 0) {
      pts_vect  <- tryCatch(suppressWarnings(terra::vect(test_pts_time)), error = function(e) NULL)
      extracted <- if (!is.null(pts_vect)) {
        tryCatch(terra::extract(binary_raster, pts_vect)[, 2],
                 error = function(e) rep(NA_real_, n_pres))
      } else rep(NA_real_, n_pres)

      extracted[is.na(extracted)] <- 0
      tp  <- sum(extracted == 1, na.rm = TRUE)
      fn  <- sum(extracted == 0, na.rm = TRUE)
      tot <- tp + fn
      sensitivity <- if (tot > 0) tp / tot else NA_real_
      cbp <- if (tot > 0 && total_cells > 0) {
        stats::dbinom(tp, size = tot, prob = suitable_cells / total_cells)
      } else NA_real_
    }
  }

  row <- data.frame(
    Fold         = fold,
    Pct_Suitable = round(pct_suitable, 4),
    N_Pres       = n_pres,
    TP           = tp,
    FN           = fn,
    Sensitivity  = round(sensitivity, 4),
    CBP          = cbp,
    N_Abs        = if (has_pa) n_abs                else NA_integer_,
    TN           = if (has_pa) tn                   else NA_integer_,
    FP           = if (has_pa) fp                   else NA_integer_,
    Specificity  = if (has_pa) round(specificity, 4) else NA_real_,
    TSS          = if (has_pa) round(tss,          4) else NA_real_,
    stringsAsFactors = FALSE
  )

  for (tc in time_cols) row[[tc]] <- time_values[[tc]]
  row
}

#' @keywords internal
#' @noRd
.resolve_raster_paths <- function(variable_patterns, dynamic_vars, static_vars,
                                  time_cols, time_values, raster_dir) {
  all_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
  paths     <- character(0)

  for (var in dynamic_vars) {
    fname <- variable_patterns[[var]]
    for (tc in time_cols) {
      if (grepl(tc, fname, ignore.case = TRUE)) {
        fname <- gsub(tc, as.character(time_values[[tc]]), fname, ignore.case = TRUE)
      }
    }

    anchored_pat <- paste0(fname, "(?=[^0-9]|$)")
    matches <- all_files[grepl(anchored_pat, basename(all_files),
                               ignore.case = TRUE, perl = TRUE)]
    if (length(matches) == 0) {
      matches <- all_files[grepl(fname, basename(all_files), ignore.case = TRUE)]
    }
    if (length(matches) == 0) {
      warning(paste0("Missing raster file for dynamic variable '", var,
                     "' matching pattern '", fname, "'."))
      return(NULL)
    }
    if (length(matches) > 1) {
      stop(paste0(
        "ERROR: Multiple raster files found for dynamic variable '", var,
        "' matching pattern '", fname, "'.",
        "\nMatching files:\n", paste(matches, collapse = "\n")
      ))
    }
    paths <- c(paths, matches)
  }

  for (var in static_vars) {
    fname   <- variable_patterns[[var]]
    matches <- all_files[grepl(fname, basename(all_files), ignore.case = TRUE)]
    if (length(matches) == 0) {
      warning(paste0("Missing raster file for static variable '", var,
                     "' matching pattern '", fname, "'."))
      return(NULL)
    }
    if (length(matches) > 1) {
      stop(paste0(
        "ERROR: Multiple raster files found for static variable '", var,
        "' matching pattern '", fname, "'.",
        "\nMatching files:\n", paste(matches, collapse = "\n")
      ))
    }
    paths <- c(paths, matches)
  }

  paths
}

#' @keywords internal
#' @noRd
.predict_hypervolume_raster <- function(model, rast_stack, fold, time_label) {
  dims     <- model@Dimensionality
  n_layers <- terra::nlyr(rast_stack)
  if (dims != n_layers) {
    warning(paste0("Fold ", fold, ": hypervolume dimensionality (", dims,
                   ") != raster layers (", n_layers, ") at ", time_label, " -- skipping."))
    return(NULL)
  }
  tryCatch({
    proj <- suppressMessages(suppressWarnings(hypervolume::hypervolume_project(
      model, rasters = rast_stack, type = "inclusion",
      fast.or.accurate = "fast", verbose = FALSE
    )))
    out <- rast_stack[[1]]
    v   <- as.integer(as.vector(proj))
    if (length(v) == terra::ncell(out)) {
      terra::values(out) <- v
    } else {
      terra::values(out) <- NA_integer_
    }
    out
  }, error = function(e) {
    warning(paste0("Fold ", fold, ": hypervolume projection failed at ",
                   time_label, ": ", e$message))
    NULL
  })
}

#' @keywords internal
#' @noRd
.bstart_match <- function(nms, pfx) {
  which(startsWith(nms, pfx) &
          (nchar(nms) == nchar(pfx) |
             !grepl("^[A-Za-z0-9]", substr(nms, nchar(pfx) + 1, nchar(pfx) + 1))))
}

#' @keywords internal
#' @noRd
.predict_parametric_raster <- function(model, rast_stack, model_vars, threshold,
                                       model_type, variable_patterns = NULL,
                                       time_cols = NULL, time_values = NULL) {
  cell_vals  <- as.data.frame(terra::values(rast_stack))
  rast_names <- names(cell_vals)

  name_map <- vapply(model_vars, function(v) {
    if (!is.null(variable_patterns) && v %in% names(variable_patterns) &&
        !is.null(time_cols) && !is.null(time_values)) {
      fname <- variable_patterns[[v]]
      for (tc in time_cols) {
        if (grepl(tc, fname, ignore.case = TRUE)) {
          fname <- gsub(tc, as.character(time_values[[tc]]), fname, ignore.case = TRUE)
        }
      }
      for (candidate in c(paste0(fname, "_Scaled"), fname)) {
        m <- which(rast_names == candidate)
        if (length(m) > 0) return(rast_names[m[1]])
      }
      m <- .bstart_match(rast_names, paste0(fname, "_Scaled"))
      if (length(m) > 0) return(rast_names[m[1]])
      m <- .bstart_match(rast_names, fname)
      if (length(m) > 0) return(rast_names[m[1]])
    }
    m <- which(rast_names == v)
    if (length(m) > 0) return(rast_names[m[1]])
    m <- .bstart_match(rast_names, v)
    if (length(m) > 0) return(rast_names[m[1]])
    NA_character_
  }, character(1))

  unmatched <- model_vars[is.na(name_map)]
  if (length(unmatched) > 0) {
    warning(paste0("Could not match model variables to raster layers: ",
                   paste(unmatched, collapse = ", "),
                   ". Available: ", paste(rast_names, collapse = ", ")))
    return(NULL)
  }

  pred_df        <- cell_vals[, name_map, drop = FALSE]
  names(pred_df) <- model_vars
  complete_idx   <- which(stats::complete.cases(pred_df))
  pred_vals      <- rep(NA_real_, nrow(pred_df))

  if (length(complete_idx) > 0) {
    pred_sub <- pred_df[complete_idx, , drop = FALSE]
    raw_pred <- tryCatch({
      if (model_type == "rf") {
        pred_sub$presence <- factor(0, levels = c(0, 1))
        stats::predict(model, newdata = pred_sub, type = "prob")[, "1"]
      } else {
        stats::predict(model, newdata = pred_sub, type = "response")
      }
    }, error = function(e) { warning(paste("Raster prediction failed:", e$message)); NULL })

    if (!is.null(raw_pred)) {
      pred_vals[complete_idx] <- as.integer(raw_pred >= threshold)
    }
  }

  out_rast <- rast_stack[[1]]
  terra::values(out_rast) <- pred_vals
  out_rast
}

#' @keywords internal
#' @noRd
.tseries <- function(x_all, y_all, folds_all, x_sum, y_sum,
                     all_folds, fold_colors, col_mean,
                     col_cbp_good = NULL, col_cbp_bad = NULL,
                     fold_alpha = 0.3, fold_lwd = 0.8,
                     sum_col = col_mean, sum_lwd = 2,
                     cbp_good_sum = NULL) {
  for (f in all_folds) {
    idx <- !is.na(folds_all) & folds_all == f
    if (sum(idx) == 0) next
    fc  <- grDevices::adjustcolor(fold_colors[as.character(f)], alpha.f = fold_alpha)
    xf  <- x_all[idx]
    yf  <- y_all[idx]
    v   <- !is.na(xf) & !is.na(yf)
    if (sum(v) > 0) {
      ord <- order(xf[v])
      graphics::lines(xf[v][ord], yf[v][ord], col = fc, lwd = fold_lwd)
      graphics::points(xf[v][ord], yf[v][ord], col = fc, pch = 16, cex = 0.5)
    }
  }
  v2 <- !is.na(x_sum) & !is.na(y_sum)
  if (sum(v2) > 0) {
    ord2 <- order(x_sum[v2])
    graphics::lines(x_sum[v2][ord2], y_sum[v2][ord2], col = sum_col, lwd = sum_lwd)
    if (!is.null(cbp_good_sum)) {
      pt_col <- ifelse(cbp_good_sum[v2], col_cbp_good, col_cbp_bad)
      graphics::points(x_sum[v2][ord2], y_sum[v2][ord2], col = pt_col, pch = 21, bg = "white",
                       cex = 1.1, lwd = 1.2)
    } else {
      graphics::points(x_sum[v2][ord2], y_sum[v2][ord2], col = sum_col, pch = 19, cex = 1.0)
    }
  }
}

#' @keywords internal
#' @noRd
.ref_lines <- function(lookup_df, val_col, all_folds, fold_colors,
                       lty = 2, lwd = 1.2, alpha = 0.65) {
  if (is.null(lookup_df) || !val_col %in% names(lookup_df)) return(invisible(NULL))
  for (f in all_folds) {
    rows <- lookup_df[!is.na(lookup_df$Fold) & lookup_df$Fold == f, ]
    if (nrow(rows) == 0) next
    val <- rows[[val_col]][1]
    if (is.na(val)) next
    graphics::abline(h = val,
                     col = grDevices::adjustcolor(fold_colors[as.character(f)], alpha.f = alpha),
                     lty = lty, lwd = lwd)
  }
}

#' @keywords internal
#' @noRd
.draw_xaxis <- function(x_tck, x_labels_attr) {
  if (!is.null(x_labels_attr)) {
    graphics::axis(1, at = x_tck,
                   labels = x_labels_attr[x_tck],
                   las = 2, cex.axis = 0.7)
  } else {
    graphics::axis(1, at = x_tck)
  }
}

#' @keywords internal
#' @noRd
.new_plot <- function(y_vals, y_ref_vals = NULL, title, ylab,
                      x_lim, x_axis_label, x_tck, x_labels_attr,
                      y_floor = NULL, y_ceil = NULL, right_mar = 12,
                      draw_fn = NULL) {
  all_y <- c(y_vals, y_ref_vals)
  y_rng <- range(all_y, na.rm = TRUE)
  y_pad <- diff(y_rng) * 0.08
  y_lo  <- if (!is.null(y_floor)) max(y_floor, y_rng[1] - y_pad) else y_rng[1] - y_pad
  y_hi  <- if (!is.null(y_ceil))  min(y_ceil,  y_rng[2] + y_pad) else y_rng[2] + y_pad
  if (is.na(y_lo) || is.na(y_hi) || y_lo >= y_hi) { y_lo <- 0; y_hi <- 1 }
  bot_mar <- if (!is.null(x_labels_attr)) {
    max(nchar(as.character(x_labels_attr[x_tck])), na.rm = TRUE) * 0.4 + 2.5
  } else {
    4
  }
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mar = c(bot_mar, 5, 3.5, right_mar), xpd = FALSE)
  if (!is.null(x_labels_attr)) {
    graphics::plot(NULL, xlim = x_lim, ylim = c(y_lo, y_hi),
                   xlab = "", ylab = ylab, main = title, las = 1, xaxt = "n")
    graphics::mtext(x_axis_label, side = 1, line = bot_mar - 1.2, cex = 0.9)
  } else {
    graphics::plot(NULL, xlim = x_lim, ylim = c(y_lo, y_hi),
                   xlab = x_axis_label, ylab = ylab, main = title, las = 1, xaxt = "n")
  }
  .draw_xaxis(x_tck, x_labels_attr)
  if (!is.null(draw_fn)) draw_fn()
}

#' @keywords internal
#' @noRd
.plot_facet_metric <- function(an, metric_col, ylab, title_prefix,
                               secondary_combos, primary_col, secondary_cols,
                               fold_colors, col_mean,
                               y_floor = NULL, y_ceil = NULL) {
  n_combos <- nrow(secondary_combos)

  all_y <- numeric(0)
  for (ci in seq_len(n_combos)) {
    combo <- secondary_combos[ci, , drop = FALSE]
    idx   <- rep(TRUE, nrow(an))
    for (sc in secondary_cols) idx <- idx & an[[sc]] == combo[[sc]]
    sub   <- an[idx, , drop = FALSE]
    if (metric_col %in% names(sub)) all_y <- c(all_y, sub[[metric_col]])
  }
  if (length(all_y) == 0 || all(is.na(all_y))) return(invisible(NULL))
  y_rng <- range(all_y, na.rm = TRUE)
  y_pad <- diff(y_rng) * 0.08
  y_lo  <- if (!is.null(y_floor)) max(y_floor, y_rng[1] - y_pad) else y_rng[1] - y_pad
  y_hi  <- if (!is.null(y_ceil))  min(y_ceil,  y_rng[2] + y_pad) else y_rng[2] + y_pad
  if (is.na(y_lo) || is.na(y_hi) || y_lo >= y_hi) { y_lo <- 0; y_hi <- 1 }

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mfrow = c(n_combos, 1), mar = c(4, 5, 2.5, 1), oma = c(0, 0, 2.5, 0))
  for (ci in seq_len(n_combos)) {
    combo       <- secondary_combos[ci, , drop = FALSE]
    combo_label <- paste(unlist(combo), collapse = "_")
    idx   <- rep(TRUE, nrow(an))
    for (sc in secondary_cols) idx <- idx & an[[sc]] == combo[[sc]]
    sub   <- an[idx, , drop = FALSE]
    if (nrow(sub) == 0 || !metric_col %in% names(sub)) next
    x_lim_f <- range(sub[[primary_col]], na.rm = TRUE)
    graphics::plot(NULL, xlim = x_lim_f, ylim = c(y_lo, y_hi),
                   xlab = primary_col, ylab = ylab,
                   main = combo_label, las = 1)
    for (fld in sort(unique(sub$Fold))) {
      fld_key <- as.character(fld)
      fd      <- sub[sub$Fold == fld, , drop = FALSE]
      fd      <- fd[order(fd[[primary_col]]), ]
      if (nrow(fd) == 0) next
      graphics::lines(fd[[primary_col]], fd[[metric_col]],
                      col = fold_colors[fld_key], lwd = 1.5)
      graphics::points(fd[[primary_col]], fd[[metric_col]],
                       col = fold_colors[fld_key], pch = 19, cex = 0.6)
    }
    mn_vals <- tapply(sub[[metric_col]], sub[[primary_col]], mean, na.rm = TRUE)
    graphics::lines(as.numeric(names(mn_vals)), as.numeric(mn_vals),
                    col = col_mean, lwd = 2.5, lty = 1)
  }
  graphics::mtext(paste(title_prefix, "by", paste(secondary_cols, collapse = "/")),
                  side = 3, outer = TRUE, line = 1, cex = 1, font = 2)
}

#' @keywords internal
#' @noRd
.line_legend <- function(lab, col, lty, lwd) {
  usr <- graphics::par("usr")
  x_leg <- usr[2] + diff(usr[1:2]) * 0.03
  y_leg <- usr[4]
  old_xpd <- graphics::par(xpd = NA)
  on.exit(graphics::par(xpd = old_xpd$xpd), add = TRUE)
  graphics::legend(
    x         = x_leg,
    y         = y_leg,
    legend    = lab,
    col       = col,
    lty       = lty,
    lwd       = lwd,
    bty       = "n",
    bg        = "transparent",
    cex       = 0.82,
    seg.len   = 1.2,
    xpd       = NA,
    yjust     = 1,
    x.intersp = 0.5
  )
}

#' @keywords internal
#' @noRd
.bar_legend <- function(bar_labs, bar_cols, pt_lab, pt_col) {
  usr   <- graphics::par("usr")
  x_leg <- usr[2] + diff(usr[1:2]) * 0.03
  y_leg <- usr[4]
  old_xpd <- graphics::par(xpd = NA)
  on.exit(graphics::par(xpd = old_xpd$xpd), add = TRUE)

  graphics::legend(
    x         = x_leg,
    y         = y_leg,
    legend    = bar_labs,
    fill      = bar_cols,
    border    = "gray30",
    bty       = "n",
    bg        = "transparent",
    cex       = 0.82,
    xpd       = NA,
    yjust     = 1,
    x.intersp = 0.5
  )

  n_bar  <- length(bar_labs)
  lh_est <- diff(usr[3:4]) * 0.065
  y_pt   <- y_leg - n_bar * lh_est
  graphics::legend(
    x         = x_leg,
    y         = y_pt,
    legend    = pt_lab,
    col       = "gray20",
    pch       = 21,
    pt.bg     = pt_col,
    bty       = "n",
    bg        = "transparent",
    cex       = 0.82,
    xpd       = NA,
    yjust     = 1,
    x.intersp = 0.5
  )
}

#' @keywords internal
#' @noRd
.format_time_estimate <- function(seconds) {
  hours   <- seconds / 3600
  minutes <- seconds / 60
  if (hours > 1) {
    paste0(round(hours,   1), " hours")
  } else if (minutes > 1) {
    paste0(round(minutes, 1), " minutes")
  } else {
    paste0(round(seconds, 1), " seconds")
  }
}
