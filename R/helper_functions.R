generate_absences <- function(partition_result,
                              reference_shapefile_path,
                              raster_dir,
                              variable_patterns,
                              method                   = "random",
                              ratio                    = 1,
                              buffer_distance          = NULL,
                              env_percentile           = 0.05,
                              time_cols                = NULL,
                              pseudoabsence_times      = NULL,
                              min_points_per_timestep  = 1,
                              user_absence_data        = NULL,
                              xcol                     = NULL,
                              ycol                     = NULL,
                              points_crs               = NULL,
                              create_plot              = TRUE,
                              plot_by_fold             = FALSE,
                              plot_palette             = "Dark 2",
                              output_file              = NULL,
                              verbose                  = TRUE) {

  if (missing(partition_result)) {
    stop(paste0("ERROR: 'partition_result' is required. Provide output from ",
                "spatiotemporal_partition() or a path to an .rds file."))
  }

  if (missing(reference_shapefile_path)) {
    stop("ERROR: 'reference_shapefile_path' is required.")
  }

  if (missing(raster_dir)) {
    stop("ERROR: 'raster_dir' is required for all methods.")
  }

  if (missing(variable_patterns)) {
    stop(paste0("ERROR: 'variable_patterns' is required. Provide a named character vector ",
                "mapping variable names to filename patterns, ",
                "e.g. c(\"Var1\" = \"Var1_YEAR\", \"StaticVar\" = \"StaticVar\")."))
  }

  .validate_variable_patterns(variable_patterns)

  if (!dir.exists(raster_dir)) {
    stop(paste0("ERROR: 'raster_dir' does not exist: ", raster_dir))
  }

  if (!is.null(output_file) && tolower(tools::file_ext(output_file)) != "rds") {
    stop("ERROR: 'output_file' must have a .rds extension.")
  }

  method <- match.arg(method, c("random", "buffer", "environmental", "user_data"))

  if (method == "user_data" && is.null(user_absence_data)) {
    stop("ERROR: 'user_absence_data' is required when method = 'user_data'.")
  }

  if (!is.numeric(ratio) || ratio < 0) {
    stop("ERROR: 'ratio' must be a non-negative number.")
  }

  if (!is.numeric(min_points_per_timestep) || min_points_per_timestep < 0) {
    stop("ERROR: 'min_points_per_timestep' must be a non-negative number.")
  }

  if (ratio == 0 && min_points_per_timestep == 0) {
    stop(paste0("ERROR: 'ratio' and 'min_points_per_timestep' cannot both be 0. ",
                "Set 'ratio' > 0 for proportional allocation, or ",
                "'min_points_per_timestep' > 0 for a fixed number of points per time step."))
  }

  if (method == "buffer" && is.null(buffer_distance)) {
    stop("ERROR: 'buffer_distance' is required when method = 'buffer'.")
  }

  if (!is.null(pseudoabsence_times)) {
    if (!is.vector(pseudoabsence_times) || is.list(pseudoabsence_times)) {
      stop("ERROR: 'pseudoabsence_times' must be a plain vector of time step values.")
    }
  }

  partition_result <- .load_partition_result(partition_result, verbose = verbose)

  use_buffer_with_env <- method == "environmental" && !is.null(buffer_distance)

  time_cols    <- if (is.null(time_cols)) character(0) else as.character(time_cols)
  var_classes  <- .classify_variable_patterns(variable_patterns, time_cols)
  dynamic_vars <- var_classes$dynamic
  static_vars  <- var_classes$static

  if (verbose) message(paste("Dynamic variables:", if (length(dynamic_vars) > 0) paste(dynamic_vars, collapse = ", ") else "none"))
  if (verbose) message(paste("Static variables:",  if (length(static_vars)  > 0) paste(static_vars,  collapse = ", ") else "none"))

  reference_shapefile <- .load_shapefile_input(reference_shapefile_path, "reference_shapefile_path")

  pts_sf <- partition_result$points_sf
  if (!"fold" %in% names(pts_sf)) {
    stop("ERROR: 'partition_result$points_sf' must contain a 'fold' column.")
  }

  reference_shapefile <- sf::st_transform(reference_shapefile, crs = sf::st_crs(pts_sf))
  ref_union           <- sf::st_union(reference_shapefile)
  all_folds           <- sort(unique(pts_sf$fold[!is.na(pts_sf$fold)]))

  has_time <- length(time_cols) > 0 && all(time_cols %in% names(pts_sf))
  if (length(time_cols) > 0 && !has_time) {
    missing_tc <- time_cols[!time_cols %in% names(pts_sf)]
    warning(paste0("time_cols column(s) not found in points_sf: ",
                   paste(missing_tc, collapse = ", "),
                   ". Pseudoabsences will not be matched to time steps."))
  }

  if (!is.null(pseudoabsence_times) && !has_time) {
    stop("ERROR: 'pseudoabsence_times' requires 'time_cols' to be specified and present in the data.")
  }

  if (!is.null(pseudoabsence_times)) {
    pseudoabsence_times <- sort(unique(pseudoabsence_times))
  }

  pts_sf$presence <- 1

  ### Helper to build a ts_key string from a combo list
  .make_ts_key <- function(combo) {
    if (has_time && length(combo) > 0) {
      paste(vapply(time_cols, function(tc) as.character(combo[[tc]]), character(1)), collapse = "_")
    } else {
      "static"
    }
  }

  ### Helper to load rasters for a time combo, returning NULL on failure
  .load_rasters_for_combo <- function(combo, warn_prefix = "") {
    time_values  <- if (has_time && length(combo) > 0) {
      as.data.frame(combo, stringsAsFactors = FALSE)
    } else {
      as.data.frame(stats::setNames(lapply(time_cols, function(tc) NA_character_), time_cols),
                    stringsAsFactors = FALSE)
    }
    raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                          time_cols, time_values, raster_dir)
    if (is.null(raster_paths)) return(NULL)
    tryCatch(
      terra::rast(raster_paths),
      error = function(e) {
        warning(paste0(warn_prefix, "Could not load rasters: ", e$message))
        NULL
      }
    )
  }

  env_candidate_coords_by_year <- list()
  env_var_names                <- NULL

  if (method == "environmental") {

    if (verbose) message("Profiling environmental space from presence points...")

    unique_time_steps <- if (has_time) {
      pts_df_all <- sf::st_drop_geometry(pts_sf)
      sub        <- unique(pts_df_all[, time_cols, drop = FALSE])
      sub        <- sub[do.call(order, sub), , drop = FALSE]
      lapply(seq_len(nrow(sub)), function(i) as.list(sub[i, , drop = FALSE]))
    } else list(list())

    lower_bounds_by_timestep <- list()
    upper_bounds_by_timestep <- list()

    for (combo in unique_time_steps) {
      ts_key <- .make_ts_key(combo)
      r_ts   <- .load_rasters_for_combo(combo, warn_prefix = paste0("Time step ", ts_key, ": "))
      if (is.null(r_ts)) {
        warning(paste0("Could not load rasters for time step ", ts_key,
                       ". Skipping environmental profiling for this time step."))
        next
      }

      ts_pts <- if (has_time) {
        pts_df <- sf::st_drop_geometry(pts_sf)
        pts_sf[Reduce(`&`, lapply(time_cols, function(tc) !is.na(pts_df[[tc]]) & pts_df[[tc]] == combo[[tc]])), ]
      } else pts_sf

      if (nrow(ts_pts) == 0) next

      presence_vect    <- suppressWarnings(terra::vect(ts_pts))
      env_at_presences <- terra::extract(r_ts, presence_vect, ID = FALSE)
      env_at_presences <- env_at_presences[stats::complete.cases(env_at_presences), , drop = FALSE]

      if (nrow(env_at_presences) == 0) next
      if (is.null(env_var_names)) env_var_names <- names(env_at_presences)

      lower_bounds_by_timestep[[ts_key]] <- vapply(env_var_names, function(v) {
        stats::quantile(env_at_presences[[v]], probs = env_percentile, na.rm = TRUE)
      }, numeric(1))
      upper_bounds_by_timestep[[ts_key]] <- vapply(env_var_names, function(v) {
        stats::quantile(env_at_presences[[v]], probs = 1 - env_percentile, na.rm = TRUE)
      }, numeric(1))

      ### Candidate cell identification done in same pass as profiling
      lower_bounds <- lower_bounds_by_timestep[[ts_key]]
      upper_bounds <- upper_bounds_by_timestep[[ts_key]]

      all_cell_idx <- seq_len(terra::ncell(r_ts))
      all_env_vals <- as.data.frame(terra::values(r_ts))
      names(all_env_vals) <- env_var_names
      cell_xy      <- terra::xyFromCell(r_ts, all_cell_idx)
      ref_vect     <- terra::vect(ref_union)
      cell_vect    <- suppressWarnings(
        terra::vect(as.data.frame(cell_xy), geom = c("x", "y"),
                    crs = terra::crs(r_ts))
      )
      inside_ref <- terra::is.related(cell_vect, ref_vect, "intersects")

      outside_envelope <- rep(FALSE, nrow(all_env_vals))
      for (v in env_var_names) {
        v_vals           <- all_env_vals[, v]
        outside_envelope <- outside_envelope |
          (!is.na(v_vals) & (v_vals < lower_bounds[v] | v_vals > upper_bounds[v]))
      }

      candidate_mask   <- inside_ref & outside_envelope &
        !apply(is.na(all_env_vals), 1, any)
      coords_candidate <- cell_xy[candidate_mask, , drop = FALSE]
      valid_coords     <- !apply(is.na(coords_candidate), 1, any)
      env_candidate_coords_by_year[[ts_key]] <- coords_candidate[valid_coords, , drop = FALSE]

      if (verbose) message(paste0("  ", ts_key, ": ",
                                  nrow(env_candidate_coords_by_year[[ts_key]]),
                                  " candidate cells outside tolerance envelope."))
    }

    if (is.null(env_var_names)) {
      stop("ERROR: No environmental values could be extracted at presence points. ",
           "Check that raster_dir contains files matching variable_patterns and ",
           "that the rasters overlap with the presence locations.")
    }

    if (verbose) message(paste0("  Tolerance envelope: ", env_percentile * 100, "th to ",
                                (1 - env_percentile) * 100, "th percentiles across ",
                                length(env_var_names), " variables."))
    if (use_buffer_with_env) {
      if (verbose) message(paste0("  Spatial buffer constraint: ", buffer_distance, " units."))
    }
  }

  if (method == "user_data") {

    if (verbose) message("Processing user-supplied absence data...")

    ### Resolve input format, build ud_df and capture crs in one pass
    if (inherits(user_absence_data, "sf")) {
      ud_df  <- sf::st_drop_geometry(user_absence_data)
      ud_crs <- sf::st_crs(user_absence_data)
    } else if (inherits(user_absence_data, "SpatialPointsDataFrame")) {
      tmp_sf <- sf::st_as_sf(user_absence_data)
      ud_df  <- sf::st_drop_geometry(tmp_sf)
      ud_crs <- sf::st_crs(tmp_sf)
    } else if (is.character(user_absence_data) && file.exists(user_absence_data)) {
      file_ext_lower <- tolower(tools::file_ext(user_absence_data))
      if (file_ext_lower == "csv") {
        if (is.null(xcol))       stop("ERROR: 'xcol' is required when 'user_absence_data' is a .csv file.")
        if (is.null(ycol))       stop("ERROR: 'ycol' is required when 'user_absence_data' is a .csv file.")
        if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when 'user_absence_data' is a .csv file.")
        ud_df  <- utils::read.csv(user_absence_data, stringsAsFactors = FALSE)
        ud_crs <- sf::st_crs(points_crs)
      } else if (file_ext_lower %in% c("shp", "geojson", "gpkg")) {
        tmp_sf <- sf::st_read(user_absence_data, quiet = TRUE)
        ud_df  <- sf::st_drop_geometry(tmp_sf)
        ud_crs <- sf::st_crs(tmp_sf)
      } else {
        stop(paste0("ERROR: Unsupported file format '.", file_ext_lower,
                    "'. 'user_absence_data' must be a .csv, .shp, .geojson, or .gpkg file."))
      }
    } else if (is.data.frame(user_absence_data)) {
      if (is.null(xcol))       stop("ERROR: 'xcol' is required when 'user_absence_data' is a data frame.")
      if (is.null(ycol))       stop("ERROR: 'ycol' is required when 'user_absence_data' is a data frame.")
      if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when 'user_absence_data' is a data frame.")
      ud_df  <- user_absence_data
      ud_crs <- sf::st_crs(points_crs)
    } else {
      stop("ERROR: 'user_absence_data' must be an sf object, SpatialPointsDataFrame, data frame, or file path to a .csv, .shp, .geojson, or .gpkg file.")
    }

    ### Resolve coordinate columns once after format detection
    xcol_ud <- if (!is.null(xcol)) xcol else if ("X" %in% names(ud_df)) "X" else if ("x" %in% names(ud_df)) "x" else NULL
    ycol_ud <- if (!is.null(ycol)) ycol else if ("Y" %in% names(ud_df)) "Y" else if ("y" %in% names(ud_df)) "y" else NULL

    if (is.null(xcol_ud) || is.null(ycol_ud)) {
      stop("ERROR: Could not determine coordinate columns for 'user_absence_data'. ",
           "Supply 'xcol' and 'ycol' explicitly.")
    }
    if (!xcol_ud %in% names(ud_df)) stop(paste0("ERROR: Column '", xcol_ud, "' not found in 'user_absence_data'."))
    if (!ycol_ud %in% names(ud_df)) stop(paste0("ERROR: Column '", ycol_ud, "' not found in 'user_absence_data'."))

    if (has_time && !all(time_cols %in% names(ud_df))) {
      missing_tc <- time_cols[!time_cols %in% names(ud_df)]
      stop(paste0("ERROR: 'user_absence_data' is missing time column(s): ",
                  paste(missing_tc, collapse = ", "), "."))
    }

    ud_sf <- sf::st_as_sf(ud_df, coords = c(xcol_ud, ycol_ud),
                          crs = ud_crs, remove = FALSE)

    if (!identical(sf::st_crs(ud_sf), sf::st_crs(pts_sf))) {
      ud_sf <- sf::st_transform(ud_sf, crs = sf::st_crs(pts_sf))
    }

    voronoi_sf    <- partition_result$voronoi_folds
    pts_df_full   <- sf::st_drop_geometry(pts_sf)

    ud_sf$fold           <- NA_integer_
    ud_sf$presence       <- 0L
    ud_sf$temporal_block <- NA_integer_

    if (!is.null(voronoi_sf) && "fold" %in% names(voronoi_sf)) {

      sentinel   <- max(all_folds) + 1L
      voronoi_t  <- sf::st_transform(voronoi_sf, crs = sf::st_crs(ud_sf))
      spatial_polys <- voronoi_t[!is.na(voronoi_t$fold) & voronoi_t$fold != sentinel, ]

      if (nrow(spatial_polys) > 0) {
        spatial_polys$voronoi_fold <- spatial_polys$fold
        spatial_polys <- spatial_polys[, "voronoi_fold"]

        ud_sf$.row_id <- seq_len(nrow(ud_sf))
        ud_joined     <- suppressMessages(suppressWarnings(
          sf::st_join(ud_sf, spatial_polys, left = TRUE)
        ))
        ud_joined     <- ud_joined[!duplicated(ud_joined$.row_id), ]
        joined_folds  <- ud_joined$voronoi_fold[match(seq_len(nrow(ud_sf)), ud_joined$.row_id)]
        ud_sf$fold    <- joined_folds
        ud_sf$.row_id <- NULL
      }
    }

    unassigned_idx <- which(is.na(ud_sf$fold))

    if (length(unassigned_idx) > 0 && has_time) {
      spatial_fold_ids  <- sort(unique(ud_sf$fold[!is.na(ud_sf$fold)]))
      temporal_fold_ids <- sort(setdiff(all_folds, spatial_fold_ids))

      if (length(temporal_fold_ids) > 0) {
        time_col1 <- time_cols[1]
        tf_ranges <- lapply(temporal_fold_ids, function(f) {
          tv <- pts_df_full[[time_col1]][!is.na(pts_df_full$fold) & pts_df_full$fold == f]
          tv <- tv[!is.na(tv)]
          if (length(tv) == 0) return(NULL)
          list(fold = f, min = min(tv), max = max(tv))
        })
        tf_ranges <- Filter(Negate(is.null), tf_ranges)

        if (length(tf_ranges) > 0) {
          ud_times    <- ud_df[[time_col1]][unassigned_idx]
          assigned_tf <- vapply(ud_times, function(t) {
            if (is.na(t)) return(NA_integer_)
            for (r in tf_ranges) if (t >= r$min && t <= r$max) return(r$fold)
            dists <- vapply(tf_ranges, function(r) min(abs(t - r$min), abs(t - r$max)), numeric(1))
            tf_ranges[[which.min(dists)]]$fold
          }, integer(1))
          ud_sf$fold[unassigned_idx] <- assigned_tf
        }
      }
    }

    n_unassigned <- sum(is.na(ud_sf$fold))
    if (n_unassigned > 0) {
      warning(paste0(n_unassigned, " point(s) in 'user_absence_data' could not be ",
                     "assigned to any fold and were dropped."))
      ud_sf <- ud_sf[!is.na(ud_sf$fold), ]
    }
    if (nrow(ud_sf) == 0) {
      stop("ERROR: No 'user_absence_data' points could be assigned to any fold. ",
           "Check that CRS, spatial extent, and time values match the partition.")
    }

    pseudoabs_sf <- ud_sf

    if (verbose) {
      message(paste0("  ", nrow(pseudoabs_sf),
                     " user-supplied absence points assigned across folds."))
      for (f in all_folds) {
        n_f <- sum(!is.na(pseudoabs_sf$fold) & pseudoabs_sf$fold == f)
        message(paste0("  Fold ", f, ": ", n_f, " absence points"))
      }
    }

  } else {

    if (verbose) message(paste0("Generating pseudoabsence points across ",
                                length(all_folds), " fold",
                                if (length(all_folds) > 1) "s" else "", "..."))

    result_list <- list()

    for (fold_id in all_folds) {

      fold_pts <- pts_sf[!is.na(pts_sf$fold) & pts_sf$fold == fold_id, ]
      if (nrow(fold_pts) == 0) next

      n_presences_fold <- nrow(fold_pts)
      if (ratio == 0) {
        n_pseudoabs_fold <- NA_integer_
      } else {
        n_pseudoabs_fold <- max(1, round(n_presences_fold * ratio))
      }

      if (has_time) {
        fold_df        <- sf::st_drop_geometry(fold_pts)
        sub_ft         <- unique(fold_df[, time_cols, drop = FALSE])
        sub_ft         <- sub_ft[do.call(order, sub_ft), , drop = FALSE]
        presence_times <- lapply(seq_len(nrow(sub_ft)),
                                 function(i) as.list(sub_ft[i, , drop = FALSE]))
      } else {
        presence_times <- list(list())
      }

      unique_times <- if (!is.null(pseudoabsence_times) && has_time) {
        other_cols <- time_cols[-1]
        if (length(other_cols) > 0) {
          all_df <- sf::st_drop_geometry(pts_sf)
          other_combos <- unique(all_df[, other_cols, drop = FALSE])
          cross <- expand.grid(
            c(list(pseudoabsence_times), as.list(other_combos)),
            stringsAsFactors = FALSE
          )
          names(cross) <- time_cols
          lapply(seq_len(nrow(cross)), function(i) as.list(cross[i, , drop = FALSE]))
        } else {
          lapply(pseudoabsence_times, function(v) stats::setNames(list(v), time_cols[1]))
        }
      } else {
        presence_times
      }

      if (length(unique_times) == 0 || (length(unique_times) == 1 && length(unique_times[[1]]) == 0)) {
        if (ratio == 0) {
          time_alloc <- as.integer(min_points_per_timestep)
        } else {
          time_alloc <- n_pseudoabs_fold
        }
        names(time_alloc) <- "all"
      } else {
        ts_keys <- vapply(unique_times, .make_ts_key, character(1))

        time_counts <- vapply(unique_times, function(combo) {
          if (!has_time || length(presence_times) == 0) return(0)
          fold_df <- sf::st_drop_geometry(fold_pts)
          sum(Reduce(`&`, lapply(time_cols, function(tc) !is.na(fold_df[[tc]]) & fold_df[[tc]] == combo[[tc]])))
        }, integer(1))
        names(time_counts) <- ts_keys

        if (ratio == 0) {
          time_alloc <- rep(as.integer(min_points_per_timestep), length(unique_times))
        } else {
          raw_alloc  <- (time_counts / max(sum(time_counts), 1)) * n_pseudoabs_fold
          time_alloc <- pmax(as.integer(round(raw_alloc)), as.integer(min_points_per_timestep))

          diff_pts <- sum(time_alloc) - n_pseudoabs_fold
          if (diff_pts > 0) {
            ord <- order(time_alloc, decreasing = TRUE)
            for (i in seq_len(min(diff_pts, length(ord)))) {
              if (time_alloc[ord[i]] > min_points_per_timestep) {
                time_alloc[ord[i]] <- time_alloc[ord[i]] - 1
              }
            }
          } else if (diff_pts < 0) {
            ord <- order(time_alloc, decreasing = FALSE)
            for (i in seq_len(min(abs(diff_pts), length(ord)))) {
              time_alloc[ord[i]] <- time_alloc[ord[i]] + 1
            }
          }
        }
        names(time_alloc) <- ts_keys
      }

      fold_buffer_region <- NULL
      if (method == "buffer" || use_buffer_with_env) {
        fold_buffer_region <- suppressWarnings(
          sf::st_intersection(
            sf::st_buffer(sf::st_union(fold_pts), dist = buffer_distance),
            ref_union
          )
        )
      }

      for (ti in seq_along(unique_times)) {

        combo  <- unique_times[[ti]]
        ts_key <- names(time_alloc)[ti]
        n_pts  <- time_alloc[ti]
        if (n_pts <= 0) next

        ts_pts <- if (has_time && length(combo) > 0) {
          fold_df <- sf::st_drop_geometry(fold_pts)
          fold_pts[Reduce(`&`, lapply(time_cols, function(tc) !is.na(fold_df[[tc]]) & fold_df[[tc]] == combo[[tc]])), ]
        } else fold_pts

        excl_buffer <- if (nrow(ts_pts) > 0) {
          suppressWarnings(sf::st_buffer(sf::st_union(ts_pts), dist = 1e-6))
        } else NULL

        if (method == "random") {
          sampling_region <- if (!is.null(excl_buffer)) {
            suppressWarnings(sf::st_difference(ref_union, excl_buffer))
          } else ref_union

          sampled <- tryCatch({
            s <- sf::st_sample(sampling_region, size = n_pts * 5, type = "random")
            s <- s[!sf::st_is_empty(s)]
            if (length(s) > n_pts) s[seq_len(n_pts)] else s
          }, error = function(e) {
            warning(paste0("Fold ", fold_id, " / time ", ts_key,
                           ": random sampling failed -- ", e$message))
            NULL
          })
        }

        if (method == "buffer") {
          sampling_region <- if (!is.null(fold_buffer_region) &&
                                 !sf::st_is_empty(fold_buffer_region)) {
            fold_buffer_region
          } else {
            warning(paste0("Fold ", fold_id, " / time ", ts_key,
                           ": buffer region is empty -- falling back to full study area."))
            ref_union
          }
          if (!is.null(excl_buffer)) {
            sampling_region <- suppressWarnings(sf::st_difference(sampling_region, excl_buffer))
          }
          sampled <- tryCatch({
            s <- sf::st_sample(sampling_region, size = n_pts * 5, type = "random")
            s <- s[!sf::st_is_empty(s)]
            if (length(s) > n_pts) s[seq_len(n_pts)] else s
          }, error = function(e) {
            warning(paste0("Fold ", fold_id, " / time ", ts_key,
                           ": buffer sampling failed -- ", e$message))
            NULL
          })
        }

        if (method == "environmental") {
          env_ts_key <- if (has_time && length(env_candidate_coords_by_year) > 0) {
            if (ts_key %in% names(env_candidate_coords_by_year)) {
              ts_key
            } else if ("static" %in% names(env_candidate_coords_by_year)) {
              "static"
            } else names(env_candidate_coords_by_year)[1]
          } else names(env_candidate_coords_by_year)[1]

          candidate_coords <- env_candidate_coords_by_year[[env_ts_key]]

          if (is.null(candidate_coords) || nrow(candidate_coords) == 0) {
            warning(paste0("Fold ", fold_id, " / time ", ts_key,
                           ": no candidate cells outside tolerance envelope. Skipping."))
            sampled <- NULL
          } else {
            if (use_buffer_with_env && !is.null(fold_buffer_region) &&
                !sf::st_is_empty(fold_buffer_region)) {
              cand_sf    <- sf::st_as_sf(as.data.frame(candidate_coords),
                                         coords = c("x", "y"), crs = sf::st_crs(pts_sf))
              inside_buf <- suppressMessages(suppressWarnings(
                lengths(sf::st_intersects(cand_sf, fold_buffer_region)) > 0
              ))
              candidate_coords <- candidate_coords[inside_buf, , drop = FALSE]
            }

            if (nrow(candidate_coords) == 0) {
              warning(paste0("Fold ", fold_id, " / time ", ts_key,
                             ": no candidate cells within buffer region. Skipping."))
              sampled <- NULL
            } else {
              n_clusters <- min(n_pts, nrow(candidate_coords))
              km <- if (n_clusters < nrow(candidate_coords)) {
                suppressWarnings(tryCatch(
                  stats::kmeans(candidate_coords, centers = n_clusters,
                                nstart = 5, iter.max = 300, algorithm = "Lloyd"),
                  error = function(e) {
                    warning(paste0("Fold ", fold_id, " / time ", ts_key,
                                   ": k-means failed -- ", e$message,
                                   ". Falling back to random sample."))
                    NULL
                  }
                ))
              } else NULL

              selected_coords <- if (!is.null(km)) {
                km$centers
              } else {
                idx <- if (nrow(candidate_coords) > n_clusters) {
                  sample(nrow(candidate_coords), n_clusters)
                } else seq_len(nrow(candidate_coords))
                candidate_coords[idx, , drop = FALSE]
              }

              sampled <- tryCatch({
                sf::st_as_sfc(
                  lapply(seq_len(nrow(selected_coords)), function(i) {
                    sf::st_point(selected_coords[i, ])
                  }),
                  crs = sf::st_crs(pts_sf)
                )
              }, error = function(e) {
                warning(paste0("Fold ", fold_id, " / time ", ts_key,
                               ": could not build geometry from cluster centroids -- ",
                               e$message))
                NULL
              })
            }
          }
        }

        if (is.null(sampled) || length(sampled) == 0) next

        sampled_sf           <- sf::st_as_sf(sampled)
        names(sampled_sf)[1] <- "geometry"
        sf::st_geometry(sampled_sf) <- "geometry"
        sampled_sf$fold      <- fold_id
        sampled_sf$presence  <- 0

        sampled_sf$temporal_block <- if ("temporal_block" %in% names(fold_pts)) {
          if (has_time && length(combo) > 0) {
            fold_df <- sf::st_drop_geometry(fold_pts)
            tb_v    <- fold_df$temporal_block[Reduce(`&`, lapply(time_cols, function(tc) !is.na(fold_df[[tc]]) & fold_df[[tc]] == combo[[tc]]))]
            if (length(tb_v) > 0) tb_v[1] else NA_integer_
          } else NA_integer_
        } else NA_integer_

        if (has_time && length(combo) > 0) {
          for (tc in time_cols) {
            sampled_sf[[tc]] <- combo[[tc]]
          }
        }

        result_list[[length(result_list) + 1]] <- sampled_sf
      }
    }

    if (length(result_list) == 0) {
      stop("ERROR: No pseudoabsence points could be generated. ",
           "Check method parameters, buffer distance, and raster extent.")
    }

    all_cols    <- unique(unlist(lapply(result_list, names)))
    result_list <- lapply(result_list, function(x) {
      for (col in setdiff(all_cols, names(x))) x[[col]] <- NA
      x[, all_cols, drop = FALSE]
    })
    pseudoabs_sf <- do.call(rbind, result_list)

  }

  if (verbose) message("Extracting time-matched environmental values at pseudoabsence points...")

  if (has_time && all(time_cols %in% names(pseudoabs_sf))) {
    pa_df     <- sf::st_drop_geometry(pseudoabs_sf)
    sub_pa    <- unique(pa_df[, time_cols, drop = FALSE])
    sub_pa    <- sub_pa[do.call(order, sub_pa), , drop = FALSE]
    pa_combos <- lapply(seq_len(nrow(sub_pa)),
                        function(i) as.list(sub_pa[i, , drop = FALSE]))

    for (combo in pa_combos) {
      pa_idx <- which(Reduce(`&`, lapply(time_cols, function(tc) !is.na(pa_df[[tc]]) & pa_df[[tc]] == combo[[tc]])))
      if (length(pa_idx) == 0) next
      r_ts <- .load_rasters_for_combo(combo, warn_prefix = paste0("Extraction combo ", .make_ts_key(combo), ": "))
      if (is.null(r_ts)) next
      pa_sub <- suppressWarnings(terra::vect(pseudoabs_sf[pa_idx, ]))
      env_ex <- terra::extract(r_ts, pa_sub, ID = FALSE)
      var_order <- c(dynamic_vars, static_vars)
      if (ncol(env_ex) == length(var_order)) names(env_ex) <- var_order
      for (v in var_order) {
        if (!v %in% names(env_ex)) next
        if (!v %in% names(pseudoabs_sf)) pseudoabs_sf[[v]] <- NA_real_
        pseudoabs_sf[[v]][pa_idx] <- env_ex[[v]]
      }
    }
  } else {
    r_static <- .load_rasters_for_combo(list(), warn_prefix = "Static extraction: ")
    if (!is.null(r_static)) {
      pa_vect <- suppressWarnings(terra::vect(pseudoabs_sf))
      env_ex  <- terra::extract(r_static, pa_vect, ID = FALSE)
      var_order <- c(dynamic_vars, static_vars)
      if (ncol(env_ex) == length(var_order)) names(env_ex) <- var_order
      for (v in var_order) {
        if (!v %in% names(env_ex)) next
        pseudoabs_sf[[v]] <- env_ex[[v]]
      }
    }
  }

  if (verbose) message(paste0("  Done. ", nrow(pseudoabs_sf), " pseudoabsence points generated."))

  fold_summary <- do.call(rbind, lapply(all_folds, function(f) {
    n_pr <- sum(!is.na(pts_sf$fold) & pts_sf$fold == f)
    n_pa <- sum(!is.na(pseudoabs_sf$fold) & pseudoabs_sf$fold == f)
    data.frame(fold             = f,
               n_presences      = n_pr,
               n_pseudoabsences = n_pa,
               ratio_achieved   = round(n_pa / max(n_pr, 1), 3),
               stringsAsFactors = FALSE)
  }))

  method_label <- switch(method,
                         random        = "Random",
                         buffer        = paste0("Buffer (", buffer_distance, " units)"),
                         environmental = if (use_buffer_with_env) {
                           paste0("Environmental + Buffer (", buffer_distance, " units)")
                         } else "Environmental",
                         user_data     = "User-supplied data"
  )

  plot_list    <- list()
  presence_col <- "steelblue"
  absence_col  <- "tomato"

  if (create_plot) {

    fold_colors <- .resolve_palette(length(all_folds), plot_palette)
    names(fold_colors) <- as.character(all_folds)

    if (has_time && all(time_cols %in% names(pseudoabs_sf))) {
      pa_df    <- sf::st_drop_geometry(pseudoabs_sf)
      pr_df    <- sf::st_drop_geometry(pts_sf)
      pa_time  <- apply(pa_df[, time_cols, drop = FALSE], 1,
                        function(r) paste(r, collapse = "_"))
      pr_time  <- apply(pr_df[, time_cols, drop = FALSE], 1,
                        function(r) paste(r, collapse = "_"))
      all_time <- sort(unique(c(pa_time, pr_time)))
      all_time <- all_time[!is.na(all_time) & all_time != "NA"]
      time_label <- paste(time_cols, collapse = " / ")

      count_matrix <- rbind(
        Presence      = vapply(all_time, function(t)
          sum(!is.na(pr_time) & pr_time == t), integer(1)),
        Pseudoabsence = vapply(all_time, function(t)
          sum(!is.na(pa_time) & pa_time == t), integer(1))
      )
      colnames(count_matrix) <- as.character(all_time)

      label_interval <- max(1, round(length(all_time) / 15))
      display_labels <- rep("", length(all_time))
      display_labels[seq(1, length(all_time), by = label_interval)] <-
        as.character(all_time)[seq(1, length(all_time), by = label_interval)]

      opar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(opar), add = TRUE)
      tryCatch({
        graphics::barplot(count_matrix, beside = TRUE,
                          col = c(presence_col, absence_col),
                          names.arg = display_labels, las = 2, border = NA,
                          main = "Temporal Distribution",
                          xlab = time_label, ylab = "Count")
        graphics::legend("topright",
                         legend = c("Presence", "Pseudoabsence"),
                         fill = c(presence_col, absence_col),
                         border = NA, cex = 0.8, bty = "n")
        plot_list$temporal_distribution <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not generate temporal distribution plot: ", e$message))
      })
    }

    if (plot_by_fold) {
      for (fi in seq_along(all_folds)) {
        fold_id <- all_folds[fi]
        fcol    <- fold_colors[as.character(fold_id)]
        fold_pa <- pseudoabs_sf[!is.na(pseudoabs_sf$fold) & pseudoabs_sf$fold == fold_id, ]
        fold_pr <- pts_sf[!is.na(pts_sf$fold) & pts_sf$fold == fold_id, ]

        opar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(opar), add = TRUE)
        tryCatch({
          graphics::plot(sf::st_geometry(reference_shapefile),
                         col = "gray97", border = "gray50",
                         main = paste0("Fold ", fold_id, " Pseudoabsences"))
          if (nrow(fold_pa) > 0) {
            pa_coords <- sf::st_coordinates(fold_pa)
            graphics::points(pa_coords[, 1], pa_coords[, 2],
                             col = grDevices::adjustcolor(absence_col, alpha.f = 0.6),
                             pch = 4, cex = 0.7)
          }
          if (nrow(fold_pr) > 0) {
            pr_coords <- sf::st_coordinates(fold_pr)
            graphics::points(pr_coords[, 1], pr_coords[, 2],
                             col = fcol, pch = 16, cex = 0.8)
          }
          graphics::legend("topright",
                           legend = c(paste0("Presence (n = ", nrow(fold_pr), ")"),
                                      paste0("Pseudoabsence (n = ", nrow(fold_pa), ")")),
                           col = c(fcol, grDevices::adjustcolor(absence_col, alpha.f = 0.6)),
                           pch = c(16, 4), cex = 0.8, bty = "n")
          plot_list[[paste0("spatial_fold_", fold_id)]] <- grDevices::recordPlot()
        }, error = function(e) {
          warning(paste0("Could not generate plot for fold ", fold_id, ": ", e$message))
        })
      }
    } else {
      opar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(opar), add = TRUE)
      tryCatch({
        graphics::plot(sf::st_geometry(reference_shapefile),
                       col = "gray97", border = "gray50",
                       main = "Pseudoabsences")
        if (nrow(pseudoabs_sf) > 0) {
          pa_coords <- sf::st_coordinates(pseudoabs_sf)
          graphics::points(pa_coords[, 1], pa_coords[, 2],
                           col = grDevices::adjustcolor(absence_col, alpha.f = 0.5),
                           pch = 4, cex = 0.6)
        }
        for (fi in seq_along(all_folds)) {
          fold_id <- all_folds[fi]
          fold_pr <- pts_sf[!is.na(pts_sf$fold) & pts_sf$fold == fold_id, ]
          if (nrow(fold_pr) == 0) next
          pr_coords <- sf::st_coordinates(fold_pr)
          graphics::points(pr_coords[, 1], pr_coords[, 2],
                           col = fold_colors[as.character(fold_id)], pch = 16, cex = 0.8)
        }
        graphics::legend("topright",
                         legend = c(paste0("Presence - Fold ", all_folds),
                                    paste0("Pseudoabsence (n = ", nrow(pseudoabs_sf), ")")),
                         col = c(fold_colors,
                                 grDevices::adjustcolor(absence_col, alpha.f = 0.5)),
                         pch = c(rep(16, length(all_folds)), 4),
                         cex = 0.75, bty = "n")
        plot_list$spatial_combined <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste0("Could not generate spatial plot: ", e$message))
      })
    }
  }

  if (verbose) message("\nPseudoabsence generation complete.")
  if (verbose) message(paste0("  Method: ", method_label))
  if (method != "user_data") {
    if (verbose) message(paste0("  Ratio: ", ratio,
                                " | Presences: ", nrow(pts_sf),
                                " | Pseudoabsences: ", nrow(pseudoabs_sf)))
  } else {
    if (verbose) message(paste0("  Presences: ", nrow(pts_sf),
                                " | User absences: ", nrow(pseudoabs_sf)))
  }
  for (f in all_folds) {
    row <- fold_summary[fold_summary$fold == f, ]
    if (verbose) message(paste0("  Fold ", f, ": ",
                                row$n_presences, " presences, ",
                                row$n_pseudoabsences, " pseudoabsences"))
  }

  if (!is.null(output_file)) {
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(pseudoabsences = pseudoabs_sf,
                 plots          = plot_list,
                 summary        = fold_summary),
            output_file)
    if (verbose) message(paste0("Results saved to: ", output_file))
  }

  invisible(list(pseudoabsences = pseudoabs_sf,
                 plots          = plot_list,
                 summary        = fold_summary))
}
