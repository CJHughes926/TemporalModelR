#' Generate Temporally Explicit Pseudoabsence Points
#'
#' Generates pseudoabsence or background points for each fold produced by
#' \code{\link{spatiotemporal_partition}}, distributed across time steps
#' proportionally to the number of presence points in each time step within each
#' fold. Three generation methods are supported: random sampling within the
#' study area, buffer-constrained sampling around presence points, and
#' environmentally biased sampling that targets areas outside the known
#' environmental tolerance of the species.
#'
#' @usage
#' generate_absences(partition_result, reference_shapefile_path, raster_dir,
#'                   variable_patterns, method = "random", ratio = 1,
#'                   buffer_distance = NULL, env_percentile = 0.05,
#'                   time_cols = NULL, pseudoabsence_times = NULL,
#'                   min_points_per_timestep = 1, create_plot = TRUE,
#'                   plot_by_fold = FALSE, plot_palette = "Dark 2",
#'                   output_file = NULL, verbose = TRUE)
#'
#' @param partition_result List or character. Output from
#'   \code{\link{spatiotemporal_partition}} or path to an \code{.rds}
#'   file containing that output.
#' @param reference_shapefile_path Character or sf object. Path to a polygon
#'   file or an \code{sf} polygon object defining the study area.
#' @param raster_dir Character. Directory containing environmental raster
#'   files (\code{.tif}), typically the output of
#'   \code{\link{raster_align}} or \code{\link{scale_rasters}}. File names
#'   must follow the patterns supplied in \code{variable_patterns}, with any
#'   time placeholder substituted for the corresponding value from
#'   \code{time_cols}. Required for all methods.
#' @param variable_patterns Named character vector mapping clean variable names
#'   to raster filename patterns. For time-varying variables include the time
#'   placeholder in the pattern (e.g. \code{"forest_cover" = "forest_cover_YEAR"});
#'   for static variables omit it (e.g. \code{"elevation" = "elevation"}). Time
#'   placeholders must match entries in \code{time_cols}.
#' @param method Character. Pseudoabsence generation method. One of
#'   \code{"random"}, \code{"buffer"}, or \code{"environmental"}.
#'   Default is \code{"random"}.
#' @param ratio Numeric. Number of pseudoabsence points to generate per
#'   presence point. Default is \code{1}. Values of 2, 10, 50, etc. are
#'   accepted. Points are always distributed proportionally across time steps
#'   within each fold. Set to \code{0} to disable proportional allocation and
#'   use a fixed number of points per time step instead, in which case
#'   \code{min_points_per_timestep} must be greater than 0. \code{ratio} and
#'   \code{min_points_per_timestep} cannot both be 0.
#' @param buffer_distance Numeric. Distance in the units of the CRS (typically
#'   meters for projected CRS) within which pseudoabsence points are sampled.
#'   Required when \code{method = "buffer"}. When \code{method =
#'   "environmental"}, supplying a value automatically applies a spatial buffer
#'   constraint before environmental profiling, following the three-step approach
#'   of Senay et al. (2013). If \code{NULL} for the environmental method, no
#'   spatial constraint is applied. Default is \code{NULL}.
#' @param env_percentile Numeric between 0 and 1. Quantile threshold used to
#'   define the boundary of the known environmental tolerance when
#'   \code{method = "environmental"}. Environmental cells within this quantile
#'   range across all variables are excluded from pseudoabsence sampling.
#'   Default is \code{0.05} (5th to 95th percentile envelope).
#' @param time_cols Character or character vector. Name of the column(s)
#'   containing the time step values. Must match \code{time_cols} used in
#'   \code{\link{spatiotemporal_partition}} and the time placeholders used
#'   in \code{variable_patterns}. Default is \code{NULL}.
#' @param pseudoabsence_times Vector. Optional vector of specific time step
#'   values (for the first time column) at which to generate pseudoabsences.
#'   When \code{NULL} (default), all time steps present in the occurrence data
#'   are used.
#' @param min_points_per_timestep Integer. Minimum number of pseudoabsence
#'   points to generate per time step per fold. Default is \code{1}. When
#'   \code{ratio = 0}, this value sets the exact (fixed) number of points
#'   generated per time step per fold, independent of the number of
#'   presence points. \code{ratio} and \code{min_points_per_timestep} cannot
#'   both be 0.
#' @param create_plot Logical. If \code{TRUE} (default), generates diagnostic
#'   plots showing the spatial and temporal distribution of generated
#'   pseudoabsence points alongside presence points.
#' @param plot_by_fold Logical. If \code{TRUE}, generates one map per fold. If
#'   \code{FALSE} (default), generates a single combined map.
#' @param plot_palette Character. Name of an HCL or RColorBrewer palette used
#'   to color folds in diagnostic plots. Accepts any HCL palette name (see
#'   \code{\link[grDevices]{hcl.pals}}) or, if \pkg{RColorBrewer} is installed,
#'   any Brewer palette name. Default is \code{"Dark 2"}.
#' @param output_file Character. Optional path to save the result as an
#'   \code{.rds} file. The parent directory will be created if it does not
#'   exist. Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes per-fold and per-time-step
#'   pseudoabsence counts.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{pseudoabsences}: An sf object of all generated pseudoabsence
#'     points with columns \code{fold}, \code{temporal_block}, \code{presence}
#'     (always 0), the time column(s) if provided, and extracted environmental
#'     variable values matched to each point's time step.
#'   \item \code{plots}: A named list of recorded plot objects when
#'     \code{create_plot = TRUE}. Contains \code{temporal_distribution} and
#'     either \code{spatial_combined} or one \code{spatial_fold_N} entry per
#'     fold. Plots can be replayed with \code{grDevices::replayPlot()}.
#'   \item \code{summary}: A data frame summarising points generated per fold
#'     with columns \code{fold}, \code{n_presences}, \code{n_pseudoabsences},
#'     and \code{ratio_achieved}.
#' }
#'
#' @details
#' Generates sets of background data based on user-specified methodology that
#' can be used as pseudoabsence data for the purposes of training
#' presence/absence models.
#'
#' The three generation methods differ in how the sampling region is defined:
#' \itemize{
#'   \item \strong{Random}: Points are sampled uniformly at random from the
#'     full study area, excluding a negligible buffer around presence locations
#'     to prevent exact overlap.
#'   \item \strong{Buffer}: Points are sampled within buffers of radius
#'     \code{buffer_distance} drawn around all fold presences, clipped to the
#'     reference shapefile boundary.
#'   \item \strong{Environmental}: Raster cells whose values fall outside the
#'     species tolerance envelope in at least one variable are identified as
#'     candidates. K-means clustering then selects a spatially
#'     representative subset. If \code{buffer_distance} is supplied the
#'     environmental filtering is applied only within that buffered region,
#'     implementing the full three-step approach of Senay et al. (2013).
#' }
#'
#' @references
#' Senay SD, Worner SP, Ikeda T (2013) Novel Three-Step Pseudo-Absence
#' Selection Technique for Improved Species Distribution Modeling.
#' PLoS ONE 8(8): e71218.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_partition}},
#'   \code{\link{temporally_explicit_extraction}}
#'
#' Modeling: \code{\link{build_temporal_hv}},
#'   \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_rf}}
#'
#' @examples
#' data(tmr_partition, package = "TemporalModelR")
#'
#' tmr_partition$points_sf <- do.call(rbind, lapply(
#'   unique(tmr_partition$points_sf$fold),
#'   function(f) {
#'     rows <- tmr_partition$points_sf[tmr_partition$points_sf$fold == f, ]
#'     rows[sample(nrow(rows), ceiling(nrow(rows) / 2)), ]
#'   }
#' ))
#'
#' scl_dir   <- system.file("extdata/rasters_scaled",
#'                          package = "TemporalModelR")
#'
#' ref_file  <- system.file("extdata/rasters_raw/elevation.tif",
#'                          package = "TemporalModelR")
#'
#' study_crs <- sf::st_crs(terra::rast(ref_file))
#'
#' study_area_sf <- sf::st_as_sf(sf::st_as_sfc(
#'   sf::st_bbox(c(xmin = 0, xmax = 3000, ymin = 0, ymax = 1500),
#'               crs = study_crs)
#' ))
#'
#' generate_absences(
#'   partition_result         = tmr_partition,
#'   reference_shapefile_path = study_area_sf,
#'   raster_dir               = scl_dir,
#'   variable_patterns        = c(
#'     "elevation"    = "elevation",
#'     "forest_cover" = "forest_cover_YEAR",
#'     "prseas"       = "prseas_YEAR_SEASON"
#'   ),
#'   method                   = "random",
#'   ratio                    = 1,
#'   time_cols                = c("year", "season"),
#'   create_plot              = FALSE,
#'   verbose                  = FALSE
#' )

#' @export
#' @importFrom sf st_as_sf st_as_sfc st_buffer st_coordinates st_crs
#' @importFrom sf st_difference st_drop_geometry st_geometry st_intersection
#' @importFrom sf st_intersects st_is_empty st_point st_read st_sample st_sf
#' @importFrom sf st_transform st_union
#' @importFrom terra rast extract values vect ncell crs xyFromCell is.related
#' @importFrom graphics barplot legend par plot points
#' @importFrom grDevices adjustcolor recordPlot
#' @importFrom stats complete.cases kmeans quantile setNames
#' @importFrom tools file_ext
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

  method <- match.arg(method, c("random", "buffer", "environmental"))

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

  env_candidate_coords_by_year <- list()
  env_var_names                <- NULL

  if (method == "environmental") {

    if (verbose) message("Profiling environmental space from presence points...")

    unique_time_steps <- if (has_time) {
      {
        pts_df_all <- sf::st_drop_geometry(pts_sf)
        sub        <- unique(pts_df_all[, time_cols, drop = FALSE])
        sub        <- sub[do.call(order, sub), , drop = FALSE]
        lapply(seq_len(nrow(sub)), function(i) as.list(sub[i, , drop = FALSE]))
      }
    } else list(list())

    lower_bounds_by_timestep <- list()
    upper_bounds_by_timestep <- list()

    for (combo in unique_time_steps) {
      ts_key      <- if (has_time) paste(vapply(time_cols, function(tc) as.character(combo[[tc]]), character(1)), collapse = "_") else "static"
      time_values <- if (has_time) as.data.frame(combo, stringsAsFactors = FALSE) else
        as.data.frame(stats::setNames(lapply(time_cols, function(tc) NA_character_), time_cols),
                      stringsAsFactors = FALSE)
      raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                            time_cols, time_values, raster_dir)
      r_ts <- if (!is.null(raster_paths)) tryCatch(
        terra::rast(raster_paths),
        error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL }
      ) else NULL
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

    for (combo in unique_time_steps) {
      ts_key       <- if (has_time) paste(vapply(time_cols, function(tc) as.character(combo[[tc]]), character(1)), collapse = "_") else "static"
      lower_bounds <- lower_bounds_by_timestep[[ts_key]]
      upper_bounds <- upper_bounds_by_timestep[[ts_key]]
      if (is.null(lower_bounds)) next

      time_values  <- if (has_time) as.data.frame(combo, stringsAsFactors = FALSE) else
        as.data.frame(stats::setNames(lapply(time_cols, function(tc) NA_character_), time_cols),
                      stringsAsFactors = FALSE)
      raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                            time_cols, time_values, raster_dir)
      r_ts <- if (!is.null(raster_paths)) tryCatch(
        terra::rast(raster_paths),
        error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL }
      ) else NULL
      if (is.null(r_ts)) next

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
  }

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
      ts_keys <- vapply(unique_times, function(combo) {
        if (has_time) paste(vapply(time_cols, function(tc) as.character(combo[[tc]]), character(1)), collapse = "_") else "all"
      }, character(1))

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

  if (verbose) message("Extracting time-matched environmental values at pseudoabsence points...")

  if (has_time && all(time_cols %in% names(pseudoabs_sf))) {
    pa_df       <- sf::st_drop_geometry(pseudoabs_sf)
    sub_pa    <- unique(pa_df[, time_cols, drop = FALSE])
    sub_pa    <- sub_pa[do.call(order, sub_pa), , drop = FALSE]
    pa_combos <- lapply(seq_len(nrow(sub_pa)),
                        function(i) as.list(sub_pa[i, , drop = FALSE]))

    for (combo in pa_combos) {
      pa_idx  <- which(Reduce(`&`, lapply(time_cols, function(tc) !is.na(pa_df[[tc]]) & pa_df[[tc]] == combo[[tc]])))
      if (length(pa_idx) == 0) next
      raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                            time_cols,
                                            as.data.frame(combo, stringsAsFactors = FALSE),
                                            raster_dir)
      r_ts <- if (!is.null(raster_paths)) tryCatch(
        terra::rast(raster_paths),
        error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL }
      ) else NULL
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
    static_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, static_vars,
                                          time_cols,
                                          as.data.frame(stats::setNames(
                                            lapply(time_cols, function(tc) NA_character_),
                                            time_cols), stringsAsFactors = FALSE),
                                          raster_dir)
    r_static <- if (!is.null(static_paths)) tryCatch(
      terra::rast(static_paths),
      error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL }
    ) else NULL
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
                         } else "Environmental"
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
  if (verbose) message(paste0("  Ratio: ", ratio,
                              " | Presences: ", nrow(pts_sf),
                              " | Pseudoabsences: ", nrow(pseudoabs_sf)))
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
