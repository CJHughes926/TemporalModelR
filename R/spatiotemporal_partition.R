#' Spatiotemporal Cross-Validation Partitioning
#'
#' Preprocesses species occurrence data by partitioning it into spatially and
#' temporally structured folds for cross-validation. Supports creation of
#' spatial-only folds, temporal-only folds, and random folds.
#'
#' Works better with smaller numbers of folds and may have difficulties
#' creating even folds for large numbers of groups or where sample sizes
#' are very small.
#'
#' @param reference_shapefile_path Character or sf object. Path to a polygon
#'   file or an \code{sf} polygon object defining the study area.
#' @param points_file_path Character, sf object, sfc object, Spatial object,
#'   or data frame. Path to occurrence data (\code{.csv}, \code{.shp},
#'   \code{.geojson}, \code{.gpkg}) or a spatial object.
#' @param time_cols Character. Name of a single column containing temporal
#'   values (e.g. year). Used to define temporal blocks. Required when using
#'   temporal folds. Must be a single column name; does not support more than one
#'   time column unlike other functions in this package. Compound time representations
#'   (e.g. year + season) should be encoded into a single ordered numeric column
#'   before partitioning, or only one (e.g. year) should be used.
#' @param xcol Character. Name of the x-coordinate column. Required when
#'   \code{points_file_path} is a CSV file or data frame.
#' @param ycol Character. Name of the y-coordinate column. Required when
#'   \code{points_file_path} is a CSV file or data frame.
#' @param points_crs Character or CRS object. CRS of the input points.
#'   Required when \code{points_file_path} is a CSV file or data frame.
#' @param n_spatial_folds Integer. Number of spatially explicit folds. Ignored
#'   when using random folds. Default is \code{0}.
#' @param n_temporal_folds Integer. Number of temporally explicit folds.
#'   When used alone (with \code{n_spatial_folds = 0}), creates temporal-only
#'   folds where each fold spans the full study area but covers a distinct
#'   slice of the time series. When combined with \code{n_spatial_folds},
#'   creates a spatiotemporal design. Ignored when using random folds.
#'   Default is \code{0}.
#' @param n_balanced_folds Integer. Reserved for future use. Default is
#'   \code{0} (disabled).
#' @param n_random_folds Integer. Number of random folds with no spatial or
#'   temporal structure. Overrides all other fold parameters. Default is
#'   \code{0}.
#' @param single_fold Logical. If \code{TRUE}, bypasses all partitioning and
#'   assigns all points to a single fold (fold 1). In this mode all points are
#'   used for both training and testing, producing a single model trained on
#'   the full dataset. All downstream functions accept the result identically
#'   to a standard multi-fold partition. Overrides all fold count parameters.
#'   Default is \code{FALSE}.
#' @param max_imbalance Numeric. Maximum allowed fold size imbalance as a
#'   proportion between 0 and 1. Default is \code{0.05}.
#' @param max_attempts Integer. Maximum number of partitioning attempts for
#'   spatiotemporal and balanced modes. Each attempt re-runs the spatial block
#'   construction; the attempt with the lowest imbalance is returned. Ignored
#'   for random and spatial-only modes. Default is \code{10}.
#' @param create_plot Logical. If \code{TRUE} (default), generates
#'   diagnostic plots showing fold distributions.
#' @param plot_palette Character. Name of an HCL or RColorBrewer palette used
#'   to color folds in diagnostic plots. Accepts any HCL palette name (see
#'   \code{\link[grDevices]{hcl.pals}}) or, if \pkg{RColorBrewer} is installed,
#'   any Brewer palette name. Default is \code{"Dark 2"}.
#' @param output_file Character. Optional path to save the result as an
#'   \code{.rds} file. The parent directory will be created if it does not
#'   exist. Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes the partition mode, fold
#'   structure, per-fold point counts, and file-save confirmation.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{folds}: Data frame of fold assignments with a \code{fold}
#'     column identifying each point's cross-validation fold.
#'   \item \code{points_sf}: sf object of occurrence points with assigned
#'     folds.
#'   \item \code{voronoi_folds}: sf object of Voronoi polygons representing
#'     the spatial fold boundaries. \code{NULL} for random folds, temporal-only
#'     folds, and single-fold mode.
#'   \item \code{summary}: Data frame of partitioning summary statistics.
#'   \item \code{plots}: Named list of recorded plot objects when
#'     \code{create_plot = TRUE}. Empty list in single-fold mode.
#' }
#'
#' @details
#' The function partitions data into folds using one of five modes:
#' \itemize{
#'   \item \strong{Single fold}: All points are assigned to fold 1 and used
#'     for both training and testing. This produces a single model trained on
#'     the full dataset with no held-out validation. Useful when sample sizes
#'     are too small for cross-validation, or as a final production model step
#'     after cross-validation has already established model quality. Set
#'     \code{single_fold = TRUE}. All downstream functions accept the result
#'     identically to standard multi-fold output.
#'
#'   \item \strong{Random}: Points are assigned to folds by random shuffling
#'     with no spatial or temporal structure. Each fold is a simple random
#'     sample of the full dataset, intended as a naive baseline that makes no
#'     attempt to reduce spatial or temporal autocorrelation between training
#'     and test sets. Use \code{n_random_folds}.
#'
#'   \item \strong{Spatial-only}: The study area is divided into \eqn{k}
#'     contiguous spatial regions using a recursive k-d tree bisection
#'     algorithm. At each step the point set is split along its longest
#'     spatial axis, recursively halving until the target number of folds
#'     is reached. A centroid reassignment pass then refines boundaries
#'     to improve balance. Each region becomes one fold, so training
#'     always occurs on data from geographically distinct areas relative
#'     to the test fold. No temporal separation is imposed, meaning that
#'     points from any time period may appear in any fold. Use
#'     \code{n_spatial_folds} alone.
#'
#'   \item \strong{Temporal-only}: Each fold covers the full spatial extent
#'     of the study area but is restricted to a distinct, non-overlapping
#'     slice of the time series. The global time series is divided into
#'     \code{n_temporal_folds} equal intervals using quantile-based breaks,
#'     and all points within each interval form one fold. This design tests
#'     model transferability across time while retaining full spatial coverage
#'     in every fold. Use \code{n_temporal_folds} alone (with
#'     \code{n_spatial_folds = 0}). Requires \code{time_cols}.
#'
#'   \item \strong{Spatiotemporal}: Folds are assigned using the same
#'     recursive k-d tree bisection as spatial-only mode, operating on
#'     the full point set to produce spatially contiguous groups. The
#'     resulting groups are then split into a spatial pool
#'     (\code{n_spatial_folds} folds drawn from geographically distinct
#'     regions) and a temporal pool (\code{n_temporal_folds} folds each
#'     restricted to a distinct slice of the time series but spanning
#'     the full study area). Together the two pools assess both geographic
#'     and temporal transferability in a single cross-validation design.
#'     Use \code{n_spatial_folds} and \code{n_temporal_folds} together.
#'     Requires \code{time_cols}.
#' }
#'
#' Fold assignment uses a recursive k-d tree bisection algorithm that splits
#' points along their longest spatial axis at each step, followed by a
#' centroid reassignment pass to improve boundary regularity and point-count
#' balance. Voronoi tessellation on fold centroids is used only for
#' visualisation of the resulting spatial boundaries. For temporal
#' mode, temporal blocks are defined by dividing the global time series
#' into equal intervals using quantile-based breaks. For spatiotemporal mode,
#' the typical spatial assignment is done, but with one larger spatial block
#' made with enough points to represent all of the temporal folds, then the
#' temporal blocking is applied to those points.
#'
#' Partitioned datasets are suitable for cross-validation in modeling
#' workflows, ensuring spatial and/or temporal independence between folds.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_rarefaction}},
#'   \code{\link{temporally_explicit_extraction}},
#'   \code{\link{generate_absences}}
#'
#' Modeling: \code{\link{build_temporal_hv}},
#'   \code{\link{build_temporal_glm}}, \code{\link{build_temporal_gam}},
#'   \code{\link{build_temporal_rf}}
#'
#' @examples
#' \donttest{
#'   pts_file <- system.file(
#'     "extdata/points/extracted_seasonal_Scaled_Values.csv",
#'     package = "TemporalModelR"
#'   )
#'
#'   ref_file <- system.file("extdata/rasters_raw/elevation.tif",
#'                           package = "TemporalModelR")
#'
#'   study_crs <- sf::st_crs(terra::rast(ref_file))
#'
#'   study_area_sf <- sf::st_as_sf(sf::st_as_sfc(
#'     sf::st_bbox(c(xmin = 0, xmax = 3000, ymin = 0, ymax = 1500),
#'                 crs = study_crs)
#'   ))
#'
#'   spatiotemporal_partition(
#'     reference_shapefile_path = study_area_sf,
#'     points_file_path         = pts_file,
#'     xcol                     = "x",
#'     ycol                     = "y",
#'     points_crs               = study_crs,
#'     time_cols                = "year",
#'     n_spatial_folds          = 2,
#'     n_temporal_folds         = 2,
#'     create_plot              = FALSE,
#'     verbose                  = FALSE
#'   )
#' }

#' @export
#' @importFrom sf st_read st_as_sf st_transform st_coordinates st_crs st_bbox
#' @importFrom sf st_drop_geometry st_intersection st_polygon st_sf st_sfc st_union
#' @importFrom deldir deldir tile.list
#' @importFrom stats complete.cases median quantile
#' @importFrom tools file_ext
#' @importFrom utils capture.output
spatiotemporal_partition <- function(reference_shapefile_path,
                                     points_file_path,
                                     time_cols         = NULL,
                                     xcol             = NULL,
                                     ycol             = NULL,
                                     points_crs       = NULL,
                                     n_spatial_folds  = 0,
                                     n_temporal_folds = 0,
                                     n_balanced_folds = 0,
                                     n_random_folds   = 0,
                                     single_fold      = FALSE,
                                     max_imbalance    = 0.05,
                                     max_attempts     = 10,
                                     create_plot   = TRUE,
                                     plot_palette     = "Dark 2",
                                     output_file      = NULL,
                                     verbose          = TRUE) {

  if (missing(reference_shapefile_path) || is.null(reference_shapefile_path)) {
    stop(paste0("ERROR: 'reference_shapefile_path' is required but was not provided. ",
                "Please provide a file path or sf object defining the study area."))
  }
  if (missing(points_file_path) || is.null(points_file_path)) {
    stop(paste0("ERROR: 'points_file_path' is required but was not provided. ",
                "Please provide a file path, sf object, or data frame of occurrence points."))
  }

  if (!is.null(output_file)) {
    if (tolower(tools::file_ext(output_file)) != "rds") {
      stop("ERROR: 'output_file' must have a '.rds' extension.")
    }
    out_dir <- dirname(output_file)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  if (isTRUE(single_fold)) {
    use_random            <- FALSE
    use_balanced          <- FALSE
    temporal_partitioning <- FALSE
    total_folds           <- 1
    n_spatial             <- 1
    n_temporal            <- 1
    partition_mode        <- "single_fold"
  } else {

    use_random   <- n_random_folds > 0
    use_balanced <- n_balanced_folds > 0

    modes_active <- sum(c(
      use_random,
      use_balanced,
      (n_spatial_folds > 0 || n_temporal_folds > 0)
    ))
    if (modes_active == 0) {
      stop(paste0("ERROR: Must specify one of n_random_folds, n_balanced_folds, ",
                  "n_spatial_folds/n_temporal_folds, or set single_fold = TRUE."))
    }
    if (modes_active > 1) {
      stop(paste0("ERROR: Only one partitioning mode is allowed at a time. ",
                  "Use n_random_folds, n_balanced_folds, or n_spatial_folds/n_temporal_folds."))
    }
    if (!use_random) {
      if (n_spatial_folds == 1) {
        stop("ERROR: n_spatial_folds must be 0 or >= 2.")
      }
      if (n_temporal_folds == 1) {
        stop("ERROR: n_temporal_folds must be 0 or >= 2.")
      }
    }

    if (use_random) {
      total_folds           <- n_random_folds
      n_spatial             <- 1
      n_temporal            <- 1
      temporal_partitioning <- FALSE
      partition_mode        <- "random"

    } else if (use_balanced) {
      if (is.null(time_cols)) {
        stop("ERROR: 'time_cols' is required for balanced folds. Please specify the column name containing time or year values.")
      }
      total_folds           <- n_balanced_folds
      n_spatial             <- n_balanced_folds * ifelse(n_balanced_folds %% 2 == 0, 4, 5)
      n_temporal            <- 2
      temporal_partitioning <- TRUE
      partition_mode        <- "balanced"
      warning("Balanced folds are more difficult to perfectly create. You may need a more flexible max_imbalance threshold.")

    } else {
      total_folds <- n_spatial_folds + n_temporal_folds
      if (n_spatial_folds == 0 && n_temporal_folds > 0) {
        if (is.null(time_cols)) {
          stop("ERROR: 'time_cols' is required for temporal folds. Please specify the column name containing time or year values.")
        }
        warning("Only temporal folds specified. Spatial structure will not be considered.")
        n_spatial             <- 1
        n_temporal            <- n_temporal_folds
        temporal_partitioning <- TRUE
        partition_mode        <- "temporal_only"
      } else if (n_temporal_folds == 0 && n_spatial_folds > 0) {
        warning("Only spatial folds specified. Temporal structure will not be considered.")
        n_spatial             <- NULL
        n_temporal            <- 1
        temporal_partitioning <- FALSE
        partition_mode        <- "spatial_only"
      } else {
        if (is.null(time_cols)) {
          stop("ERROR: 'time_cols' is required when temporal partitioning is enabled. Please specify the column name containing time or year values.")
        }
        n_spatial             <- n_spatial_folds * 2 * n_temporal_folds
        n_temporal            <- n_temporal_folds
        temporal_partitioning <- TRUE
        partition_mode        <- "spatiotemporal"
      }
    }

  }

  if (is.character(reference_shapefile_path) && verbose) {
    message(paste("Reading shapefile:", basename(reference_shapefile_path)))
  }
  reference_shapefile <- .load_shapefile_input(reference_shapefile_path, "reference_shapefile_path")

  pts        <- .load_points_data(points_file_path, xcol, ycol, points_crs,
                                  verbose = verbose)
  xcol       <- attr(pts, "xcol")
  ycol       <- attr(pts, "ycol")
  points_crs <- attr(pts, "crs")

  n_original <- nrow(pts)
  pts        <- pts[stats::complete.cases(pts), ]
  n_removed  <- n_original - nrow(pts)
  if (n_removed > 0) {
    warning(paste0("Removed ", n_removed, " incomplete rows (",
                   round(n_removed / n_original * 100, 2), "%)"))
  }
  if (nrow(pts) == 0) {
    stop("ERROR: No complete rows remaining after removing incomplete data. Please check your input data.")
  }
  if (!xcol %in% names(pts)) {
    stop(paste0("ERROR: Column '", xcol, "' not found in the input data."))
  }
  if (!ycol %in% names(pts)) {
    stop(paste0("ERROR: Column '", ycol, "' not found in the input data."))
  }
  if (temporal_partitioning && !time_cols %in% names(pts)) {
    stop(paste0("ERROR: Column '", time_cols, "' not found in the input data."))
  }

  pts_sf       <- sf::st_as_sf(pts, coords = c(xcol, ycol), crs = points_crs)
  pts_sf       <- sf::st_transform(pts_sf, crs = sf::st_crs(reference_shapefile))
  total_points <- nrow(pts_sf)

  if (isTRUE(single_fold)) {
    if (verbose) message("\nPartition mode: SINGLE FOLD (all points used for training and testing)")

    pts_sf$fold <- 1

    summary_stats <- data.frame(
      parameter = c("total_folds", "n_spatial_folds", "n_temporal_folds",
                    "n_balanced_folds", "n_random_folds", "partition_mode",
                    "total_points", "points_removed", "pct_rows_removed",
                    "final_imbalance_pct"),
      value = c(1, 0, 0, 0, 0, "single_fold",
                total_points, n_removed,
                ifelse(n_removed > 0, round(n_removed / n_original * 100, 2), 0),
                0),
      stringsAsFactors = FALSE
    )
    if (verbose) message(paste0("  1 fold | ", total_points, " points (100%)"))

    result <- list(
      folds        = sf::st_drop_geometry(pts_sf)[, "fold", drop = FALSE],
      points_sf    = pts_sf,
      voronoi_folds = NULL,
      summary      = summary_stats,
      plots        = list()
    )

    if (!is.null(output_file)) {
      out_dir <- dirname(output_file)
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(result, output_file)
      if (verbose) message(paste0("Results saved to: ", output_file))
    }

    return(invisible(result))
  }

  if (total_points < total_folds) {
    stop(paste0("ERROR: Not enough points (", total_points, ") for ", total_folds,
                " folds. Please reduce the number of folds or provide more occurrence data."))
  }

  if (partition_mode == "spatial_only") {
    n_unique_coords <- nrow(unique(sf::st_coordinates(pts_sf)))
    n_spatial <- min(
      floor(total_points / n_spatial_folds),
      n_unique_coords - 1
    )
    n_spatial <- max(n_spatial, n_spatial_folds * 2)
  }


  if (use_random) {
    if (verbose) message(paste0("\nPartition mode: RANDOM (", n_random_folds, " folds)"))

    shuffled   <- sample(total_points)
    fold_sizes <- diff(round(seq(0, total_points, length.out = n_random_folds + 1)))
    pts_sf$fold          <- NA_integer_
    pts_sf$fold[shuffled] <- rep(seq_len(n_random_folds), times = fold_sizes)
    pts_sf$spatial_block  <- NA_integer_
    pts_sf$temporal_block <- NA_integer_
    pts_sf$block_type     <- "random"

    final_fold_counts <- table(factor(pts_sf$fold, levels = seq_len(n_random_folds)))
    mean_per_fold     <- mean(final_fold_counts)
    imbalance         <- max(abs(final_fold_counts - mean_per_fold)) / mean_per_fold

    plot_list <- if (create_plot) {
      .plot_partitions_base(pts_sf, reference_shapefile, final_fold_counts, mean_per_fold,
                            n_random_folds, "random", time_cols, !is.null(time_cols), 1, NULL,
                            plot_palette)
    } else {
      list()
    }

    result <- .build_partition_result(pts_sf, NULL, final_fold_counts, mean_per_fold, imbalance,
                                      n_spatial, n_temporal, total_folds, partition_mode,
                                      temporal_partitioning, FALSE, TRUE,
                                      0, 0, 0, n_random_folds,
                                      n_removed, n_original, plot_list, output_file,
                                      verbose = verbose)
    if (verbose) message(paste(utils::capture.output(
      print(result$summary, row.names = FALSE)), collapse = "\n"))
    return(invisible(result))
  }

  if (n_spatial >= total_points) {
    stop(paste0("ERROR: Spatial blocks (", n_spatial, ") >= points (", total_points,
                "). Please reduce the number of folds or provide more occurrence data."))
  }

  if (partition_mode == "spatial_only") {
    if (verbose) message(paste0("\nPartition mode: SPATIAL_ONLY (", n_spatial_folds, " spatial folds)"))

    pts_work               <- pts_sf
    pts_work$spatial_block <- NA_integer_
    pts_work$temporal_block <- 1
    pts_work$fold          <- NA_integer_
    pts_work$block_type    <- NA_character_

    spatial_info           <- .build_spatial_blocks(pts_work, reference_shapefile, n_spatial)
    pts_work$spatial_block <- spatial_info$clusters

    pts_work <- .assign_spatiotemporal_folds(pts_work, n_spatial_folds, n_temporal_folds,
                                             total_folds, n_spatial, n_temporal,
                                             spatial_info$centers, spatial_info$dist_matrix,
                                             spatial_info$adjacency, temporal_partitioning,
                                             max_imbalance, use_balanced)

    pts_work$spatial_block <- pts_work$fold

    fold_coords <- sf::st_coordinates(pts_work)
    fold_cx <- vapply(seq_len(total_folds), function(f) {
      mean(fold_coords[pts_work$fold == f, 1])
    }, numeric(1))
    fold_cy <- vapply(seq_len(total_folds), function(f) {
      mean(fold_coords[pts_work$fold == f, 2])
    }, numeric(1))
    bbox  <- sf::st_bbox(reference_shapefile)
    vor   <- deldir::deldir(fold_cx, fold_cy,
                            rw = c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]))
    tiles <- deldir::tile.list(vor)
    polys <- lapply(seq_along(tiles), function(i) {
      tile <- tiles[[i]]
      sf::st_polygon(list(cbind(c(tile$x, tile$x[1]), c(tile$y, tile$y[1]))))
    })
    voronoi_sf <- sf::st_sf(fold_group = seq_len(total_folds),
                            geometry = sf::st_sfc(polys, crs = sf::st_crs(pts_work)))
    voronoi_sf <- suppressWarnings(sf::st_intersection(voronoi_sf, sf::st_union(reference_shapefile)))

    fc  <- table(factor(pts_work$fold, levels = seq_len(total_folds)))
    mpf <- mean(fc)
    imb <- max(abs(fc - mpf)) / mpf

    if (imb > max_imbalance) {
      warning(paste0("Spatial-only folds are imbalanced (", round(imb * 100, 2),
                     "%) due to uneven point density. This is expected for spatial CV."))
    }

    plot_list <- if (create_plot) {
      .plot_partitions_base(pts_work, reference_shapefile, fc, mpf,
                            total_folds, partition_mode, time_cols, temporal_partitioning,
                            n_temporal, voronoi_sf, plot_palette)
    } else {
      list()
    }

    result <- .build_partition_result(pts_work, voronoi_sf, fc, mpf, imb,
                                      n_spatial, n_temporal, total_folds, partition_mode,
                                      temporal_partitioning, use_balanced, FALSE,
                                      n_spatial_folds, n_temporal_folds, n_balanced_folds, 0,
                                      n_removed, n_original, plot_list, output_file,
                                      verbose = verbose)
    if (verbose) message(paste(utils::capture.output(
      print(result$summary, row.names = FALSE)), collapse = "\n"))
    return(invisible(result))
  }

  best_imbalance       <- Inf
  best_results         <- NULL
  best_voronoi         <- NULL

  for (attempt in seq_len(max_attempts)) {
    pts_work               <- pts_sf
    pts_work$spatial_block  <- NA_integer_
    pts_work$temporal_block <- NA_integer_
    pts_work$fold           <- NA_integer_
    pts_work$block_type     <- NA_character_

    if (attempt == 1) {
      if (partition_mode == "spatiotemporal") {
        if (verbose) message(paste0("\nPartition mode: SPATIOTEMPORAL (",
                                    n_spatial_folds, " spatial folds, ",
                                    n_temporal_folds, " temporal folds, ",
                                    total_folds, " total folds)"))
      } else {
        if (verbose) message(paste0("\nPartition mode: ", toupper(partition_mode), " (",
                                    n_spatial, " spatial blocks, ",
                                    n_temporal, " temporal blocks, ",
                                    total_folds, " folds)"))
      }
    }

    spatial_info            <- .build_spatial_blocks(pts_work, reference_shapefile, n_spatial)
    pts_work$spatial_block  <- spatial_info$clusters
    voronoi_sf              <- spatial_info$voronoi_sf

    if (temporal_partitioning) {
      tv <- pts_work[[time_cols]]

      if (n_temporal == 1) {
        pts_work$temporal_block <- 1
      } else if (n_temporal == 2) {
        global_mid <- stats::median(tv, na.rm = TRUE)
        pts_work$temporal_block <- ifelse(tv <= global_mid, 1, 2)
      } else {
        breaks <- stats::quantile(tv, probs = seq(0, 1, length.out = n_temporal + 1),
                                  na.rm = TRUE)
        if (any(duplicated(breaks))) {
          breaks <- unique(breaks)
          warning("Duplicate temporal breaks detected.")
        }
        pts_work$temporal_block <- as.integer(cut(tv, breaks = breaks,
                                                  labels = FALSE, include.lowest = TRUE))
      }

      if (n_temporal == 2 && attempt > 1) {
        tb_counts <- table(factor(pts_work$temporal_block, levels = 1:2))
        tb_imb    <- max(tb_counts) / sum(tb_counts)
        if (tb_imb > 0.65) {
          pool_idx <- which(!is.na(pts_work$spatial_block))
          tv_pool  <- tv[pool_idx]
          alt_mid  <- stats::median(tv_pool, na.rm = TRUE)
          if (!is.na(alt_mid) && alt_mid != global_mid) {
            pts_work$temporal_block[pool_idx] <- ifelse(tv_pool <= alt_mid, 1, 2)
          }
        }
      }

    } else {
      pts_work$temporal_block <- 1
    }

    if (use_balanced) {
      pts_work <- .assign_balanced_folds(pts_work, n_balanced_folds, n_spatial,
                                         spatial_info$centers, spatial_info$dist_matrix,
                                         spatial_info$adjacency, max_imbalance, total_folds)
    } else {
      pts_work <- .assign_spatiotemporal_folds(pts_work, n_spatial_folds, n_temporal_folds,
                                               total_folds, n_spatial, n_temporal,
                                               spatial_info$centers, spatial_info$dist_matrix,
                                               spatial_info$adjacency, temporal_partitioning,
                                               max_imbalance, use_balanced)
    }

    unassigned <- which(is.na(pts_work$fold))
    for (idx in unassigned) {
      counts <- table(factor(pts_work$fold[!is.na(pts_work$fold)],
                             levels = seq_len(total_folds)))
      pts_work$fold[idx]       <- as.integer(names(which.min(counts))[1])
      pts_work$block_type[idx] <- "remainder"
    }

    fc  <- table(factor(pts_work$fold, levels = seq_len(total_folds)))
    mpf <- mean(fc)
    imb <- max(abs(fc - mpf)) / mpf

    if (imb < best_imbalance) {
      best_imbalance <- imb
      best_results   <- pts_work
      best_voronoi   <- voronoi_sf
      best_counts    <- fc
      best_mean      <- mpf
    }
    if (imb <= max_imbalance) break
  }

  if (best_imbalance > max_imbalance) {
    warning(paste0("Could not achieve target balance within ",
                   max_attempts, " attempts. ",
                   "Final imbalance: ", round(best_imbalance * 100, 2), "%. ",
                   "Returning best result achieved. Try increasing max_imbalance or adjusting the fold configuration."))
  }

  pts_sf            <- best_results
  voronoi_sf        <- best_voronoi
  final_fold_counts <- best_counts
  mean_per_fold     <- best_mean
  imbalance         <- best_imbalance

  plot_list <- if (create_plot) {
    .plot_partitions_base(pts_sf, reference_shapefile, final_fold_counts, mean_per_fold,
                          total_folds, partition_mode, time_cols, temporal_partitioning,
                          n_temporal, voronoi_sf, plot_palette)
  } else {
    list()
  }

  result <- .build_partition_result(pts_sf, voronoi_sf, final_fold_counts, mean_per_fold,
                                    imbalance, n_spatial, n_temporal, total_folds, partition_mode,
                                    temporal_partitioning, use_balanced, FALSE,
                                    n_spatial_folds, n_temporal_folds, n_balanced_folds, 0,
                                    n_removed, n_original, plot_list, output_file,
                                    verbose = verbose)
  if (verbose) message(paste(utils::capture.output(
    print(result$summary, row.names = FALSE)), collapse = "\n"))
  invisible(result)
}
