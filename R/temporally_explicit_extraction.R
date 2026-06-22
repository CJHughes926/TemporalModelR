#' Extract Time-Aligned Environmental Values at Species Occurrences
#'
#' Preprocessing function that extracts raster values to species occurrence
#' records based on temporal components. Matches environmental layers to
#' occurrence timestamps and optionally computes scaling parameters for standardization.
#'
#' @usage
#' temporally_explicit_extraction(points_sp, raster_dir, variable_patterns,
#'                                time_cols, xcol = NULL, ycol = NULL,
#'                                points_crs = NULL, output_dir,
#'                                output_prefix = "temp_explicit_df",
#'                                save_raw = TRUE, save_scaled = TRUE,
#'                                save_scaling_params = TRUE,
#'                                verbose = TRUE)
#'
#' @param points_sp sf object, SpatialPointsDataFrame, file path to
#'   .csv/.shp/.geojson/.gpkg, or data frame with coordinate columns.
#' @param raster_dir Character. Directory containing environmental raster
#'   files (\code{.tif}), typically the output of \code{\link{raster_align}}.
#'   File names must follow the patterns supplied in \code{variable_patterns}, with any
#'   time placeholder substituted for the corresponding value from
#'   \code{time_cols}.
#' @param variable_patterns Named character vector mapping clean variable names
#'   to raster filename patterns. For time-varying variables include the time
#'   placeholder in the pattern (e.g. \code{"forest_cover" = "forest_cover_YEAR"});
#'   for static variables omit it (e.g. \code{"elevation" = "elevation"}). Time
#'   placeholders must match entries in \code{time_cols}.
#' @param time_cols Character vector of time column names present in the point
#'   data (e.g., c("YEAR"), c("YEAR", "MONTH")).
#' @param xcol Character. Name of the x-coordinate column. Required when
#'   \code{points_sp} is a CSV file or data frame.
#' @param ycol Character. Name of the y-coordinate column. Required when
#'   \code{points_sp} is a CSV file or data frame.
#' @param points_crs Character or CRS object. CRS of the input points.
#'   Required when \code{points_sp} is a CSV file or data frame.
#' @param output_dir Character. Directory to write output files.
#' @param output_prefix Character. Prefix for output filenames. Default is
#'   "temp_explicit_df".
#' @param save_raw Logical. If \code{TRUE} (default), writes raw extracted values
#'   CSV. If \code{FALSE}, skips raw values output.
#' @param save_scaled Logical. If \code{TRUE} (default), writes z-scaled values
#'   CSV. If \code{FALSE}, skips scaled values output.
#' @param save_scaling_params Logical. If \code{TRUE} (default), writes CSV of
#'   per-variable means and standard deviations. If \code{FALSE}, skips scaling
#'   parameters output.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes file loading, extraction
#'   progress, and file-save confirmation.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{raw_values}: Data frame of raw extracted values at each
#'     occurrence record (when \code{save_raw = TRUE}; \code{NULL} otherwise).
#'   \item \code{scaled_values}: Data frame of z-scaled extracted values
#'     (when \code{save_scaled = TRUE}; \code{NULL} otherwise).
#'   \item \code{scaling_params}: Data frame of per-variable means and standard
#'     deviations used for scaling (when \code{save_scaling_params = TRUE};
#'     \code{NULL} otherwise). Pass this to \code{\link{scale_rasters}}.
#'   \item \code{files_created}: Named list of file paths written, with elements
#'     \code{raw}, \code{scaled}, and \code{scaling_params} (each \code{NULL}
#'     when the corresponding save flag is \code{FALSE}).
#' }
#'
#' @details
#' Extracts raster values to species occurrence records based on matched temporal components.
#' Matches environmental layers to occurrence timestamps and optionally computes scaling
#' parameters for standardization.
#'
#' Output CSV files are written to output_dir containing raw values, scaled
#' values, and scaling parameters.
#'
#' Scaling parameters (mean and standard deviation) are optionally computed across
#' all occurrence records for each variable. These parameters should be used with
#' \code{\link{scale_rasters}} to standardize prediction layers.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_rarefaction}},
#'   \code{\link{scale_rasters}}, \code{\link{spatiotemporal_partition}}
#'
#' @examples
#' pts_file <- system.file(
#'   "extdata/points/synthetic_occurrence_points.csv",
#'   package = "TemporalModelR"
#' )
#'
#' aln_dir  <- system.file("extdata/rasters_aligned",
#'                         package = "TemporalModelR")
#'
#' ref_file <- system.file("extdata/rasters_raw/elevation.tif",
#'                         package = "TemporalModelR")
#'
#' out_dir  <- file.path(tempdir(), "extracted")
#'
#' temporally_explicit_extraction(
#'   points_sp           = pts_file,
#'   raster_dir          = aln_dir,
#'   variable_patterns   = c(
#'     "elevation"    = "elevation",
#'     "forest_cover" = "forest_cover_YEAR",
#'     "prseas"       = "prseas_YEAR_SEASON"
#'   ),
#'   time_cols           = c("year", "season"),
#'   xcol                = "x",
#'   ycol                = "y",
#'   points_crs          = terra::crs(terra::rast(ref_file)),
#'   output_dir          = out_dir,
#'   save_raw            = TRUE,
#'   save_scaled         = FALSE,
#'   save_scaling_params = TRUE,
#'   verbose             = FALSE
#' )

#' @export
#' @importFrom terra rast extract vect
#' @importFrom sf st_as_sf st_drop_geometry st_coordinates
#' @importFrom stats sd setNames
#' @importFrom utils write.csv
temporally_explicit_extraction <- function(points_sp,
                                           raster_dir,
                                           variable_patterns,
                                           time_cols,
                                           xcol = NULL,
                                           ycol = NULL,
                                           points_crs = NULL,
                                           output_dir,
                                           output_prefix       = "temp_explicit_df",
                                           save_raw            = TRUE,
                                           save_scaled         = TRUE,
                                           save_scaling_params = TRUE,
                                           verbose             = TRUE) {

  pts_df    <- .load_points_data(points_sp, xcol, ycol, points_crs, verbose = verbose)
  xcol      <- attr(pts_df, "xcol")
  ycol      <- attr(pts_df, "ycol")
  pts_crs   <- attr(pts_df, "crs")
  points_sp <- sf::st_as_sf(pts_df, coords = c(xcol, ycol), crs = pts_crs)

  .validate_variable_patterns(variable_patterns)

  if (missing(time_cols) || !is.character(time_cols) || length(time_cols) == 0) {
    stop("ERROR: 'time_cols' must be a character vector with at least one column name.")
  }
  missing_cols <- setdiff(time_cols, names(points_sp))
  if (length(missing_cols) > 0) {
    stop(paste0("ERROR: The following time_cols are missing from data: ",
                paste(missing_cols, collapse = ", ")))
  }

  if (!dir.exists(raster_dir)) stop(paste0("ERROR: 'raster_dir' does not exist: ", raster_dir))
  if (verbose) message(paste("Found",
                             length(list.files(raster_dir, "\\.tif$", recursive = TRUE)),
                             "raster files"))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  if (verbose) message(paste("Processing", nrow(points_sp), "points"))

  var_classes  <- .classify_variable_patterns(variable_patterns, time_cols)
  dynamic_vars <- var_classes$dynamic
  static_vars  <- var_classes$static

  if (verbose) message(paste("Dynamic variables:",
                             if (length(dynamic_vars) > 0) paste(dynamic_vars, collapse = ", ") else "none"))
  if (verbose) message(paste("Static variables:",
                             if (length(static_vars) > 0) paste(static_vars, collapse = ", ") else "none"))

  for (var in names(variable_patterns)) points_sp[[var]] <- NA

  if (length(static_vars) > 0) {
    if (verbose) message("Extracting static variables...")
    static_paths <- .resolve_raster_paths(variable_patterns, character(0), static_vars,
                                          time_cols,
                                          as.data.frame(stats::setNames(
                                            lapply(time_cols, function(tc) NA_character_),
                                            time_cols), stringsAsFactors = FALSE),
                                          raster_dir)
    if (!is.null(static_paths)) {
      r_static <- tryCatch(
        terra::rast(static_paths),
        error = function(e) { warning(paste("Could not load static rasters:", e$message)); NULL }
      )
      if (!is.null(r_static)) {
        env_ex <- terra::extract(r_static, terra::vect(points_sp), ID = FALSE)
        if (ncol(env_ex) == length(static_vars)) names(env_ex) <- static_vars
        for (v in static_vars) {
          if (v %in% names(env_ex)) points_sp[[v]] <- env_ex[[v]]
        }
      }
    }
  }

  if (length(dynamic_vars) > 0) {
    if (verbose) message("Extracting dynamic variables...")

    pts_data <- sf::st_drop_geometry(points_sp)
    time_combinations <- unique(pts_data[, time_cols, drop = FALSE])
    time_combinations <- time_combinations[do.call(order, time_combinations), , drop = FALSE]

    if (verbose) message(paste("Extracting values for", nrow(time_combinations), "time periods"))
    for (i in seq_len(nrow(time_combinations))) {
      time_values <- time_combinations[i, , drop = FALSE]
      time_filter <- Reduce(`&`, lapply(time_cols,
                                        function(tc) points_sp[[tc]] == time_values[[tc]]))
      points_subset <- points_sp[time_filter, ]
      if (nrow(points_subset) == 0) next

      raster_paths <- .resolve_raster_paths(variable_patterns, dynamic_vars, character(0),
                                            time_cols, time_values, raster_dir)
      if (is.null(raster_paths)) next

      r_ts <- tryCatch(
        terra::rast(raster_paths),
        error = function(e) { warning(paste("Could not load rasters:", e$message)); NULL }
      )
      if (is.null(r_ts)) next

      env_ex <- terra::extract(r_ts, terra::vect(points_subset), ID = FALSE)
      if (ncol(env_ex) == length(dynamic_vars)) names(env_ex) <- dynamic_vars
      for (v in dynamic_vars) {
        if (v %in% names(env_ex)) points_sp[[v]][time_filter] <- env_ex[[v]]
      }
    }
  }

  coords_mat         <- sf::st_coordinates(points_sp)
  coords_df          <- data.frame(x = coords_mat[, 1], y = coords_mat[, 2])
  points_with_coords <- cbind(sf::st_drop_geometry(points_sp), coords_df)

  raw_output_file    <- NULL
  params_file        <- NULL
  scaled_output_file <- NULL
  scaled_data        <- NULL

  if (save_raw) {
    raw_output_file <- file.path(output_dir, paste0(output_prefix, "_Raw_Values.csv"))
    utils::write.csv(points_with_coords, raw_output_file, row.names = FALSE)
    if (verbose) message(paste("Raw values saved to:", basename(raw_output_file)))
  }

  if (verbose) message("Calculating scaling parameters...")
  scaling_params <- data.frame(variable = character(), mean = numeric(), sd = numeric(),
                               stringsAsFactors = FALSE)
  for (var_name in names(variable_patterns)) {
    values <- points_sp[[var_name]]
    values <- values[!is.na(values)]
    if (length(values) > 0) {
      scaling_params <- rbind(scaling_params,
                              data.frame(variable = var_name, mean = mean(values),
                                         sd = stats::sd(values), stringsAsFactors = FALSE))
    } else {
      warning("No valid values found for ", var_name)
    }
  }

  if (save_scaling_params) {
    params_file <- file.path(output_dir, paste0(output_prefix, "_Scaling_Parameters.csv"))
    utils::write.csv(scaling_params, params_file, row.names = FALSE)
    if (verbose) message(paste("Scaling parameters saved to:", basename(params_file)))
  }

  if (save_scaled) {
    if (verbose) message("Applying scaling...")
    scaled_data <- sf::st_drop_geometry(points_sp)
    for (var_name in names(variable_patterns)) {
      var_params <- scaling_params[scaling_params$variable == var_name, ]
      if (nrow(var_params) > 0) {
        scaled_data[[var_name]] <- (scaled_data[[var_name]] - var_params$mean) / var_params$sd
      } else {
        warning("No scaling parameters found for ", var_name)
      }
    }
    scaled_data <- cbind(scaled_data, coords_df)
    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    utils::write.csv(scaled_data, scaled_output_file, row.names = FALSE)
    if (verbose) message(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  if (verbose) message("Processing complete.")

  invisible(list(
    raw_values     = if (save_raw)            points_with_coords else NULL,
    scaled_values  = if (save_scaled)         scaled_data        else NULL,
    scaling_params = if (save_scaling_params) scaling_params     else NULL,
    files_created  = list(
      raw            = if (save_raw)            raw_output_file else NULL,
      scaled         = if (save_scaled)         scaled_output_file else NULL,
      scaling_params = if (save_scaling_params) params_file else NULL
    )
  ))
}
