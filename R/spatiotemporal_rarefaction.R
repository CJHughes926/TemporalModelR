#' Spatiotemporal Rarefaction of Species Occurrence Data
#'
#' Preprocessing function that rarefies species occurrence data to one point per
#' raster pixel, optionally accounting for temporal components. Reduces sampling
#' bias and spatial autocorrelation in occurrence datasets.
#'
#' @usage
#' spatiotemporal_rarefaction(points_sp, output_dir, reference_raster,
#'                            time_cols = NULL, xcol = NULL, ycol = NULL,
#'                            points_crs = NULL, output_prefix = "Pts_Database",
#'                            verbose = TRUE)
#'
#' @param points_sp Input point data. Accepts an sf object, data frame,
#'   SpatialPointsDataFrame, or file path to a \code{.csv}, \code{.shp},
#'   \code{.geojson}, or \code{.gpkg} file.
#' @param output_dir Character. Directory where output CSV files will be saved.
#' @param reference_raster Character, SpatRaster, or RasterLayer. Raster used
#'   to define pixel boundaries for rarefaction. Accepts a file path,
#'   RasterLayer, or SpatRaster.
#' @param time_cols Character vector. Column names in \code{points_sp} defining
#'   temporal grouping for spatiotemporal rarefaction. When \code{NULL}
#'   (default), only spatial rarefaction is performed.
#' @param xcol Character. Name of the x-coordinate column. Required when
#'   \code{points_sp} is a CSV file or data frame.
#' @param ycol Character. Name of the y-coordinate column. Required when
#'   \code{points_sp} is a CSV file or data frame.
#' @param points_crs Character or CRS object. CRS of the input points.
#'   Required when \code{points_sp} is a CSV file or data frame.
#' @param output_prefix Character. Prefix for output file names. Default is
#'   \code{"Pts_Database"}. Output files are named
#'   \code{<output_prefix>_OnePerPix.csv} and, when \code{time_cols} are
#'   provided, \code{<output_prefix>_OnePerPixPerTimeStep.csv}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes file loading, missing-data removal
#'   counts, rarefaction summaries, and a final comparison of spatial
#'   vs spatiotemporal point counts.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{input_points}: Integer. Number of input points after CRS
#'     alignment and removal of rows with missing time column values.
#'   \item \code{spatial_points}: Integer. Number of points retained after
#'     spatial-only rarefaction (one per pixel).
#'   \item \code{spatiotemporal_points}: Integer. Number of points retained
#'     after spatiotemporal rarefaction (one per pixel per time combination).
#'     \code{NA} when \code{time_cols} is not provided.
#'   \item \code{time_cols_used}: Character vector of time columns used.
#'     \code{NULL} when \code{time_cols} is not provided.
#'   \item \code{spatial_table}: Data frame of spatially rarefied points with
#'     columns \code{pixel_id}, \code{X}, \code{Y}, and any \code{time_cols}.
#'   \item \code{spatiotemporal_table}: Data frame of spatiotemporally rarefied
#'     points with columns \code{pixel_id}, \code{X}, \code{Y}, and
#'     \code{time_cols}. \code{NULL} when \code{time_cols} is not provided.
#'   \item \code{files_created}: Named list of file paths written. Always
#'     contains \code{$spatial}; additionally contains
#'     \code{$spatiotemporal} when \code{time_cols} are provided.
#' }
#'
#' @details
#' The function assigns each point to a raster pixel using the resolution and
#' extent of \code{reference_raster}, then performs:
#' \itemize{
#'   \item Spatial rarefaction: retains one point per pixel, written to
#'     \code{<output_prefix>_OnePerPix.csv}.
#'   \item Spatiotemporal rarefaction (when \code{time_cols} are provided):
#'     retains one point per pixel per unique combination of time column values,
#'     written to \code{<output_prefix>_OnePerPixPerTimeStep.csv}.
#' }
#'
#' Output CSV files are suitable as direct input to
#' \code{\link{temporally_explicit_extraction}}.
#'
#' @seealso
#' Preprocessing: \code{\link{temporally_explicit_extraction}},
#'   \code{\link{spatiotemporal_partition}}
#'
#' @examples
#' pts_file <- system.file(
#'   "extdata/points/synthetic_occurrence_points.csv",
#'   package = "TemporalModelR"
#' )
#'
#' ref_file <- system.file("extdata/rasters_raw/elevation.tif",
#'                         package = "TemporalModelR")
#'
#' out_dir  <- file.path(tempdir(), "rarefied")
#'
#' spatiotemporal_rarefaction(
#'   points_sp        = pts_file,
#'   output_dir       = out_dir,
#'   reference_raster = ref_file,
#'   time_cols        = c("year", "season"),
#'   xcol             = "x",
#'   ycol             = "y",
#'   points_crs       = terra::crs(terra::rast(ref_file)),
#'   verbose          = FALSE
#' )

#' @export
#' @importFrom terra rast res ext ncell crs extract vect
#' @importFrom sf st_as_sf st_transform st_coordinates st_read st_drop_geometry st_crs
#' @importFrom tools file_ext
#' @importFrom utils read.csv write.csv
#' @importFrom stats complete.cases
spatiotemporal_rarefaction <- function(points_sp,
                                       output_dir,
                                       reference_raster,
                                       time_cols     = NULL,
                                       xcol          = NULL,
                                       ycol          = NULL,
                                       points_crs    = NULL,
                                       output_prefix = "Pts_Database",
                                       verbose       = TRUE) {


  if (missing(points_sp)) {
    stop(paste0("ERROR: 'points_sp' is required but was not provided. Please provide ",
                "an sf object, data frame, SpatialPointsDataFrame, or file path."))
  }
  if (missing(output_dir) || is.null(output_dir) || output_dir == "") {
    stop(paste0("ERROR: 'output_dir' is required but was not provided. ",
                "Please specify the directory where output files should be saved."))
  }
  if (missing(reference_raster)) {
    stop(paste0("ERROR: 'reference_raster' is required but was not provided. Please specify ",
                "a raster file path or raster object used to define pixel boundaries."))
  }

  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) {
      stop(paste0("ERROR: File does not exist: ", points_sp))
    }

    file_ext_lower <- tolower(tools::file_ext(points_sp))

    if (file_ext_lower == "csv") {
      if (is.null(xcol))       stop("ERROR: 'xcol' is required when reading CSV files.")
      if (is.null(ycol))       stop("ERROR: 'ycol' is required when reading CSV files.")
      if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when reading CSV files.")

      if (verbose) message(paste("Reading CSV file:", basename(points_sp)))
      points_data <- utils::read.csv(points_sp, stringsAsFactors = FALSE)

      if (!xcol %in% names(points_data)) stop(paste0("ERROR: Column '", xcol, "' not found in CSV."))
      if (!ycol %in% names(points_data)) stop(paste0("ERROR: Column '", ycol, "' not found in CSV."))

      points_sp <- sf::st_as_sf(points_data, coords = c(xcol, ycol), crs = points_crs)

    } else if (file_ext_lower %in% c("shp", "geojson", "gpkg")) {
      if (verbose) message(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- sf::st_read(points_sp, quiet = TRUE)
    } else {
      stop(paste0("ERROR: Unsupported file format: '", file_ext_lower,
                  "'. Supported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (inherits(points_sp, "sf")) {
    if (verbose) message("Using provided sf object...")

  } else if (inherits(points_sp, "SpatialPointsDataFrame")) {
    warning("Converting 'SpatialPointsDataFrame' to sf object for terra compatibility.")
    points_sp <- sf::st_as_sf(points_sp)

  } else if (is.data.frame(points_sp)) {
    if (is.null(xcol))       stop("ERROR: 'xcol' is required when providing a data frame.")
    if (is.null(ycol))       stop("ERROR: 'ycol' is required when providing a data frame.")
    if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when providing a data frame.")

    if (verbose) message("Converting data frame to sf object...")
    points_sp <- sf::st_as_sf(points_sp, coords = c(xcol, ycol), crs = points_crs)

  } else {
    stop("ERROR: 'points_sp' must be an sf object, data frame, SpatialPointsDataFrame, or file path.")
  }

  reference_raster <- .load_raster_input(reference_raster, "reference_raster")
  if (is.na(terra::crs(reference_raster)) || nchar(terra::crs(reference_raster)) == 0) {
    stop("ERROR: Reference raster has no defined CRS. Please assign a coordinate reference system before processing.")
  }

  if (sf::st_crs(points_sp) != sf::st_crs(reference_raster)) {
    if (verbose) message("Reprojecting points to match reference raster CRS...")
    points_sp <- sf::st_transform(points_sp, terra::crs(reference_raster))
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(time_cols) || length(time_cols) == 0) {
    if (verbose) message("No time columns provided. Performing spatial-only rarefaction.")
    perform_spatiotemporal <- FALSE
  } else {
    if (!is.character(time_cols)) {
      stop("ERROR: 'time_cols' must be a character vector. Please provide column names as a character vector, e.g. c(\"Year\", \"Month\").")
    }

    missing_cols <- setdiff(time_cols, names(points_sp))
    if (length(missing_cols) > 0) {
      stop(paste0("ERROR: The following time_cols are missing from the input data: ",
                  paste(missing_cols, collapse = ", "),
                  ". Available columns: ",
                  paste(names(points_sp)[names(points_sp) != "geometry"], collapse = ", ")))
    }

    n_original <- nrow(points_sp)
    pts_data      <- sf::st_drop_geometry(points_sp)
    complete_rows <- stats::complete.cases(pts_data[, time_cols, drop = FALSE])
    points_sp     <- points_sp[complete_rows, ]

    n_removed <- n_original - nrow(points_sp)
    if (n_removed > 0) {
      pct_removed <- round(n_removed / n_original * 100, 2)
      if (verbose) message(paste0("Removed ", n_removed, " rows (", pct_removed,
                                  "%) with missing values in time columns"))
    }

    perform_spatiotemporal <- TRUE
  }

  res_x           <- terra::res(reference_raster)[1]
  res_y           <- terra::res(reference_raster)[2]
  pixel_id_raster <- terra::rast(terra::ext(reference_raster), res = c(res_x, res_y))
  terra::crs(pixel_id_raster) <- terra::crs(reference_raster)
  pixel_id_raster[] <- seq_len(terra::ncell(pixel_id_raster))

  points_sp$pixel_id <- terra::extract(pixel_id_raster, terra::vect(points_sp))[, 2]

  coords        <- sf::st_coordinates(points_sp)
  points_sp$X  <- coords[, 1]
  points_sp$Y  <- coords[, 2]

  freq_table           <- as.data.frame(table(points_sp$pixel_id))
  colnames(freq_table) <- c("pixel_id", "Freq")

  if (perform_spatiotemporal) {
    pts_df     <- sf::st_drop_geometry(points_sp)
    group_cols <- c("pixel_id", time_cols)
    pts_df     <- pts_df[!duplicated(pts_df[, group_cols, drop = FALSE]), ]
    points_subset_st <- merge(pts_df, freq_table, by = "pixel_id")

    n_spatiotemporal    <- nrow(points_subset_st)
    cols_to_save        <- c(time_cols, "pixel_id", "X", "Y")
    spatiotemporal_file <- file.path(output_dir,
                                     paste0(output_prefix, "_OnePerPixPerTimeStep.csv"))
    utils::write.csv(points_subset_st[, cols_to_save], spatiotemporal_file, row.names = FALSE)

    time_combinations_n <- nrow(unique(pts_df[, time_cols, drop = FALSE]))
    if (verbose) message(paste("Spatiotemporal file saved:", basename(spatiotemporal_file)))
    if (verbose) message(paste("Retained 1 point per pixel across",
                               time_combinations_n, "unique time combinations"))
  }

  pts_df_spatial        <- sf::st_drop_geometry(points_sp)
  pts_df_spatial        <- pts_df_spatial[!duplicated(pts_df_spatial$pixel_id), ]
  points_subset_spatial <- merge(pts_df_spatial, freq_table, by = "pixel_id")

  n_spatial    <- nrow(points_subset_spatial)
  cols_to_save <- if (perform_spatiotemporal) {
    c(time_cols, "pixel_id", "X", "Y")
  } else {
    c("pixel_id", "X", "Y")
  }

  spatial_file <- file.path(output_dir, paste0(output_prefix, "_OnePerPix.csv"))
  utils::write.csv(points_subset_spatial[, cols_to_save], spatial_file, row.names = FALSE)

  if (verbose) message(paste("Spatial file saved:", basename(spatial_file)))

  if (perform_spatiotemporal) {
    additional_points <- n_spatiotemporal - n_spatial
    pct_additional    <- round((additional_points / n_spatial) * 100, 2)
    if (verbose) message(paste0("Spatial: ", n_spatial, " points | Spatiotemporal: ", n_spatiotemporal,
                                " points | Additional retained: ", additional_points,
                                " (", pct_additional, "% increase)"))
  }

  if (verbose) message("Processing complete.")

  invisible(list(
    input_points          = nrow(points_sp),
    spatial_points        = n_spatial,
    spatiotemporal_points = if (perform_spatiotemporal) n_spatiotemporal else NA,
    time_cols_used        = if (perform_spatiotemporal) time_cols else NULL,
    spatial_table         = points_subset_spatial[, cols_to_save],
    spatiotemporal_table  = if (perform_spatiotemporal) points_subset_st[, cols_to_save] else NULL,
    files_created         = if (perform_spatiotemporal) {
      list(spatiotemporal = spatiotemporal_file, spatial = spatial_file)
    } else {
      list(spatial = spatial_file)
    }
  ))
}
