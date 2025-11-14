#' Spatiotemporal Rarefication of Species Occurrence Data
#'
#' Preprocessing function that rarefies species occurrence data to one point per
#' raster pixel, optionally accounting for temporal components. Reduces sampling
#' bias and spatial autocorrelation in occurrence datasets.
#'
#' @param points_sp Input point data. May be an sf object, data frame,
#'   SpatialPointsDataFrame, or file path to .csv, .shp, .geojson, or .gpkg.
#' @param output_dir Character. Directory where output files will be saved.
#' @param reference_raster Character or raster object. Raster used to define
#'   pixel boundaries for rarefication. May be a file path, RasterLayer, or
#'   SpatRaster.
#' @param time_cols Character vector. Column names defining temporal grouping
#'   for spatiotemporal rarefication. Default is NULL for spatial-only
#'   rarefication.
#' @param xcol Character. Name of x-coordinate column when input is a data frame
#'   or CSV. Required for tabular inputs.
#' @param ycol Character. Name of y-coordinate column when input is a data frame
#'   or CSV. Required for tabular inputs.
#' @param points_crs Character or CRS object. CRS of the input points. Required
#'   when coordinates are supplied in tabular form.
#' @param output_prefix Character. Prefix for output file names. Default is
#'   "Pts_Database".
#'
#' @return Invisibly returns a list with summary information. CSV files are
#'   written to output_dir containing rarefied occurrence data.
#'
#' @details
#' The function assigns each point to a raster pixel and performs:
#' \itemize{
#'   \item Spatial rarefication: retains one point per pixel
#'   \item Spatiotemporal rarefication (if time_cols provided): retains one
#'     point per pixel per time combination
#' }
#'
#' Output files include spatial-only and optionally spatiotemporal rarefied
#' datasets suitable for \code{\link{temporally_explicit_extraction}}.
#'
#' @seealso
#' Preprocessing: \code{\link{temporally_explicit_extraction}},
#' \code{\link{spatiotemporal_partition}}
#'
#' @examples
#' \dontrun{
#' spatiotemporal_rarefication(
#'   points_sp = "occurrences.csv",
#'   output_dir = "rarefied_data/",
#'   reference_raster = "reference.tif",
#'   time_cols = c("Year", "Month"),
#'   xcol = "longitude",
#'   ycol = "latitude",
#'   points_crs = "EPSG:4326"
#' )
#' }
#'
#' @export
#' @importFrom terra rast res extent ncell crs extract vect
#' @importFrom sf st_as_sf st_transform st_coordinates st_read st_drop_geometry
#' @importFrom dplyr group_by across all_of slice ungroup select distinct
#' @importFrom utils write.csv
spatiotemporal_rarefication <- function(points_sp,
                                        output_dir,
                                        reference_raster,
                                        time_cols = NULL,
                                        xcol = NULL,
                                        ycol = NULL,
                                        points_crs = NULL,
                                        output_prefix = "Pts_Database") {

  require(terra)
  require(sf)
  require(dplyr)

  ### Input validation and conversion - points
  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) {
      stop(paste0("ERROR: File does not exist: ", points_sp))
    }

    file_ext <- tolower(tools::file_ext(points_sp))

    if (file_ext == "csv") {
      if (is.null(xcol)) stop("ERROR: 'xcol' is required when reading CSV files.")
      if (is.null(ycol)) stop("ERROR: 'ycol' is required when reading CSV files.")
      if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when reading CSV files.")

      print(paste("Reading CSV file:", basename(points_sp)))
      points_data <- read.csv(points_sp, stringsAsFactors = FALSE)

      if (!xcol %in% names(points_data)) stop(paste0("ERROR: Column '", xcol, "' not found in CSV."))
      if (!ycol %in% names(points_data)) stop(paste0("ERROR: Column '", ycol, "' not found in CSV."))

      points_sp <- st_as_sf(points_data, coords = c(xcol, ycol), crs = points_crs)

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- st_read(points_sp, quiet = TRUE)
    } else {
      stop(paste0("ERROR: Unsupported file format: ", file_ext, ". Supported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (is.data.frame(points_sp)) {
    if (is.null(xcol)) stop("ERROR: 'xcol' is required when providing a data frame.")
    if (is.null(ycol)) stop("ERROR: 'ycol' is required when providing a data frame.")
    if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when providing a data frame.")

    print("Converting data frame to sf object...")
    points_sp <- st_as_sf(points_sp, coords = c(xcol, ycol), crs = points_crs)

  } else if (inherits(points_sp, "SpatialPointsDataFrame")) {
    warning("Converting 'SpatialPointsDataFrame' to sf object for terra compatibility.")
    points_sp <- st_as_sf(points_sp)
  } else if (!inherits(points_sp, "sf")) {
    stop("ERROR: points_sp must be an sf object, data frame, SpatialPointsDataFrame, or file path.")
  }

  ### Input validation and conversion - raster
  if (is.character(reference_raster)) {
    if (!file.exists(reference_raster)) {
      stop(paste0("ERROR: Reference raster file not found at path: ", reference_raster))
    }
    binary_mask <- rast(reference_raster)
  } else if (inherits(reference_raster, "RasterLayer")) {
    warning("Converting 'RasterLayer' to 'SpatRaster' for terra compatibility.")
    binary_mask <- rast(reference_raster)
  } else if (inherits(reference_raster, "SpatRaster")) {
    binary_mask <- reference_raster
  } else {
    stop("ERROR: reference_raster must be a file path, RasterLayer, or SpatRaster object.")
  }

  ### Check CRS of raster
  if (is.na(terra::crs(binary_mask))) {
    stop("ERROR: Reference raster has no defined CRS. Please assign a coordinate reference system before processing.")
  }

  ### Reproject points to match raster CRS if needed
  if (st_crs(points_sp) != st_crs(binary_mask)) {
    print("Reprojecting points to match reference raster CRS...")
    points_sp <- st_transform(points_sp, crs(binary_mask))
  }

  ### Validate output directory
  if (is.null(output_dir) || output_dir == "") {
    stop("ERROR: output_dir must be provided and cannot be empty.")
  }

  ### Validate time columns
  if (is.null(time_cols) || length(time_cols) == 0) {
    print("No time columns provided. Performing spatial-only rarefaction.")
    perform_spatiotemporal <- FALSE
  } else {
    if (!is.character(time_cols)) stop("ERROR: time_cols must be a character vector.")

    missing_cols <- setdiff(time_cols, names(points_sp))
    if (length(missing_cols) > 0) {
      stop(paste0("ERROR: The following time_cols are missing from the input data: ",
                  paste(missing_cols, collapse = ", "),
                  " Available columns: ",
                  paste(names(points_sp), collapse = ", ")))
    }

    n_original <- nrow(points_sp)
    for (tc in time_cols) {
      n_missing <- sum(is.na(points_sp[[tc]]))
      if (n_missing > 0) {
        pct_missing <- round(n_missing / n_original * 100, 2)
        warning(paste0("WARNING: ", n_missing, " rows (", pct_missing,
                       "%) have missing values in time column '", tc, "'"))
      }
    }

    perform_spatiotemporal <- TRUE
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ### Process coordinates and pixel IDs
  resolution_x <- res(binary_mask)[1]
  resolution_y <- res(binary_mask)[2]
  pixel_id_raster <- rast(ext(binary_mask), res = c(resolution_x, resolution_y))
  pixel_id_raster[] <- 1:ncell(pixel_id_raster)

  # Extract pixel IDs using terra::extract (sf compatible)
  points_sp$pixel_id <- terra::extract(pixel_id_raster, vect(points_sp))[, 2]

  coords <- st_coordinates(points_sp)
  points_sp$X <- coords[, 1]
  points_sp$Y <- coords[, 2]

  Freq_Table <- as.data.frame(table(points_sp$pixel_id))
  colnames(Freq_Table) <- c("pixel_id", "Freq")

  ### Spatiotemporal rarefaction
  if (perform_spatiotemporal) {
    points_subset <- points_sp %>%
      group_by(pixel_id, across(all_of(time_cols))) %>%
      slice(1) %>%
      ungroup() %>%
      merge(Freq_Table, by = "pixel_id")

    n_spatiotemporal <- nrow(points_subset)
    cols_to_save <- c(time_cols, "pixel_id", "X", "Y")

    spatiotemporal_file <- file.path(output_dir, paste0(output_prefix, "_OnePerPixPerTimeStep.csv"))
    write.csv(st_drop_geometry(points_subset)[, cols_to_save], spatiotemporal_file, row.names = FALSE)
    print(paste("Spatiotemporal file saved:", basename(spatiotemporal_file)))

    time_combinations <- points_subset %>%
      select(all_of(time_cols)) %>%
      distinct() %>%
      nrow()
    print(paste("  Retained 1 point per pixel across", time_combinations, "unique time combinations"))
  }

  ### Spatial-only rarefaction
  points_subset_perpixel <- points_sp %>%
    group_by(pixel_id) %>%
    slice(1) %>%
    ungroup() %>%
    merge(Freq_Table, by = "pixel_id")

  n_spatial <- nrow(points_subset_perpixel)
  cols_to_save <- if (perform_spatiotemporal) c(time_cols, "pixel_id", "X", "Y") else c("pixel_id", "X", "Y")

  spatial_file <- file.path(output_dir, paste0(output_prefix, "_OnePerPix.csv"))
  write.csv(st_drop_geometry(points_subset_perpixel)[, cols_to_save], spatial_file, row.names = FALSE)
  print(paste("Spatial file saved:", basename(spatial_file)))

  if (perform_spatiotemporal) {
    additional_points <- n_spatiotemporal - n_spatial
    pct_additional <- round((additional_points / n_spatial) * 100, 2)
    print(paste0("Spatial: ", n_spatial, " points | Spatiotemporal: ", n_spatiotemporal,
                 " points | Additional retained: ", additional_points,
                 " (", pct_additional, "% increase)"))
  }

  ### Summary
  summary_info <- list(
    input_points = nrow(points_sp),
    spatiotemporal_points = if (perform_spatiotemporal) n_spatiotemporal else NA,
    time_cols_used = if (perform_spatiotemporal) time_cols else NULL,
    files_created = list()
  )

  if (perform_spatiotemporal) summary_info$files_created$spatiotemporal <- spatiotemporal_file
  summary_info$files_created$spatial <- spatial_file

  print("Processing complete!")
  invisible(summary_info)
}
