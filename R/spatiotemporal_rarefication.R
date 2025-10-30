#' Spatiotemporal rarefaction
#' @export
#' @importFrom raster raster res extent ncell crs extract
#' @importFrom sp SpatialPointsDataFrame spTransform coordinates
#' @importFrom dplyr group_by across all_of slice ungroup
#' @importFrom utils write.csv
spatiotemporal_rarefication <- function(points_sp,
                                        output_dir,
                                        reference_raster,
                                        time_cols = NULL,
                                        xcol = NULL,
                                        ycol = NULL,
                                        points_crs = NULL,
                                        create_spatial_only = FALSE,
                                        output_prefix = "Pts_Database") {
  require(raster)
  require(sp)
  require(dplyr)

  ### INPUT VALIDATION AND CONVERSION - POINTS

  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) {
      stop(paste("Error: File does not exist:", points_sp))
    }

    file_ext <- tolower(tools::file_ext(points_sp))

    if (file_ext == "csv") {

      # Require xcol, ycol, and points_crs for CSV files
      if (is.null(xcol)) {
        stop("Error: 'xcol' is required when reading CSV files.")
      }
      if (is.null(ycol)) {
        stop("Error: 'ycol' is required when reading CSV files.")
      }
      if (is.null(points_crs)) {
        stop("Error: 'points_crs' is required when reading CSV files.")
      }

      print(paste("Reading CSV file:", basename(points_sp)))
      points_data <- read.csv(points_sp, stringsAsFactors = FALSE)

      if (!xcol %in% names(points_data)) {
        stop(paste0("Error: Column '", xcol, "' not found in CSV. Available columns: ",
                    paste(names(points_data), collapse = ", ")))
      }
      if (!ycol %in% names(points_data)) {
        stop(paste0("Error: Column '", ycol, "' not found in CSV. Available columns: ",
                    paste(names(points_data), collapse = ", ")))
      }

      coordinates(points_data) <- c(xcol, ycol)
      proj4string(points_data) <- CRS(points_crs)
      points_sp <- points_data

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- rgdal::readOGR(points_sp, verbose = FALSE)
    } else {
      stop(paste("Error: Unsupported file format:", file_ext,
                 "\nSupported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (is.data.frame(points_sp)) {

    # Require xcol, ycol, and points_crs for data frames
    if (is.null(xcol)) {
      stop("Error: 'xcol' is required when providing a data frame.")
    }
    if (is.null(ycol)) {
      stop("Error: 'ycol' is required when providing a data frame.")
    }
    if (is.null(points_crs)) {
      stop("Error: 'points_crs' is required when providing a data frame.")
    }

    print("Converting data frame to SpatialPointsDataFrame...")

    if (!xcol %in% names(points_sp)) {
      stop(paste0("Error: Column '", xcol, "' not found in data frame. Available columns: ",
                  paste(names(points_sp), collapse = ", ")))
    }
    if (!ycol %in% names(points_sp)) {
      stop(paste0("Error: Column '", ycol, "' not found in data frame. Available columns: ",
                  paste(names(points_sp), collapse = ", ")))
    }

    coordinates(points_sp) <- c(xcol, ycol)
    proj4string(points_sp) <- CRS(points_crs)

  } else if (!inherits(points_sp, "SpatialPointsDataFrame")) {
    stop("Error: points_sp must be a SpatialPointsDataFrame, data frame with x/y columns, or file path")
  }

  ### INPUT VALIDATION AND CONVERSION - RASTER

  if (inherits(reference_raster, "RasterLayer")) {
    binary_mask <- reference_raster
  } else {
    if (!file.exists(reference_raster)) {
      stop(paste("Error: Reference raster file not found at path:", reference_raster))
    }
    binary_mask <- raster(reference_raster)
  }

  ### VALIDATE OUTPUT DIRECTORY

  if (is.null(output_dir) || output_dir == "") {
    stop("Error: output_dir must be provided and cannot be empty")
  }

  ### VALIDATE TIME COLUMNS

  if (is.null(time_cols) || length(time_cols) == 0) {
    warning("No time columns provided. Defaulting to spatial-only rarefaction.")
    create_spatial_only <- TRUE
    perform_spatiotemporal <- FALSE
  } else {
    missing_cols <- setdiff(time_cols, names(points_sp@data))
    if (length(missing_cols) > 0) {
      stop(paste("Error: The following time columns were not found in the data:",
                 paste(missing_cols, collapse = ", ")))
    }

    # Check for missing values in time columns
    n_original <- nrow(points_sp@data)
    for (tc in time_cols) {
      n_missing <- sum(is.na(points_sp@data[[tc]]))
      if (n_missing > 0) {
        pct_missing <- round(n_missing / n_original * 100, 2)
        warning(paste0("Warning: ", n_missing, " rows (", pct_missing,
                       "%) have missing values in time column '", tc, "'"))
      }
    }

    perform_spatiotemporal <- TRUE
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ### PROCESS COORDINATES AND PIXEL IDs

  coords <- coordinates(points_sp)
  coord_names <- colnames(coords)
  coord_x <- coord_names[1]
  coord_y <- coord_names[2]

  resolution_x <- res(binary_mask)[1]
  resolution_y <- res(binary_mask)[2]
  pixel_id_raster <- raster(extent(binary_mask), res = c(resolution_x, resolution_y))
  pixel_id_raster[] <- 1:ncell(pixel_id_raster)

  points_sp <- spTransform(points_sp, crs(binary_mask))
  points_sp$pixel_id <- raster::extract(pixel_id_raster, points_sp)

  coords_transformed <- coordinates(points_sp)
  points_sp@data[[coord_x]] <- coords_transformed[, 1]
  points_sp@data[[coord_y]] <- coords_transformed[, 2]

  Freq_Table <- as.data.frame(table(points_sp$pixel_id))
  colnames(Freq_Table) <- c("pixel_id", "Freq")

  ### SPATIOTEMPORAL RAREFACTION

  if (perform_spatiotemporal) {
    points_subset <- points_sp@data %>%
      group_by(pixel_id, across(all_of(time_cols))) %>%
      slice(1) %>%
      ungroup() %>%
      merge(Freq_Table, by = "pixel_id")

    points_sp_subset <- SpatialPointsDataFrame(
      points_subset[, c(coord_x, coord_y)],
      points_subset,
      proj4string = crs(binary_mask)
    )

    n_spatiotemporal <- nrow(points_subset)

    cols_to_save <- c(time_cols, "pixel_id", coord_x, coord_y)

    spatiotemporal_file <- paste0(output_dir, "/", output_prefix, "_OnePerPixPerTimeStep.csv")
    write.csv(
      points_sp_subset@data[, cols_to_save],
      spatiotemporal_file,
      row.names = FALSE
    )

    print(paste("Spatiotemporal file saved:", basename(spatiotemporal_file)))
  }

  ### SPATIAL-ONLY RAREFACTION

  if (create_spatial_only || !perform_spatiotemporal) {
    points_subset_perpixel <- points_sp@data %>%
      group_by(pixel_id) %>%
      slice(1) %>%
      ungroup() %>%
      merge(Freq_Table, by = "pixel_id")

    n_spatial <- nrow(points_subset_perpixel)

    if (perform_spatiotemporal) {
      cols_to_save <- c(time_cols, "pixel_id", coord_x, coord_y)
    } else {
      cols_to_save <- c("pixel_id", coord_x, coord_y)
    }

    spatial_file <- paste0(output_dir, "/", output_prefix, "_OnePerPix.csv")
    write.csv(
      points_subset_perpixel[, cols_to_save],
      spatial_file,
      row.names = FALSE
    )

    print(paste("Spatial file saved:", basename(spatial_file)))

    if (perform_spatiotemporal) {
      additional_points <- n_spatiotemporal - n_spatial
      pct_additional <- round((additional_points / n_spatial) * 100, 2)
      print(paste0("Spatial: ", n_spatial, " points | Spatiotemporal: ", n_spatiotemporal,
                   " points | Additional retained: ", additional_points,
                   " (", pct_additional, "% increase)"))
    }
  }

  ### SUMMARY

  summary_info <- list(
    input_points = nrow(points_sp@data),
    spatial_points = if (create_spatial_only || !perform_spatiotemporal) n_spatial else NA,
    spatiotemporal_points = if (perform_spatiotemporal) n_spatiotemporal else NA,
    files_created = list()
  )

  if (perform_spatiotemporal) {
    summary_info$files_created$spatiotemporal <- spatiotemporal_file
  }
  if (create_spatial_only || !perform_spatiotemporal) {
    summary_info$files_created$spatial <- spatial_file
  }

  print("Processing complete!")

  invisible(summary_info)
}
