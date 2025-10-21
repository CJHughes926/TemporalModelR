#' Spatiotemporal rarefaction
#' @export
#' @importFrom raster raster res extent ncell crs extract
#' @importFrom sp SpatialPointsDataFrame spTransform coordinates
#' @importFrom dplyr group_by across all_of slice ungroup
#' @importFrom utils write.csv
spatiotemporal_rarefication <- function(points_sp,
                                        output_dir,
                                        reference_raster,
                                        time_cols,
                                        create_spatial_only = FALSE,
                                        output_prefix = "Pts_Database") {

  # Load required libraries
  require(raster)
  require(sp)
  require(dplyr)

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load reference raster
  binary_mask <- raster(reference_raster)

  print(paste("Total points in input:", nrow(points_sp@data)))

  # Get coordinate column names
  coords <- coordinates(points_sp)
  coord_names <- colnames(coords)
  coord_x <- coord_names[1]
  coord_y <- coord_names[2]

  # Create pixel ID raster
  resolution_x <- res(binary_mask)[1]
  resolution_y <- res(binary_mask)[2]
  pixel_id_raster <- raster(extent(binary_mask), res = c(resolution_x, resolution_y))
  pixel_id_raster[] <- 1:ncell(pixel_id_raster)

  # Transform and extract pixel IDs
  points_sp <- spTransform(points_sp, crs(binary_mask))
  points_sp$pixel_id <- raster::extract(pixel_id_raster, points_sp)

  # Add coordinates to data slot for easier processing
  coords_transformed <- coordinates(points_sp)
  points_sp@data[[coord_x]] <- coords_transformed[, 1]
  points_sp@data[[coord_y]] <- coords_transformed[, 2]

  # Generate frequency table
  Freq_Table <- as.data.frame(table(points_sp$pixel_id))
  colnames(Freq_Table) <- c("pixel_id", "Freq")

  # Keep one point per pixel per time step(s)
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

  # Print summary
  print(paste("Unfiltered points:", nrow(points_sp@data)))
  print(paste("Filtered points (one per pixel per time step):", nrow(points_subset)))

  # Prepare columns to save (time cols + pixel_id + coordinates)
  cols_to_save <- c(time_cols, "pixel_id", coord_x, coord_y)

  # Save results - one per pixel per time step
  write.csv(
    points_sp_subset@data[, cols_to_save],
    paste0(output_dir, "/", output_prefix, "_OnePerPixPerTimeStep.csv"),
    row.names = FALSE
  )
  print(paste("Saved:", paste0(output_dir, "/", output_prefix, "_OnePerPixPerTimeStep.csv")))

  # Create one-per-pixel subset (ignoring time) if requested
  if (create_spatial_only) {
    points_subset_perpixel <- points_sp@data %>%
      group_by(pixel_id) %>%
      slice(1) %>%
      ungroup() %>%
      merge(Freq_Table, by = "pixel_id")

    write.csv(
      points_subset_perpixel[, cols_to_save],
      paste0(output_dir, "/", output_prefix, "_OnePerPix.csv"),
      row.names = FALSE
    )
    print(paste("Saved:", paste0(output_dir, "/", output_prefix, "_OnePerPix.csv")))
  } else {
    print("Skipped creating spatial-only dataset (create_spatial_only = FALSE)")
  }
}
