#' Align and Mask Raster Files to a Reference Raster
#'
#' This utility function aligns a batch of raster files within a directory to a
#' given reference raster. It performs reprojection, resampling, and masking
#' operations so that all rasters share the same CRS, resolution, and spatial
#' extent. The processed rasters are saved to an output directory with a
#' specified suffix.
#'
#' @param input_dir Character. Path to the directory containing input raster files.
#' @param output_dir Character. Path to the output directory where processed rasters will be saved.
#' @param reference_raster Character or `RasterLayer`. File path or object of the reference raster
#'   used for reprojection, resampling, and masking.
#' @param output_suffix Character. Suffix appended to output raster filenames
#'   (default = "_Masked_Updated").
#' @param pattern Character. Regex pattern to match raster files in `input_dir`
#'   (default = ".*\\.tif$").
#' @param overwrite Logical. If `TRUE`, existing output files will be overwritten.
#'   If `FALSE`, already processed files will be skipped (default = FALSE).
#'
#' @details
#' The function iterates through all raster files matching the pattern in
#' `input_dir`. Each raster is reprojected to the CRS of the reference raster,
#' resampled using bilinear interpolation, and masked using the reference
#' raster’s non-NA area. The final rasters are written to `output_dir` in
#' GeoTIFF format.
#'
#' This is especially useful when working with multi-source raster datasets that
#' must be aligned before stacking or extracting values.
#'
#' @return Invisibly returns `NULL`. Output rasters are saved to disk.
#'
#' @importFrom raster raster projectRaster resample mask writeRaster crs
#' @export
raster_align <- function(input_dir,
                         output_dir,
                         reference_raster,
                         output_suffix = "_Masked_Updated",
                         pattern = ".*\\.tif$",
                         overwrite = FALSE) {
  # Load required library
  require(raster)

  # Check for required inputs
  if (missing(input_dir)) {
    stop("ERROR: 'input_dir' is required but was not provided. ",
         "Please specify the directory containing input raster files.")
  }

  if (missing(output_dir)) {
    stop("ERROR: 'output_dir' is required but was not provided. ",
         "Please specify the directory where processed rasters should be saved.")
  }

  if (missing(reference_raster)) {
    stop("ERROR: 'reference_raster' is required but was not provided. ",
         "Please specify the path to the reference raster file that will be used ",
         "for alignment (projection, extent, and resolution).")
  }

  # Validate that paths exist and are accessible
  if (!dir.exists(input_dir)) {
    stop("ERROR: Input directory does not exist: '", input_dir, "'. ",
         "Please check the path and try again.")
  }

  if (!file.exists(reference_raster)) {
    stop("ERROR: Reference raster file does not exist: '", reference_raster, "'. ",
         "Please check the file path and try again.")
  }

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load reference raster
  reference_raster <- raster(reference_raster)

  # List all raster files
  all_tif_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  total_files <- length(all_tif_files)

  if (total_files == 0) {
    stop("ERROR: No raster files found in input directory matching pattern '", pattern, "'. ",
         "Please check that the directory contains .tif files or adjust the pattern parameter.")
  }

  print(paste("Found", total_files, "total raster files in directory"))

  # Determine which files to process based on overwrite setting
  if (overwrite) {
    print("Overwrite mode: ON - Will process and overwrite all files")
    tif_files <- all_tif_files
    files_to_process <- length(tif_files)
  } else {
    print("Overwrite mode: OFF - Will skip files that already exist")

    # Check which files already have outputs
    output_exists <- sapply(all_tif_files, function(f) {
      original_file_name <- basename(f)
      updated_file_name <- sub("\\.tif$", paste0(output_suffix, ".tif"), original_file_name)
      output_path <- file.path(output_dir, updated_file_name)
      file.exists(output_path)
    })

    already_processed <- sum(output_exists)
    percent_processed <- round((already_processed / total_files) * 100, 1)
    print(paste(already_processed, " files already exist in output directory (", percent_processed, "%)", sep = ""))

    # Filter to only files that need processing
    tif_files <- all_tif_files[!output_exists]
    files_to_process <- length(tif_files)

    if (files_to_process == 0) {
      print("All rasters have already been processed. No new files to process.")
      return(invisible(NULL))
    }
  }

  print(paste("Processing", files_to_process, "files"))

  # Process each raster
  for (i in 1:length(tif_files)) {
    original_file_name <- basename(tif_files[i])
    updated_file_name <- sub("\\.tif$", paste0(output_suffix, ".tif"), original_file_name)
    output_path <- file.path(output_dir, updated_file_name)

    print(paste0("Processing raster ", i, " of ", files_to_process, ": ", original_file_name))

    r <- raster(tif_files[i])

    if (is.na(crs(r))) {
      stop(paste0("ERROR: Raster '", original_file_name, "' has no defined CRS. ",
                  "Please assign a coordinate reference system to this raster before processing. ",
                  "You can use raster::crs() to set the CRS if you know what it should be."))
    }

    r <- projectRaster(r, crs = crs(reference_raster))
    print("  - Reprojected")

    r <- resample(r, reference_raster, method = "bilinear")
    print("  - Resampled")

    r <- mask(r, reference_raster, maskvalue = 0)
    print("  - Masked")

    writeRaster(r, output_path, format = "GTiff", overwrite = TRUE)
    print("  - Saved")
  }

  print("Processing complete!")
  return(invisible(NULL))
}
