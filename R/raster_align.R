#' Align and Mask Raster Files to a Reference Raster
#'
#' Aligns a batch of raster files to a specified reference raster by performing
#' reprojection, resampling, and masking operations. All output rasters are
#' written to an output directory with a user-defined filename suffix.
#'
#' @description
#' This utility function standardizes multiple raster files so that they share
#' the same CRS, resolution, and spatial extent as a reference raster. Useful
#' for ensuring that all rasters are properly aligned and don't cause errors later
#' on in modeling or analyses.
#'
#' @param input_dir Character. Directory containing the input raster files.
#' @param output_dir Character. Directory where processed rasters will be saved.
#' @param reference_raster Character or `SpatRaster` or `RasterLayer`.
#'   File path or raster object used as the alignment reference.
#' @param output_suffix Character. Suffix appended to output filenames.
#'   Default is `"_Masked_Updated"`.
#' @param pattern Character. Character pattern used to match raster files within
#'   `input_dir`. Default is `"\\.tif$"`.
#' @param overwrite Logical. If `TRUE`, all files are creating overwriting files
#' with the same name in a given directory. If `FALSE`,
#'   existing processed files are skipped. Default is `FALSE`.
#'
#' @details
#' For each raster in `input_dir` matching `pattern`, the function:
#' \enumerate{
#'   \item Reprojects it to the CRS of the reference raster.
#'   \item Resamples it (bilinear) to match the reference resolution.
#'   \item Masks values outside the reference raster's non-`NA` area.
#'   \item Saves the result as a GeoTIFF to `output_dir`.
#' }
#'
#' @importFrom terra rast project resample mask writeRaster crs ext
#' @importFrom methods is
#' @export

raster_align <- function(input_dir,
                         output_dir,
                         reference_raster,
                         output_suffix = "_Masked_Updated",
                         pattern = ".*\\.tif$",
                         overwrite = FALSE) {

  require(terra)

  ### Validate required inputs

  if (missing(input_dir)) {
    stop("ERROR: 'input_dir' is required but was not provided. Please specify the directory containing input raster files.")
  }

  if (missing(output_dir)) {
    stop("ERROR: 'output_dir' is required but was not provided. Please specify the directory where processed rasters should be saved.")
  }

  if (missing(reference_raster)) {
    stop("ERROR: 'reference_raster' is required but was not provided. Please specify the path to the reference raster file or a raster object that will be used for alignment (projection, extent, and resolution).")
  }

  ### Validate that paths exist and are accessible

  if (!dir.exists(input_dir)) {
    stop(paste0("ERROR: Input directory does not exist: '", input_dir, "'. Please check the path and try again."))
  }

  ### Handle reference_raster as either file path, terra object, or raster object

  if (is.character(reference_raster)) {
    if (!file.exists(reference_raster)) {
      stop(paste0("ERROR: Reference raster file does not exist: '", reference_raster, "'. Please check the file path and try again."))
    }
    reference_raster <- terra::rast(reference_raster)
  } else if (inherits(reference_raster, "RasterLayer")) {
    warning(paste0("WARNING: 'reference_raster' is a 'RasterLayer' object from the raster package.\n",
                   "         Converting to 'SpatRaster' for terra compatibility."))
    reference_raster <- terra::rast(reference_raster)
  } else if (!inherits(reference_raster, "SpatRaster")) {
    stop(paste0("ERROR: 'reference_raster' must be a file path, 'RasterLayer', or 'SpatRaster'. Provided object is of class: ", class(reference_raster)[1]))
  }
  if (is.na(terra::crs(reference_raster))) {
    stop("ERROR: The reference raster has no defined CRS. Please assign a coordinate reference system before processing.")
  }

  ### Replace NA values with zero in reference (for mask logic consistency)
  reference_raster[is.na(reference_raster)] <- 0

  ### Create output directory

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ### List all raster files

  all_tif_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  total_files <- length(all_tif_files)

  if (total_files == 0) {
    stop(paste0("ERROR: No raster files found in input directory matching pattern '", pattern, "'. Please check that the directory contains .tif files or adjust the pattern parameter."))
  }

  print(paste("Found", total_files, "total raster files in directory"))

  ### Determine which files to process based on overwrite setting

  if (overwrite) {
    print("Overwrite mode: ON - Will process and overwrite all files")
    tif_files <- all_tif_files
    files_to_process <- length(tif_files)
  } else {
    print("Overwrite mode: OFF - Will skip files that already exist")

    output_exists <- sapply(all_tif_files, function(f) {
      original_file_name <- basename(f)
      updated_file_name <- sub("\\.tif$", paste0(output_suffix, ".tif"), original_file_name)
      output_path <- file.path(output_dir, updated_file_name)
      file.exists(output_path)
    })

    already_processed <- sum(output_exists)
    percent_processed <- round((already_processed / total_files) * 100, 1)
    print(paste0(already_processed, " files already exist in output directory (", percent_processed, "%)"))

    tif_files <- all_tif_files[!output_exists]
    files_to_process <- length(tif_files)

    if (files_to_process == 0) {
      print("All rasters have already been processed. No new files to process.")
      return(invisible(NULL))
    }
  }

  print(paste("Processing", files_to_process, "files"))

  ### Process each raster

  for (i in seq_along(tif_files)) {
    original_file_name <- basename(tif_files[i])
    updated_file_name <- sub("\\.tif$", paste0(output_suffix, ".tif"), original_file_name)
    output_path <- file.path(output_dir, updated_file_name)

    print(paste0("Processing raster ", i, " of ", files_to_process, ": ", original_file_name))

    r <- terra::rast(tif_files[i])

    if (is.na(terra::crs(r))) {
      stop(paste0("ERROR: Raster '", original_file_name, "' has no defined CRS. Please assign a coordinate reference system to this raster before processing."))
    }

    r <- terra::project(r, terra::crs(reference_raster))
    print("  - Reprojected")

    r <- terra::resample(r, reference_raster, method = "bilinear")
    print("  - Resampled")

    r <- terra::mask(r, reference_raster, maskvalue = 0)
    print("  - Masked")

    terra::writeRaster(r, output_path, overwrite = TRUE)
    print("  - Saved")
  }

  print("Processing complete!")
  return(invisible(NULL))
}

