#' Align and Standardize Raster Files to a Reference Raster
#'
#' Preprocessing function that aligns a batch of raster files to a specified
#' reference raster by performing reprojection, resampling, and masking
#' operations. Ensures all rasters share identical CRS, resolution, and spatial
#' extent for later analyses.
#'
#' @param input_dir Character. Directory containing the input raster files.
#' @param output_dir Character. Directory where processed rasters will be saved.
#' @param reference_raster Character, SpatRaster, or RasterLayer. File path or
#'   raster object used as the alignment reference.
#' @param output_suffix Character. Suffix appended to output filenames. Default
#'   is \code{"_Masked_Updated"}. Output files are named
#'   \code{<original_name><output_suffix>.tif}.
#' @param pattern Character. Regular expression used to match raster files
#'   within \code{input_dir}. Default is \code{".*\\.tif$"}.
#' @param resample_method Character. Resampling method passed to
#'   \code{\link[terra]{resample}}. Default is \code{"bilinear"}, which is
#'   appropriate for continuous variables. Use \code{"near"} for categorical
#'   rasters (e.g. land cover classes) to avoid interpolating between class
#'   codes. Other accepted values include \code{"cubic"} and \code{"lanczos"}.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing output files.
#'   If \code{FALSE} (default), existing files are skipped.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes file counts, overwrite mode,
#'   per-file progress, and a completion summary.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{output_files}: Character vector of file paths written to
#'     \code{output_dir}.
#'   \item \code{n_processed}: Integer. Number of rasters successfully processed
#'     in this call (excludes skipped files when \code{overwrite = FALSE}).
#' }
#'
#' @details
#' For each raster in \code{input_dir} matching \code{pattern}, the function:
#' \enumerate{
#'   \item Reprojects to the CRS of the reference raster
#'   \item Resamples to match the reference resolution using \code{resample_method}
#'   \item Masks values outside the reference raster's non-NA extent
#'   \item Saves the result as a GeoTIFF to \code{output_dir}
#' }
#'
#' This preprocessing step ensures spatial alignment before applying other downstream
#' analyses such as \code{\link{temporally_explicit_extraction}} or \code{\link{scale_rasters}}.
#'
#' @seealso
#' Preprocessing: \code{\link{scale_rasters}},
#'   \code{\link{temporally_explicit_extraction}}
#'
#' @examples
#' \donttest{
#'   raw_dir <- system.file("extdata/rasters_raw", package = "TemporalModelR")
#'
#'   ref     <- file.path(raw_dir, "elevation.tif")
#'
#'   out_dir <- file.path(tempdir(), "aligned")
#'
#'   raster_align(
#'     input_dir        = raw_dir,
#'     output_dir       = out_dir,
#'     reference_raster = ref,
#'     overwrite        = TRUE,
#'     verbose          = FALSE
#'   )
#' }

#' @export
#' @importFrom terra rast project resample mask writeRaster crs classify
raster_align <- function(input_dir,
                         output_dir,
                         reference_raster,
                         output_suffix   = "_Masked_Updated",
                         pattern         = ".*\\.tif$",
                         resample_method = "bilinear",
                         overwrite       = FALSE,
                         verbose         = TRUE) {

  if (missing(input_dir)) {
    stop(paste0("ERROR: 'input_dir' is required but was not provided. ",
                "Please specify the directory containing input raster files."))
  }
  if (missing(output_dir)) {
    stop(paste0("ERROR: 'output_dir' is required but was not provided. ",
                "Please specify the directory where processed rasters should be saved."))
  }
  if (missing(reference_raster)) {
    stop(paste0("ERROR: 'reference_raster' is required but was not provided. Please specify ",
                "the path to the reference raster file or a raster object used for alignment."))
  }
  if (!dir.exists(input_dir)) {
    stop(paste0("ERROR: Input directory does not exist: '", input_dir,
                "'. Please check the path and try again."))
  }
  if (!resample_method %in% c("bilinear", "near", "cubic", "lanczos",
                              "cubicspline", "average", "mode")) {
    stop(paste0("ERROR: 'resample_method' must be one of 'bilinear', 'near', 'cubic', 'lanczos', ",
                "'cubicspline', 'average', or 'mode'. Please adjust resample_method and try again."))
  }

  reference_raster <- .load_raster_input(reference_raster, "reference_raster")
  if (is.na(terra::crs(reference_raster)) || nchar(terra::crs(reference_raster)) == 0) {
    stop("ERROR: Reference raster has no defined CRS. Please assign a coordinate reference system before processing.")
  }

  reference_raster <- terra::classify(reference_raster, cbind(NA, 0))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  all_tif_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  total_files   <- length(all_tif_files)

  if (total_files == 0) {
    stop(paste0("ERROR: No raster files found in '", input_dir,
                "' matching pattern '", pattern, "'. Please check that the directory contains .tif files or adjust the pattern parameter."))
  }

  if (verbose) message(paste("Found", total_files, "raster files in directory"))

  if (overwrite) {
    if (verbose) message("Overwrite mode: ON - will process and overwrite all files")
    tif_files        <- all_tif_files
    files_to_process <- total_files
  } else {
    if (verbose) message("Overwrite mode: OFF - will skip files that already exist")

    output_exists <- vapply(all_tif_files, function(f) {
      out_name <- sub("\\.tif$", paste0(output_suffix, ".tif"), basename(f))
      file.exists(file.path(output_dir, out_name))
    }, logical(1))

    already_processed <- sum(output_exists)
    percent_processed <- round((already_processed / total_files) * 100, 1)
    if (verbose) message(paste0(already_processed, " of ", total_files,
                                " files already exist in output directory (",
                                percent_processed, "%)"))

    tif_files        <- all_tif_files[!output_exists]
    files_to_process <- length(tif_files)

    if (files_to_process == 0) {
      if (verbose) message("All rasters already processed. No new files to process.")
      return(invisible(list(output_files = character(0), n_processed = 0)))
    }
  }

  if (verbose) message(paste("Processing", files_to_process, "files"))

  output_paths <- character(files_to_process)

  for (i in seq_along(tif_files)) {
    original_name <- basename(tif_files[i])
    output_name   <- sub("\\.tif$", paste0(output_suffix, ".tif"), original_name)
    output_path   <- file.path(output_dir, output_name)

    if (verbose) message(paste0("  [", i, "/", files_to_process, "] ", original_name))

    r <- terra::rast(tif_files[i])

    if (is.na(terra::crs(r))) {
      stop(paste0("ERROR: Raster '", original_name, "' has no defined CRS. Please assign a coordinate reference system to this raster before processing."))
    }

    r <- terra::project(r, terra::crs(reference_raster))
    r <- terra::resample(r, reference_raster, method = resample_method)
    r <- terra::mask(r, reference_raster, maskvalue = 0)

    terra::writeRaster(r, output_path, overwrite = TRUE)
    output_paths[i] <- output_path
  }

  if (verbose) message("Processing complete.")
  invisible(list(
    output_files = output_paths,
    n_processed  = files_to_process
  ))
}
