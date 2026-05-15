#' Summarize Prediction Rasters into Consensus Outputs
#'
#' Postprocessing function that synthesizes per-time-step fold-vote rasters
#' (output from \code{\link{generate_spatiotemporal_predictions}}) into binary
#' consensus predictions and a temporal frequency summary. Each input raster
#' contains integer vote counts per pixel the number of cross-validation
#' folds that classified that pixel as suitable. The \code{consensus} threshold
#' controls how many folds must agree for a pixel to be classified as suitable
#' in the output binary rasters.
#'
#' @param predictions_dir Character. Directory containing prediction raster
#'   files, typically the \code{output_dir} used in
#'   \code{\link{generate_spatiotemporal_predictions}}.
#' @param output_dir Character. Directory for output files. Defaults to
#'   \code{predictions_dir} if \code{NULL}.
#' @param consensus Integer. Minimum number of folds that must agree on
#'   suitability for a pixel to be classified as suitable in the binary output.
#'   For example, \code{consensus = 1} marks any pixel suitable if at least one
#'   fold predicts it suitable; \code{consensus = 4} requires at least four folds to
#'   agree. Must be between 1 and the maximum vote count in the rasters (number of
#'   folds). Default is \code{1}.
#' @param file_pattern Character. Regular expression to match prediction raster
#'   files. Default is \code{"Prediction_.*\\.tif$"}.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing output files.
#'   If \code{FALSE} (default), existing files are skipped.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{consensus_stack}: \code{SpatRaster} stack of per-time-step
#'     binary consensus rasters (one layer per input file).
#'   \item \code{frequency_raster}: \code{SpatRaster} showing the proportion of
#'     time steps during which each pixel met the consensus threshold. Values
#'     range from 0 (never suitable) to 1 (always suitable).
#'   \item \code{consensus}: the consensus threshold used.
#'   \item \code{n_timesteps}: number of time steps processed.
#'   \item \code{consensus_dir}: path to the directory containing per-time-step
#'     binary consensus rasters.
#'   \item \code{frequency_file}: path to the written frequency raster file.
#' }
#'
#' @details
#' Input rasters are fold-vote-count rasters produced by
#' \code{\link{generate_spatiotemporal_predictions}}, where each pixel value is
#' the number of cross-validation fold models that predicted suitability at that
#' location for that time step. The \code{consensus} threshold is applied as:
#' \code{binary = as.integer(vote_count >= consensus)}.
#'
#' The frequency raster (proportion of time steps suitable under the chosen
#' consensus) serves as input to \code{\link[TemporalModelR]{analyze_temporal_patterns}} for
#' identifying long-term trends in suitability.
#'
#' @seealso
#' Upstream: \code{\link{generate_spatiotemporal_predictions}}
#'
#' Downstream: \code{\link[TemporalModelR]{analyze_temporal_patterns}}
#'
#' @examples
#' \donttest{
#'   pred_dir <- system.file("extdata/predictions",
#'                           package = "TemporalModelR")
#'
#'   summarize_raster_outputs(
#'     predictions_dir = pred_dir,
#'     output_dir      = tempdir(),
#'     consensus       = 3,
#'     overwrite       = TRUE,
#'     verbose         = FALSE
#'   )
#' }

#' @export
#' @importFrom terra rast nlyr app writeRaster values global
summarize_raster_outputs <- function(predictions_dir,
                                     output_dir   = NULL,
                                     consensus    = 1,
                                     file_pattern = "Prediction_.*\\.tif$",
                                     overwrite    = FALSE,
                                     verbose      = TRUE) {

  if (missing(predictions_dir) || is.null(predictions_dir) || predictions_dir == "") {
    stop("ERROR: 'predictions_dir' is required and cannot be NULL or empty.")
  }
  if (!dir.exists(predictions_dir)) {
    stop(paste0("ERROR: predictions_dir does not exist: ", predictions_dir))
  }
  if (!is.numeric(consensus) || length(consensus) != 1 || consensus < 1) {
    stop("ERROR: 'consensus' must be a single positive integer >= 1.")
  }
  consensus <- as.integer(consensus)

  prediction_files <- list.files(
    path       = predictions_dir,
    pattern    = file_pattern,
    full.names = TRUE
  )

  if (length(prediction_files) == 0) {
    stop(paste0("ERROR: No files matching '", file_pattern,
                "' found in: ", predictions_dir))
  }

  prediction_files <- .natural_sort(prediction_files)
  if (verbose) message(paste("Found", length(prediction_files), "prediction rasters."))

  if (is.null(output_dir)) output_dir <- predictions_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE,
                                          showWarnings = FALSE)

  consensus_dir  <- file.path(output_dir,
                              paste0("Consensus_", consensus, "_Rasters"))
  frequency_file <- file.path(output_dir,
                              paste0("Frequency_Consensus_", consensus, ".tif"))

  if (dir.exists(consensus_dir) && file.exists(frequency_file) && !overwrite) {
    if (verbose) message("Outputs already exist and overwrite = FALSE. Reloading from disk.")
    con_files <- .natural_sort(list.files(consensus_dir, pattern = "\\.tif$",
                                          full.names = TRUE))
    con_stack      <- terra::rast(con_files)
    frequency_rast <- terra::rast(frequency_file)
    return(invisible(list(
      consensus_stack  = con_stack,
      frequency_raster = frequency_rast,
      consensus        = consensus,
      n_timesteps      = terra::nlyr(con_stack),
      consensus_dir    = consensus_dir,
      frequency_file   = frequency_file
    )))
  }

  if (verbose) message("Loading prediction rasters...")
  vote_stack  <- terra::rast(prediction_files)
  clean_names <- gsub("\\.tif$", "", basename(prediction_files))
  names(vote_stack) <- clean_names

  max_votes <- max(vapply(seq_len(terra::nlyr(vote_stack)), function(i) {
    v <- terra::global(vote_stack[[i]], "max", na.rm = TRUE)[1, 1]
    if (is.na(v)) 0 else as.integer(v)
  }, integer(1)), na.rm = TRUE)

  if (verbose) message(paste0("Maximum fold vote count detected: ", max_votes,
                              " | Consensus threshold: ", consensus))

  if (consensus > max_votes) {
    stop(paste0(
      "ERROR: consensus = ", consensus, " exceeds the maximum vote count (",
      max_votes, ") in the prediction rasters. ",
      "Reduce consensus to <= ", max_votes, "."
    ))
  }

  if (!dir.exists(consensus_dir)) dir.create(consensus_dir, recursive = TRUE,
                                             showWarnings = FALSE)

  if (verbose) message(paste0("Applying consensus threshold (>= ", consensus,
                              " folds) to ", terra::nlyr(vote_stack), " time steps..."))

  consensus_files <- character(terra::nlyr(vote_stack))

  for (i in seq_len(terra::nlyr(vote_stack))) {
    layer <- vote_stack[[i]]
    con_layer <- terra::app(layer, function(x) {
      ifelse(is.na(x), NA_integer_, as.integer(x >= consensus))
    })

    con_file  <- file.path(consensus_dir,
                           paste0(clean_names[i], "_consensus", consensus, ".tif"))
    terra::writeRaster(con_layer, con_file, overwrite = TRUE,
                       datatype = "INT1U", NAflag = 255)
    consensus_files[i] <- con_file
    if (verbose) message(paste0("  [", i, "/", terra::nlyr(vote_stack), "] ",
                                clean_names[i], " done."))
  }

  con_stack <- terra::rast(consensus_files)
  names(con_stack) <- clean_names

  if (verbose) message("Calculating frequency raster (proportion of time steps suitable)...")
  frequency_rast <- terra::app(con_stack, mean, na.rm = TRUE)

  terra::writeRaster(frequency_rast, frequency_file, overwrite = TRUE)
  if (verbose) message(paste("Saved frequency raster:", basename(frequency_file)))

  n_suitable_always <- sum(terra::values(frequency_rast) == 1, na.rm = TRUE)
  n_suitable_ever   <- sum(terra::values(frequency_rast) >  0, na.rm = TRUE)
  n_total           <- sum(!is.na(terra::values(frequency_rast)))

  if (verbose) message("--- Summarize Raster Outputs Complete ---")
  if (verbose) message(paste0("  Consensus threshold : ", consensus, " of ", max_votes, " folds"))
  if (verbose) message(paste0("  Time steps processed: ", terra::nlyr(vote_stack)))
  if (verbose) message(paste0("  Ever suitable       : ", n_suitable_ever, " / ", n_total,
                              " pixels (", round(n_suitable_ever / n_total * 100, 1), "%)"))
  if (verbose) message(paste0("  Always suitable     : ", n_suitable_always, " / ", n_total,
                              " pixels (", round(n_suitable_always / n_total * 100, 1), "%)"))
  if (verbose) message(paste0("  Output dir          : ", output_dir))

  invisible(list(
    consensus_stack    = con_stack,
    frequency_raster   = frequency_rast,
    consensus          = consensus,
    n_timesteps        = terra::nlyr(vote_stack),
    consensus_dir      = consensus_dir,
    frequency_file     = frequency_file
  ))
}
