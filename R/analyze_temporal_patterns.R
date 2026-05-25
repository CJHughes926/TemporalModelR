#' Analyze Temporal Patterns in Binary Raster Time Series
#'
#' Postprocessing function that applies changepoint detection methods to
#' identify temporal trends in habitat suitability across consecutive
#' predictions. Classifies pixels as stable, increasing in quality, or
#' decreasing in quality, and identifies time periods of significant change.
#'
#' @usage
#' analyze_temporal_patterns(binary_stack, summary_raster, time_steps,
#'                           fastcpd_params = list(), output_dir = NULL,
#'                           n_tiles_x = 1, n_tiles_y = 1, alpha = 0.05,
#'                           spatial_autocorrelation = TRUE, verbose = TRUE,
#'                           estimate_time = TRUE, overwrite = FALSE)
#'
#' @param binary_stack RasterStack, RasterBrick, or character. Stack of binary
#'   raster layers across time, or path to directory containing binary rasters.
#'   Typically from \code{\link{summarize_raster_outputs}}.
#' @param summary_raster RasterLayer. Per-pixel proportion of time periods where
#'  a pixel is suitable. From \code{\link{summarize_raster_outputs}}.
#' @param time_steps Integer vector. Time labels corresponding to raster layers
#'   (same length as number of layers).
#' @param fastcpd_params List. Named list of parameters passed to fastcpd
#'   changepoint detection function. Default is empty list. Supports parameterization
#'   from \code{\link[fastcpd]{fastcpd_binomial}}.
#' @param output_dir Character. Optional. Output directory for pattern rasters.
#'   When \code{NULL} (default), rasters are written to temporary files and not
#'   saved persistently. Provide a path to write named output files to disk.
#' @param n_tiles_x Integer. Number of tiles in the x direction for tiled
#'   processing. Default is \code{1}. Increase to reduce peak memory use for
#'   large rasters and prevent crashes.
#' @param n_tiles_y Integer. Number of tiles in the y direction for tiled
#'   processing. Default is \code{1}. Increase to reduce peak memory use for
#'   large rasters and prevent crashes.
#' @param alpha Numeric. Significance level for changepoint detection. Default
#'   is \code{0.05}.
#' @param spatial_autocorrelation Logical. If \code{TRUE} (default), includes a
#'   neighbor variable in the changepoint analysis to account for spatial
#'   autocorrelation.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing.
#' @param estimate_time Logical. If \code{TRUE} (default), estimates runtime from
#'   a sample of pixels before full processing begins. If \code{FALSE}, proceeds
#'   directly to processing without a time estimate.
#' @param overwrite Logical. If \code{TRUE}, overwrites existing output files.
#'   If \code{FALSE} (default), existing files are skipped.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{pattern}: \code{SpatRaster} classifying pixels as integer
#'     values 1-6, corresponding to "Never Suitable", "Always Suitable",
#'     "No Pattern", "Increasing Suitability", "Decreasing Suitability",
#'     or "Fluctuating".
#'   \item \code{time_decrease}: \code{SpatRaster} showing the time step of
#'     first significant decrease for pixels classified as decreasing.
#'   \item \code{time_increase}: \code{SpatRaster} showing the time step of
#'     first significant increase for pixels classified as increasing.
#' }
#'
#' @details
#' Applies changepoint detection using fastcpd to identify significant temporal
#' shifts in spatial suitability. Accounts for spatial and temporal
#' autocorrelation when \code{spatial_autocorrelation = TRUE}. The
#' \code{fastcpd_params} list allows customization of the changepoint detection
#' algorithm.
#'
#' Pattern classifications enable identification of expanding, contracting, or
#' stable g-space distributions over time or site level assessments of directional
#' change in suitability.
#'
#' Classification assumes consecutive rasters. Time periods
#' shorter than ~15 time steps may be too short to classify increases or decreases.
#'
#' @seealso
#' Postprocessing: \code{\link{summarize_raster_outputs}}
#'
#' External: \code{\link[fastcpd]{fastcpd}}
#'
#' @examples
#' \donttest{
#'   con_file <- system.file("extdata/binary/consensus_stack.tif",
#'         package = "TemporalModelR")
#'
#'   frq_file <- system.file("extdata/binary/frequency_raster.tif",
#'         package = "TemporalModelR")
#'
#'   binary_stack   <- terra::rast(con_file)
#'
#'   summary_raster <- terra::rast(frq_file)
#'
#'   time_steps <- expand.grid(
#'     year    = 1:15,
#'     season  = "Spring",
#'     stringsAsFactors = FALSE
#'   )
#'
#'   analyze_temporal_patterns(
#'     binary_stack   = binary_stack,
#'     summary_raster = summary_raster,
#'     time_steps     = time_steps,
#'     output_dir     = tempdir(),
#'     spatial_autocorrelation = FALSE,
#'     overwrite      = TRUE,
#'     estimate_time  = FALSE,
#'     verbose        = FALSE
#'   )
#' }

#' @export
#' @importFrom terra rast ext res crop nlyr app focal values ncell freq plot writeRaster mosaic
#' @importFrom graphics par legend
#' @importFrom grDevices heat.colors terrain.colors
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @importFrom tools file_path_sans_ext
analyze_temporal_patterns <- function(binary_stack,
           summary_raster,
           time_steps,
           fastcpd_params = list(),
           output_dir     = NULL,
           n_tiles_x      = 1,
           n_tiles_y      = 1,
           alpha          = 0.05,
           spatial_autocorrelation = TRUE,
           verbose        = TRUE,
           estimate_time  = TRUE,
           overwrite      = FALSE) {


  if (!requireNamespace("fastcpd", quietly = TRUE)) {
    stop(paste0(
      "ERROR: The 'fastcpd' package is required for analyze_temporal_patterns(). ",
      "Install with: install.packages('fastcpd')"
    ))
  }

  if (missing(binary_stack)) {
    stop("ERROR: 'binary_stack' is required")
  }
  if (missing(summary_raster)) {
    stop("ERROR: 'summary_raster' is required")
  }
  if (missing(time_steps)) {
    stop("ERROR: 'time_steps' is required")
  }

  ts_norm  <- .normalize_time_steps(time_steps, "analyze_temporal_patterns")
  time_steps        <- ts_norm$time_steps
  secondary_filters <- ts_norm$secondary_filters

  if (is.character(binary_stack) && dir.exists(binary_stack)) {
    if (verbose) message(paste0("Loading binary rasters from directory: ", binary_stack))
    binary_files <- list.files(binary_stack, pattern = "\\.tif$", full.names = TRUE)
    if (length(binary_files) == 0) {
      stop(paste0("ERROR: No .tif files found in provided binary_stack directory: ", binary_stack))
    }
    binary_files <- .natural_sort(binary_files)
    binary_stack <- terra::rast(binary_files)
    names(binary_stack) <- tools::file_path_sans_ext(basename(binary_files))
    if (verbose) message(paste("Loaded", terra::nlyr(binary_stack), "binary raster layers"))
  } else if (is.character(binary_stack) && file.exists(binary_stack)) {
    if (verbose) message(paste0("Loading binary raster from file: ", binary_stack))
    binary_stack <- terra::rast(binary_stack)
  } else if (inherits(binary_stack, "SpatRaster")) {
    if (verbose) message("Using provided raster object")
  } else {
    stop("ERROR: 'binary_stack' must be a directory path, file path, or SpatRaster object.")
  }

  summary_raster <- .load_raster_input(summary_raster, "summary_raster")

  if (length(secondary_filters) > 0) {
    lyr_names <- names(binary_stack)
    keep_idx  <- seq_len(terra::nlyr(binary_stack))
    for (sc in names(secondary_filters)) {
      val      <- as.character(secondary_filters[[sc]])
      keep_idx <- keep_idx[grepl(val, lyr_names[keep_idx], ignore.case = TRUE)]
    }
    if (length(keep_idx) == 0) {
      stop(paste0(
        "ERROR: No binary raster layers matched the secondary filter(s): ",
        paste(mapply(function(k, v) paste0(k, "=", v),
   names(secondary_filters), secondary_filters), collapse = ", "),
        ". Layer names in binary_stack: ",
        paste(head(lyr_names, 5), collapse = ", "),
        if (length(lyr_names) > 5) paste0(" ... (", length(lyr_names), " total)") else ""
      ))
    }
    if (length(keep_idx) != length(time_steps)) {
      stop(paste0(
        "ERROR: After filtering by secondary column(s), found ", length(keep_idx),
        " matching raster layer(s) but time_steps has ", length(time_steps), " value(s). ",
        "Ensure time_steps values match the primary time steps in your predictions."
      ))
    }
    binary_stack <- binary_stack[[keep_idx]]
    if (verbose) message(paste0(
      "Filtered binary_stack to ", terra::nlyr(binary_stack), " layer(s) matching: ",
      paste(mapply(function(k, v) paste0(k, "=", v),
          names(secondary_filters), secondary_filters), collapse = ", ")
    ))
  } else {
    lyr_names <- names(binary_stack)
    if (terra::nlyr(binary_stack) > length(time_steps)) {
      match_counts <- vapply(as.character(time_steps), function(ts) {
        sum(grepl(paste0("(^|[^0-9])", ts, "([^0-9]|$)"), lyr_names))
      }, integer(1))
      ambiguous <- time_steps[match_counts > 1]
      if (length(ambiguous) > 0) {
        stop(paste0(
 "ERROR: The binary stack has ", terra::nlyr(binary_stack),
 " layers but time_steps only has ", length(time_steps), " value(s). ",
 "Time step value(s) ", paste(ambiguous, collapse = ", "),
 " each match more than one raster layer, suggesting the stack contains ",
 "predictions for multiple secondary time values (e.g. seasons). ",
 "Supply time_steps as a multi-column data frame specifying the secondary ",
 "column value to use. For example: ",
 "time_steps <- expand.grid(time_step = 1:15, season = \"Spring\", stringsAsFactors = FALSE)"
        ))
      }
    }
  }

  if (!is.list(fastcpd_params)) {
    stop("ERROR: 'fastcpd_params' must be a named list")
  }

  save_output <- !is.null(output_dir)
  if (save_output && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  ts_range      <- paste0(min(time_steps), "_", max(time_steps))
  pattern_file  <- if (save_output) file.path(output_dir,
          paste0("pattern_raster_",     ts_range, ".tif")) else tempfile(fileext = ".tif")
  decrease_file <- if (save_output) file.path(output_dir,
          paste0("time_first_decrease_", ts_range, ".tif")) else tempfile(fileext = ".tif")
  increase_file <- if (save_output) file.path(output_dir,
          paste0("time_first_increase_", ts_range, ".tif")) else tempfile(fileext = ".tif")

  if (file.exists(pattern_file) && file.exists(decrease_file) && file.exists(increase_file) && !overwrite) {
    if (verbose) message("Output rasters exist. Set overwrite = TRUE to rerun.")
    return(invisible(list(
      pattern = terra::rast(pattern_file),
      time_decrease = terra::rast(decrease_file),
      time_increase = terra::rast(increase_file)
    )))
  }

  if (spatial_autocorrelation) {
    if (verbose) message("Spatial autocorrelation: ENABLED (neighbor variable included)")
  } else {
    if (verbose) message("Spatial autocorrelation: DISABLED (neighbor variable excluded)")
  }

  skip_tiling <- (is.null(n_tiles_x) || n_tiles_x == 1) && (is.null(n_tiles_y) || n_tiles_y == 1)

  if (skip_tiling) {
    if (verbose) message("Processing entire raster without tiling...")

    if (estimate_time) {
      if (verbose) message("Estimating processing time...")

      n_times <- terra::nlyr(binary_stack)
      n_middle <- n_times - 2

      summary_vals <- terra::values(summary_raster, mat = FALSE)
      valid_indices <- which(!is.na(summary_vals))

      if (length(valid_indices) > 0) {
        mean_vals <- summary_vals[valid_indices]
        n_quick <- sum(mean_vals < 0.01 | mean_vals > 0.99)
        n_complex <- sum(mean_vals >= 0.01 & mean_vals <= 0.99)

        if (verbose) message(paste0("Quick pixels (always absent/present): ", format(n_quick, big.mark = ",")))
        if (verbose) message(paste0("Complex pixels (changepoint analysis): ", format(n_complex, big.mark = ",")))

        if (n_complex > 0) {
 if (verbose) message("Timing sample pixels...")

 middle_times <- binary_stack[[2:(n_times - 1)]]
 lag_stack <- binary_stack[[1:(n_times - 2)]]

 if (spatial_autocorrelation) {
   middle_neighbor <- terra::rast(lapply(seq_len(terra::nlyr(middle_times)), function(i) {
     terra::focal(middle_times[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
   }))
   predictor_stack <- c(middle_times, lag_stack, middle_neighbor, summary_raster)
 } else {
   predictor_stack <- c(middle_times, lag_stack, summary_raster)
 }

 pred_vals_all <- terra::values(predictor_stack, mat = TRUE)
 valid_pred <- which(!is.na(pred_vals_all[, 1]))

 if (spatial_autocorrelation) {
   mean_pred <- pred_vals_all[valid_pred, (3 * n_middle + 1)]
 } else {
   mean_pred <- pred_vals_all[valid_pred, (2 * n_middle + 1)]
 }

 complex_indices <- valid_pred[mean_pred >= 0.01 & mean_pred <= 0.99]

 sample_size <- min(100, length(complex_indices))
 sample_indices <- sample(complex_indices, size = sample_size, replace = FALSE)

 start_time <- Sys.time()
 for (idx in sample_indices) {
   result <- .classify_pixel_with_times(pred_vals_all[idx, ], n_middle,
             time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
 }
 end_time <- Sys.time()

 time_per_pixel <- as.numeric(difftime(end_time, start_time, units = "secs")) / sample_size
 if (verbose) message(paste("Average time per complex pixel:", round(time_per_pixel, 4), "seconds"))

 base_seconds <- time_per_pixel * n_complex
 overhead <- 30

 lower_seconds <- (base_seconds * 0.8) + (overhead * 0.5)
 upper_seconds <- (base_seconds * 1.2) + (overhead * 1.5)

 if (verbose) message(paste("Estimated processing time:",
          .format_time_estimate(lower_seconds), "to",
          .format_time_estimate(upper_seconds)))

 rm(middle_times, lag_stack, predictor_stack, pred_vals_all)
 if (spatial_autocorrelation) rm(middle_neighbor)
 gc(verbose = FALSE)
        }
      }
    }

    if (verbose) message("Processing raster...")

    n_times <- terra::nlyr(binary_stack)
    n_middle <- n_times - 2

    middle_times <- binary_stack[[2:(n_times - 1)]]
    lag_stack <- binary_stack[[1:(n_times - 2)]]

    if (spatial_autocorrelation) {
      middle_neighbor <- terra::rast(lapply(seq_len(terra::nlyr(middle_times)), function(i) {
        terra::focal(middle_times[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
      }))
      predictor_stack <- c(middle_times, lag_stack, middle_neighbor, summary_raster)
    } else {
      predictor_stack <- c(middle_times, lag_stack, summary_raster)
    }

    if (verbose) {
      n_cells <- terra::ncell(predictor_stack)
      pb <- utils::txtProgressBar(min = 0, max = n_cells, style = 3, width = 50)

      pred_vals <- terra::values(predictor_stack, mat = TRUE)
      pattern_vals <- numeric(n_cells)
      decrease_vals <- numeric(n_cells)
      increase_vals <- numeric(n_cells)

      for (cell_i in seq_len(n_cells)) {
        if (!any(is.na(pred_vals[cell_i, ]))) {
 result <- .classify_pixel_with_times(pred_vals[cell_i, ], n_middle,
           time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
 pattern_vals[cell_i] <- result[1]
 decrease_vals[cell_i] <- result[2]
 increase_vals[cell_i] <- result[3]
        } else {
 pattern_vals[cell_i] <- NA
 decrease_vals[cell_i] <- NA
 increase_vals[cell_i] <- NA
        }
        if (cell_i %% 10 == 0) utils::setTxtProgressBar(pb, cell_i)
      }
      close(pb)
      if (verbose) message("")

      pattern_raster <- terra::rast(predictor_stack, nlyr = 1)
      decrease_raster <- terra::rast(predictor_stack, nlyr = 1)
      increase_raster <- terra::rast(predictor_stack, nlyr = 1)

      terra::values(pattern_raster) <- pattern_vals
      terra::values(decrease_raster) <- decrease_vals
      terra::values(increase_raster) <- increase_vals

    } else {
      result_matrix <- terra::app(predictor_stack,
       fun = function(x) .classify_pixel_with_times(x, n_middle, time_steps,
                fastcpd_params, alpha, use_neighbor = spatial_autocorrelation))

      pattern_raster <- result_matrix[[1]]
      decrease_raster <- result_matrix[[2]]
      increase_raster <- result_matrix[[3]]
    }

  } else {

    tiles_dir <- if (save_output) file.path(output_dir, "tiles") else file.path(tempdir(), "tiles")
    if (!dir.exists(tiles_dir)) dir.create(tiles_dir, recursive = TRUE, showWarnings = FALSE)

    if (verbose) message("Calculating tile extents...")

    full_ext <- terra::ext(binary_stack)
    x_min <- full_ext[1]
    x_max <- full_ext[2]
    y_min <- full_ext[3]
    y_max <- full_ext[4]

    res_vals <- terra::res(binary_stack)
    res_x <- res_vals[1]
    res_y <- res_vals[2]

    x_range <- x_max - x_min
    y_range <- y_max - y_min

    tile_width <- x_range / n_tiles_x
    tile_height <- y_range / n_tiles_y

    tile_extents <- list()
    tile_idx <- 1

    for (i in seq_len(n_tiles_y)) {
      for (j in seq_len(n_tiles_x)) {
        tile_x_min <- x_min + (j - 1) * tile_width
        tile_x_max <- x_min + j * tile_width
        tile_y_min <- y_min + (i - 1) * tile_height
        tile_y_max <- y_min + i * tile_height

        tile_extents[[tile_idx]] <- terra::ext(tile_x_min, tile_x_max, tile_y_min, tile_y_max)
        tile_idx <- tile_idx + 1
      }
    }

    n_tiles <- length(tile_extents)
    if (verbose) message(paste("Created", n_tiles, "tiles"))

    if (estimate_time) {
      if (verbose) message("Estimating processing time...")

      n_times <- terra::nlyr(binary_stack)
      n_middle <- n_times - 2

      total_complex <- 0
      total_quick <- 0
      time_per_pixel <- NULL

      for (tile_i in seq_len(n_tiles)) {
        if (verbose) message(paste0("Scanning tile ", tile_i, "/", n_tiles, "..."))

        tile_ext <- tile_extents[[tile_i]]

        tryCatch({
 suppressWarnings({
   tile_binary <- terra::crop(binary_stack, tile_ext)
   tile_summary <- terra::crop(summary_raster, tile_ext)
 })

 summary_vals <- terra::values(tile_summary, mat = FALSE)
 valid_indices <- which(!is.na(summary_vals))

 if (length(valid_indices) > 0) {
   mean_vals <- summary_vals[valid_indices]
   n_quick <- sum(mean_vals < 0.01 | mean_vals > 0.99)
   n_complex <- sum(mean_vals >= 0.01 & mean_vals <= 0.99)

   total_quick <- total_quick + n_quick
   total_complex <- total_complex + n_complex

   if (is.null(time_per_pixel) && n_complex > 0) {
     if (verbose) message("Timing sample pixels...")

     middle_times <- tile_binary[[2:(n_times - 1)]]
     lag_stack <- tile_binary[[1:(n_times - 2)]]

     if (spatial_autocorrelation) {
       middle_neighbor <- terra::rast(lapply(seq_len(terra::nlyr(middle_times)), function(i) {
         terra::focal(middle_times[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
       }))
       predictor_stack <- c(middle_times, lag_stack, middle_neighbor, tile_summary)
     } else {
       predictor_stack <- c(middle_times, lag_stack, tile_summary)
     }

     pred_vals_all <- terra::values(predictor_stack, mat = TRUE)
     valid_pred <- which(!is.na(pred_vals_all[, 1]))

     if (spatial_autocorrelation) {
       mean_pred <- pred_vals_all[valid_pred, (3 * n_middle + 1)]
     } else {
       mean_pred <- pred_vals_all[valid_pred, (2 * n_middle + 1)]
     }

     complex_indices <- valid_pred[mean_pred >= 0.01 & mean_pred <= 0.99]

     sample_size <- min(100, length(complex_indices))
     sample_indices <- sample(complex_indices, size = sample_size, replace = FALSE)

     start_time <- Sys.time()
     for (idx in sample_indices) {
       result <- .classify_pixel_with_times(pred_vals_all[idx, ], n_middle,
        time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
     }
     end_time <- Sys.time()

     time_per_pixel <- as.numeric(difftime(end_time, start_time, units = "secs")) / sample_size
     if (verbose) message(paste("Average time per complex pixel:", round(time_per_pixel, 4), "seconds"))
   }
 }
        }, error = function(e) {
 if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
   stop("ERROR: Memory error. Increase n_tiles_x and n_tiles_y to use smaller tiles.")
 } else {
   stop(e)
 }
        })
      }

      if (verbose) message(paste0("Quick pixels (always absent/present): ", format(total_quick, big.mark = ",")))
      if (verbose) message(paste0("Complex pixels (changepoint analysis): ", format(total_complex, big.mark = ",")))

      if (!is.null(time_per_pixel) && total_complex > 0) {
        base_seconds <- time_per_pixel * total_complex
        overhead_per_tile <- 15
        total_overhead <- (overhead_per_tile * n_tiles) + 30

        lower_seconds <- (base_seconds * 0.8) + (total_overhead * 0.5)
        upper_seconds <- (base_seconds * 1.2) + (total_overhead * 1.5)

        if (verbose) message(paste("Estimated processing time:",
        .format_time_estimate(lower_seconds), "to",
        .format_time_estimate(upper_seconds)))
      }
    }

    if (verbose) message("Processing tiles...")

    tile_files_pattern <- character(n_tiles)
    tile_files_decrease <- character(n_tiles)
    tile_files_increase <- character(n_tiles)

    for (tile_i in seq_len(n_tiles)) {
      if (verbose) message(paste("Processing tile", tile_i, "of", n_tiles))

      tile_file_pattern <- file.path(tiles_dir, paste0("pattern_tile_", tile_i, ".tif"))
      tile_file_decrease <- file.path(tiles_dir, paste0("decrease_tile_", tile_i, ".tif"))
      tile_file_increase <- file.path(tiles_dir, paste0("increase_tile_", tile_i, ".tif"))

      tile_files_pattern[tile_i] <- tile_file_pattern
      tile_files_decrease[tile_i] <- tile_file_decrease
      tile_files_increase[tile_i] <- tile_file_increase

      tile_ext <- tile_extents[[tile_i]]

      tryCatch({
        suppressWarnings({
 tile_binary <- terra::crop(binary_stack, tile_ext)
 tile_summary <- terra::crop(summary_raster, tile_ext)
        })

        n_times <- terra::nlyr(tile_binary)
        n_middle <- n_times - 2

        middle_times <- tile_binary[[2:(n_times - 1)]]
        lag_stack <- tile_binary[[1:(n_times - 2)]]

        if (spatial_autocorrelation) {
 middle_neighbor <- terra::rast(lapply(seq_len(terra::nlyr(middle_times)), function(i) {
   terra::focal(middle_times[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
 }))
 predictor_stack <- c(middle_times, lag_stack, middle_neighbor, tile_summary)
        } else {
 predictor_stack <- c(middle_times, lag_stack, tile_summary)
        }

        if (verbose) {
 n_cells <- terra::ncell(predictor_stack)
 pb <- utils::txtProgressBar(min = 0, max = n_cells, style = 3, width = 50)

 pred_vals <- terra::values(predictor_stack, mat = TRUE)
 pattern_vals <- numeric(n_cells)
 decrease_vals <- numeric(n_cells)
 increase_vals <- numeric(n_cells)

 for (cell_i in seq_len(n_cells)) {
   if (!any(is.na(pred_vals[cell_i, ]))) {
     result <- .classify_pixel_with_times(pred_vals[cell_i, ], n_middle,
      time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
     pattern_vals[cell_i] <- result[1]
     decrease_vals[cell_i] <- result[2]
     increase_vals[cell_i] <- result[3]
   } else {
     pattern_vals[cell_i] <- NA
     decrease_vals[cell_i] <- NA
     increase_vals[cell_i] <- NA
   }
   if (cell_i %% 10 == 0) utils::setTxtProgressBar(pb, cell_i)
 }
 close(pb)
 if (verbose) message("")

 tile_pattern <- terra::rast(predictor_stack, nlyr = 1)
 tile_decrease <- terra::rast(predictor_stack, nlyr = 1)
 tile_increase <- terra::rast(predictor_stack, nlyr = 1)

 terra::values(tile_pattern) <- pattern_vals
 terra::values(tile_decrease) <- decrease_vals
 terra::values(tile_increase) <- increase_vals

        } else {
 result_matrix <- terra::app(predictor_stack,
           fun = function(x) .classify_pixel_with_times(x, n_middle, time_steps,
           fastcpd_params, alpha, use_neighbor = spatial_autocorrelation))

 tile_pattern <- result_matrix[[1]]
 tile_decrease <- result_matrix[[2]]
 tile_increase <- result_matrix[[3]]
        }

        terra::writeRaster(tile_pattern, tile_file_pattern, overwrite = TRUE,
         datatype = "INT1U", gdal = c("COMPRESS=LZW"))
        terra::writeRaster(tile_decrease, tile_file_decrease, overwrite = TRUE,
         datatype = "INT2S", gdal = c("COMPRESS=LZW"))
        terra::writeRaster(tile_increase, tile_file_increase, overwrite = TRUE,
         datatype = "INT2S", gdal = c("COMPRESS=LZW"))

        rm(tile_binary, tile_summary, middle_times, lag_stack,
  predictor_stack, tile_pattern, tile_decrease, tile_increase)
        if (spatial_autocorrelation) rm(middle_neighbor)
        if (exists("pred_vals")) rm(pred_vals, pattern_vals, decrease_vals, increase_vals)
        gc(verbose = FALSE)

      }, error = function(e) {
        if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
 stop(paste0("ERROR: Memory error on tile ", tile_i, ". Increase n_tiles_x and n_tiles_y."))
        } else {
 stop(e)
        }
      })
    }

    if (n_tiles > 1) {
      if (verbose) message("Merging tiles...")

      tryCatch({
        tile_rasters_pattern <- lapply(tile_files_pattern, rast)
        pattern_raster <- do.call(mosaic, c(tile_rasters_pattern, fun = "mean"))

        tile_rasters_decrease <- lapply(tile_files_decrease, rast)
        decrease_raster <- do.call(mosaic, c(tile_rasters_decrease, fun = "mean"))

        tile_rasters_increase <- lapply(tile_files_increase, rast)
        increase_raster <- do.call(mosaic, c(tile_rasters_increase, fun = "mean"))

      }, error = function(e) {
        if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
 stop("ERROR: Memory error merging tiles. Increase n_tiles_x and n_tiles_y.")
        } else {
 stop(e)
        }
      })
    } else {
      if (verbose) message("Single tile detected, skipping merge...")

      pattern_raster <- terra::rast(tile_files_pattern[1])
      decrease_raster <- terra::rast(tile_files_decrease[1])
      increase_raster <- terra::rast(tile_files_increase[1])
    }

    if (verbose) message("Cleaning up tile files...")
    unlink(tiles_dir, recursive = TRUE)
    if (verbose) message("Tiles removed")
  }

  if (verbose) message("Saving results...")
  terra::writeRaster(pattern_raster, pattern_file, overwrite = TRUE,
   datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  terra::writeRaster(decrease_raster, decrease_file, overwrite = TRUE,
   datatype = "INT2S", gdal = c("COMPRESS=LZW"))
  terra::writeRaster(increase_raster, increase_file, overwrite = TRUE,
   datatype = "INT2S", gdal = c("COMPRESS=LZW"))

  if (verbose) message("Generating plots...")

  pat_cols   <- c("#730000", "#267300", "#B2B2B2", "#A3FF73", "#FF7F7F", "#A900E6", "#eed202")
  pat_labels <- c("Always Absent", "Always Present", "No Pattern",
         "Increasing", "Decreasing", "Fluctuating", "Failed")

  opar1 <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(6.5, 4, 3, 2))
  terra::plot(pattern_raster,
     col    = pat_cols,
     main   = paste("Pattern Classification\n", min(time_steps), "-", max(time_steps)),
     breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5),
     legend = FALSE)
  graphics::legend(
    x      = "bottom",
    legend = pat_labels,
    fill   = pat_cols,
    cex    = 0.7,
    bty    = "n",
    horiz  = FALSE,
    ncol   = 4,
    xpd    = TRUE,
    inset  = c(0, -0.42)
  )
  graphics::par(opar1)

  dec_vals <- terra::values(decrease_raster, na.rm = TRUE)
  if (length(dec_vals) > 0 && !all(is.na(dec_vals))) {
    opar2 <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(4, 4, 3, 2))
    terra::plot(decrease_raster,
       main = paste("Time of First Decrease\n", min(time_steps), "-", max(time_steps)),
       col  = rev(heat.colors(50)))
    graphics::par(opar2)
  } else {
    if (verbose) message("Skipping decrease raster plot: no pixels with a significant decrease detected.")
  }

  inc_vals <- terra::values(increase_raster, na.rm = TRUE)
  if (length(inc_vals) > 0 && !all(is.na(inc_vals))) {
    opar3 <- graphics::par(no.readonly = TRUE)
    graphics::par(mar = c(4, 4, 3, 2))
    terra::plot(increase_raster,
       main = paste("Time of First Increase\n", min(time_steps), "-", max(time_steps)),
       col  = terrain.colors(50))
    graphics::par(opar3)
  } else {
    if (verbose) message("Skipping increase raster plot: no pixels with a significant increase detected.")
  }

  if (verbose) message("Analysis complete")
  if (verbose) message(paste("Period:", min(time_steps), "-", max(time_steps)))
  if (!skip_tiling) {
    if (verbose) message(paste("Tiles:", n_tiles))
  }
  if (verbose) message(paste("Spatial autocorrelation:", ifelse(spatial_autocorrelation, "ENABLED", "DISABLED")))

  pattern_freq <- terra::freq(pattern_raster)
  if (!is.null(pattern_freq) && nrow(pattern_freq) > 0) {
    pattern_freq <- as.data.frame(pattern_freq)
    pattern_freq$proportion <- round(pattern_freq$count / sum(pattern_freq$count), 3)
    pattern_names <- c("Always Absent", "Always Present", "No Pattern",
     "Increasing", "Decreasing", "Fluctuating", "Failed")
    pattern_freq$pattern <- pattern_names[pattern_freq$value]

    if (verbose) message("Pattern Classifications:")
    if (verbose) message(pattern_freq[, c("value", "pattern", "count", "proportion")])
  }

  dec_freq <- as.data.frame(terra::freq(decrease_raster))
  dec_freq <- dec_freq[!is.na(dec_freq$value), ]
  if (nrow(dec_freq) > 0) {
    n_decreasing <- sum(dec_freq$count)
    if (verbose) message(paste("Decreasing pixels:", format(n_decreasing, big.mark = ",")))
    if (verbose) message(paste("Time range:", min(dec_freq$value), "-", max(dec_freq$value)))
    if (verbose) message(paste("Most common:", dec_freq$value[which.max(dec_freq$count)],
    paste0("(", format(max(dec_freq$count), big.mark = ","), " pixels)")))
  }

  inc_freq <- as.data.frame(terra::freq(increase_raster))
  inc_freq <- inc_freq[!is.na(inc_freq$value), ]
  if (nrow(inc_freq) > 0) {
    n_increasing <- sum(inc_freq$count)
    if (verbose) message(paste("Increasing pixels:", format(n_increasing, big.mark = ",")))
    if (verbose) message(paste("Time range:", min(inc_freq$value), "-", max(inc_freq$value)))
    if (verbose) message(paste("Most common:", inc_freq$value[which.max(inc_freq$count)],
    paste0("(", format(max(inc_freq$count), big.mark = ","), " pixels)")))
  }

  invisible(list(
    pattern = pattern_raster,
    time_decrease = decrease_raster,
    time_increase = increase_raster
  ))
}
