#' Analyze Temporal Patterns in Binary Raster Time Series
#'
#' Postprocessing function that applies changepoint detection methods to
#' identify temporal trends in habitat suitability across consecutive
#' predictions. Classifies pixels as stable, increasing in quality, or
#' decreasing in quality, and identifies time periods of significant change.
#'
#' @param binary_stack RasterStack, RasterBrick, or character. Stack of binary
#'   raster layers across time, or path to directory containing binary rasters.
#'   Typically from \code{\link{summarize_raster_outputs}}.
#' @param summary_raster RasterLayer. Per-pixel summary statistic (e.g., mean
#'   suitability across time). From \code{\link{summarize_raster_outputs}}.
#' @param time_steps Integer vector. Time labels corresponding to raster layers
#'   (same length as number of layers).
#' @param fastcpd_params List. Named list of parameters passed to fastcpd
#'   changepoint detection function. Default is empty list.
#' @param output_dir Character. Output directory for pattern rasters. Default is
#'   "output".
#' @param n_tiles_x Integer. Number of tiles in x direction for parallel
#'   processing. Default is 1.
#' @param n_tiles_y Integer. Number of tiles in y direction for parallel
#'   processing. Default is 1.
#' @param alpha Numeric. Significance level for changepoint detection. Default
#'   is 0.05.
#' @param spatial_autocorrelation Logical. If TRUE, includes neighbor variable in
#'   analysis to account for spatial autocorrelation. If FALSE, spatial
#'   autocorrelation is not considered. Default is TRUE.
#' @param show_progress Logical. If TRUE, displays progress bar during
#'   processing. If FALSE, runs silently. Default is TRUE.
#' @param estimate_time Logical. If TRUE, estimates runtime from sampled pixels
#'   before processing. If FALSE, proceeds directly to processing. Default is
#'   TRUE.
#' @param overwrite Logical. If TRUE, overwrites existing output files. If
#'   FALSE, skips files that already exist. Default is FALSE.
#'
#' @return A list containing:
#' \itemize{
#'   \item pattern: RasterLayer classifying pixels as "Never Suitable", "Always
#'     Suitable", "No Pattern", "Increasing Suitability", "Decreasing
#'     Suitability", or "Fluctuating"
#'   \item time_decrease: RasterLayer showing time period of first significant
#'     decrease (for decreasing pixels)
#'   \item time_increase: RasterLayer showing time period of first significant
#'     increase (for increasing pixels)
#' }
#'
#' @details
#' Applies changepoint detection using fastcpd to identify significant temporal
#' shifts in habitat suitability. Accounts for spatial and temporal
#' autocorrelation when spatial_autocorrelation is TRUE. The fastcpd_params list
#' allows customization of the changepoint detection algorithm.
#'
#' Pattern classifications enable identification of expanding, contracting, or
#' stable habitat distributions over time.
#'
#' @seealso
#' Postprocessing: \code{\link{summarize_raster_outputs}}
#'
#' External: \code{\link[fastcpd]{fastcpd}}
#'
#' @examples
#' \dontrun{
#' pattern_results <- analyze_temporal_patterns(
#'   binary_stack = consensus_stack,
#'   summary_raster = summary_raster,
#'   time_steps = 2000:2020,
#'   fastcpd_params = list(),
#'   output_dir = "temporal_patterns/"
#' )
#' }
#'
#' @export
#' @importFrom terra rast ext res crop nlyr subset focal values ncell mosaic
#'   writeRaster freq plot
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom graphics par legend
#' @importFrom grDevices heat.colors terrain.colors
analyze_temporal_patterns <- function(binary_stack,
                                      summary_raster,
                                      time_steps,
                                      fastcpd_params = list(),
                                      output_dir = "output",
                                      n_tiles_x = 1,
                                      n_tiles_y = 1,
                                      alpha = 0.05,
                                      spatial_autocorrelation = TRUE,
                                      show_progress = TRUE,
                                      estimate_time = TRUE,
                                      overwrite = FALSE) {

  require(terra)
  require(fastcpd)

  ### Input validation

  if (missing(binary_stack)) {
    stop("ERROR: binary_stack is required")
  }
  if (missing(summary_raster)) {
    stop("ERROR: summary_raster is required")
  }
  if (missing(time_steps)) {
    stop("ERROR: time_steps is required")
  }

  ### Handle binary_stack as directory or SpatRaster

  if (is.character(binary_stack) && dir.exists(binary_stack)) {
    print(paste0("Loading binary rasters from directory: ", binary_stack))
    binary_files <- list.files(binary_stack, pattern = "\\.tif$", full.names = TRUE)
    if (length(binary_files) == 0) {
      stop(paste0("ERROR: No .tif files found in provided binary_stack directory: ", binary_stack))
    }
    binary_stack <- rast(binary_files)
    print(paste("Loaded", nlyr(binary_stack), "binary raster layers"))
  } else if (is.character(binary_stack) && file.exists(binary_stack)) {
    print(paste0("Loading binary raster from file: ", binary_stack))
    binary_stack <- rast(binary_stack)
  } else if (inherits(binary_stack, "SpatRaster")) {
    print("Using provided raster object")
  } else {
    stop("ERROR: binary_stack must be a directory path, file path, or SpatRaster object.")
  }

  ### Handle summary_raster as file or SpatRaster

  if (is.character(summary_raster)) {
    if (!file.exists(summary_raster)) {
      stop(paste0("ERROR: summary_raster file does not exist: ", summary_raster))
    }
    summary_raster <- rast(summary_raster)
  } else if (!inherits(summary_raster, "SpatRaster")) {
    stop("ERROR: summary_raster must be a file path or SpatRaster object")
  }

  ### Validate fastcpd_params

  if (!is.list(fastcpd_params)) {
    stop("ERROR: fastcpd_params must be a named list")
  }

  ### Setup output directories

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  tiles_dir <- file.path(output_dir, "tiles")
  if (!dir.exists(tiles_dir)) dir.create(tiles_dir, recursive = TRUE)

  ### Define output files

  pattern_file <- file.path(output_dir, paste0("pattern_raster_", min(time_steps), "_", max(time_steps), ".tif"))
  decrease_file <- file.path(output_dir, paste0("year_first_decrease_", min(time_steps), "_", max(time_steps), ".tif"))
  increase_file <- file.path(output_dir, paste0("year_first_increase_", min(time_steps), "_", max(time_steps), ".tif"))

  if (file.exists(pattern_file) && file.exists(decrease_file) && file.exists(increase_file) && !overwrite) {
    print("Output rasters exist. Set overwrite = TRUE to rerun.")
    return(list(
      pattern = rast(pattern_file),
      year_decrease = rast(decrease_file),
      year_increase = rast(increase_file)
    ))
  }

  ### Spatial autocorrelation message

  if (spatial_autocorrelation) {
    print("Spatial autocorrelation: ENABLED (neighbor variable included)")
  } else {
    print("Spatial autocorrelation: DISABLED (neighbor variable excluded)")
  }

  ### Calculate tile extents

  print("Calculating tile extents...")

  full_ext <- ext(binary_stack)
  x_min <- full_ext[1]
  x_max <- full_ext[2]
  y_min <- full_ext[3]
  y_max <- full_ext[4]

  res_vals <- res(binary_stack)
  res_x <- res_vals[1]
  res_y <- res_vals[2]

  x_range <- x_max - x_min
  y_range <- y_max - y_min

  tile_width <- x_range / n_tiles_x
  tile_height <- y_range / n_tiles_y

  tile_extents <- list()
  tile_idx <- 1

  for (i in 1:n_tiles_y) {
    for (j in 1:n_tiles_x) {
      tile_x_min <- x_min + (j - 1) * tile_width
      tile_x_max <- x_min + j * tile_width
      tile_y_min <- y_min + (i - 1) * tile_height
      tile_y_max <- y_min + i * tile_height

      tile_extents[[tile_idx]] <- ext(tile_x_min, tile_x_max, tile_y_min, tile_y_max)
      tile_idx <- tile_idx + 1
    }
  }

  n_tiles <- length(tile_extents)
  print(paste("Created", n_tiles, "tiles"))

  ### Time estimation

  if (estimate_time) {
    print("Estimating processing time...")

    n_years <- nlyr(binary_stack)
    n_middle <- n_years - 2

    total_complex <- 0
    total_quick <- 0
    time_per_pixel <- NULL

    for (tile_i in 1:n_tiles) {
      cat(paste0("Scanning tile ", tile_i, "/", n_tiles, "...\r"))

      tile_ext <- tile_extents[[tile_i]]

      tryCatch({
        suppressWarnings({
          tile_binary <- crop(binary_stack, tile_ext)
          tile_summary <- crop(summary_raster, tile_ext)
        })

        summary_vals <- values(tile_summary, mat = FALSE)
        valid_indices <- which(!is.na(summary_vals))

        if (length(valid_indices) > 0) {
          mean_vals <- summary_vals[valid_indices]
          n_quick <- sum(mean_vals < 0.01 | mean_vals > 0.99)
          n_complex <- sum(mean_vals >= 0.01 & mean_vals <= 0.99)

          total_quick <- total_quick + n_quick
          total_complex <- total_complex + n_complex

          if (is.null(time_per_pixel) && n_complex > 0) {
            print("Timing sample pixels...")

            middle_years <- tile_binary[[2:(n_years - 1)]]
            lag_stack <- tile_binary[[1:(n_years - 2)]]

            if (spatial_autocorrelation) {
              middle_neighbor <- rast(lapply(1:nlyr(middle_years), function(i) {
                focal(middle_years[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
              }))
              predictor_stack <- c(middle_years, lag_stack, middle_neighbor, tile_summary)
            } else {
              predictor_stack <- c(middle_years, lag_stack, tile_summary)
            }

            pred_vals_all <- values(predictor_stack, mat = TRUE)
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
              result <- classify_pixel_with_years(pred_vals_all[idx, ], n_middle,
                                                  time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
            }
            end_time <- Sys.time()

            time_per_pixel <- as.numeric(difftime(end_time, start_time, units = "secs")) / sample_size
            print(paste("Average time per complex pixel:", round(time_per_pixel, 4), "seconds"))
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

    print(paste0("Quick pixels (always absent/present): ", format(total_quick, big.mark = ",")))
    print(paste0("Complex pixels (changepoint analysis): ", format(total_complex, big.mark = ",")))

    if (!is.null(time_per_pixel) && total_complex > 0) {
      base_seconds <- time_per_pixel * total_complex
      overhead_per_tile <- 15
      total_overhead <- (overhead_per_tile * n_tiles) + 30

      lower_seconds <- (base_seconds * 0.8) + (total_overhead * 0.5)
      upper_seconds <- (base_seconds * 1.2) + (total_overhead * 1.5)

      format_time <- function(seconds) {
        hours <- seconds / 3600
        minutes <- seconds / 60
        if (hours > 1) {
          return(paste0(round(hours, 1), " hours"))
        } else if (minutes > 1) {
          return(paste0(round(minutes, 1), " minutes"))
        } else {
          return(paste0(round(seconds, 1), " seconds"))
        }
      }

      print(paste("Estimated processing time:", format_time(lower_seconds), "to", format_time(upper_seconds)))
    }
  }

  ### Process tiles

  print("Processing tiles...")

  tile_files_pattern <- character(n_tiles)
  tile_files_decrease <- character(n_tiles)
  tile_files_increase <- character(n_tiles)

  for (tile_i in 1:n_tiles) {
    print(paste("Processing tile", tile_i, "of", n_tiles))

    tile_file_pattern <- file.path(tiles_dir, paste0("pattern_tile_", tile_i, ".tif"))
    tile_file_decrease <- file.path(tiles_dir, paste0("decrease_tile_", tile_i, ".tif"))
    tile_file_increase <- file.path(tiles_dir, paste0("increase_tile_", tile_i, ".tif"))

    tile_files_pattern[tile_i] <- tile_file_pattern
    tile_files_decrease[tile_i] <- tile_file_decrease
    tile_files_increase[tile_i] <- tile_file_increase

    tile_ext <- tile_extents[[tile_i]]

    tryCatch({
      suppressWarnings({
        tile_binary <- crop(binary_stack, tile_ext)
        tile_summary <- crop(summary_raster, tile_ext)
      })

      n_years <- nlyr(tile_binary)
      n_middle <- n_years - 2

      middle_years <- tile_binary[[2:(n_years - 1)]]
      lag_stack <- tile_binary[[1:(n_years - 2)]]

      if (spatial_autocorrelation) {
        middle_neighbor <- rast(lapply(1:nlyr(middle_years), function(i) {
          focal(middle_years[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
        }))
        predictor_stack <- c(middle_years, lag_stack, middle_neighbor, tile_summary)
      } else {
        predictor_stack <- c(middle_years, lag_stack, tile_summary)
      }

      if (show_progress) {
        n_cells <- ncell(predictor_stack)
        pb <- txtProgressBar(min = 0, max = n_cells, style = 3, width = 50)

        pred_vals <- values(predictor_stack, mat = TRUE)
        pattern_vals <- numeric(n_cells)
        decrease_vals <- numeric(n_cells)
        increase_vals <- numeric(n_cells)

        for (cell_i in 1:n_cells) {
          if (!any(is.na(pred_vals[cell_i, ]))) {
            result <- classify_pixel_with_years(pred_vals[cell_i, ], n_middle,
                                                time_steps, fastcpd_params, alpha, use_neighbor = spatial_autocorrelation)
            pattern_vals[cell_i] <- result[1]
            decrease_vals[cell_i] <- result[2]
            increase_vals[cell_i] <- result[3]
          } else {
            pattern_vals[cell_i] <- NA
            decrease_vals[cell_i] <- NA
            increase_vals[cell_i] <- NA
          }
          if (cell_i %% 10 == 0) setTxtProgressBar(pb, cell_i)
        }
        close(pb)
        print("")

        tile_pattern <- rast(predictor_stack, nlyr = 1)
        tile_decrease <- rast(predictor_stack, nlyr = 1)
        tile_increase <- rast(predictor_stack, nlyr = 1)

        values(tile_pattern) <- pattern_vals
        values(tile_decrease) <- decrease_vals
        values(tile_increase) <- increase_vals

      } else {
        result_matrix <- app(predictor_stack,
                             fun = function(x) classify_pixel_with_years(x, n_middle, time_steps,
                                                                         fastcpd_params, alpha, use_neighbor = spatial_autocorrelation))

        tile_pattern <- result_matrix[[1]]
        tile_decrease <- result_matrix[[2]]
        tile_increase <- result_matrix[[3]]
      }

      writeRaster(tile_pattern, tile_file_pattern, overwrite = TRUE,
                  datatype = "INT1U", gdal = c("COMPRESS=LZW"))
      writeRaster(tile_decrease, tile_file_decrease, overwrite = TRUE,
                  datatype = "INT2S", gdal = c("COMPRESS=LZW"))
      writeRaster(tile_increase, tile_file_increase, overwrite = TRUE,
                  datatype = "INT2S", gdal = c("COMPRESS=LZW"))

      rm(tile_binary, tile_summary, middle_years, lag_stack,
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

  ### Merge tiles

  print("Merging tiles...")

  if (n_tiles == 1) {
    print("Single tile detected, skipping merge...")
    pattern_raster <- rast(tile_files_pattern[1])
    decrease_raster <- rast(tile_files_decrease[1])
    increase_raster <- rast(tile_files_increase[1])

  } else {
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

  ### Save results

  print("Saving results...")
  writeRaster(pattern_raster, pattern_file, overwrite = TRUE,
              datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  writeRaster(decrease_raster, decrease_file, overwrite = TRUE,
              datatype = "INT2S", gdal = c("COMPRESS=LZW"))
  writeRaster(increase_raster, increase_file, overwrite = TRUE,
              datatype = "INT2S", gdal = c("COMPRESS=LZW"))

  ### Visualize

  print("Generating plots...")

  par(mfrow = c(2, 2))

  terra::plot(pattern_raster,
       col = c("#730000", "#267300", "#B2B2B2", "#A3FF73", "#FF7F7F", "#A900E6", "#eed202"),
       main = paste("Pattern Classification\n", min(time_steps), "-", max(time_steps)),
       breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5),
       legend = FALSE)
  legend("topright",
         legend = c("Always Absent", "Always Present", "No Pattern",
                    "Increasing", "Decreasing", "Fluctuating", "Failed"),
         fill = c("#730000", "#267300", "#B2B2B2", "#A3FF73", "#FF7F7F", "#A900E6", "#eed202"),
         cex = 0.6)

  terra::plot(decrease_raster,
       main = paste("Year of First Decrease\n", min(time_steps), "-", max(time_steps)),
       col = rev(heat.colors(50)))

  terra::plot(increase_raster,
       main = paste("Year of First Increase\n", min(time_steps), "-", max(time_steps)),
       col = terrain.colors(50))

  par(mfrow = c(1, 1))

  ### Summary

  print("Analysis complete")
  print(paste("Period:", min(time_steps), "-", max(time_steps)))
  print(paste("Tiles:", n_tiles))
  print(paste("Spatial autocorrelation:", ifelse(spatial_autocorrelation, "ENABLED", "DISABLED")))

  pattern_freq <- freq(pattern_raster)
  if (!is.null(pattern_freq) && nrow(pattern_freq) > 0) {
    pattern_freq <- as.data.frame(pattern_freq)
    pattern_freq$proportion <- round(pattern_freq$count / sum(pattern_freq$count), 3)
    pattern_names <- c("Always Absent", "Always Present", "No Pattern",
                       "Increasing", "Decreasing", "Fluctuating", "Failed")
    pattern_freq$pattern <- pattern_names[pattern_freq$value]

    print("Pattern Classifications:")
    print(pattern_freq[, c("value", "pattern", "count", "proportion")])
  }

  dec_freq <- as.data.frame(freq(decrease_raster))
  dec_freq <- dec_freq[!is.na(dec_freq$value), ]
  if (nrow(dec_freq) > 0) {
    n_decreasing <- sum(dec_freq$count)
    print(paste("Decreasing pixels:", format(n_decreasing, big.mark = ",")))
    print(paste("Year range:", min(dec_freq$value), "-", max(dec_freq$value)))
    print(paste("Most common:", dec_freq$value[which.max(dec_freq$count)],
                paste0("(", format(max(dec_freq$count), big.mark = ","), " pixels)")))
  }

  inc_freq <- as.data.frame(freq(increase_raster))
  inc_freq <- inc_freq[!is.na(inc_freq$value), ]
  if (nrow(inc_freq) > 0) {
    n_increasing <- sum(inc_freq$count)
    print(paste("Increasing pixels:", format(n_increasing, big.mark = ",")))
    print(paste("Year range:", min(inc_freq$value), "-", max(inc_freq$value)))
    print(paste("Most common:", inc_freq$value[which.max(inc_freq$count)],
                paste0("(", format(max(inc_freq$count), big.mark = ","), " pixels)")))
  }

  ### Auto-cleanup tiles

  print("Cleaning up tile files...")
  unlink(tiles_dir, recursive = TRUE)
  print("Tiles removed")

  return(list(
    pattern = pattern_raster,
    year_decrease = decrease_raster,
    year_increase = increase_raster
  ))
}
