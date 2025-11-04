#' Analyze temporal patterns in binary raster time series
#'
#' Splits a binary raster stack into tiles, classifies per-pixel temporal patterns (e.g., increase/decrease) using a chosen method, and writes pattern and change-year rasters.
#'
#' @param binary_stack RasterStack/Brick of binary layers across time.
#' @param summary_raster RasterLayer with per-pixel summary (e.g., mean/persistence).
#' @param time_steps Integer vector of time labels (same length as layers).
#' @param method Character; changepoint/selection criterion (e.g., "BIC", "MBIC", "MDL").
#' @param output_dir Output directory.
#' @param n_tiles_x,n_tiles_y Number of tiles in x/y.
#' @param overlap Overlap (in cells) between tiles.
#' @param alpha Significance level.
#' @param spatial_autocorrelation Logical; if TRUE (default), includes neighbor variable in analysis. If FALSE, spatial autocorrelation is not considered.
#' @param show_progress Logical; show progress bar.
#' @param estimate_time Logical; quick runtime estimate from sampled pixels.
#' @param overwrite Logical; overwrite existing outputs.
#'
#' @return A list with rasters: \code{pattern}, \code{year_decrease}, \code{year_increase}.
#'
#' @details Requires an auxiliary \code{classify_pixel_with_years()} function available in scope.
#'
#' @seealso \code{classify_pixel_with_years}
#'
#' @export
#' @importFrom raster raster extent res crop nlayers subset stack focal getValues values ncell mosaic writeRaster freq
#' @importFrom utils txtProgressBar setTxtProgressBar write.csv
#' @importFrom graphics par plot legend
#' @importFrom grDevices heat.colors terrain.colors

analyze_temporal_patterns <- function(
    binary_stack,
    summary_raster,
    time_steps,
    method,
    output_dir = "output",
    n_tiles_x = 1,
    n_tiles_y = 1,
    overlap = 10,
    alpha = 0.05,
    spatial_autocorrelation = TRUE,
    show_progress = TRUE,
    estimate_time = TRUE,
    overwrite = FALSE) {

  ### Check required parameters
  if (missing(binary_stack)) stop("binary_stack is required", call. = FALSE)
  if (missing(summary_raster)) stop("summary_raster is required", call. = FALSE)
  if (missing(time_steps)) stop("time_steps is required", call. = FALSE)
  if (missing(method)) stop("method is required (e.g., 'BIC', 'MBIC', 'MDL')", call. = FALSE)

  ### Setup directories
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  tiles_dir <- file.path(output_dir, "tiles")
  if (!dir.exists(tiles_dir)) dir.create(tiles_dir, recursive = TRUE)

  ### Define output files
  pattern_file <- file.path(output_dir, paste0("pattern_raster_", min(time_steps), "_", max(time_steps), ".tif"))
  decrease_file <- file.path(output_dir, paste0("year_first_decrease_", min(time_steps), "_", max(time_steps), ".tif"))
  increase_file <- file.path(output_dir, paste0("year_first_increase_", min(time_steps), "_", max(time_steps), ".tif"))

  if (file.exists(pattern_file) && file.exists(decrease_file) && file.exists(increase_file) && !overwrite) {
    cat("Output rasters exist. Set overwrite = TRUE to rerun.\n")
    return(list(
      pattern = raster(pattern_file),
      year_decrease = raster(decrease_file),
      year_increase = raster(increase_file)
    ))
  }

  ### Display spatial autocorrelation setting
  if (spatial_autocorrelation) {
    cat("\nSpatial autocorrelation: ENABLED (neighbor variable included)\n")
  } else {
    cat("\nSpatial autocorrelation: DISABLED (neighbor variable excluded)\n")
  }

  ### Calculate tile extents
  cat("\nCalculating tile extents...\n")

  full_ext <- extent(binary_stack)
  x_min <- full_ext@xmin
  x_max <- full_ext@xmax
  y_min <- full_ext@ymin
  y_max <- full_ext@ymax

  res_x <- res(binary_stack)[1]
  res_y <- res(binary_stack)[2]

  x_range <- x_max - x_min
  y_range <- y_max - y_min

  tile_width <- x_range / n_tiles_x
  tile_height <- y_range / n_tiles_y

  overlap_x <- overlap * res_x
  overlap_y <- overlap * res_y

  tile_extents <- list()
  tile_idx <- 1

  for (i in 1:n_tiles_y) {
    for (j in 1:n_tiles_x) {
      tile_x_min <- x_min + (j - 1) * tile_width
      tile_x_max <- x_min + j * tile_width
      tile_y_min <- y_min + (i - 1) * tile_height
      tile_y_max <- y_min + i * tile_height

      if (j > 1) tile_x_min <- tile_x_min - overlap_x
      if (j < n_tiles_x) tile_x_max <- tile_x_max + overlap_x
      if (i > 1) tile_y_min <- tile_y_min - overlap_y
      if (i < n_tiles_y) tile_y_max <- tile_y_max + overlap_y

      tile_extents[[tile_idx]] <- extent(tile_x_min, tile_x_max, tile_y_min, tile_y_max)
      tile_idx <- tile_idx + 1
    }
  }

  n_tiles <- length(tile_extents)
  cat(paste("Created", n_tiles, "tiles\n"))

  ### Estimate processing time
  if (estimate_time) {
    cat("\nEstimating processing time...\n")

    n_years <- nlayers(binary_stack)
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

        summary_vals <- getValues(tile_summary)
        valid_indices <- which(!is.na(summary_vals))

        if (length(valid_indices) > 0) {
          mean_vals <- summary_vals[valid_indices]
          n_quick <- sum(mean_vals < 0.01 | mean_vals > 0.99)
          n_complex <- sum(mean_vals >= 0.01 & mean_vals <= 0.99)

          total_quick <- total_quick + n_quick
          total_complex <- total_complex + n_complex

          if (is.null(time_per_pixel) && n_complex > 0) {
            cat("\nTiming sample pixels...       \n")

            middle_years <- raster::subset(tile_binary, 2:(n_years - 1))
            lag_stack <- raster::subset(tile_binary, 1:(n_years - 2))

            if (spatial_autocorrelation) {
              middle_neighbor <- stack(lapply(1:nlayers(middle_years), function(i) {
                focal(middle_years[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
              }))
              predictor_stack <- stack(middle_years, lag_stack, middle_neighbor, tile_summary)
            } else {
              predictor_stack <- stack(middle_years, lag_stack, tile_summary)
            }

            pred_vals_all <- getValues(predictor_stack)
            valid_pred <- which(!is.na(pred_vals_all[, 1]))

            if (spatial_autocorrelation) {
              mean_pred <- pred_vals_all[valid_pred, (3 * n_middle + 1)]
            } else {
              mean_pred <- pred_vals_all[valid_pred, (2 * n_middle + 1)]
            }

            complex_indices <- valid_pred[mean_pred >= 0.01 & mean_pred <= 0.99]

            sample_size <- min(100, length(complex_indices))
            sample_indices <- sample(complex_indices, size = sample_size, replace = FALSE)
            pred_vals_sample <- pred_vals_all[sample_indices, ]

            start_time <- Sys.time()
            for (i in 1:sample_size) {
              classify_pixel_with_years(pred_vals_sample[i, ], n_middle,
                                        time_steps, method, alpha, use_neighbor = spatial_autocorrelation)
            }
            end_time <- Sys.time()

            time_per_pixel <- as.numeric(difftime(end_time, start_time, units = "secs")) / sample_size

            rm(middle_years, lag_stack, predictor_stack,
               pred_vals_all, valid_pred, mean_pred, complex_indices, pred_vals_sample)
            if (spatial_autocorrelation) rm(middle_neighbor)
            gc(verbose = FALSE)
          }
        }

        rm(tile_binary, tile_summary, summary_vals, valid_indices)
        gc(verbose = FALSE)

      }, error = function(e) {
        if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
          stop("Memory error. Increase n_tiles_x and n_tiles_y to use smaller tiles.", call. = FALSE)
        } else {
          stop(e)
        }
      })
    }

    cat("\n")
    cat(paste0("Quick pixels (always absent/present): ", format(total_quick, big.mark = ","), "\n"))
    cat(paste0("Complex pixels (changepoint analysis): ", format(total_complex, big.mark = ","), "\n"))

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

      cat(paste0("\nEstimated time: ", format_time(lower_seconds), " - ", format_time(upper_seconds), "\n"))
      cat(paste0("(", round(time_per_pixel, 3), " sec/pixel + overhead)\n"))
    }
  } else {
    cat("\nSkipping time estimation. This may take a long time.\n")
    cat("Set estimate_time = TRUE to get time estimate.\n")
  }

  ### Process tiles
  cat("\nProcessing tiles...\n")

  tile_files_pattern <- character(n_tiles)
  tile_files_decrease <- character(n_tiles)
  tile_files_increase <- character(n_tiles)

  for (tile_i in 1:n_tiles) {
    cat(sprintf("\nTile %d/%d\n", tile_i, n_tiles))

    tile_ext <- tile_extents[[tile_i]]
    tile_file_pattern <- file.path(tiles_dir, paste0("tile_pattern_", tile_i, ".tif"))
    tile_file_decrease <- file.path(tiles_dir, paste0("tile_decrease_", tile_i, ".tif"))
    tile_file_increase <- file.path(tiles_dir, paste0("tile_increase_", tile_i, ".tif"))

    tile_files_pattern[tile_i] <- tile_file_pattern
    tile_files_decrease[tile_i] <- tile_file_decrease
    tile_files_increase[tile_i] <- tile_file_increase

    if (file.exists(tile_file_pattern) && file.exists(tile_file_decrease) &&
        file.exists(tile_file_increase) && !overwrite) {
      cat("Already processed, skipping\n")
      next
    }

    tryCatch({
      suppressWarnings({
        tile_binary <- crop(binary_stack, tile_ext)
        tile_summary <- crop(summary_raster, tile_ext)
      })

      n_years <- nlayers(tile_binary)
      n_middle <- n_years - 2

      middle_years <- raster::subset(tile_binary, 2:(n_years - 1))
      lag_stack <- raster::subset(tile_binary, 1:(n_years - 2))

      if (spatial_autocorrelation) {
        middle_neighbor <- stack(lapply(1:nlayers(middle_years), function(i) {
          focal(middle_years[[i]], w = matrix(1/9, 3, 3), fun = mean, na.rm = TRUE)
        }))
        predictor_stack <- stack(middle_years, lag_stack, middle_neighbor, tile_summary)
      } else {
        predictor_stack <- stack(middle_years, lag_stack, tile_summary)
      }

      if (show_progress) {
        n_cells <- ncell(predictor_stack)
        pb <- txtProgressBar(min = 0, max = n_cells, style = 3, width = 50)

        pred_vals <- getValues(predictor_stack)
        pattern_vals <- numeric(n_cells)
        decrease_vals <- numeric(n_cells)
        increase_vals <- numeric(n_cells)

        for (cell_i in 1:n_cells) {
          if (!any(is.na(pred_vals[cell_i, ]))) {
            result <- classify_pixel_with_years(pred_vals[cell_i, ], n_middle,
                                                time_steps, method, alpha, use_neighbor = spatial_autocorrelation)
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
        cat("\n")

        tile_pattern <- raster(predictor_stack, layer = 1)
        tile_decrease <- raster(predictor_stack, layer = 1)
        tile_increase <- raster(predictor_stack, layer = 1)

        values(tile_pattern) <- pattern_vals
        values(tile_decrease) <- decrease_vals
        values(tile_increase) <- increase_vals

      } else {
        result_matrix <- calc(predictor_stack,
                              fun = function(x) classify_pixel_with_years(x, n_middle, time_steps,
                                                                          method, alpha, use_neighbor = spatial_autocorrelation))

        tile_pattern <- raster::subset(result_matrix, 1)
        tile_decrease <- raster::subset(result_matrix, 2)
        tile_increase <- raster::subset(result_matrix, 3)
      }

      writeRaster(tile_pattern, tile_file_pattern, overwrite = TRUE,
                  datatype = "INT1U", options = c("COMPRESS=LZW"))
      writeRaster(tile_decrease, tile_file_decrease, overwrite = TRUE,
                  datatype = "INT2S", options = c("COMPRESS=LZW"))
      writeRaster(tile_increase, tile_file_increase, overwrite = TRUE,
                  datatype = "INT2S", options = c("COMPRESS=LZW"))

      rm(tile_binary, tile_summary, middle_years, lag_stack,
         predictor_stack, tile_pattern, tile_decrease, tile_increase)
      if (spatial_autocorrelation) rm(middle_neighbor)
      if (exists("pred_vals")) rm(pred_vals, pattern_vals, decrease_vals, increase_vals)
      gc(verbose = FALSE)

    }, error = function(e) {
      if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
        stop("Memory error on tile ", tile_i, ". Increase n_tiles_x and n_tiles_y.", call. = FALSE)
      } else {
        stop(e)
      }
    })
  }

  ### Merge tiles
  cat("\nMerging tiles...\n")

  tryCatch({
    tile_rasters_pattern <- lapply(tile_files_pattern, raster)
    tile_rasters_pattern$fun <- mean
    pattern_raster <- do.call(mosaic, tile_rasters_pattern)

    tile_rasters_decrease <- lapply(tile_files_decrease, raster)
    tile_rasters_decrease$fun <- mean
    decrease_raster <- do.call(mosaic, tile_rasters_decrease)

    tile_rasters_increase <- lapply(tile_files_increase, raster)
    tile_rasters_increase$fun <- mean
    increase_raster <- do.call(mosaic, tile_rasters_increase)

  }, error = function(e) {
    if (grepl("cannot allocate vector", e$message, ignore.case = TRUE)) {
      stop("Memory error merging tiles. Increase n_tiles_x and n_tiles_y.", call. = FALSE)
    } else {
      stop(e)
    }
  })

  ### Save results
  cat("\nSaving results...\n")
  writeRaster(pattern_raster, pattern_file, overwrite = TRUE,
              datatype = "INT1U", options = c("COMPRESS=LZW"))
  writeRaster(decrease_raster, decrease_file, overwrite = TRUE,
              datatype = "INT2S", options = c("COMPRESS=LZW"))
  writeRaster(increase_raster, increase_file, overwrite = TRUE,
              datatype = "INT2S", options = c("COMPRESS=LZW"))

  ### Visualize
  cat("\nGenerating plots...\n")

  par(mfrow = c(2, 2))

  plot(pattern_raster,
       col = c("darkgray", "lightgray", "gray", "green4", "red3", "purple", "orange"),
       main = paste("Pattern Classification\n", min(time_steps), "-", max(time_steps)),
       breaks = c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5),
       legend = FALSE)
  legend("topright",
         legend = c("Always Absent", "Always Present", "No Pattern",
                    "Increasing", "Decreasing", "Fluctuating", "Failed"),
         fill = c("darkgray", "lightgray", "gray", "green4", "red3", "purple", "orange"),
         cex = 0.6)

  plot(decrease_raster,
       main = paste("Year of First Decrease\n", min(time_steps), "-", max(time_steps)),
       col = rev(heat.colors(50)))

  plot(increase_raster,
       main = paste("Year of First Increase\n", min(time_steps), "-", max(time_steps)),
       col = terrain.colors(50))

  par(mfrow = c(1, 1))

  ### Summary
  cat("\n========================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("========================================\n")
  cat(paste("Period:", min(time_steps), "-", max(time_steps), "\n"))
  cat(paste("Tiles:", n_tiles, "\n"))
  cat(paste("Spatial autocorrelation:", ifelse(spatial_autocorrelation, "ENABLED", "DISABLED"), "\n\n"))

  pattern_freq <- freq(pattern_raster)
  if (!is.null(pattern_freq) && nrow(pattern_freq) > 0) {
    pattern_freq <- as.data.frame(pattern_freq)
    pattern_freq$proportion <- round(pattern_freq$count / sum(pattern_freq$count), 3)
    pattern_names <- c("Always Absent", "Always Present", "No Pattern",
                       "Increasing", "Decreasing", "Fluctuating", "Failed")
    pattern_freq$pattern <- pattern_names[pattern_freq$value]

    cat("Pattern Classifications:\n")
    print(pattern_freq[, c("value", "pattern", "count", "proportion")])
  }

  cat("\n")
  dec_freq <- freq(decrease_raster, useNA = "no")
  if (!is.null(dec_freq) && nrow(dec_freq) > 0) {
    dec_freq <- as.data.frame(dec_freq)
    n_decreasing <- sum(dec_freq$count)
    cat(paste("Decreasing pixels:", format(n_decreasing, big.mark = ","), "\n"))
    cat(paste("Year range:", min(dec_freq$value), "-", max(dec_freq$value), "\n"))
    cat(paste("Most common:", dec_freq$value[which.max(dec_freq$count)],
              "(", format(max(dec_freq$count), big.mark = ","), "pixels)\n"))
  }

  cat("\n")
  inc_freq <- freq(increase_raster, useNA = "no")
  if (!is.null(inc_freq) && nrow(inc_freq) > 0) {
    inc_freq <- as.data.frame(inc_freq)
    n_increasing <- sum(inc_freq$count)
    cat(paste("Increasing pixels:", format(n_increasing, big.mark = ","), "\n"))
    cat(paste("Year range:", min(inc_freq$value), "-", max(inc_freq$value), "\n"))
    cat(paste("Most common:", inc_freq$value[which.max(inc_freq$count)],
              "(", format(max(inc_freq$count), big.mark = ","), "pixels)\n"))
  }

  ### Cleanup option
  cat("\nDelete tile files? (y/n): ")
  response <- readLines(con = "stdin", n = 1)
  if (tolower(response) == "y") {
    unlink(tiles_dir, recursive = TRUE)
    cat("Tiles removed\n")
  } else {
    cat(paste("Tiles kept:", tiles_dir, "\n"))
  }

  return(list(
    pattern = pattern_raster,
    year_decrease = decrease_raster,
    year_increase = increase_raster
  ))
}
