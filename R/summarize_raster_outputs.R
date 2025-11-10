#' Summarize prediction rasters into a binary stack and mean
#' @export
#' @param predictions_dir Directory with prediction rasters.
#' @param output_dir Optional output directory (defaults to \code{predictions_dir}).
#' @param file_pattern Regex to match prediction files.
#' @param overwrite Overwrite existing outputs.
#' @return A list with \code{binary_stack} (RasterStack) and \code{summary_raster} (RasterLayer).
#' @importFrom raster stack raster calc cellStats nlayers writeRaster
summarize_raster_outputs <- function(predictions_dir,
                                      output_dir = NULL,
                                      file_pattern = "Prediction_.*\\.tif$",
                                      overwrite = TRUE) {

  require(raster)

  ### Validate required parameters

  if (missing(predictions_dir)) {
    stop("ERROR: 'predictions_dir' is required but was not provided.")
  }

  if (is.null(predictions_dir) || predictions_dir == "") {
    stop("ERROR: 'predictions_dir' cannot be NULL or empty. Please provide a valid directory path containing prediction rasters.")
  }

  if (!dir.exists(predictions_dir)) {
    stop(paste0("ERROR: The specified predictions_dir does not exist: ", predictions_dir))
  }

  ### Check for raster files

  all_raster_files <- list.files(
    path = predictions_dir,
    pattern = "\\.(tif|tiff|grd|img|asc)$",
    full.names = FALSE,
    ignore.case = TRUE
  )

  if (length(all_raster_files) == 0) {
    stop(paste0("ERROR: No raster files found in the predictions_dir: ", predictions_dir))
  }

  print("Creating binary stack and mean summary raster...")

  ### Set output directory

  if (is.null(output_dir)) output_dir <- predictions_dir
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ### Define output paths

  binary_dir <- file.path(output_dir, "Binary_Rasters")
  summary_file <- file.path(output_dir, "Summary_Binary_Mean.tif")

  ### Skip existing outputs

  if (dir.exists(binary_dir) && file.exists(summary_file) && !overwrite) {
    print("Binary rasters and summary already exist. Reloading from disk...")
    binary_files <- list.files(binary_dir, pattern = "\\.tif$", full.names = TRUE)
    binary_stack <- stack(binary_files)
    summary_raster <- raster(summary_file)
    return(list(binary_stack = binary_stack, summary_raster = summary_raster))
  }

  ### Find prediction rasters

  prediction_files <- list.files(
    path = predictions_dir,
    pattern = file_pattern,
    full.names = TRUE
  )

  prediction_names <- list.files(
    path = predictions_dir,
    pattern = file_pattern,
    full.names = FALSE
  )

  if (length(prediction_files) == 0) {
    stop(paste0("ERROR: No prediction files found matching pattern '", file_pattern, "' in directory ", predictions_dir))
  }

  print(paste("Found", length(prediction_files), "prediction rasters"))

  ### Stack rasters

  print("Stacking rasters...")
  prediction_stack <- tryCatch(
    stack(prediction_files),
    error = function(e) {
      stop(paste0("ERROR: Failed to stack rasters. Original error: ", e$message))
    }
  )

  clean_names <- gsub("\\.tif$", "", prediction_names)
  names(prediction_stack) <- clean_names

  ### Create binary raster directory

  if (!dir.exists(binary_dir)) dir.create(binary_dir, recursive = TRUE)

  ### Binary reclassification and save

  print("Performing binary reclassification (max value = 1, else = 0)...")
  binary_files <- c()

  for (i in 1:nlayers(prediction_stack)) {
    layer <- prediction_stack[[i]]
    max_val <- cellStats(layer, max, na.rm = TRUE)
    binary_raster <- calc(layer, function(x) ifelse(x == max_val, 1, 0))

    binary_file <- file.path(binary_dir, paste0(clean_names[i], "_binary.tif"))
    writeRaster(binary_raster, binary_file, format = "GTiff", overwrite = TRUE)
    binary_files <- c(binary_files, binary_file)
  }

  ### Build stack from saved binary rasters

  binary_stack <- stack(binary_files)
  names(binary_stack) <- clean_names

  ### Calculate mean summary

  print("Calculating mean across all time steps...")
  summary_raster <- calc(binary_stack, mean, na.rm = TRUE)

  ### Save mean summary raster

  print(paste("Saving mean summary raster to:", summary_file))
  writeRaster(summary_raster, summary_file, format = "GTiff", overwrite = TRUE)

  print(paste("Complete! Created", nlayers(binary_stack), "layer binary stack in:"))
  print(binary_dir)

  return(list(
    binary_stack = binary_stack,
    summary_raster = summary_raster
  ))
}
