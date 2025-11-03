summarize_raster_outputs <- function(predictions_dir,
                                        output_dir = NULL,
                                        file_pattern = "Prediction_.*\\.tif$",
                                        overwrite = TRUE) {

  require(raster)

  ### Validate required parameters
  if (missing(predictions_dir)) {
    stop("Error: 'predictions_dir' is required but was not provided.\n",
         "Usage: create_binary_stack_summary(predictions_dir = 'path/to/predictions')")
  }

  if (is.null(predictions_dir) || predictions_dir == "") {
    stop("Error: 'predictions_dir' cannot be NULL or empty.\n",
         "Please provide a valid directory path containing prediction rasters.")
  }

  if (!dir.exists(predictions_dir)) {
    stop("Error: The specified predictions_dir does not exist:\n  ",
         predictions_dir, "\n",
         "Please check the path and try again.")
  }

  ### Check if directory contains any raster files at all
  all_raster_files <- list.files(
    path = predictions_dir,
    pattern = "\\.(tif|tiff|grd|img|asc)$",
    full.names = FALSE,
    ignore.case = TRUE
  )

  if (length(all_raster_files) == 0) {
    stop("Error: No raster files found in the predictions_dir:\n  ",
         predictions_dir, "\n",
         "  The directory exists but contains no raster files (.tif, .tiff, .grd, .img, .asc)\n",
         "  Please verify you're pointing to the correct directory.")
  }

  print("Creating binary stack and mean summary raster...")

  ### Set output directory (default to predictions_dir)
  if (is.null(output_dir)) {
    output_dir <- predictions_dir
  }

  ### Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  ### Define output file paths
  stack_file <- file.path(output_dir, "Binary_Stack.RData")
  summary_file <- file.path(output_dir, "Summary_Binary_Mean.tif")

  ### Check if outputs already exist
  if (file.exists(stack_file) && file.exists(summary_file) && !overwrite) {
    print("Binary stack and summary already exist. Loading from files...")
    load(stack_file)
    summary_raster <- raster(summary_file)

    return(list(
      binary_stack = binary_stack,
      summary_raster = summary_raster
    ))
  }

  ### Find prediction files
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

  print(paste("Found", length(prediction_files), "prediction rasters"))

  if (length(prediction_files) == 0) {
    stop("Error: No prediction files found matching pattern '", file_pattern, "'\n",
         "  Directory searched: ", predictions_dir, "\n",
         "  Raster files found in directory: ", length(all_raster_files), "\n",
         "  Please verify:\n",
         "    1. The file_pattern is correct: ", file_pattern, "\n",
         "    2. Files match the expected naming convention\n",
         "  Available raster files:\n    ",
         paste(head(all_raster_files, 10), collapse = "\n    "))
  }

  ### Stack rasters
  print("Stacking rasters...")
  prediction_stack <- tryCatch(
    stack(prediction_files),
    error = function(e) {
      stop("Error: Failed to stack rasters.\n",
           "  Original error: ", e$message, "\n",
           "  This may indicate corrupted or incompatible raster files.")
    }
  )

  clean_names <- gsub("\\.tif$", "", prediction_names)
  names(prediction_stack) <- clean_names

  ### Binary reclassification (max value -> 1, else -> 0)
  print("Performing binary reclassification (max value = 1, else = 0)...")
  binary_layers <- list()

  for (i in 1:nlayers(prediction_stack)) {
    layer <- prediction_stack[[i]]
    max_val <- cellStats(layer, max, na.rm = TRUE)
    binary_layers[[i]] <- calc(layer, function(x) ifelse(x == max_val, 1, 0))
  }

  binary_stack <- stack(binary_layers)
  names(binary_stack) <- clean_names

  ### Calculate mean across all time steps
  print("Calculating mean across all time steps...")
  summary_raster <- calc(binary_stack, mean, na.rm = TRUE)

  ### Save binary stack as RData file
  print(paste("Saving binary stack to:", stack_file))
  save(binary_stack, file = stack_file)

  ### Save mean summary raster as GeoTIFF
  print(paste("Saving mean summary raster to:", summary_file))
  writeRaster(summary_raster, summary_file, format = "GTiff", overwrite = TRUE)

  print(paste("Complete! Created", nlayers(binary_stack), "layer binary stack"))

  return(list(
    binary_stack = binary_stack,
    summary_raster = summary_raster
  ))
}
