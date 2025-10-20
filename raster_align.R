raster_align <- function(input_dir, 
                         output_dir, 
                         reference_raster,
                         output_suffix = "_Masked_Updated",
                         pattern = ".*\\.tif$",
                         overwrite = FALSE) {
  
  # Load required library
  require(raster)
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load reference raster
  reference_raster <- raster(reference_raster)
  
  # List all raster files
  all_tif_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  total_files <- length(all_tif_files)
  
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