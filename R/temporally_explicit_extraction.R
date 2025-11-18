#' Extract Time-Aligned Environmental Values at Species Occurrences
#'
#' Preprocessing function that extracts raster values to species occurrence
#' records based on temporal components. Matches environmental layers to
#' occurrence timestamps and computes scaling parameters for standardization.
#'
#' @param points_sp sf object, SpatialPointsDataFrame, file path to
#'   .csv/.shp/.geojson/.gpkg, or data frame with coordinate columns.
#' @param raster_dir Character. Directory containing source raster files (.tif).
#' @param variable_patterns Named character vector mapping variable names to
#'   filename patterns. Time placeholders (e.g., YEAR, MONTH) must match
#'   time_cols. Static variables have no placeholders.
#' @param time_cols Character vector of time column names present in the point
#'   data (e.g., c("YEAR"), c("YEAR", "MONTH")).
#' @param xcol Character. Coordinate column name when reading CSV or data frame
#'   inputs. Optional.
#' @param ycol Character. Coordinate column name when reading CSV or data frame
#'   inputs. Optional.
#' @param points_crs Character or CRS object. CRS when reading CSV or data frame
#'   inputs. Optional.
#' @param output_dir Character. Directory to write output files.
#' @param output_prefix Character. Prefix for output filenames. Default is
#'   "temp_explicit_df".
#' @param save_raw Logical. If TRUE, writes raw extracted values CSV. If FALSE,
#'   skips raw values output. Default is TRUE.
#' @param save_scaled Logical. If TRUE, writes z-scaled values CSV. If FALSE,
#'   skips scaled values output. Default is TRUE.
#' @param save_scaling_params Logical. If TRUE, writes CSV of per-variable means
#'   and standard deviations. If FALSE, skips scaling parameters output. Default
#'   is TRUE.
#'
#' @details
#' Dynamic variables are detected when patterns contain placeholders matching
#' time_cols. Static variables are extracted once. Raster files are matched by
#' substituting time values into pattern placeholders.
#'
#' Output CSV files are written to output_dir containing raw values, scaled
#' values, and scaling parameters.
#'
#' Scaling parameters (mean and standard deviation) are computed across all
#' occurrence records for each variable. These parameters should be used with
#' \code{\link{scale_rasters}} to standardize prediction layers.
#'
#' @seealso
#' Preprocessing: \code{\link{spatiotemporal_rarefication}},
#' \code{\link{scale_rasters}}, \code{\link{spatiotemporal_partition}}
#'
#' @examples
#' \dontrun{
#' variable_patterns <- c(
#'   "bio1" = "bio1_YEAR",
#'   "temp" = "temp_YEAR_MONTH",
#'   "elevation" = "elevation"
#' )
#'
#' temporally_explicit_extraction(
#'   points_sp = "occurrences.csv",
#'   raster_dir = "environmental_layers/",
#'   variable_patterns = variable_patterns,
#'   time_cols = c("YEAR", "MONTH"),
#'   xcol = "longitude",
#'   ycol = "latitude",
#'   points_crs = "EPSG:4326",
#'   output_dir = "extracted_values/"
#' )
#' }
#'
#' @export
#' @importFrom terra rast extract crs vect
#' @importFrom sf st_as_sf st_transform st_read st_drop_geometry st_coordinates
#' @importFrom dplyr select distinct arrange across all_of
#' @importFrom utils read.csv write.csv
temporally_explicit_extraction <- function(points_sp,
                                           raster_dir,
                                           variable_patterns,
                                           time_cols,
                                           xcol = NULL,
                                           ycol = NULL,
                                           points_crs = NULL,
                                           output_dir,
                                           output_prefix = "temp_explicit_df",
                                           save_raw = TRUE,
                                           save_scaled = TRUE,
                                           save_scaling_params = TRUE) {

  require(terra)
  require(sf)
  require(dplyr)
  require(readr)

  ### Input validation and conversion - points
  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) stop(paste0("ERROR: File does not exist: ", points_sp))
    file_ext <- tolower(tools::file_ext(points_sp))

    if (file_ext == "csv") {
      if (is.null(xcol)) stop("ERROR: 'xcol' is required when reading CSV files.")
      if (is.null(ycol)) stop("ERROR: 'ycol' is required when reading CSV files.")
      if (is.null(points_crs)) stop("ERROR: 'points_crs' is required when reading CSV files.")

      print(paste("Reading CSV file:", basename(points_sp)))
      points_data <- read.csv(points_sp, stringsAsFactors = FALSE)
      if (!xcol %in% names(points_data)) stop(paste0("ERROR: Column '", xcol, "' not found in CSV."))
      if (!ycol %in% names(points_data)) stop(paste0("ERROR: Column '", ycol, "' not found in CSV."))

      points_sp <- st_as_sf(points_data, coords = c(xcol, ycol), crs = points_crs)

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- st_read(points_sp, quiet = TRUE)
    } else {
      stop(paste0("ERROR: Unsupported file format: ", file_ext))
    }

  } else if (is.data.frame(points_sp)) {
    if (is.null(xcol) || is.null(ycol) || is.null(points_crs)) stop("ERROR: xcol, ycol, points_crs required for data frame")
    points_sp <- st_as_sf(points_sp, coords = c(xcol, ycol), crs = points_crs)
  } else if (inherits(points_sp, "SpatialPointsDataFrame")) {
    warning("Converting SpatialPointsDataFrame to sf object")
    points_sp <- st_as_sf(points_sp)
  } else if (!inherits(points_sp, "sf")) {
    stop("ERROR: points_sp must be an sf object, data.frame, SpatialPointsDataFrame, or file path")
  }

  ### Validate variable_patterns
  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) stop("ERROR: variable_patterns must be a named vector")
  if (any(names(variable_patterns) == "")) stop("ERROR: All elements in variable_patterns must be named")

  ### Validate time_cols
  if (missing(time_cols) || !is.character(time_cols) || length(time_cols) == 0) stop("ERROR: time_cols must be a character vector with at least one column")
  missing_cols <- setdiff(time_cols, names(points_sp))
  if (length(missing_cols) > 0) stop(paste0("ERROR: time_cols missing from data: ", paste(missing_cols, collapse=", ")))

  ### Validate raster directory
  if (!dir.exists(raster_dir)) stop(paste0("ERROR: raster_dir does not exist: ", raster_dir))
  all_files <- list.files(path = raster_dir, pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) stop(paste0("ERROR: No .tif files found in raster_dir: ", raster_dir))
  print(paste("Found", length(all_files), "raster files"))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  print(paste("Processing", nrow(points_sp), "points"))

  ### Determine static vs dynamic variables
  static_vars <- c(); dynamic_vars <- c()
  for (var_name in names(variable_patterns)) {
    pattern <- variable_patterns[var_name]
    has_time_component <- any(sapply(time_cols, function(tc) grepl(tc, pattern, ignore.case=TRUE)))
    if (has_time_component) dynamic_vars <- c(dynamic_vars, var_name) else static_vars <- c(static_vars, var_name)
  }

  if (length(dynamic_vars) > 0) {
    pattern_time_placeholders <- unique(unlist(lapply(dynamic_vars, function(var_name) {
      parts <- strsplit(variable_patterns[var_name], "_")[[1]]
      parts[toupper(parts) %in% toupper(time_cols)]
    })))
    missing_in_patterns <- toupper(time_cols)[!toupper(time_cols) %in% pattern_time_placeholders]
    if (length(missing_in_patterns) > 0) warning("time_cols includes columns not found in variable_patterns: ", paste(missing_in_patterns, collapse=", "))
  }

  ### Initialize variable columns
  for (var in names(variable_patterns)) points_sp[[var]] <- NA

  ### Extract static variables
  if (length(static_vars) > 0) {
    print("Extracting static variables...")
    for (var_name in static_vars) {
      pattern_parts <- strsplit(variable_patterns[var_name], "_")[[1]]
      matching_files <- all_files[sapply(all_files, function(f) all(sapply(pattern_parts, function(p) grepl(p, basename(f), ignore.case=TRUE))))]

      if (length(matching_files) > 1) {
        stop(paste0("ERROR: Multiple raster files found for static variable '", var_name,
                    "' matching pattern '", variable_patterns[var_name], "'.\nMatching files:\n",
                    paste(matching_files, collapse = "\n")))
      }

      if (length(matching_files) == 1) {
        raster_layer <- rast(matching_files[1])
        points_sp[[var_name]] <- terra::extract(raster_layer, vect(points_sp))[,2]
      } else {
        warning("No file found for ", var_name, " matching pattern: ", variable_patterns[var_name])
      }
    }
  }

  ### Extract dynamic variables
  if (length(dynamic_vars) > 0) {
    print("Extracting dynamic variables...")
    time_combinations <- points_sp %>% st_drop_geometry() %>% dplyr::select(all_of(time_cols)) %>% distinct() %>% arrange(across(all_of(time_cols)))
    print(paste("Extracting values for", nrow(time_combinations), "time periods"))
    for (i in 1:nrow(time_combinations)) {
      time_values <- time_combinations[i, , drop=FALSE]
      time_filter <- Reduce(`&`, lapply(time_cols, function(tc) points_sp[[tc]] == time_values[[tc]]))
      points_subset <- points_sp[time_filter, ]
      for (var_name in dynamic_vars) {
        search_pattern <- variable_patterns[var_name]
        for (tc in time_cols) search_pattern <- gsub(tc, as.character(time_values[[tc]]), search_pattern, ignore.case=TRUE)
        pattern_parts <- strsplit(search_pattern, "_")[[1]]
        matching_files <- all_files[sapply(all_files, function(f) all(sapply(pattern_parts, function(p) grepl(p, basename(f), ignore.case=TRUE))))]

        if (length(matching_files) > 1) {
          stop(paste0("ERROR: Multiple raster files found for dynamic variable '", var_name,
                      "' matching pattern '", search_pattern, "'.\nMatching files:\n",
                      paste(matching_files, collapse = "\n")))
        }

        if (length(matching_files) == 1) {
          raster_layer <- rast(matching_files[1])
          points_sp[[var_name]][time_filter] <- terra::extract(raster_layer, vect(points_subset))[,2]
        } else {
          warning("No file found for ", var_name, " matching pattern: ", search_pattern)
        }
      }
    }
  }

  ### Add coordinates back
  coords_df <- st_coordinates(points_sp) %>% as.data.frame()
  names(coords_df) <- c("x", "y")
  points_with_coords <- cbind(st_drop_geometry(points_sp), coords_df)

  ### Save raw values
  if (save_raw) {
    raw_output_file <- file.path(output_dir, paste0(output_prefix, "_Raw_Values.csv"))
    write.csv(points_with_coords, raw_output_file, row.names=FALSE)
    print(paste("Raw values saved to:", basename(raw_output_file)))
  }

  ### Calculate scaling parameters
  print("Calculating scaling parameters...")
  scaling_params <- data.frame(variable=character(), mean=numeric(), sd=numeric(), stringsAsFactors=FALSE)
  for (var_name in names(variable_patterns)) {
    values <- points_sp[[var_name]]
    values <- values[!is.na(values)]
    if (length(values) > 0) {
      scaling_params <- rbind(scaling_params, data.frame(variable=var_name, mean=mean(values), sd=sd(values), stringsAsFactors=FALSE))
    } else warning("No valid values found for ", var_name)
  }

  ### Save scaling parameters
  if (save_scaling_params) {
    params_file <- file.path(output_dir, paste0(output_prefix, "_Scaling_Parameters.csv"))
    write.csv(scaling_params, params_file, row.names=FALSE)
    print(paste("Scaling parameters saved to:", basename(params_file)))
  }

  ### Apply scaling and save
  if (save_scaled) {
    print("Applying scaling...")
    scaled_data <- st_drop_geometry(points_sp)
    for (var_name in names(variable_patterns)) {
      var_params <- scaling_params[scaling_params$variable==var_name, ]
      if (nrow(var_params) > 0) scaled_data[[var_name]] <- (scaled_data[[var_name]] - var_params$mean) / var_params$sd
      else warning("No scaling parameters found for ", var_name)
    }
    ### Add coordinates to scaled output
    scaled_data <- cbind(scaled_data, coords_df)
    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    write.csv(scaled_data, scaled_output_file, row.names=FALSE)
    print(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  print("Processing complete!")
  return(invisible(NULL))
}
