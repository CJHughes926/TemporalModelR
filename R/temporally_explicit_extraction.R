#' Extract time-aligned raster values at points and compute scaling
#' @export
#' @description For each unique combination in \code{time_cols}, matches rasters by
#'   \code{variable_patterns}, extracts values at \code{points_sp}, writes raw values,
#'   per-variable scaling parameters (mean, sd), and optionally a z-scaled table.
#'
#' @param points_sp A \code{SpatialPointsDataFrame}, a path to a CSV/SHP/GEOJSON/GPKG,
#'   or a data.frame with coordinate columns (see \code{xcol}, \code{ycol}, \code{points_crs}).
#' @param raster_dir Directory containing source rasters (.tif).
#' @param variable_patterns Named character vector mapping variable names to filename
#'   patterns. Time placeholders (e.g., \code{YEAR}, \code{MONTH}) must match \code{time_cols}.
#'   Static variables have no placeholders (e.g., \code{"elevation" = "elevation"}).
#' @param time_cols Character vector of time column names present in \code{points_sp} data.
#' @param xcol,ycol Optional. Names of longitude/latitude columns when \code{points_sp} is a CSV or data.frame.
#' @param points_crs Optional. CRS string (e.g., \code{"+proj=longlat +datum=WGS84"}) when
#'   \code{points_sp} is a CSV or data.frame.
#' @param output_dir Directory to write outputs.
#' @param output_prefix Prefix for output filenames.
#' @param save_raw Logical; write raw extracted values CSV.
#' @param save_scaled Logical; write z-scaled values CSV.
#' @param save_scaling_params Logical; write per-variable means/SDs CSV.
#'
#' @returns Invisibly returns \code{NULL}. Side effects: writes
#'   \code{*_Raw_Values.csv}, \code{*_Scaling_Parameters.csv}, and optionally \code{*_Scaled_Values.csv}.
#'
#' @details
#' Dynamic variables are detected when their pattern contains any of \code{time_cols}
#' placeholders; static variables are extracted once. Filenames are matched by
#' variable prefix and placeholder substitution (bounded by underscores or extension).
#'
#' @seealso \code{\link{scale_rasters}} for applying saved scaling to rasters.
#'
#' @importFrom terra rast extract crs vect
#' @importFrom sf st_as_sf st_transform
#' @importFrom readr read_csv
#' @importFrom dplyr select distinct arrange across all_of
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
    time_combinations <- points_sp %>% st_drop_geometry() %>% select(all_of(time_cols)) %>% distinct() %>% arrange(across(all_of(time_cols)))
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

  ### Save raw values
  if (save_raw) {
    raw_output_file <- file.path(output_dir, paste0(output_prefix, "_Raw_Values.csv"))
    write.csv(st_drop_geometry(points_sp), raw_output_file, row.names=FALSE)
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
    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    write.csv(scaled_data, scaled_output_file, row.names=FALSE)
    print(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  print("Processing complete!")
  return(invisible(NULL))
}
