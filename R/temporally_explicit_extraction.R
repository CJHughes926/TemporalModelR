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
#' @importFrom raster raster extract
#' @importFrom sp coordinates proj4string CRS
#' @importFrom rgdal readOGR
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

  require(raster)
  require(sp)
  require(readr)
  require(dplyr)

  ### Input validation and conversion - points

  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) {
      stop(paste0("ERROR: File does not exist: ", points_sp))
    }

    file_ext <- tolower(tools::file_ext(points_sp))

    if (file_ext == "csv") {
      if (is.null(xcol)) {
        stop("ERROR: 'xcol' is required when reading CSV files.")
      }
      if (is.null(ycol)) {
        stop("ERROR: 'ycol' is required when reading CSV files.")
      }
      if (is.null(points_crs)) {
        stop("ERROR: 'points_crs' is required when reading CSV files.")
      }

      print(paste("Reading CSV file:", basename(points_sp)))
      points_data <- read.csv(points_sp, stringsAsFactors = FALSE)

      if (!xcol %in% names(points_data)) {
        stop(paste0("ERROR: Column '", xcol, "' not found in CSV. Available columns: ", paste(names(points_data), collapse = ", ")))
      }
      if (!ycol %in% names(points_data)) {
        stop(paste0("ERROR: Column '", ycol, "' not found in CSV. Available columns: ", paste(names(points_data), collapse = ", ")))
      }

      coords_matrix <- as.matrix(points_data[, c(xcol, ycol)])
      data_without_coords <- points_data

      points_sp_temp <- SpatialPointsDataFrame(
        coords = coords_matrix,
        data = data_without_coords,
        proj4string = if (is.numeric(points_crs)) {
          CRS(paste0("+init=epsg:", points_crs))
        } else if (inherits(points_crs, "CRS")) {
          points_crs
        } else if (is.character(points_crs)) {
          CRS(points_crs)
        } else {
          stop("ERROR: points_crs must be a character string (proj4 or WKT), CRS object, or numeric EPSG code")
        }
      )

      points_sp <- points_sp_temp

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- rgdal::readOGR(points_sp, verbose = FALSE)
    } else {
      stop(paste0("ERROR: Unsupported file format: ", file_ext, " Supported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (is.data.frame(points_sp)) {
    if (is.null(xcol)) {
      stop("ERROR: 'xcol' is required when providing a data frame.")
    }
    if (is.null(ycol)) {
      stop("ERROR: 'ycol' is required when providing a data frame.")
    }
    if (is.null(points_crs)) {
      stop("ERROR: 'points_crs' is required when providing a data frame.")
    }

    print("Converting data frame to SpatialPointsDataFrame...")

    if (!xcol %in% names(points_sp)) {
      stop(paste0("ERROR: Column '", xcol, "' not found in data frame. Available columns: ", paste(names(points_sp), collapse = ", ")))
    }
    if (!ycol %in% names(points_sp)) {
      stop(paste0("ERROR: Column '", ycol, "' not found in data frame. Available columns: ", paste(names(points_sp), collapse = ", ")))
    }

    coords_matrix <- as.matrix(points_sp[, c(xcol, ycol)])
    data_without_coords <- points_sp

    points_sp_temp <- SpatialPointsDataFrame(
      coords = coords_matrix,
      data = data_without_coords,
      proj4string = if (is.numeric(points_crs)) {
        CRS(paste0("+init=epsg:", points_crs))
      } else if (inherits(points_crs, "CRS")) {
        points_crs
      } else if (is.character(points_crs)) {
        CRS(points_crs)
      } else {
        stop("ERROR: points_crs must be a character string (proj4 or WKT), CRS object, or numeric EPSG code")
      }
    )

    points_sp <- points_sp_temp

  } else if (!inherits(points_sp, "SpatialPointsDataFrame")) {
    stop("ERROR: points_sp must be a SpatialPointsDataFrame, data frame with x/y columns, or file path")
  }

  ### Validate variable_patterns format

  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop(paste(
      "ERROR: variable_patterns must be a named vector.",
      "",
      "### Example format:",
      "#",
      "# variable_patterns <- c(",
      "#   \"VariableName1\" = \"VariableName1_YEAR\",",
      "#   \"VariableName2\" = \"VariableName2_MONTH_YEAR\",",
      "#   \"StaticVariable\" = \"StaticVariable\"",
      "# )",
      "#",
      "# - NAME (left): output column name",
      "# - VALUE (right): pattern to match in raster filenames",
      "# - Time-varying: use placeholders (YEAR, MONTH, etc.) matching time_cols",
      "# - Static: use simple pattern with no time placeholders",
      sep = "\n"
    ))
  }

  if (any(names(variable_patterns) == "")) {
    stop("ERROR: All elements in variable_patterns must be named")
  }

  ### Validate time_cols

  if (missing(time_cols)) {
    stop(paste(
      "ERROR: time_cols is required. Provide a character vector of time column names.",
      "",
      "### Examples:",
      "# time_cols = \"Year\"",
      "# time_cols = c(\"Year\", \"Month\")",
      "#",
      "# Must match: (1) column names in points data, (2) placeholders in variable_patterns",
      sep = "\n"
    ))
  }

  if (!is.character(time_cols) || length(time_cols) == 0) {
    stop("ERROR: time_cols must be a character vector with at least one time column name")
  }

  missing_cols <- setdiff(time_cols, names(points_sp@data))
  if (length(missing_cols) > 0) {
    stop(paste0("ERROR: The following time_cols are missing from the input data: ",
                paste(missing_cols, collapse = ", "),
                " Available columns: ",
                paste(names(points_sp@data), collapse = ", ")))
  }

  ### Validate raster directory

  if (!dir.exists(raster_dir)) {
    stop(paste0("ERROR: raster_dir does not exist: ", raster_dir))
  }

  all_files <- list.files(path = raster_dir, pattern = "tif",
                          recursive = TRUE, full.names = TRUE)

  if (length(all_files) == 0) {
    stop(paste0("ERROR: No .tif files found in raster_dir: ", raster_dir))
  }

  print(paste("Found", length(all_files), "raster files"))

  ### Create output directory

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print(paste("Processing", nrow(points_sp@data), "points"))

  ### Determine static vs dynamic variables

  static_vars <- c()
  dynamic_vars <- c()

  for (var_name in names(variable_patterns)) {
    pattern <- variable_patterns[var_name]
    has_time_component <- any(sapply(time_cols, function(tc) {
      grepl(tc, pattern, ignore.case = TRUE)
    }))

    if (has_time_component) {
      dynamic_vars <- c(dynamic_vars, var_name)
    } else {
      static_vars <- c(static_vars, var_name)
    }
  }

  ### Validate time_cols match variable_patterns

  if (length(dynamic_vars) > 0) {
    pattern_time_placeholders <- c()
    for (var_name in dynamic_vars) {
      pattern <- variable_patterns[var_name]
      pattern_parts <- strsplit(pattern, "_")[[1]]

      for (part in pattern_parts) {
        if (toupper(part) %in% toupper(time_cols)) {
          pattern_time_placeholders <- c(pattern_time_placeholders, toupper(part))
        }
      }
    }

    pattern_time_placeholders <- unique(pattern_time_placeholders)

    time_cols_upper <- toupper(time_cols)
    missing_in_patterns <- time_cols_upper[!time_cols_upper %in% pattern_time_placeholders]
    extra_in_patterns <- pattern_time_placeholders[!pattern_time_placeholders %in% time_cols_upper]

    if (length(missing_in_patterns) > 0 && length(pattern_time_placeholders) > 0) {
      warning(paste(
        "WARNING: time_cols includes columns not found as placeholders in variable_patterns:",
        paste(missing_in_patterns, collapse = ", "),
        "This may indicate a mismatch between your time_cols and variable_patterns.",
        "",
        "Your time_cols:", paste(time_cols, collapse = ", "),
        "Placeholders found in variable_patterns:", paste(pattern_time_placeholders, collapse = ", ")
      ))
    }

    if (length(extra_in_patterns) > 0) {
      stop(paste(
        "ERROR: variable_patterns contain time placeholders not in time_cols:",
        paste(extra_in_patterns, collapse = ", "),
        "",
        "Your time_cols:", paste(time_cols, collapse = ", "),
        "Placeholders in patterns:", paste(pattern_time_placeholders, collapse = ", "),
        "",
        "### Fix by adding missing columns to time_cols:",
        "# time_cols = c(\"", paste(c(time_cols, tolower(extra_in_patterns)), collapse = "\", \""), "\")",
        sep = "\n"
      ))
    }
  }

  ### Check for missing values in time columns

  n_original <- nrow(points_sp@data)
  for (tc in time_cols) {
    n_missing <- sum(is.na(points_sp@data[[tc]]))
    if (n_missing > 0) {
      pct_missing <- round(n_missing / n_original * 100, 2)
      warning(paste0("WARNING: ", n_missing, " rows (", pct_missing,
                     "%) have missing values in time column '", tc, "'"))
    }
  }

  print(paste("Dynamic variables:", ifelse(length(dynamic_vars) > 0,
                                           paste(dynamic_vars, collapse = ", "),
                                           "none")))
  print(paste("Static variables:", ifelse(length(static_vars) > 0,
                                          paste(static_vars, collapse = ", "),
                                          "none")))

  ### Initialize variable columns

  for (var in names(variable_patterns)) {
    points_sp@data[[var]] <- NA
  }

  ### Extract static variables

  static_vars_extracted <- 0

  if (length(static_vars) > 0) {
    print("Extracting static variables...")
    static_rasters_found <- c()

    for (var_name in static_vars) {
      search_pattern <- variable_patterns[var_name]

      first_part <- var_name
      pattern_parts <- strsplit(search_pattern, "_")[[1]]

      matching_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)

        starts_correctly <- grepl(paste0("^", first_part), filename, ignore.case = TRUE)

        all_parts_present <- all(sapply(pattern_parts, function(part) {
          grepl(part, filename, ignore.case = TRUE)
        }))

        starts_correctly && all_parts_present
      })]

      if (length(matching_files) > 0) {
        raster_file <- matching_files[1]
        static_rasters_found <- c(static_rasters_found, basename(raster_file))
        raster_layer <- raster(raster_file)
        values <- raster::extract(raster_layer, points_sp)
        points_sp@data[[var_name]] <- values
        static_vars_extracted <- static_vars_extracted + 1
      } else {
        print(paste("WARNING: No file found for", var_name, "matching pattern:", search_pattern))
      }
    }

    if (length(static_rasters_found) > 0) {
      print(paste("Found", length(static_rasters_found), "static raster(s):"))
      for (rf in static_rasters_found) {
        print(paste("  -", rf))
      }
    }
  }

  ### Extract dynamic variables

  dynamic_vars_extracted <- 0

  if (length(dynamic_vars) > 0) {
    print("Extracting dynamic variables...")

    invisible(time_combinations <- points_sp@data %>%
                select(all_of(time_cols)) %>%
                distinct() %>%
                arrange(across(all_of(time_cols))))

    print(paste("Extracting values for", nrow(time_combinations), "time periods..."))

    for (i in 1:nrow(time_combinations)) {
      time_values <- time_combinations[i, , drop = FALSE]

      time_filter <- rep(TRUE, nrow(points_sp@data))
      for (tc in time_cols) {
        time_filter <- time_filter & (points_sp@data[[tc]] == time_values[[1, tc]])
      }
      points_subset <- points_sp[time_filter, ]

      time_label <- paste(sapply(time_cols, function(tc) {
        paste0(tc, "=", time_values[[1, tc]])
      }), collapse = ", ")
      n_points <- sum(time_filter)
      print(paste0("Processing time period ", i, "/", nrow(time_combinations), ": ", time_label, " (", n_points, " points)"))

      rasters_found <- c()

      for (var_name in dynamic_vars) {
        search_pattern <- variable_patterns[var_name]

        time_value_list <- list()
        for (tc in time_cols) {
          time_value_list[[tc]] <- as.character(time_values[[1, tc]])
          search_pattern <- gsub(tc, time_value_list[[tc]],
                                 search_pattern, ignore.case = TRUE)
        }

        pattern_parts <- strsplit(search_pattern, "_")[[1]]
        first_part <- var_name

        matching_files <- all_files[sapply(all_files, function(f) {
          filename <- basename(f)

          starts_correctly <- grepl(paste0("^", first_part), filename, ignore.case = TRUE)

          all_parts_present <- all(sapply(pattern_parts, function(part) {
            grepl(part, filename, ignore.case = TRUE)
          }))

          time_values_bounded <- all(sapply(names(time_value_list), function(tc) {
            time_val <- time_value_list[[tc]]
            grepl(paste0("(^|_)", time_val, "(_|\\.)"), filename, ignore.case = TRUE)
          }))

          starts_correctly && all_parts_present && time_values_bounded
        })]

        if (length(matching_files) > 0) {
          raster_file <- matching_files[1]
          rasters_found <- c(rasters_found, basename(raster_file))
          raster_layer <- raster(raster_file)
          values <- raster::extract(raster_layer, points_subset)
          points_sp@data[time_filter, var_name] <- values
          if (i == 1) {
            dynamic_vars_extracted <- dynamic_vars_extracted + 1
          }
        } else {
          print(paste("WARNING: No file found for", var_name, "matching pattern:", search_pattern))
        }
      }

      if (length(rasters_found) > 0) {
        print(paste("Found", length(rasters_found), "raster(s):"))
        for (rf in rasters_found) {
          print(paste("  -", rf))
        }
      }
    }
  }

  ### Check if no variables were extracted

  total_vars_extracted <- static_vars_extracted + dynamic_vars_extracted
  if (total_vars_extracted == 0) {
    stop(paste(
      "ERROR: No raster files could be matched for any variables.",
      "",
      "### Example variable_patterns:",
      "# variable_patterns <- c(",
      "#   \"VariableName1\" = \"VariableName1_YEAR\",",
      "#   \"VariableName2\" = \"VariableName2_MONTH_YEAR\",",
      "#   \"StaticVariable\" = \"StaticVariable\"",
      "# )",
      "#",
      "# Ensure patterns match your actual raster filenames",
      sep = "\n"
    ))
  }

  if (static_vars_extracted > 0 && dynamic_vars_extracted == 0 && length(dynamic_vars) > 0) {
    stop(paste(
      "ERROR: Only static variables matched. No raster files found for dynamic variables.",
      "",
      "### Check that:",
      "# 1. Raster filenames contain time values (e.g., 2020, 2021)",
      "# 2. variable_patterns use correct placeholders matching time_cols",
      "# 3. Time values in filenames match values in your points data",
      sep = "\n"
    ))
  }

  print("Value extraction complete")

  ### Save raw values

  if (save_raw) {
    raw_output_file <- file.path(output_dir, paste0(output_prefix, "_Raw_Values.csv"))
    write.csv(points_sp@data, raw_output_file, row.names = FALSE)
    print(paste("Raw values saved to:", basename(raw_output_file)))
  }

  ### Calculate scaling parameters

  print("Calculating scaling parameters...")
  scaling_params <- data.frame(
    variable = character(),
    mean = numeric(),
    sd = numeric(),
    stringsAsFactors = FALSE
  )

  all_vars <- names(variable_patterns)
  for (var_name in all_vars) {
    values <- points_sp@data[[var_name]]
    values <- values[!is.na(values)]

    if (length(values) > 0) {
      var_mean <- mean(values)
      var_sd <- sd(values)

      scaling_params <- rbind(scaling_params, data.frame(
        variable = var_name,
        mean = var_mean,
        sd = var_sd,
        stringsAsFactors = FALSE
      ))
    } else {
      print(paste("WARNING: No valid values found for", var_name))
    }
  }

  ### Save scaling parameters

  if (save_scaling_params) {
    params_file <- file.path(output_dir, paste0(output_prefix, "_Scaling_Parameters.csv"))
    write.csv(scaling_params, params_file, row.names = FALSE)
    print(paste("Scaling parameters saved to:", basename(params_file)))
  }

  ### Apply scaling and save

  if (save_scaled) {
    print("Applying scaling...")
    scaled_data <- points_sp@data

    for (var_name in all_vars) {
      var_params <- scaling_params[scaling_params$variable == var_name, ]

      if (nrow(var_params) > 0) {
        var_mean <- var_params$mean
        var_sd <- var_params$sd
        scaled_data[[var_name]] <- (scaled_data[[var_name]] - var_mean) / var_sd
      } else {
        print(paste("WARNING: No scaling parameters found for", var_name))
      }
    }

    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    write.csv(scaled_data, scaled_output_file, row.names = FALSE)
    print(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  print("Processing complete!")
  return(invisible(NULL))
}
