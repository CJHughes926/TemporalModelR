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

  ### Input validation and conversion

  # Handle points_sp input
  if (is.character(points_sp)) {
    if (!file.exists(points_sp)) {
      stop(paste("Error: File does not exist:", points_sp))
    }

    file_ext <- tolower(tools::file_ext(points_sp))

    if (file_ext == "csv") {

      # Require xcol, ycol, and points_crs for CSV files
      if (is.null(xcol)) {
        stop("Error: 'xcol' is required when reading CSV files.")
      }
      if (is.null(ycol)) {
        stop("Error: 'ycol' is required when reading CSV files.")
      }
      if (is.null(points_crs)) {
        stop("Error: 'points_crs' is required when reading CSV files.")
      }

      print(paste("Reading CSV file:", basename(points_sp)))
      points_data <- read.csv(points_sp, stringsAsFactors = FALSE)

      if (!xcol %in% names(points_data)) {
        stop(paste0("Error: Column '", xcol, "' not found in CSV. Available columns: ",
                    paste(names(points_data), collapse = ", ")))
      }
      if (!ycol %in% names(points_data)) {
        stop(paste0("Error: Column '", ycol, "' not found in CSV. Available columns: ",
                    paste(names(points_data), collapse = ", ")))
      }

      coordinates(points_data) <- c(xcol, ycol)
      proj4string(points_data) <- CRS(points_crs)
      points_sp <- points_data

    } else if (file_ext %in% c("shp", "geojson", "gpkg")) {
      print(paste("Reading spatial file:", basename(points_sp)))
      points_sp <- rgdal::readOGR(points_sp)
    } else {
      stop(paste("Error: Unsupported file format:", file_ext,
                 "\nSupported formats: .csv, .shp, .geojson, .gpkg"))
    }

  } else if (is.data.frame(points_sp)) {

    # Require xcol, ycol, and points_crs for data frames
    if (is.null(xcol)) {
      stop("Error: 'xcol' is required when providing a data frame.")
    }
    if (is.null(ycol)) {
      stop("Error: 'ycol' is required when providing a data frame.")
    }
    if (is.null(points_crs)) {
      stop("Error: 'points_crs' is required when providing a data frame.")
    }

    print("Converting data frame to SpatialPointsDataFrame...")

    if (!xcol %in% names(points_sp)) {
      stop(paste0("Error: Column '", xcol, "' not found in data frame. Available columns: ",
                  paste(names(points_sp), collapse = ", ")))
    }
    if (!ycol %in% names(points_sp)) {
      stop(paste0("Error: Column '", ycol, "' not found in data frame. Available columns: ",
                  paste(names(points_sp), collapse = ", ")))
    }

    coordinates(points_sp) <- c(xcol, ycol)
    proj4string(points_sp) <- CRS(points_crs)

  } else if (!inherits(points_sp, "SpatialPointsDataFrame")) {
    stop("Error: points_sp must be a SpatialPointsDataFrame, data frame with x/y columns, or file path")
  }

  # Validate variable_patterns format
  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop(paste(
      "Error: variable_patterns must be a named vector.",
      "\n",
      "\n### Define variable patterns as follows in a named vector:",
      "\n#",
      "\n# my_variable_patterns <- c(",
      "\n#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "\n#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "\n#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "\n#   \"elevation\" = \"elevation\"",
      "\n# )",
      "\n#",
      "\n# Where:",
      "\n#   - The NAME (left side) is the variable name for the output column",
      "\n#   - The VALUE (right side) is the pattern to match in raster filenames",
      "\n#",
      "\n# For time-varying variables:",
      "\n#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "\n#   - These correspond to column names in your time_cols parameter",
      "\n#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "\n#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "\n#",
      "\n# For static variables:",
      "\n#   - Use a simple pattern with no time placeholders",
      "\n#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  if (any(names(variable_patterns) == "")) {
    stop("Error: All elements in variable_patterns must be named")
  }

  # Validate time_cols exist in data
  missing_cols <- setdiff(time_cols, names(points_sp@data))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following time_cols are missing from the input data:",
               paste(missing_cols, collapse = ", ")))
  }

  # Check for missing values in time columns
  n_original <- nrow(points_sp@data)
  for (tc in time_cols) {
    n_missing <- sum(is.na(points_sp@data[[tc]]))
    if (n_missing > 0) {
      pct_missing <- round(n_missing / n_original * 100, 2)
      warning(paste0("Warning: ", n_missing, " rows (", pct_missing,
                     "%) have missing values in time column '", tc, "'"))
    }
  }

  # Validate raster_dir
  if (!dir.exists(raster_dir)) {
    stop(paste("Error: Raster directory does not exist:", raster_dir))
  }

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print(paste("Processing", nrow(points_sp@data), "points"))

  ### Initialize columns for each variable
  for (var_name in names(variable_patterns)) {
    points_sp[[var_name]] <- NA
  }

  ### Get all files in raster directory
  all_files <- list.files(path = raster_dir, pattern = "tif",
                          recursive = TRUE, full.names = TRUE)
  print(paste("Found", length(all_files), "raster files"))

  if (length(all_files) == 0) {
    stop(paste("Error: No .tif files found in raster directory:", raster_dir))
  }

  ### Determine which variables are static vs dynamic
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

  print(paste("Dynamic variables:", ifelse(length(dynamic_vars) > 0,
                                           paste(dynamic_vars, collapse = ", "),
                                           "none")))
  print(paste("Static variables:", ifelse(length(static_vars) > 0,
                                          paste(static_vars, collapse = ", "),
                                          "none")))

  ### Extract static variables first
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
        print(paste("Warning: No file found for", var_name, "matching pattern:", search_pattern))
      }
    }

    ### Print confirmation of static rasters found
    if (length(static_rasters_found) > 0) {
      print(paste("Found", length(static_rasters_found), "static raster(s):"))
      for (rf in static_rasters_found) {
        print(paste("  -", rf))
      }
    }
  }

  ### Get unique time combinations and sort them
  time_combinations <- points_sp@data %>%
    dplyr::select(all_of(time_cols)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(across(all_of(time_cols)))

  print(paste("Extracting values for", nrow(time_combinations), "time periods..."))

  ### Extract values for each time combination
  dynamic_vars_extracted <- 0
  for (i in 1:nrow(time_combinations)) {
    time_values <- time_combinations[i, ]

    time_filter <- rep(TRUE, nrow(points_sp@data))
    for (tc in time_cols) {
      time_filter <- time_filter & (points_sp@data[[tc]] == time_values[[tc]])
    }
    points_subset <- points_sp[time_filter, ]

    ### Print time period information
    time_label <- paste(sapply(time_cols, function(tc) {
      paste0(tc, "=", time_values[[tc]])
    }), collapse = ", ")
    print(paste0("Processing time period ", i, "/", nrow(time_combinations), ": ", time_label))

    ### Collect rasters found for this time period
    rasters_found <- c()

    ### Extract dynamic variables
    for (var_name in dynamic_vars) {
      search_pattern <- variable_patterns[var_name]

      time_value_list <- list()
      for (tc in time_cols) {
        time_value_list[[tc]] <- as.character(time_values[[tc]])
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
        print(paste("Warning: No file found for", var_name, "matching pattern:", search_pattern))
      }
    }

    ### Print confirmation of rasters found
    if (length(rasters_found) > 0) {
      print(paste("Found", length(rasters_found), "raster(s):"))
      for (rf in rasters_found) {
        print(paste("  -", rf))
      }
    }
  }

  ### Check if no variables were extracted
  total_vars_extracted <- static_vars_extracted + dynamic_vars_extracted
  if (total_vars_extracted == 0) {
    stop(paste(
      "Error: No raster files could be matched for any of the specified variables.",
      "\n",
      "\n### Define variable patterns as follows in a named vector:",
      "\n#",
      "\n# my_variable_patterns <- c(",
      "\n#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "\n#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "\n#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "\n#   \"elevation\" = \"elevation\"",
      "\n# )",
      "\n#",
      "\n# Where:",
      "\n#   - The NAME (left side) is the variable name for the output column",
      "\n#   - The VALUE (right side) is the pattern to match in raster filenames",
      "\n#",
      "\n# For time-varying variables:",
      "\n#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "\n#   - These correspond to column names in your time_cols parameter",
      "\n#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "\n#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "\n#",
      "\n# For static variables:",
      "\n#   - Use a simple pattern with no time placeholders",
      "\n#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  ### Check if only static variables were extracted
  if (static_vars_extracted > 0 && dynamic_vars_extracted == 0 && length(dynamic_vars) > 0) {
    stop(paste(
      "Error: Only static variables could be matched. No raster files found for dynamic variables.",
      "\n",
      "\n### Define variable patterns as follows in a named vector:",
      "\n#",
      "\n# my_variable_patterns <- c(",
      "\n#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "\n#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "\n#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "\n#   \"elevation\" = \"elevation\"",
      "\n# )",
      "\n#",
      "\n# Where:",
      "\n#   - The NAME (left side) is the variable name for the output column",
      "\n#   - The VALUE (right side) is the pattern to match in raster filenames",
      "\n#",
      "\n# For time-varying variables:",
      "\n#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "\n#   - These correspond to column names in your time_cols parameter",
      "\n#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "\n#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "\n#",
      "\n# For static variables:",
      "\n#   - Use a simple pattern with no time placeholders",
      "\n#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  print("Value extraction complete")

  ### Save raw values if requested
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
      print(paste("Warning: No valid values found for", var_name))
    }
  }

  ### Save scaling parameters if requested
  if (save_scaling_params) {
    params_file <- file.path(output_dir, paste0(output_prefix, "_Scaling_Parameters.csv"))
    write.csv(scaling_params, params_file, row.names = FALSE)
    print(paste("Scaling parameters saved to:", basename(params_file)))
  }

  ### Apply scaling and save if requested
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
        print(paste("Warning: No scaling parameters found for", var_name))
      }
    }

    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    write.csv(scaled_data, scaled_output_file, row.names = FALSE)
    print(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  print("Processing complete!")
  return(invisible(NULL))
}
