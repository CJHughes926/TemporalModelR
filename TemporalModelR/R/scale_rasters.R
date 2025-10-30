scale_rasters <- function(input_dir,
                          output_dir,
                          scaling_params_file,
                          variable_patterns,
                          time_cols = NULL,
                          output_suffix = "_Scaled",
                          overwrite = TRUE) {

  # Load required libraries
  require(raster)
  require(readr)

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

  # Validate input directory
  if (!dir.exists(input_dir)) {
    stop(paste("Error: Input directory does not exist:", input_dir))
  }

  # Validate scaling parameters file
  if (!file.exists(scaling_params_file)) {
    stop(paste("Error: Scaling parameters file does not exist:", scaling_params_file))
  }

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Load scaling parameters
  scaling_params <- read_csv(scaling_params_file, show_col_types = FALSE)
  print(paste("Loaded scaling parameters for", nrow(scaling_params), "variables"))

  # Check if all variables in variable_patterns have scaling parameters
  missing_params <- setdiff(names(variable_patterns), scaling_params$variable)
  if (length(missing_params) > 0) {
    warning(paste("Warning: The following variables lack scaling parameters:",
                  paste(missing_params, collapse = ", ")))
  }

  # Get all files
  all_files <- list.files(path = input_dir, pattern = "tif",
                          recursive = TRUE, full.names = TRUE)
  print(paste("Found", length(all_files), "raster files"))

  if (length(all_files) == 0) {
    stop(paste("Error: No .tif files found in input directory:", input_dir))
  }

  # Print overwrite mode
  if (overwrite) {
    print("Overwrite mode: ON - Will process and overwrite all files")
  } else {
    print("Overwrite mode: OFF - Will skip files that already exist")
  }

  # Determine which variables are static vs dynamic
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

  # Track successful scaling
  static_vars_scaled <- 0
  dynamic_vars_scaled <- 0

  # Process static variables
  if (length(static_vars) > 0) {
    print("Scaling static variables...")

    for (var_name in static_vars) {
      output_file <- file.path(output_dir, paste0(var_name, output_suffix, ".tif"))

      if (file.exists(output_file) && !overwrite) {
        print(paste("Skipping", var_name, "- already scaled"))
        next
      }

      print(paste("Scaling static variable:", var_name))

      # Get scaling parameters
      var_params <- scaling_params[scaling_params$variable == var_name, ]

      if (nrow(var_params) == 0) {
        print(paste("Warning: No scaling parameters found for", var_name))
        next
      }

      var_mean <- var_params$mean
      var_sd <- var_params$sd

      # Find matching file
      search_pattern <- variable_patterns[var_name]
      pattern_parts <- strsplit(search_pattern, "_")[[1]]

      matching_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)

        # Check if filename starts with variable name
        starts_correctly <- grepl(paste0("^", var_name), filename, ignore.case = TRUE)

        # Check if all parts exist
        all_parts_present <- all(sapply(pattern_parts, function(part) {
          grepl(part, filename, ignore.case = TRUE)
        }))

        starts_correctly && all_parts_present
      })]

      if (length(matching_files) > 0) {
        print(paste("  Using file:", basename(matching_files[1])))

        # Load, scale, and save
        raster_layer <- raster(matching_files[1])
        raster_scaled <- (raster_layer - var_mean) / var_sd
        writeRaster(raster_scaled, output_file, format = "GTiff", overwrite = TRUE)
        print(paste("  Saved to:", basename(output_file)))
        static_vars_scaled <- static_vars_scaled + 1
      } else {
        print(paste("Warning: No file found for", var_name, "matching pattern:", search_pattern))
      }
    }
  }

  # Process dynamic variables
  if (length(dynamic_vars) > 0) {
    print("Scaling dynamic variables...")

    for (var_name in dynamic_vars) {
      print(paste("Processing dynamic variable:", var_name))

      # Get scaling parameters
      var_params <- scaling_params[scaling_params$variable == var_name, ]

      if (nrow(var_params) == 0) {
        print(paste("Warning: No scaling parameters found for", var_name))
        next
      }

      var_mean <- var_params$mean
      var_sd <- var_params$sd

      # Find all matching files for this variable
      matching_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)
        grepl(paste0("^", var_name), filename, ignore.case = TRUE)
      })]

      if (length(matching_files) == 0) {
        print(paste("Warning: No files found for", var_name))
        next
      }

      print(paste("  Found", length(matching_files), "files"))

      # Extract time values from filenames
      time_file_pairs <- list()

      for (file in matching_files) {
        filename <- basename(file)

        # Extract time values for each time column
        time_vals <- list()
        all_found <- TRUE

        for (tc in time_cols) {
          tc_lower <- tolower(tc)

          # Extract based on column type
          if (tc_lower == "year") {
            # Extract 4-digit year
            year_match <- regmatches(filename, gregexpr("\\d{4}", filename))[[1]]
            if (length(year_match) > 0) {
              time_vals[[tc]] <- year_match[1]
            } else {
              all_found <- FALSE
              break
            }
          } else if (tc_lower == "month") {
            # Extract 2-digit month
            month_match <- regmatches(filename, gregexpr("_(\\d{2})_", filename))[[1]]
            if (length(month_match) > 0) {
              time_vals[[tc]] <- gsub("_", "", month_match[1])
            } else {
              all_found <- FALSE
              break
            }
          } else if (tc_lower == "doy") {
            # Extract 3-digit day of year
            doy_match <- regmatches(filename, gregexpr("\\d{3}", filename))[[1]]
            if (length(doy_match) > 0) {
              time_vals[[tc]] <- doy_match[1]
            } else {
              all_found <- FALSE
              break
            }
          } else {
            # Generic extraction for custom columns (e.g., Season)
            # Look for pattern: columnname_value_ or columnname_value.
            # First try to find the column name in the filename
            pattern_with_col <- paste0(tc, "_([^_\\.]+)")
            col_match <- regmatches(filename, regexpr(pattern_with_col, filename, ignore.case = TRUE))

            if (length(col_match) > 0) {
              # Extract the value after the column name
              value_match <- regmatches(col_match, regexpr("_(.+)$", col_match))
              if (length(value_match) > 0) {
                time_vals[[tc]] <- gsub("^_", "", value_match[1])
              } else {
                all_found <- FALSE
                break
              }
            } else {
              # Fallback: try to extract any alphanumeric value
              # Look for values between underscores that might match
              parts <- strsplit(filename, "_")[[1]]
              found_value <- FALSE

              for (part in parts) {
                # Check if this part contains alphanumeric values
                clean_part <- gsub("\\.(tif|TIF)$", "", part)
                if (nchar(clean_part) > 0 && clean_part != var_name) {
                  time_vals[[tc]] <- clean_part
                  found_value <- TRUE
                  break
                }
              }

              if (!found_value) {
                all_found <- FALSE
                break
              }
            }
          }
        }

        if (all_found) {
          time_file_pairs[[file]] <- time_vals
        }
      }

      if (length(time_file_pairs) == 0) {
        print(paste("  Warning: No valid time information found in filenames for", var_name))
        next
      }

      print(paste("  Processing", length(time_file_pairs), "time periods"))

      # Process each file
      processed_count <- 0
      skipped_count <- 0

      for (file in names(time_file_pairs)) {
        time_vals <- time_file_pairs[[file]]

        # Build output filename
        time_string <- paste(sapply(time_cols, function(tc) time_vals[[tc]]), collapse = "_")
        output_file <- file.path(output_dir, paste0(var_name, "_", time_string, output_suffix, ".tif"))

        if (file.exists(output_file) && !overwrite) {
          skipped_count <- skipped_count + 1
          next
        }

        # Load, scale, and save
        raster_layer <- raster(file)
        raster_scaled <- (raster_layer - var_mean) / var_sd
        writeRaster(raster_scaled, output_file, format = "GTiff", overwrite = TRUE)
        processed_count <- processed_count + 1
      }

      if (processed_count > 0) {
        dynamic_vars_scaled <- dynamic_vars_scaled + 1
      }

      print(paste("  Processed:", processed_count, "| Skipped:", skipped_count))
    }
  }

  # Check if no variables were scaled
  total_vars_scaled <- static_vars_scaled + dynamic_vars_scaled
  if (total_vars_scaled == 0) {
    stop(paste(
      "Error: No raster files could be matched and scaled for any of the specified variables.",
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

  # Check if only static variables were scaled
  if (static_vars_scaled > 0 && dynamic_vars_scaled == 0 && length(dynamic_vars) > 0) {
    stop(paste(
      "Error: Only static variables could be matched and scaled. No raster files found for dynamic variables.",
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

  print("All scaling complete!")
  print(paste("Successfully scaled", total_vars_scaled, "variables"))
  print(paste("  Static:", static_vars_scaled, "| Dynamic:", dynamic_vars_scaled))

  return(invisible(NULL))
}
