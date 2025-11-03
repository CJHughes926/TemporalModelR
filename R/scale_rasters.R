#' Scale rasters using precomputed mean/sd
#' @export
#' @param input_dir Input directory with .tif files.
#' @param output_dir Output directory for scaled rasters.
#' @param scaling_params_file CSV with columns `variable,mean,sd`.
#' @param variable_patterns Named vector mapping variable -> filename pattern.
#' @param time_cols Optional character vector of time placeholders.
#' @param output_suffix Suffix for output files.
#' @param overwrite Overwrite existing outputs.
#' @importFrom raster raster writeRaster
#' @importFrom readr read_csv
scale_rasters <- function(input_dir,
                          output_dir,
                          scaling_params_file,
                          variable_patterns,
                          time_cols = NULL,
                          output_suffix = "_Scaled",
                          overwrite = TRUE) {

  require(raster)
  require(readr)

  ### Validate variable_patterns format
  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop(paste(
      "Error: variable_patterns must be a named vector.",
      "",
      "### Define variable patterns as follows in a named vector:",
      "#",
      "# my_variable_patterns <- c(",
      "#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "#   \"elevation\" = \"elevation\"",
      "# )",
      "#",
      "# Where:",
      "#   - The NAME (left side) is the variable name for the output column",
      "#   - The VALUE (right side) is the pattern to match in raster filenames",
      "#",
      "# For time-varying variables:",
      "#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "#   - These correspond to column names in your time_cols parameter",
      "#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "#",
      "# For static variables:",
      "#   - Use a simple pattern with no time placeholders",
      "#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  if (any(names(variable_patterns) == "")) {
    stop("Error: All elements in variable_patterns must be named")
  }

  ### Validate input directory
  if (!dir.exists(input_dir)) {
    stop(paste("Error: Input directory does not exist:", input_dir))
  }

  ### Validate scaling parameters file
  if (!file.exists(scaling_params_file)) {
    stop(paste("Error: Scaling parameters file does not exist:", scaling_params_file))
  }

  ### Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ### Load scaling parameters
  scaling_params <- read_csv(scaling_params_file, show_col_types = FALSE)
  print(paste("Loaded scaling parameters for", nrow(scaling_params), "variables"))

  ### Check if all variables in variable_patterns have scaling parameters
  missing_params <- setdiff(names(variable_patterns), scaling_params$variable)
  if (length(missing_params) > 0) {
    warning(paste("Warning: The following variables lack scaling parameters:",
                  paste(missing_params, collapse = ", ")))
  }

  ### Get all files
  all_files <- list.files(path = input_dir, pattern = "tif",
                          recursive = TRUE, full.names = TRUE)
  print(paste("Found", length(all_files), "raster files"))

  if (length(all_files) == 0) {
    stop(paste("Error: No .tif files found in input directory:", input_dir))
  }

  ### Print overwrite mode
  if (overwrite) {
    print("Overwrite mode: ON - Will process and overwrite all files")
  } else {
    print("Overwrite mode: OFF - Will skip files that already exist")
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

  ### Validate time_cols match variable_patterns

  if (length(dynamic_vars) > 0) {
    if (is.null(time_cols) || length(time_cols) == 0) {
      stop(paste(
        "Error: time_cols is required when variable_patterns contain time placeholders.",
        "",
        "Your variable_patterns include dynamic variables:", paste(dynamic_vars, collapse = ", "),
        "",
        "### You must specify time_cols:",
        "#",
        "# For single time dimension:",
        "# time_cols = \"Year\"",
        "#",
        "# For multiple time dimensions:",
        "# time_cols = c(\"Year\", \"Month\")",
        "# time_cols = c(\"Year\", \"DOY\")",
        "#",
        "# The time_cols must match the placeholders used in your variable_patterns.",
        sep = ""
      ))
    }

    ### Extract time placeholders from patterns
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
        "Warning: time_cols includes columns not found as placeholders in variable_patterns:",
        paste(missing_in_patterns, collapse = ", "),
        "This may indicate a mismatch between your time_cols and variable_patterns.",
        "",
        "Your time_cols:", paste(time_cols, collapse = ", "),
        "Placeholders found in variable_patterns:", paste(pattern_time_placeholders, collapse = ", ")
      ))
    }

    if (length(extra_in_patterns) > 0) {
      stop(paste(
        "Error: variable_patterns contain time placeholders not specified in time_cols:",
        paste(extra_in_patterns, collapse = ", "),
        "",
        "Your time_cols:", ifelse(is.null(time_cols), "NULL", paste(time_cols, collapse = ", ")),
        "Placeholders found in variable_patterns:", paste(pattern_time_placeholders, collapse = ", "),
        "",
        "### To fix this issue:",
        "#",
        "# Option 1: Add missing time columns to time_cols:",
        "# time_cols = c(\"", paste(c(time_cols, tolower(extra_in_patterns)), collapse = "\", \""), "\")",
        "#",
        "# Option 2: Update your variable_patterns to match time_cols:",
        "# Ensure all placeholders in patterns (YEAR, MONTH, etc.) are in time_cols",
        "#",
        "# Example of matching time_cols and patterns:",
        "# time_cols = c(\"Year\", \"Month\")",
        "# variable_patterns = c(",
        "#   \"Temperature\" = \"Temperature_YEAR_MONTH\",",
        "#   \"Precipitation\" = \"Precipitation_YEAR_MONTH\"",
        "# )",
        sep = ""
      ))
    }
  }

  if (length(dynamic_vars) > 0 && (is.null(time_cols) || length(time_cols) == 0)) {
    stop(paste(
      "Error: Dynamic variables detected in variable_patterns, but time_cols is NULL or empty.",
      "",
      "Dynamic variables found:", paste(dynamic_vars, collapse = ", "),
      "",
      "### To fix this issue:",
      "#",
      "# Provide time_cols parameter matching the placeholders in your patterns:",
      "# time_cols = c(\"Year\")  # for patterns like \"Variable_YEAR\"",
      "# time_cols = c(\"Year\", \"Month\")  # for patterns like \"Variable_YEAR_MONTH\"",
      sep = ""
    ))
  }

  ### Track successful scaling
  static_vars_scaled <- 0
  dynamic_vars_scaled <- 0

  ### Process static variables
  if (length(static_vars) > 0) {
    print("Scaling static variables...")

    for (var_name in static_vars) {
      output_file <- file.path(output_dir, paste0(var_name, output_suffix, ".tif"))

      if (file.exists(output_file) && !overwrite) {
        print(paste("Skipping", var_name, "- already scaled"))
        next
      }

      print(paste("Scaling static variable:", var_name))

      ### Get scaling parameters
      var_params <- scaling_params[scaling_params$variable == var_name, ]

      if (nrow(var_params) == 0) {
        print(paste("Warning: No scaling parameters found for", var_name))
        next
      }

      var_mean <- var_params$mean
      var_sd <- var_params$sd

      ### Find matching file
      search_pattern <- variable_patterns[var_name]
      pattern_parts <- strsplit(search_pattern, "_")[[1]]
      first_part <- var_name

      matching_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)

        starts_correctly <- grepl(paste0("^", first_part), filename, ignore.case = TRUE)

        all_parts_present <- all(sapply(pattern_parts, function(part) {
          grepl(part, filename, ignore.case = TRUE)
        }))

        starts_correctly && all_parts_present
      })]

      if (length(matching_files) > 0) {
        print(paste("  Using file:", basename(matching_files[1])))

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

  ### Process dynamic variables
  if (length(dynamic_vars) > 0) {
    print("Scaling dynamic variables...")

    for (var_name in dynamic_vars) {
      print(paste("Processing dynamic variable:", var_name))

      ### Get scaling parameters
      var_params <- scaling_params[scaling_params$variable == var_name, ]

      if (nrow(var_params) == 0) {
        print(paste("Warning: No scaling parameters found for", var_name))
        next
      }

      var_mean <- var_params$mean
      var_sd <- var_params$sd

      ### Get the pattern for this variable
      search_pattern <- variable_patterns[var_name]

      ### Find all matching files for this variable
      candidate_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)
        grepl(paste0("^", var_name), filename, ignore.case = TRUE)
      })]

      if (length(candidate_files) == 0) {
        print(paste("Warning: No files found for", var_name))
        next
      }

      print(paste("  Found", length(candidate_files), "candidate files"))

      ### Extract time values from filenames using pattern structure
      time_file_pairs <- list()

      for (file in candidate_files) {
        filename <- basename(file)

        ### Try to extract time values for each time column
        time_vals <- list()
        all_found <- TRUE

        ### For each time column, look for numeric values in appropriate positions
        ### We'll extract all numbers and match them by position in the pattern

        ### Split filename into parts
        filename_no_ext <- gsub("\\.(tif|TIF)$", "", filename)
        filename_parts <- strsplit(filename_no_ext, "_")[[1]]

        ### Get pattern parts
        pattern_parts <- strsplit(search_pattern, "_")[[1]]

        ### Try to match pattern structure to filename structure
        for (tc in time_cols) {
          ### Find position of time column placeholder in pattern
          tc_positions <- which(toupper(pattern_parts) == toupper(tc))

          if (length(tc_positions) > 0) {
            ### Use the first occurrence
            tc_pos <- tc_positions[1]

            ### Try to get the corresponding part from filename
            if (tc_pos <= length(filename_parts)) {
              ### Extract numeric value from this position
              value <- filename_parts[tc_pos]

              ### Clean the value (remove non-numeric characters if it's a number)
              if (grepl("^\\d+$", value)) {
                time_vals[[tc]] <- value
              } else {
                ### Try to extract numbers from the part
                num_match <- regmatches(value, gregexpr("\\d+", value))[[1]]
                if (length(num_match) > 0) {
                  time_vals[[tc]] <- num_match[1]
                } else {
                  all_found <- FALSE
                  break
                }
              }
            } else {
              all_found <- FALSE
              break
            }
          } else {
            all_found <- FALSE
            break
          }
        }

        if (all_found && length(time_vals) == length(time_cols)) {
          ### Verify this is a valid match by reconstructing the pattern
          test_pattern <- search_pattern
          for (tc in time_cols) {
            test_pattern <- gsub(tc, time_vals[[tc]], test_pattern, ignore.case = TRUE)
          }

          ### Check if reconstructed pattern matches major parts of filename
          test_parts <- strsplit(test_pattern, "_")[[1]]
          match_count <- sum(sapply(test_parts, function(part) {
            grepl(part, filename, ignore.case = TRUE)
          }))

          ### If most parts match, consider it valid
          if (match_count >= length(test_parts) * 0.8) {
            time_file_pairs[[file]] <- time_vals
          }
        }
      }

      if (length(time_file_pairs) == 0) {
        print(paste("  Warning: No valid time information found in filenames for", var_name))
        next
      }

      print(paste("  Processing", length(time_file_pairs), "time periods"))

      ### Process each file
      processed_count <- 0
      skipped_count <- 0

      for (file in names(time_file_pairs)) {
        time_vals <- time_file_pairs[[file]]

        ### Build output filename
        time_string <- paste(sapply(time_cols, function(tc) time_vals[[tc]]), collapse = "_")
        output_file <- file.path(output_dir, paste0(var_name, "_", time_string, output_suffix, ".tif"))

        if (file.exists(output_file) && !overwrite) {
          skipped_count <- skipped_count + 1
          next
        }

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

  ### Check if no variables were scaled
  total_vars_scaled <- static_vars_scaled + dynamic_vars_scaled
  if (total_vars_scaled == 0) {
    stop(paste(
      "Error: No raster files could be matched and scaled for any of the specified variables.",
      "",
      "### Define variable patterns as follows in a named vector:",
      "#",
      "# my_variable_patterns <- c(",
      "#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "#   \"elevation\" = \"elevation\"",
      "# )",
      "#",
      "# Where:",
      "#   - The NAME (left side) is the variable name for the output column",
      "#   - The VALUE (right side) is the pattern to match in raster filenames",
      "#",
      "# For time-varying variables:",
      "#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "#   - These correspond to column names in your time_cols parameter",
      "#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "#",
      "# For static variables:",
      "#   - Use a simple pattern with no time placeholders",
      "#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  ### Check if only static variables were scaled
  if (static_vars_scaled > 0 && dynamic_vars_scaled == 0 && length(dynamic_vars) > 0) {
    stop(paste(
      "Error: Only static variables could be matched and scaled. No raster files found for dynamic variables.",
      "",
      "### Define variable patterns as follows in a named vector:",
      "#",
      "# my_variable_patterns <- c(",
      "#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "#   \"Open_Percentage2\" = \"Open_Percentage2_YEAR\",",
      "#   \"Forest_Percentage2\" = \"Forest_Percentage2_YEAR\",",
      "#   \"elevation\" = \"elevation\"",
      "# )",
      "#",
      "# Where:",
      "#   - The NAME (left side) is the variable name for the output column",
      "#   - The VALUE (right side) is the pattern to match in raster filenames",
      "#",
      "# For time-varying variables:",
      "#   - Use placeholders like YEAR, MONTH, DAY in the pattern",
      "#   - These correspond to column names in your time_cols parameter",
      "#   - Example: \"Forest_YEAR\" with time_cols = \"YEAR\"",
      "#   - Example: \"Temp_MONTH_YEAR\" with time_cols = c(\"MONTH\", \"YEAR\")",
      "#",
      "# For static variables:",
      "#   - Use a simple pattern with no time placeholders",
      "#   - Example: \"elevation\" = \"elevation\"",
      sep = ""
    ))
  }

  print("All scaling complete!")
  print(paste("Successfully scaled", total_vars_scaled, "variables"))
  print(paste("  Static:", static_vars_scaled, "| Dynamic:", dynamic_vars_scaled))

  return(invisible(NULL))
}
