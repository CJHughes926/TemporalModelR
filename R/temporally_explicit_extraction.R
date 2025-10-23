#' Temporally explicit extraction
#' @export
#' @importFrom raster raster extract
#' @importFrom dplyr select distinct all_of
#' @importFrom utils write.csv
temporally_explicit_extraction <- function(points_sp,
                                           raster_dir,
                                           variable_patterns, ### # Define variable patterns as follows in a c() list
                                           #my_variable_patterns <- c(
                                          #   "Developed_Percentage2" = "Developed_Percentage2_YEAR",
                                          #   "Open_Percentage2" = "Open_Percentage2_YEAR",
                                          #   "Forest_Percentage2" = "Forest_Percentage2_YEAR",
                                          #   "elevation" = "elevation"
                                           #) Where the first component is the Variable name without time, and the second component is the variable name with the pattern stand
                                          #in which corrosponds to a given time collumn, for example "YEAR", "MONTH", "DAY"... Should be able to accept multiple columns/patterns
                                          # for example Variable_MONTH_YEAR, in which case "time cols" should be a concatinated list c("Month", "Year")
                                          #Can also handle varuavkes with no time component- like elevation.
                                           time_cols,
                                           output_dir,
                                           output_prefix = "temp_explicit_df",
                                           save_raw = TRUE,
                                           save_scaled = TRUE,
                                           save_scaling_params = TRUE) {

  # Load required libraries
  require(raster)
  require(sp)
  require(readr)
  require(dplyr)

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  print(paste("Processing", nrow(points_sp@data), "points"))

  # Initialize columns for each variable
  for (var_name in names(variable_patterns)) {
    points_sp[[var_name]] <- NA
  }

  # Get all files in raster directory
  all_files <- list.files(path = raster_dir, pattern = "tif",
                          recursive = TRUE, full.names = TRUE)
  print(paste("Found", length(all_files), "raster files"))

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

  # Get unique time combinations
  time_combinations <- points_sp@data %>%
    dplyr::select(all_of(time_cols)) %>%
    dplyr::distinct()

  print(paste("Extracting values for", nrow(time_combinations), "time periods..."))

  # Extract values for each time combination
  for (i in 1:nrow(time_combinations)) {
    time_values <- time_combinations[i, ]

    # Create filter for this time combination
    time_filter <- rep(TRUE, nrow(points_sp@data))
    for (tc in time_cols) {
      time_filter <- time_filter & (points_sp@data[[tc]] == time_values[[tc]])
    }
    points_subset <- points_sp[time_filter, ]

    # Extract dynamic variables
    for (var_name in dynamic_vars) {
      # Build search pattern by replacing time placeholders
      search_pattern <- variable_patterns[var_name]

      # Store the actual time values that were substituted
      time_value_list <- list()
      for (tc in time_cols) {
        time_value_list[[tc]] <- as.character(time_values[[tc]])
        search_pattern <- gsub(tc, time_value_list[[tc]],
                               search_pattern, ignore.case = TRUE)
      }

      # Split pattern into parts
      pattern_parts <- strsplit(search_pattern, "_")[[1]]

      # Use the variable name itself as the first part to check
      first_part <- var_name

      matching_files <- all_files[sapply(all_files, function(f) {
        filename <- basename(f)

        # Check if filename starts with the variable name (case-insensitive)
        starts_correctly <- grepl(paste0("^", first_part), filename, ignore.case = TRUE)

        # Check if all parts exist in the filename
        all_parts_present <- all(sapply(pattern_parts, function(part) {
          grepl(part, filename, ignore.case = TRUE)
        }))

        # Check that each time value appears as a complete segment (bounded by _ or .)
        time_values_bounded <- all(sapply(names(time_value_list), function(tc) {
          time_val <- time_value_list[[tc]]
          # Match time value only when bounded by underscore or dot on both sides
          grepl(paste0("(^|_)", time_val, "(_|\\.)"), filename, ignore.case = TRUE)
        }))

        # All three conditions must be true
        starts_correctly && all_parts_present && time_values_bounded
      })]

      if (length(matching_files) > 0) {
        raster_file <- matching_files[1]
        raster_layer <- raster(raster_file)
        values <- raster::extract(raster_layer, points_subset)
        points_sp@data[time_filter, var_name] <- values
      } else {
        print(paste("  Warning: No file found for", var_name, "matching pattern:", search_pattern))
      }
    }

    # Extract static variables (only once, on first iteration)
    if (i == 1) {
      for (var_name in static_vars) {
        search_pattern <- variable_patterns[var_name]

        # Use the variable name itself as the first part to check
        first_part <- var_name
        pattern_parts <- strsplit(search_pattern, "_")[[1]]

        matching_files <- all_files[sapply(all_files, function(f) {
          filename <- basename(f)

          # Check if filename starts with the variable name
          starts_correctly <- grepl(paste0("^", first_part), filename, ignore.case = TRUE)

          # Check if all parts exist in the filename
          all_parts_present <- all(sapply(pattern_parts, function(part) {
            grepl(part, filename, ignore.case = TRUE)
          }))

          starts_correctly && all_parts_present
        })]

        if (length(matching_files) > 0) {
          raster_file <- matching_files[1]
          raster_layer <- raster(raster_file)
          values <- raster::extract(raster_layer, points_sp)
          points_sp@data[[var_name]] <- values
        } else {
          print(paste("  Warning: No file found for", var_name, "matching pattern:", search_pattern))
        }
      }
    }
  }

  print("Value extraction complete")

  # Save raw values if requested
  if (save_raw) {
    raw_output_file <- file.path(output_dir, paste0(output_prefix, "_Raw_Values.csv"))
    write.csv(points_sp@data, raw_output_file, row.names = FALSE)
    print(paste("Raw values saved to:", basename(raw_output_file)))
  }

  # Calculate scaling parameters
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
      print(paste("  Warning: No valid values found for", var_name))
    }
  }

  # Save scaling parameters if requested
  if (save_scaling_params) {
    params_file <- file.path(output_dir, paste0(output_prefix, "_Scaling_Parameters.csv"))
    write.csv(scaling_params, params_file, row.names = FALSE)
    print(paste("Scaling parameters saved to:", basename(params_file)))
  }

  # Apply scaling and save if requested
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
        print(paste("  Warning: No scaling parameters found for", var_name))
      }
    }

    scaled_output_file <- file.path(output_dir, paste0(output_prefix, "_Scaled_Values.csv"))
    write.csv(scaled_data, scaled_output_file, row.names = FALSE)
    print(paste("Scaled values saved to:", basename(scaled_output_file)))
  }

  print("Processing complete!")
  return(invisible(NULL))
}
