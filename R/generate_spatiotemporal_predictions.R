test_cp_likelihood <- function(data, cp, alpha = 0.05) {
  n <- nrow(data)
  seg1_data <- data[1:cp, ]
  seg2_data <- data[(cp + 1):n, ]

  if (nrow(seg1_data) < 3 || nrow(seg2_data) < 3) {
    return(list(p_value = NA, significant = FALSE, test_statistic = NA))
  }

  tryCatch({
    model_full <- glm(y ~ lag1 + neighbor, data = data, family = binomial())
    loglik_full <- logLik(model_full)[1]

    model_seg1 <- glm(y ~ lag1 + neighbor, data = seg1_data, family = binomial())
    model_seg2 <- glm(y ~ lag1 + neighbor, data = seg2_data, family = binomial())
    loglik_seg <- logLik(model_seg1)[1] + logLik(model_seg2)[1]

    lr_stat <- 2 * (loglik_seg - loglik_full)
    df <- length(coef(model_seg1)) + length(coef(model_seg2)) - length(coef(model_full))
    p_value <- 1 - pchisq(lr_stat, df = max(df, 1))

    list(p_value = p_value, significant = p_value < alpha, test_statistic = lr_stat, df = df)
  }, error = function(e) {
    list(p_value = NA, significant = FALSE, test_statistic = NA)
  })
}

test_cp_permutation <- function(data, cp, n_perm = 1000, alpha = 0.05) {
  n <- nrow(data)
  seg1_mean <- mean(data$y[1:cp])
  seg2_mean <- mean(data$y[(cp + 1):n])
  obs_stat <- abs(seg1_mean - seg2_mean)

  perm_stats <- replicate(n_perm, {
    perm_y <- sample(data$y)
    abs(mean(perm_y[1:cp]) - mean(perm_y[(cp + 1):n]))
  })

  p_value <- mean(perm_stats >= obs_stat)
  list(p_value = p_value, significant = p_value < alpha,
       test_statistic = obs_stat, effect_size = obs_stat)
}

test_cp_proportion <- function(data, cp, alpha = 0.05) {
  n <- nrow(data)
  seg1_y <- data$y[1:cp]
  seg2_y <- data$y[(cp + 1):n]

  n1 <- length(seg1_y)
  n2 <- length(seg2_y)
  p1 <- mean(seg1_y)
  p2 <- mean(seg2_y)

  p_pooled <- (sum(seg1_y) + sum(seg2_y)) / (n1 + n2)

  if (p_pooled == 0 || p_pooled == 1 || p1 == p2) {
    return(list(p_value = 1, significant = FALSE, test_statistic = 0, effect_size = abs(p1 - p2)))
  }

  se <- sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))
  if (se == 0) {
    return(list(p_value = 1, significant = FALSE, test_statistic = 0, effect_size = abs(p1 - p2)))
  }

  z_stat <- (p1 - p2) / se
  p_value <- 2 * (1 - pnorm(abs(z_stat)))

  list(p_value = p_value, significant = p_value < alpha,
       test_statistic = z_stat, effect_size = abs(p1 - p2))
}

test_cp_chisquare <- function(data, cp, alpha = 0.05) {
  n <- nrow(data)
  seg1_y <- data$y[1:cp]
  seg2_y <- data$y[(cp + 1):n]

  cont_table <- matrix(c(
    sum(seg1_y == 1), sum(seg1_y == 0),
    sum(seg2_y == 1), sum(seg2_y == 0)
  ), nrow = 2, byrow = TRUE)

  test_result <- chisq.test(cont_table, correct = FALSE)

  list(p_value = test_result$p.value, significant = test_result$p.value < alpha,
       test_statistic = test_result$statistic, effect_size = abs(mean(seg1_y) - mean(seg2_y)))
}

assess_changepoint_significance <- function(data, cp_set, alpha = 0.05, n_perm = 1000) {
  if (length(cp_set) == 0) return(data.frame())

  cp_set <- sort(cp_set)
  significant_cps <- c()
  results_list <- list()

  for (i in seq_along(cp_set)) {
    cp <- cp_set[i]

    if (length(significant_cps) == 0) {
      baseline_start <- 1
    } else {
      baseline_start <- max(significant_cps) + 1
    }

    if (i < length(cp_set)) {
      current_end <- cp_set[i + 1]
    } else {
      current_end <- nrow(data)
    }

    baseline_segment <- baseline_start:cp
    current_segment <- (cp + 1):current_end

    if (length(baseline_segment) < 3 || length(current_segment) < 3) {
      results_list[[i]] <- data.frame(
        ChangePoint = cp,
        LR_PValue = NA, LR_Significant = FALSE,
        Perm_PValue = NA, Perm_Significant = FALSE, Perm_EffectSize = NA,
        Prop_PValue = NA, Prop_Significant = FALSE,
        Chi_PValue = NA, Chi_Significant = FALSE,
        Seg1_Proportion = NA, Seg2_Proportion = NA,
        Overall_Significant = FALSE
      )
      next
    }

    baseline_prop <- mean(data$y[baseline_segment])
    current_prop <- mean(data$y[current_segment])

    temp_combined <- data[c(baseline_segment, current_segment), ]
    cp_temp <- length(baseline_segment)

    lr_test <- suppressWarnings(test_cp_likelihood(temp_combined, cp_temp, alpha))
    lr_p_value <- ifelse(is.null(lr_test$p.value), lr_test$p_value, lr_test$p.value)
    lr_significant <- ifelse(is.na(lr_p_value), FALSE, lr_p_value < alpha)

    perm_test <- suppressWarnings(test_cp_permutation(temp_combined, cp_temp, n_perm, alpha))
    perm_p_value <- perm_test$p_value
    perm_significant <- perm_test$significant
    perm_effect <- perm_test$effect_size

    prop_test <- suppressWarnings(test_cp_proportion(temp_combined, cp_temp, alpha))
    prop_p_value <- prop_test$p_value
    prop_significant <- prop_test$significant

    chi_test <- suppressWarnings(tryCatch({
      test_cp_chisquare(temp_combined, cp_temp, alpha)
    }, error = function(e) {
      list(p_value = NA, significant = FALSE, test_statistic = NA)
    }))
    chi_p_value <- chi_test$p_value
    chi_significant <- ifelse(is.na(chi_p_value), FALSE, chi_p_value < alpha)

    sig_votes <- sum(c(lr_significant, perm_significant, prop_significant, chi_significant), na.rm = TRUE)
    overall_significant <- sig_votes >= 2

    if (overall_significant) {
      significant_cps <- c(significant_cps, cp)
    }

    results_list[[i]] <- data.frame(
      ChangePoint = cp,
      LR_PValue = lr_p_value, LR_Significant = lr_significant,
      Perm_PValue = perm_p_value, Perm_Significant = perm_significant,
      Perm_EffectSize = perm_effect,
      Prop_PValue = prop_p_value, Prop_Significant = prop_significant,
      Chi_PValue = chi_p_value, Chi_Significant = chi_significant,
      Seg1_Proportion = baseline_prop, Seg2_Proportion = current_prop,
      Overall_Significant = overall_significant
    )
  }

  do.call(rbind, results_list)
}

# FIXED classify_pattern function with better logic
classify_pattern <- function(sig_results) {
  # Filter to significant changepoints only
  sig <- sig_results[sig_results$Overall_Significant == TRUE, ]

  # If no significant changepoints, it's truly "No Pattern"
  if (nrow(sig) == 0) {
    return("No Pattern")
  }

  # Use a small tolerance for floating point comparisons
  tolerance <- 1e-10

  # Calculate differences with tolerance
  prop_diffs <- sig$Seg2_Proportion - sig$Seg1_Proportion

  # Check for increases and decreases
  has_inc <- any(prop_diffs > tolerance)
  has_dec <- any(prop_diffs < -tolerance)

  # Classify based on patterns
  if (has_inc & has_dec) {
    return("Fluctuating/Intermittent")
  } else if (has_inc) {
    return("Increasing")
  } else if (has_dec) {
    return("Decreasing")
  } else {
    # This catches the edge case where proportions are essentially equal
    # Could indicate data issues or very subtle changes
    return("Failed Classification")
  }
}

# ============================================================================
# FIXED classify_pixel_cpd - USES MATRIX FORMAT WITH WARNINGS SUPPRESSED
# ============================================================================
classify_pixel_cpd <- function(pixel_vals, n_middle, method = "MBIC", alpha = 0.05, n_perm = 1000) {
  y <- pixel_vals[1:n_middle]
  lag <- pixel_vals[(n_middle + 1):(2 * n_middle)]
  neighbor <- pixel_vals[(2 * n_middle + 1):(3 * n_middle)]
  mean_val <- pixel_vals[3 * n_middle + 1]

  # Handle always absent/present cases
  if (mean_val < 0.01) return(1)  # Always Absent
  if (mean_val > 0.99) return(2)  # Always Present
  if (any(is.na(c(y, lag, neighbor)))) return(NA)  # Missing data

  # FIXED: fastcpd.binomial expects cbind(y, x1, x2, ...) format
  # Create matrix with y in first column, predictors in remaining columns
  data_matrix <- cbind(y, lag, neighbor)

  # Try to run changepoint detection (suppress GLM warnings)
  cp_result <- tryCatch({
    suppressWarnings({
      fastcpd.binomial(
        data = data_matrix,
        r.progress = FALSE
      )
    })
  }, error = function(e) {
    return(NULL)
  })

  # If changepoint detection failed, mark as "Failed to Classify"
  if (is.null(cp_result)) return(7)  # Failed - CPD error

  # Extract changepoints
  cp_set <- cp_result@cp_set

  # If no changepoints detected, classify as "No Pattern"
  if (length(cp_set) == 0) return(3)  # No Pattern

  # Create data frame for significance testing
  data_df <- data.frame(y = y, lag1 = lag, neighbor = neighbor)

  # Assess significance of detected changepoints (suppress warnings)
  sig_results <- tryCatch({
    suppressWarnings({
      assess_changepoint_significance(data_df, cp_set, alpha, n_perm)
    })
  }, error = function(e) {
    return(NULL)
  })

  # If significance testing failed
  if (is.null(sig_results) || nrow(sig_results) == 0) return(7)  # Failed - sig test error

  # Classify the pattern based on significance results
  pattern <- classify_pattern(sig_results)

  # Return numeric code
  pattern_codes <- c(
    "Always Absent" = 1,
    "Always Present" = 2,
    "No Pattern" = 3,
    "Increasing" = 4,
    "Decreasing" = 5,
    "Fluctuating/Intermittent" = 6,
    "Failed Classification" = 7
  )

  return(pattern_codes[pattern])
}

library(raster)
library(fastcpd)
library(dplyr)

# ============================================================================
# INTEGRATED TILED ANALYSIS WITH YEAR-OF-CHANGE RASTERS
# ============================================================================
# This script includes all the functions from FIXED_classification_code.R
# but with an enhanced analyze_temporal_patterns_tiles_raster() that creates
# THREE rasters at once:
# 1. Pattern classification raster
# 2. Year of first decrease raster
# 3. Year of first increase raster
# ============================================================================

# Source the complete FIXED_classification_code.R first to get all the helper functions
# Then use the new function below

# ============================================================================
# ENHANCED PIXEL CLASSIFICATION WITH YEAR EXTRACTION
# ============================================================================

classify_pixel_with_years <- function(pixel_vals, n_middle, time_steps,
                                      method = method, alpha = 0.05, n_perm = 1000) {
  # Returns a vector of length 3: [classification, year_decrease, year_increase]

  y <- pixel_vals[1:n_middle]
  lag <- pixel_vals[(n_middle + 1):(2 * n_middle)]
  neighbor <- pixel_vals[(2 * n_middle + 1):(3 * n_middle)]
  mean_val <- pixel_vals[3 * n_middle + 1]

  # Handle always absent/present cases
  if (mean_val < 0.01) return(c(1, NA, NA))  # Always Absent
  if (mean_val > 0.99) return(c(2, NA, NA))  # Always Present
  if (any(is.na(c(y, lag, neighbor)))) return(c(NA, NA, NA))  # Missing data

  # Create matrix for changepoint detection
  data_matrix <- cbind(y, lag, neighbor)

  # Run changepoint detection (suppress warnings)
  cp_result <- tryCatch({
    suppressWarnings({
      fastcpd.binomial(
        data = data_matrix,
        r.progress = FALSE
      )
    })
  }, error = function(e) {
    return(NULL)
  })

  # If changepoint detection failed
  if (is.null(cp_result)) return(c(7, NA, NA))  # Failed - CPD error

  # Extract changepoints
  cp_set <- cp_result@cp_set

  # If no changepoints detected
  if (length(cp_set) == 0) return(c(3, NA, NA))  # No Pattern

  # Create data frame for significance testing
  data_df <- data.frame(y = y, lag1 = lag, neighbor = neighbor)

  # Assess significance (suppress warnings)
  sig_results <- tryCatch({
    suppressWarnings({
      assess_changepoint_significance(data_df, cp_set, alpha, n_perm)
    })
  }, error = function(e) {
    return(NULL)
  })

  # If significance testing failed
  if (is.null(sig_results) || nrow(sig_results) == 0) return(c(7, NA, NA))

  # Classify the pattern
  pattern <- classify_pattern(sig_results)

  pattern_codes <- c(
    "Always Absent" = 1,
    "Always Present" = 2,
    "No Pattern" = 3,
    "Increasing" = 4,
    "Decreasing" = 5,
    "Fluctuating/Intermittent" = 6,
    "Failed Classification" = 7
  )

  classification_code <- pattern_codes[pattern]

  # IMPORTANT: Only extract years for SIGNIFICANT changepoints
  # And only for pixels classified as "Increasing" or "Decreasing"
  year_decrease <- NA
  year_increase <- NA

  # Filter to ONLY significant changepoints
  sig_cps <- sig_results[sig_results$Overall_Significant == TRUE, ]

  if (nrow(sig_cps) > 0) {
    # Only extract year of decrease if classified as "Decreasing"
    if (pattern == "Decreasing") {
      decreasing_cps <- sig_cps[sig_cps$Seg2_Proportion < sig_cps$Seg1_Proportion, ]
      if (nrow(decreasing_cps) > 0 && length(decreasing_cps$ChangePoint) > 0) {
        first_cp_dec <- min(decreasing_cps$ChangePoint)
        year_decrease <- time_steps[first_cp_dec + 1]
      }
    }

    # Only extract year of increase if classified as "Increasing"
    if (pattern == "Increasing") {
      increasing_cps <- sig_cps[sig_cps$Seg2_Proportion > sig_cps$Seg1_Proportion, ]
      if (nrow(increasing_cps) > 0 && length(increasing_cps$ChangePoint) > 0) {
        first_cp_inc <- min(increasing_cps$ChangePoint)
        year_increase <- time_steps[first_cp_inc + 1]
      }
    }
  }

  return(c(classification_code, year_decrease, year_increase))
}

generate_spatiotemporal_predictions <- function(partition_results,
                                                hypervolume_results,
                                                model_vars,
                                                time_cols,
                                                time_steps,
                                                variable_patterns,
                                                raster_dir,
                                                mask = NULL,
                                                output_dir = "./Predictions",
                                                combine_folds = TRUE,
                                                overwrite = FALSE) {

  require(hypervolume)
  require(raster)
  require(terra)
  require(sp)
  require(dplyr)
  require(sf)

  print("Generating spatiotemporal predictions")

  ### VALIDATE AND LOAD partition_results INPUT

  if (missing(partition_results)) {
    stop("partition_results is required. Provide output from spatiotemporal_partition() as list or .rds file path")
  }

  if (is.character(partition_results)) {
    if (!file.exists(partition_results)) {
      stop(paste("partition_results file not found:", partition_results, "Working directory:", getwd()))
    }

    file_ext <- tolower(tools::file_ext(partition_results))
    if (file_ext != "rds") {
      stop(paste("partition_results must be .rds format, got:", file_ext))
    }

    print(paste("Reading partition results from:", basename(partition_results)))
    tryCatch({
      partition_data <- readRDS(partition_results)
    }, error = function(e) {
      stop(paste("Error reading partition_results .rds file:", e$message))
    })

  } else if (is.list(partition_results)) {
    print("Using provided partition results list...")
    partition_data <- partition_results
  } else {
    stop(paste("partition_results must be file path or list, got:", class(partition_results)[1]))
  }

  if (!"points_sf" %in% names(partition_data)) {
    stop(paste("partition_results missing 'points_sf' element. Available:",
               paste(names(partition_data), collapse = ", ")))
  }

  occurrence_points <- partition_data$points_sf

  if (is.null(occurrence_points) || nrow(occurrence_points) == 0) {
    stop("partition_results$points_sf is empty")
  }

  if (!"fold" %in% names(occurrence_points)) {
    stop(paste("Missing 'fold' column in partition_results$points_sf. Available columns:",
               paste(names(occurrence_points)[names(occurrence_points) != "geometry"], collapse = ", ")))
  }

  if (inherits(occurrence_points, "sf")) {
    folds <- st_drop_geometry(occurrence_points)$fold
  } else {
    folds <- occurrence_points$fold
  }

  if (all(is.na(folds))) {
    stop("All fold values are NA in partition_results")
  }

  unique_folds <- sort(unique(folds[!is.na(folds)]))
  print(paste("Loaded", nrow(occurrence_points), "points with", length(unique_folds), "folds"))

  ### VALIDATE AND LOAD hypervolume_results INPUT

  if (missing(hypervolume_results)) {
    stop("hypervolume_results is required. Provide output from build_hypervolume_models() as list or .rds file path")
  }

  if (is.character(hypervolume_results)) {
    if (!file.exists(hypervolume_results)) {
      stop(paste("hypervolume_results file not found:", hypervolume_results, "Working directory:", getwd()))
    }

    file_ext <- tolower(tools::file_ext(hypervolume_results))
    if (file_ext != "rds") {
      stop(paste("hypervolume_results must be .rds format, got:", file_ext))
    }

    print(paste("Reading hypervolume results from:", basename(hypervolume_results)))
    tryCatch({
      hv_data <- readRDS(hypervolume_results)
    }, error = function(e) {
      stop(paste("Error reading hypervolume_results .rds file:", e$message))
    })

  } else if (is.list(hypervolume_results)) {
    print("Using provided hypervolume results list...")
    hv_data <- hypervolume_results
  } else {
    stop(paste("hypervolume_results must be file path or list, got:", class(hypervolume_results)[1]))
  }

  if ("hypervolumes" %in% names(hv_data)) {
    hypervolume_models <- hv_data$hypervolumes
    print(paste("Extracted hypervolumes from hypervolume results"))
  } else if (all(sapply(hv_data, function(x) inherits(x, "Hypervolume")))) {
    hypervolume_models <- hv_data
    print("Using hypervolume list directly")
  } else {
    stop("hypervolume_results must be output from build_hypervolume_models() or a list of Hypervolume objects")
  }

  if (!is.list(hypervolume_models)) {
    stop("hypervolume_models must be a list of hypervolume objects")
  }

  if (!all(sapply(hypervolume_models, function(x) inherits(x, "Hypervolume")))) {
    stop("All elements in hypervolume_models must be Hypervolume objects")
  }

  n_folds <- length(hypervolume_models)

  if (n_folds != length(unique_folds)) {
    stop(paste("Number of hypervolume models (", n_folds,
               ") does not match number of unique folds (",
               length(unique_folds), ")", sep = ""))
  }

  print(paste("Loaded", n_folds, "hypervolume models"))

  ### VALIDATE model_vars INPUT

  if (missing(model_vars)) {
    stop("model_vars is required. Provide a character vector of model variable names")
  }

  if (!is.character(model_vars) || length(model_vars) == 0) {
    stop("model_vars must be character vector with at least one variable")
  }

  print(paste("Model variables:", paste(model_vars, collapse = ", ")))

  clean_model_vars <- gsub("^X", "", model_vars)

  missing_vars <- clean_model_vars[!clean_model_vars %in% names(occurrence_points)]
  if (length(missing_vars) > 0) {
    available_vars <- names(occurrence_points)
    available_vars <- available_vars[available_vars != "geometry"]
    warning(paste("Some model_vars not found in occurrence_points:",
                  paste(missing_vars, collapse = ", "),
                  "Available:", paste(available_vars, collapse = ", ")))
  }

  ### VALIDATE time_cols

  if (missing(time_cols)) {
    stop(paste(
      "Error: time_cols is required. Provide a character vector of time column names.",
      "",
      "### Examples:",
      "#",
      "# For single time dimension:",
      "# time_cols = \"Year\"",
      "#",
      "# For multiple time dimensions:",
      "# time_cols = c(\"Year\", \"Month\")",
      "# time_cols = c(\"Year\", \"DOY\")",
      "#",
      "# The time_cols must:",
      "#   1. Match column names in your partition_results data",
      "#   2. Match placeholders used in your variable_patterns",
      "#   3. Match column names in your time_steps data frame",
      sep = ""
    ))
  }

  if (!is.character(time_cols) || length(time_cols) == 0) {
    stop("time_cols must be a character vector with at least one time column name")
  }

  missing_time_cols <- setdiff(time_cols, names(occurrence_points))
  if (length(missing_time_cols) > 0) {
    stop(paste("time_cols not found in occurrence_points:",
               paste(missing_time_cols, collapse = ", "),
               "Available columns:",
               paste(names(occurrence_points)[names(occurrence_points) != "geometry"], collapse = ", ")))
  }

  print(paste("Time columns:", paste(time_cols, collapse = ", ")))

  ### VALIDATE time_cols match variable_patterns

  pattern_time_placeholders <- c()
  for (var_name in names(variable_patterns)) {
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
      "Your time_cols:", paste(time_cols, collapse = ", "),
      "Placeholders found in variable_patterns:", paste(pattern_time_placeholders, collapse = ", "),
      "",
      "### To fix this issue:",
      "#",
      "# Option 1: Add missing time columns to time_cols:",
      "# time_cols = c(\"", paste(c(time_cols, tolower(extra_in_patterns)), collapse = "\", \""), "\")",
      "#",
      "# Option 2: Update your variable_patterns to match time_cols:",
      "# Ensure all placeholders in patterns (YEAR, MONTH, etc.) are in time_cols",
      sep = ""
    ))
  }

  ### VALIDATE time_steps

  if (missing(time_steps)) {
    stop("time_steps is required. Provide a data frame or matrix with columns matching time_cols")
  }

  if (is.vector(time_steps) && length(time_cols) == 1) {
    time_steps_df <- data.frame(x = time_steps)
    names(time_steps_df) <- time_cols
  } else if (is.data.frame(time_steps) || is.matrix(time_steps)) {
    time_steps_df <- as.data.frame(time_steps)

    missing_cols <- setdiff(time_cols, names(time_steps_df))
    if (length(missing_cols) > 0) {
      stop(paste("time_steps missing columns:", paste(missing_cols, collapse = ", "),
                 "Required columns:", paste(time_cols, collapse = ", "),
                 "Provided columns:", paste(names(time_steps_df), collapse = ", ")))
    }
  } else {
    stop("time_steps must be a vector (for single time column), data frame, or matrix")
  }

  print(paste("Processing", nrow(time_steps_df), "time steps"))

  ### VALIDATE variable_patterns

  if (missing(variable_patterns)) {
    stop(paste(
      "Error: variable_patterns is required. Provide a named vector mapping model variables to raster file patterns.",
      "",
      "### Define variable patterns as follows in a named vector:",
      "#",
      "# my_variable_patterns <- c(",
      "#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "#   \"Temperature\" = \"Temperature_YEAR_MONTH\",",
      "#   \"elevation\" = \"elevation\"",
      "# )",
      "#",
      "# Where:",
      "#   - The NAME (left side) is the variable name matching your model variables",
      "#   - The VALUE (right side) is the pattern to match scaled raster filenames",
      "#",
      "# For time-varying variables:",
      "#   - Use placeholders like YEAR, MONTH, DAY matching your time_cols",
      "#   - Example: \"Temperature_YEAR_MONTH\" with time_cols = c(\"Year\", \"Month\")",
      "#",
      "# For static variables:",
      "#   - Use a simple pattern with no time placeholders",
      "#   - Example: \"elevation\"",
      sep = ""
    ))
  }

  if (!is.vector(variable_patterns) || is.null(names(variable_patterns))) {
    stop(paste(
      "Error: variable_patterns must be a named vector.",
      "",
      "### Define variable patterns as follows in a named vector:",
      "#",
      "# my_variable_patterns <- c(",
      "#   \"Developed_Percentage2\" = \"Developed_Percentage2_YEAR\",",
      "#   \"Temperature\" = \"Temperature_YEAR_MONTH\",",
      "#   \"elevation\" = \"elevation\"",
      "# )",
      sep = ""
    ))
  }

  if (any(names(variable_patterns) == "")) {
    stop("Error: All elements in variable_patterns must be named")
  }

  if (!all(clean_model_vars %in% names(variable_patterns))) {
    missing_patterns <- clean_model_vars[!clean_model_vars %in% names(variable_patterns)]
    stop(paste(
      "Error: Missing variable_patterns for:", paste(missing_patterns, collapse = ", "),
      "Provided patterns for:", paste(names(variable_patterns), collapse = ", ")
    ))
  }

  ### VALIDATE raster_dir

  if (missing(raster_dir)) {
    stop("raster_dir is required. Provide directory containing scaled raster files")
  }

  if (!dir.exists(raster_dir)) {
    stop(paste("raster_dir does not exist:", raster_dir, "Working directory:", getwd()))
  }

  print(paste("Reading rasters from:", raster_dir))

  ### CREATE OUTPUT DIRECTORY

  tryCatch({
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      print(paste("Created output directory:", output_dir))
    }
  }, error = function(e) {
    stop(paste("Could not create output directory:", output_dir, "-", e$message))
  })

  ### CONVERT OCCURRENCE POINTS

  if (inherits(occurrence_points, "sf")) {
    coords <- st_coordinates(occurrence_points)
    occ_df <- st_drop_geometry(occurrence_points)
    points_sp <- SpatialPointsDataFrame(
      coords, occ_df,
      proj4string = CRS("+proj=longlat +datum=WGS84")
    )
  } else if (inherits(occurrence_points, "SpatialPointsDataFrame")) {
    points_sp <- occurrence_points
  } else {
    stop("occurrence_points must be an sf object or SpatialPointsDataFrame")
  }

  points_sp$fold <- folds

  ### LOAD MASK

  if (!is.null(mask)) {
    if (is.character(mask)) {
      if (!file.exists(mask)) {
        stop(paste("Mask file not found:", mask))
      }
      tryCatch({
        mask_raster <- raster(mask)
        print(paste("Loaded mask from:", basename(mask)))
      }, error = function(e) {
        stop(paste("Error loading mask raster:", e$message))
      })
    } else if (inherits(mask, "RasterLayer")) {
      mask_raster <- mask
      print("Using provided mask raster")
    } else {
      stop("mask must be file path or RasterLayer object")
    }

    tryCatch({
      points_sp <- spTransform(points_sp, crs(mask_raster))
      print("Transformed occurrence points to mask CRS")
    }, error = function(e) {
      warning(paste("Could not transform points to mask CRS:", e$message))
    })
  } else {
    mask_raster <- NULL
  }

  ### IDENTIFY VARIABLE TYPES

  dynamic_vars <- c()
  static_vars <- c()

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

  ### INITIALIZE RESULTS

  all_model_data <- data.frame()

  ### PROCESS EACH TIME STEP

  for (i in 1:nrow(time_steps_df)) {

    time_values <- time_steps_df[i, , drop = FALSE]

    ### Create time label for output filename
    time_label <- paste(sapply(time_cols, function(tc) {
      as.character(time_values[[tc]])
    }), collapse = "_")

    output_path <- file.path(output_dir, paste0("Prediction_", time_label, ".tif"))

    if (file.exists(output_path) && !overwrite) {
      print(paste("Skipping time step", time_label, "- output already exists"))
      next
    }

    print(paste("========================================"))
    print(paste("Processing time step:", time_label))
    print(paste("========================================"))

    ### LOAD AND STACK RASTERS

    raster_files <- c()

    for (var in dynamic_vars) {
      pattern <- variable_patterns[var]

      ### Replace each time column placeholder with its value
      file_name <- pattern
      for (tc in time_cols) {
        file_name <- gsub(tc, as.character(time_values[[tc]]),
                          file_name, ignore.case = TRUE)
      }

      file_name <- paste0(file_name, "_Scaled.tif")
      raster_files <- c(raster_files, file_name)
    }

    for (var in static_vars) {
      pattern <- variable_patterns[var]
      file_name <- paste0(pattern, "_Scaled.tif")
      raster_files <- c(raster_files, file_name)
    }

    print(paste("  Loading", length(raster_files), "raster files..."))

    raster_paths <- file.path(raster_dir, raster_files)

    missing_files <- raster_paths[!file.exists(raster_paths)]
    if (length(missing_files) > 0) {
      warning(paste(
        "Missing raster files for time step", time_label, ":",
        paste(basename(missing_files), collapse = ", "),
        "  Check that your variable_patterns match the actual scaled raster filenames."
      ))
      next
    }

    tryCatch({
      s <- stack(raster_paths)
    }, error = function(e) {
      warning(paste("Error stacking rasters for time step", time_label, ":", e$message))
      next
    })

    if (!is.null(mask_raster)) {
      tryCatch({
        s_masked <- terra::mask(s, mask = sum(s))
        s2 <- stack(rast(s_masked))
      }, error = function(e) {
        warning(paste("Error applying mask for time step", time_label, ":", e$message))
        s2 <- s
      })
    } else {
      s2 <- s
    }

    s2 <- s2[[sort(names(s2))]]

    print(paste("  Raster stack created with", nlayers(s2), "layers"))

    ### PROCESS EACH FOLD

    fold_projections <- list()

    for (j in 1:n_folds) {
      current_fold <- unique_folds[j]
      print(paste("  Processing Fold", current_fold, "..."))

      hypervolume_model <- hypervolume_models[[j]]

      test_points_all <- subset(points_sp, points_sp$fold == current_fold)

      ### Filter test points for current time values
      time_filter <- rep(TRUE, nrow(points_sp@data))
      for (tc in time_cols) {
        time_filter <- time_filter & (points_sp@data[[tc]] == time_values[[tc]])
      }
      test_points_year <- subset(points_sp,
                                 points_sp$fold == current_fold & time_filter)

      tryCatch({
        hv_map <- hypervolume_project(
          hypervolume_model,
          rasters = s2,
          type = 'inclusion',
          fast.or.accurate = 'fast',
          verbose = FALSE
        )

        fold_projections[[j]] <- hv_map

      }, error = function(e) {
        warning(paste("Error projecting hypervolume for Fold", current_fold,
                      "at time step", time_label, ":", e$message))
        fold_projections[[j]] <- NULL
        next
      })

      if (is.null(fold_projections[[j]])) next

      tryCatch({
        metrics <- model_assessment_metrics(
          hypervolume_model = hypervolume_model,
          projected_raster = hv_map,
          test_points_current_year = test_points_year,
          test_points_all_years = test_points_all,
          model_vars = model_vars
        )

        ### Build time columns as a named list
        time_data <- lapply(time_cols, function(tc) time_values[[tc]])
        names(time_data) <- time_cols

        ### Build model data with time columns first, then metrics
        model_data <- data.frame(
          time_data,
          Fold = current_fold,
          G_volume = round(metrics$G_volume, 3),
          E_volume = round(metrics$E_volume, 3),
          CBP_test_G = metrics$CBP_test_G,
          TP_test_G = metrics$TP_test_G,
          FN_test_G = metrics$FN_test_G,
          sensitivity_test_G = round(metrics$sensitivity_test_G, 3),
          omission_test_G = round(metrics$omission_test_G, 3),
          CBP_test_E = metrics$CBP_test_E,
          TP_test_E = metrics$TP_test_E,
          FN_test_E = metrics$FN_test_E,
          sensitivity_test_E = round(metrics$sensitivity_test_E, 3),
          omission_test_E = round(metrics$omission_test_E, 3),
          stringsAsFactors = FALSE
        )

        all_model_data <- rbind(all_model_data, model_data)

      }, error = function(e) {
        warning(paste("Error calculating metrics for Fold", current_fold,
                      "at time step", time_label, ":", e$message))
      })
    }

    ### COMBINE FOLDS AND SAVE

    valid_projections <- fold_projections[!sapply(fold_projections, is.null)]

    if (length(valid_projections) == 0) {
      warning(paste("No valid projections for time step", time_label, "- skipping"))
      next
    }

    if (combine_folds) {
      print(paste("  Combining", length(valid_projections), "fold projections..."))

      result_raster <- valid_projections[[1]]
      if (length(valid_projections) > 1) {
        for (k in 2:length(valid_projections)) {
          result_raster <- result_raster + valid_projections[[k]]
        }
      }

      tryCatch({
        writeRaster(result_raster, output_path, overwrite = TRUE)
        print(paste("  Saved combined prediction to:", basename(output_path)))
      }, error = function(e) {
        warning(paste("Error saving combined prediction for time step", time_label, ":", e$message))
      })

    } else {
      for (k in 1:n_folds) {
        if (!is.null(fold_projections[[k]])) {
          current_fold <- unique_folds[k]
          output_path_fold <- file.path(output_dir,
                                        paste0("Prediction_Fold", current_fold,
                                               "_", time_label, ".tif"))
          tryCatch({
            writeRaster(fold_projections[[k]], output_path_fold, overwrite = TRUE)
            print(paste("  Saved Fold", current_fold, "prediction to:", basename(output_path_fold)))
          }, error = function(e) {
            warning(paste("Error saving Fold", current_fold, "prediction:", e$message))
          })
        }
      }
    }

    metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")
    tryCatch({
      write.csv(all_model_data, metrics_file, row.names = FALSE)
    }, error = function(e) {
      warning(paste("Error saving metrics:", e$message))
    })

    print(paste("  Completed time step:", time_label))
  }

  ### FINAL SAVE OF METRICS

  metrics_file <- file.path(output_dir, "Model_Assessment_Metrics.csv")
  tryCatch({
    write.csv(all_model_data, metrics_file, row.names = FALSE)
    print("========================================")
    print("All predictions complete!")
    print(paste("Metrics saved to:", metrics_file))
    print("========================================")
  }, error = function(e) {
    warning(paste("Error saving final metrics:", e$message))
  })

  return(all_model_data)
}
