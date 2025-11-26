#' Generate time spesific model assessment metrics
#' @keywords internal
#' @importFrom hypervolume hypervolume_inclusion_test get_volume
#' @importFrom terra extract values crs ncell vect
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom stats dbinom
model_assessment_metrics <- function(hypervolume_model,
                                     projected_raster,
                                     test_points_current_time_step,
                                     test_points_all_time_steps,
                                     variable_patterns) {

  require(hypervolume)
  require(terra)
  require(sf)
  require(dplyr)

  ### Derive model_vars from variable_patterns
  model_vars <- names(variable_patterns)
  clean_model_vars <- gsub("^X", "", model_vars)

  ### Validate inputs
  if (missing(hypervolume_model) || !inherits(hypervolume_model, "Hypervolume")) {
    stop("ERROR: hypervolume_model must be a Hypervolume object")
  }

  if (missing(projected_raster) || !inherits(projected_raster, "SpatRaster")) {
    stop("ERROR: projected_raster must be a SpatRaster object")
  }

  if (missing(test_points_current_time_step)) {
    stop("ERROR: test_points_current_time_step is required")
  }

  if (missing(test_points_all_time_steps)) {
    stop("ERROR: test_points_all_time_steps is required")
  }

  ### Get hypervolume name for messages
  hv_name <- if (!is.null(hypervolume_model@Name) && hypervolume_model@Name != "untitled") {
    hypervolume_model@Name
  } else {
    "Unknown"
  }

  ### Geographic space metrics (G)
  if (inherits(test_points_current_time_step, "sf")) {
    test_points_proj <- suppressWarnings(terra::vect(test_points_current_time_step))
    terra::crs(test_points_proj) <- terra::crs(projected_raster)
  } else if (inherits(test_points_current_time_step, "SpatVector")) {
    test_points_proj <- test_points_current_time_step
    terra::crs(test_points_proj) <- terra::crs(projected_raster)
  } else {
    stop("test_points_current_time_step must be an sf or SpatVector object")
  }

  n_test_points_G <- nrow(test_points_proj)
  if (n_test_points_G == 0) {
    warning(paste0("", hv_name, " has 0 time step-specific test points - G-space metrics will not be calculated"))
    TP_test_G <- 0; FN_test_G <- 0; sensitivity_test_G <- NA; omission_test_G <- NA; CBP_test_G <- NA
  } else {
    if (n_test_points_G < 10) {
      warning(paste0("", hv_name, " has only ", n_test_points_G,
                     " time step-specific test points in geographic space - results may be unreliable"))
    } else {
      print(paste("  ", hv_name, "using", n_test_points_G, "time step-specific test points in geographic space"))
    }

    hv_projected_values <- terra::extract(projected_raster, test_points_proj)[,2]
    hv_projected_values[is.na(hv_projected_values)] <- 0

    TP_test_G <- sum(hv_projected_values == 1, na.rm = TRUE)
    FN_test_G <- sum(hv_projected_values == 0, na.rm = TRUE)

    if ((TP_test_G + FN_test_G) == 0) {
      warning(paste0("No valid test point extractions for ", hv_name, " in geographic space"))
      sensitivity_test_G <- NA; omission_test_G <- NA; CBP_test_G <- NA
    } else {
      sensitivity_test_G <- TP_test_G / (TP_test_G + FN_test_G)
      omission_test_G <- 1 - sensitivity_test_G

      total_area <- terra::ncell(projected_raster) - sum(is.na(values(projected_raster)))
      total_suitable_area <- sum(values(projected_raster) == 1, na.rm = TRUE)

      CBP_test_G <- if (total_area > 0) {
        dbinom(TP_test_G, size = TP_test_G + FN_test_G, prob = total_suitable_area / total_area)
      } else {
        warning(paste0("No valid raster cells for CBP calculation in ", hv_name))
        NA
      }
    }
  }

  ### Environmental space metrics (E)
  if (inherits(test_points_all_time_steps, "sf")) {
    test_data <- st_drop_geometry(test_points_all_time_steps)
  } else if (inherits(test_points_all_time_steps, "SpatVector")) {
    test_data <- as.data.frame(test_points_all_time_steps)
  } else if (is.data.frame(test_points_all_time_steps)) {
    test_data <- test_points_all_time_steps
  } else {
    stop(paste0("ERROR: test_points_all_time_steps must be sf, SpatVector, or data.frame"))
  }

  missing_vars <- clean_model_vars[!clean_model_vars %in% names(test_data)]
  if (length(missing_vars) > 0) {
    stop(paste0("ERROR: Missing variables in test_points_all_time_steps: ",
                paste(missing_vars, collapse = ", "),
                " Available: ", paste(names(test_data), collapse = ", ")))
  }

  test_env <- dplyr::select(test_data, all_of(clean_model_vars))
  n_before <- nrow(test_env)
  test_env <- test_env[complete.cases(test_env), , drop = FALSE]
  n_after <- nrow(test_env)

  if (n_before > n_after) {
    warning(paste0("Removed ", n_before - n_after, " test points with missing environmental data for ", hv_name))
  }

  if (nrow(test_env) == 0) {
    warning(paste0("No complete test points for ", hv_name, " in environmental space"))
    TP_test_E <- 0; FN_test_E <- 0; sensitivity_test_E <- NA; omission_test_E <- NA; CBP_test_E <- NA; volume_env <- NA
  } else {
    if (nrow(test_env) < 5) {
      warning(paste0("", hv_name, " has only ", nrow(test_env),
                     " test points in environmental space - results may be unreliable"))
    }

    hv_inclusion <- suppressMessages(suppressWarnings({
      hypervolume_inclusion_test(hypervolume_model, test_env)
    }))

    TP_test_E <- sum(hv_inclusion == TRUE)
    FN_test_E <- sum(hv_inclusion == FALSE)

    if ((TP_test_E + FN_test_E) == 0) {
      warning(paste0("No valid inclusion tests for ", hv_name, " in environmental space"))
      sensitivity_test_E <- NA; omission_test_E <- NA; CBP_test_E <- NA
    } else {
      sensitivity_test_E <- TP_test_E / (TP_test_E + FN_test_E)
      omission_test_E <- 1 - sensitivity_test_E
      total_area <- terra::ncell(projected_raster) - sum(is.na(values(projected_raster)))
      total_suitable_area <- sum(values(projected_raster) == 1, na.rm = TRUE)
      CBP_test_E <- if (total_area > 0) {
        dbinom(TP_test_E, size = TP_test_E + FN_test_E, prob = total_suitable_area / total_area)
      } else { NA }
    }

    tryCatch({
      volume_env <- suppressMessages(suppressWarnings({
        get_volume(hypervolume_model)
      }))
    }, error = function(e) {
      stop(paste0("Error calculating environmental volume for ", hv_name, ": ", e$message))
      volume_env <- NA
    })
  }

  ### Volume metrics (G)
  tryCatch({
    total_cells <- terra::ncell(projected_raster) - sum(is.na(values(projected_raster)))
    suitable_cells <- sum(values(projected_raster) == 1, na.rm = TRUE)
    volume_geo <- if (total_cells > 0) suitable_cells / total_cells else NA
  }, error = function(e) {
    stop(paste0("Error calculating geographic volume for ", hv_name, ": ", e$message))
    volume_geo <- NA
  })

  return(list(
    G_volume = volume_geo,
    TP_test_G = TP_test_G,
    FN_test_G = FN_test_G,
    sensitivity_test_G = sensitivity_test_G,
    omission_test_G = omission_test_G,
    CBP_test_G = CBP_test_G,
    n_test_points_G = n_test_points_G,

    E_volume = volume_env,
    TP_test_E = TP_test_E,
    FN_test_E = FN_test_E,
    sensitivity_test_E = sensitivity_test_E,
    omission_test_E = omission_test_E,
    CBP_test_E = CBP_test_E,
    n_test_points_E = n_after
  ))
}
#' Likelihood-ratio changepoint test
#' @keywords internal
#' @importFrom stats glm binomial logLik pchisq coef
test_cp_likelihood <- function(data, cp, alpha = 0.05, use_neighbor = TRUE) {
  require(stats)

  n <- nrow(data)
  seg1_data <- data[1:cp, ]
  seg2_data <- data[(cp + 1):n, ]

  if (nrow(seg1_data) < 3 || nrow(seg2_data) < 3) {
    return(list(p_value = NA, significant = FALSE, test_statistic = NA))
  }

  tryCatch({
    if (use_neighbor) {
      model_full <- glm(y ~ lag1 + neighbor, data = data, family = binomial())
      model_seg1 <- glm(y ~ lag1 + neighbor, data = seg1_data, family = binomial())
      model_seg2 <- glm(y ~ lag1 + neighbor, data = seg2_data, family = binomial())
    } else {
      model_full <- glm(y ~ lag1, data = data, family = binomial())
      model_seg1 <- glm(y ~ lag1, data = seg1_data, family = binomial())
      model_seg2 <- glm(y ~ lag1, data = seg2_data, family = binomial())
    }

    loglik_full <- logLik(model_full)[1]
    loglik_seg <- logLik(model_seg1)[1] + logLik(model_seg2)[1]

    lr_stat <- 2 * (loglik_seg - loglik_full)
    df <- length(coef(model_seg1)) + length(coef(model_seg2)) - length(coef(model_full))
    p_value <- 1 - pchisq(lr_stat, df = max(df, 1))

    list(p_value = p_value, significant = p_value < alpha, test_statistic = lr_stat, df = df)
  }, error = function(e) {
    list(p_value = NA, significant = FALSE, test_statistic = NA)
  })
}

#' Permutation changepoint test (difference in means)
#' @keywords internal
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

#' Two-proportion z test at a changepoint
#' @keywords internal
#' @importFrom stats pnorm
test_cp_proportion <- function(data, cp, alpha = 0.05) {
  require(stats)

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


#' Chi-square test at a changepoint
#' @keywords internal
#' @importFrom stats chisq.test
test_cp_chisquare <- function(data, cp, alpha = 0.05) {
  require(stats)

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

#' Assess significance of candidate changepoints
#' @keywords internal
assess_changepoint_significance <- function(data, cp_set, alpha = 0.05, n_perm = 1000, use_neighbor = TRUE) {
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

    lr_test <- suppressWarnings(test_cp_likelihood(temp_combined, cp_temp, alpha, use_neighbor))
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

#' Classify temporal pattern from significant changepoints
#' @keywords internal
classify_pattern <- function(sig_results) {
  sig <- sig_results[sig_results$Overall_Significant == TRUE, ]

  if (nrow(sig) == 0) {
    return("No Pattern")
  }

  tolerance <- 1e-10

  prop_diffs <- sig$Seg2_Proportion - sig$Seg1_Proportion

  has_inc <- any(prop_diffs > tolerance)
  has_dec <- any(prop_diffs < -tolerance)

  if (has_inc & has_dec) {
    return("Fluctuating/Intermittent")
  } else if (has_inc) {
    return("Increasing")
  } else if (has_dec) {
    return("Decreasing")
  } else {
    return("Failed Classification")
  }
}

#' Classify a pixel with changepoint detection
#' @keywords internal
#' @importFrom fastcpd fastcpd.binomial
classify_pixel_cpd <- function(pixel_vals, n_middle, method = "MBIC", alpha = 0.05, n_perm = 1000, use_neighbor = TRUE) {
  require(fastcpd)

  y <- pixel_vals[1:n_middle]
  lag <- pixel_vals[(n_middle + 1):(2 * n_middle)]

  if (use_neighbor) {
    neighbor <- pixel_vals[(2 * n_middle + 1):(3 * n_middle)]
    mean_val <- pixel_vals[3 * n_middle + 1]
  } else {
    neighbor <- NULL
    mean_val <- pixel_vals[2 * n_middle + 1]
  }

  if (mean_val < 0.01) return(1)
  if (mean_val > 0.99) return(2)
  if (any(is.na(c(y, lag)))) return(NA)
  if (use_neighbor && any(is.na(neighbor))) return(NA)

  if (use_neighbor) {
    data_matrix <- cbind(y, lag, neighbor)
  } else {
    data_matrix <- cbind(y, lag)
  }

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

  if (is.null(cp_result)) return(7)

  cp_set <- cp_result@cp_set

  if (length(cp_set) == 0) return(3)

  if (use_neighbor) {
    data_df <- data.frame(y = y, lag1 = lag, neighbor = neighbor)
  } else {
    data_df <- data.frame(y = y, lag1 = lag)
  }

  sig_results <- tryCatch({
    suppressWarnings({
      assess_changepoint_significance(data_df, cp_set, alpha, n_perm, use_neighbor)
    })
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(sig_results) || nrow(sig_results) == 0) return(7)

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

  return(pattern_codes[pattern])
}

#' Classify a pixel and extract first increase/decrease years
#' @keywords internal
#' @importFrom fastcpd fastcpd.binomial
classify_pixel_with_years <- function(pixel_vals, n_middle, time_steps,
                                      method = method, alpha = 0.05, n_perm = 1000, use_neighbor = TRUE) {

  require(fastcpd)

  y <- pixel_vals[1:n_middle]
  lag <- pixel_vals[(n_middle + 1):(2 * n_middle)]

  if (use_neighbor) {
    neighbor <- pixel_vals[(2 * n_middle + 1):(3 * n_middle)]
    mean_val <- pixel_vals[3 * n_middle + 1]
  } else {
    neighbor <- NULL
    mean_val <- pixel_vals[2 * n_middle + 1]
  }

  if (mean_val < 0.01) return(c(1, NA, NA))
  if (mean_val > 0.99) return(c(2, NA, NA))
  if (any(is.na(c(y, lag)))) return(c(NA, NA, NA))
  if (use_neighbor && any(is.na(neighbor))) return(c(NA, NA, NA))

  if (use_neighbor) {
    data_matrix <- cbind(y, lag, neighbor)
  } else {
    data_matrix <- cbind(y, lag)
  }

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

  if (is.null(cp_result)) return(c(7, NA, NA))

  cp_set <- cp_result@cp_set

  if (length(cp_set) == 0) return(c(3, NA, NA))

  if (use_neighbor) {
    data_df <- data.frame(y = y, lag1 = lag, neighbor = neighbor)
  } else {
    data_df <- data.frame(y = y, lag1 = lag)
  }

  sig_results <- tryCatch({
    suppressWarnings({
      assess_changepoint_significance(data_df, cp_set, alpha, n_perm, use_neighbor)
    })
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(sig_results) || nrow(sig_results) == 0) return(c(7, NA, NA))

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

  year_decrease <- NA
  year_increase <- NA

  sig_cps <- sig_results[sig_results$Overall_Significant == TRUE, ]

  if (nrow(sig_cps) > 0) {
    if (pattern == "Decreasing") {
      decreasing_cps <- sig_cps[sig_cps$Seg2_Proportion < sig_cps$Seg1_Proportion, ]
      if (nrow(decreasing_cps) > 0 && length(decreasing_cps$ChangePoint) > 0) {
        first_cp_dec <- min(decreasing_cps$ChangePoint)
        year_decrease <- time_steps[first_cp_dec + 1]
      }
    }

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
