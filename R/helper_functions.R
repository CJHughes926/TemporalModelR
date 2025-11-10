#' Likelihood-ratio changepoint test
#' @keywords internal
#' @importFrom stats glm binomial logLik pchisq coef
test_cp_likelihood <- function(data, cp, alpha = 0.05, use_neighbor = TRUE) {
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

library(raster)
library(fastcpd)
library(dplyr)

#' Classify a pixel and extract first increase/decrease years
#' @keywords internal
#' @importFrom fastcpd fastcpd.binomial
classify_pixel_with_years <- function(pixel_vals, n_middle, time_steps,
                                      method = method, alpha = 0.05, n_perm = 1000, use_neighbor = TRUE) {

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
