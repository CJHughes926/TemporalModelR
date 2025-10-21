model_assessment_metrics <- function(hypervolume_model,
                                     projected_raster,
                                     test_points_current_year,
                                     test_points_all_years,
                                     model_vars) {

  require(hypervolume)
  require(raster)
  require(sp)
  require(dplyr)

  # Clean variable names
  clean_model_vars <- gsub("^X", "", model_vars)

  # =============================
  # GEOGRAPHIC SPACE METRICS (G)
  # Using only current year test points
  # =============================
  if (nrow(test_points_current_year) > 0) {
    # Ensure CRS match
    test_points_proj <- spTransform(test_points_current_year, crs(projected_raster))

    # Extract hypervolume values at test point locations
    hv_projected_values <- raster::extract(projected_raster, test_points_proj)
    hv_projected_values[is.na(hv_projected_values)] <- 0

    # Calculate confusion matrix elements
    TP_test_G <- sum(hv_projected_values == 1, na.rm = TRUE)
    FN_test_G <- sum(hv_projected_values == 0, na.rm = TRUE)

    # Calculate sensitivity and omission
    sensitivity_test_G <- TP_test_G / (TP_test_G + FN_test_G)
    omission_test_G <- 1 - sensitivity_test_G

  } else {
    TP_test_G <- 0
    FN_test_G <- 0
    sensitivity_test_G <- NA
    omission_test_G <- NA
  }

  # =============================
  # ENVIRONMENTAL SPACE METRICS (E)
  # Using all years of test points
  # =============================
  if (inherits(test_points_all_years, "SpatialPointsDataFrame")) {
    test_env <- dplyr::select(test_points_all_years@data, all_of(clean_model_vars))
  } else {
    test_env <- dplyr::select(as.data.frame(test_points_all_years), all_of(clean_model_vars))
  }

  test_env <- test_env[complete.cases(test_env), , drop = FALSE]

  if (nrow(test_env) > 0) {
    # Test inclusion in hypervolume
    hv_inclusion <- hypervolume_inclusion_test(hypervolume_model, test_env)

    # Calculate confusion matrix elements
    TP_test_E <- sum(hv_inclusion == TRUE)
    FN_test_E <- sum(hv_inclusion == FALSE)

    # Calculate sensitivity and omission
    sensitivity_test_E <- TP_test_E / (TP_test_E + FN_test_E)
    omission_test_E <- 1 - sensitivity_test_E

  } else {
    TP_test_E <- 0
    FN_test_E <- 0
    sensitivity_test_E <- NA
    omission_test_E <- NA
  }

  # =============================
  # VOLUME METRICS
  # =============================
  # Geographic volume (proportion of area predicted suitable)
  volume_geo <- sum(values(projected_raster) == 1, na.rm = TRUE) /
    sum(!is.na(values(projected_raster)))

  # Environmental volume
  volume_env <- get_volume(hypervolume_model)

  # =============================
  # CONTINUOUS BOYCE INDEX (CBP)
  # =============================
  # Calculate hypervolume-based probability
  total_area <- length(projected_raster[!is.na(projected_raster)])
  total_suitable_area <- length(projected_raster[projected_raster == 1])
  hypervolume_prob <- total_suitable_area / total_area

  # Calculate CBP for geographic space
  if (!is.na(sensitivity_test_G) && (TP_test_G + FN_test_G) > 0) {
    CBP_test_G <- dbinom(TP_test_G, size = TP_test_G + FN_test_G,
                         prob = hypervolume_prob)
  } else {
    CBP_test_G <- NA
  }

  # Calculate CBP for environmental space
  if (!is.na(sensitivity_test_E) && (TP_test_E + FN_test_E) > 0) {
    CBP_test_E <- dbinom(TP_test_E, size = TP_test_E + FN_test_E,
                         prob = hypervolume_prob)
  } else {
    CBP_test_E <- NA
  }

  # =============================
  # RETURN METRICS
  # =============================
  return(list(
    # Geographic space metrics
    G_volume = volume_geo,
    TP_test_G = TP_test_G,
    FN_test_G = FN_test_G,
    sensitivity_test_G = sensitivity_test_G,
    omission_test_G = omission_test_G,
    CBP_test_G = CBP_test_G,

    # Environmental space metrics
    E_volume = volume_env,
    TP_test_E = TP_test_E,
    FN_test_E = FN_test_E,
    sensitivity_test_E = sensitivity_test_E,
    omission_test_E = omission_test_E,
    CBP_test_E = CBP_test_E
  ))
}
