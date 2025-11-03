#' Model assessment metrics (G- and E-space)
#' @export
#' @importFrom hypervolume hypervolume_inclusion_test get_volume
#' @importFrom raster extract values crs
#' @importFrom sp spTransform
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom stats dbinom
model_assessment_metrics <- function(hypervolume_model,
                                     projected_raster,
                                     test_points_current_year,
                                     test_points_all_years,
                                     model_vars) {

  require(hypervolume)
  require(raster)
  require(sp)
  require(dplyr)

  ### VALIDATE INPUTS

  if (missing(hypervolume_model) || !inherits(hypervolume_model, "Hypervolume")) {
    stop("hypervolume_model must be a Hypervolume object")
  }

  if (missing(projected_raster) || !inherits(projected_raster, "RasterLayer")) {
    stop("projected_raster must be a RasterLayer object")
  }

  if (missing(test_points_current_year)) {
    stop("test_points_current_year is required")
  }

  if (missing(test_points_all_years)) {
    stop("test_points_all_years is required")
  }

  if (missing(model_vars)) {
    stop("model_vars is required. Provide character vector of model variable names ")
  }

  if (!is.character(model_vars) || length(model_vars) == 0) {
    stop("model_vars must be a character vector with at least one variable ")
  }

  clean_model_vars <- gsub("^X", "", model_vars)

  ### GET HYPERVOLUME NAME FOR ERROR MESSAGES

  hv_name <- if (!is.null(hypervolume_model@Name) && hypervolume_model@Name != "untitled") {
    hypervolume_model@Name
  } else {
    "Unknown"
  }

  ### GEOGRAPHIC SPACE METRICS (G)

  n_test_points_G <- nrow(test_points_current_year)

  if (n_test_points_G == 0) {
    warning(paste(hv_name, "has 0 year-specific test points - G-space metrics will not be calculated "))

    TP_test_G <- 0
    FN_test_G <- 0
    sensitivity_test_G <- NA
    omission_test_G <- NA
    CBP_test_G <- NA

  } else {

    if (n_test_points_G < 10) {
      warning(paste(hv_name, "has only", n_test_points_G,
                    "year-specific test points in geographic space - results may be unreliable "))
    } else {
      print(paste("  ", hv_name, "using", n_test_points_G, "year-specific test points in geographic space "))
    }

    tryCatch({
      test_points_proj <- spTransform(test_points_current_year, crs(projected_raster))
    }, error = function(e) {
      stop(paste("Error transforming test points to raster CRS:", e$message))
    })

    tryCatch({
      hv_projected_values <- raster::extract(projected_raster, test_points_proj)
    }, error = function(e) {
      stop(paste("Error extracting raster values at test point locations:", e$message))
    })

    hv_projected_values[is.na(hv_projected_values)] <- 0

    TP_test_G <- sum(hv_projected_values == 1, na.rm = TRUE)
    FN_test_G <- sum(hv_projected_values == 0, na.rm = TRUE)

    if ((TP_test_G + FN_test_G) == 0) {
      warning(paste("No valid test point extractions for", hv_name, "in geographic space "))
      sensitivity_test_G <- NA
      omission_test_G <- NA
      CBP_test_G <- NA
    } else {
      sensitivity_test_G <- TP_test_G / (TP_test_G + FN_test_G)
      omission_test_G <- 1 - sensitivity_test_G

      total_area <- length(projected_raster[!is.na(projected_raster)])
      total_suitable_area <- length(projected_raster[projected_raster == 1])

      if (total_area > 0) {
        hypervolume_prob <- total_suitable_area / total_area
        CBP_test_G <- dbinom(TP_test_G, size = TP_test_G + FN_test_G,
                             prob = hypervolume_prob)
      } else {
        warning(paste("No valid raster cells for CBP calculation in", hv_name))
        CBP_test_G <- NA
      }
    }
  }

  ### ENVIRONMENTAL SPACE METRICS (E)

  if (inherits(test_points_all_years, "SpatialPointsDataFrame")) {
    test_data <- test_points_all_years@data
  } else if (inherits(test_points_all_years, "sf")) {
    test_data <- st_drop_geometry(test_points_all_years)
  } else if (is.data.frame(test_points_all_years)) {
    test_data <- test_points_all_years
  } else {
    stop(paste("test_points_all_years must be SpatialPointsDataFrame, sf, or data.frame, got:",
               class(test_points_all_years)[1]))
  }

  missing_vars <- clean_model_vars[!clean_model_vars %in% names(test_data)]
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables in test_points_all_years:",
               paste(missing_vars, collapse = ", "),
               "Available:", paste(names(test_data), collapse = ", ")))
  }

  tryCatch({
    test_env <- dplyr::select(test_data, all_of(clean_model_vars))
  }, error = function(e) {
    stop(paste("Error selecting model variables from test points:", e$message))
  })

  n_before <- nrow(test_env)
  test_env <- test_env[complete.cases(test_env), , drop = FALSE]
  n_after <- nrow(test_env)

  if (n_before > n_after) {
    warning(paste("Removed", n_before - n_after, "test points with missing environmental data for", hv_name))
  }

  if (nrow(test_env) == 0) {
    warning(paste("No complete test points for", hv_name, "in environmental space"))

    TP_test_E <- 0
    FN_test_E <- 0
    sensitivity_test_E <- NA
    omission_test_E <- NA
    CBP_test_E <- NA
    volume_env <- NA

  } else {

    n_test_points_E <- nrow(test_env)
    if (n_test_points_E < 10) {
      warning(paste(hv_name, "has only", n_test_points_E,
                    "test points in environmental space - results may be unreliable"))
    }

    tryCatch({
      hv_inclusion <- hypervolume_inclusion_test(hypervolume_model, test_env)
    }, error = function(e) {
      stop(paste("Error testing hypervolume inclusion for", hv_name, ":", e$message))
    })

    TP_test_E <- sum(hv_inclusion == TRUE)
    FN_test_E <- sum(hv_inclusion == FALSE)

    if ((TP_test_E + FN_test_E) == 0) {
      warning(paste("No valid inclusion tests for", hv_name, "in environmental space"))
      sensitivity_test_E <- NA
      omission_test_E <- NA
      CBP_test_E <- NA
    } else {
      sensitivity_test_E <- TP_test_E / (TP_test_E + FN_test_E)
      omission_test_E <- 1 - sensitivity_test_E

      total_area <- length(projected_raster[!is.na(projected_raster)])
      total_suitable_area <- length(projected_raster[projected_raster == 1])

      if (total_area > 0) {
        hypervolume_prob <- total_suitable_area / total_area
        CBP_test_E <- dbinom(TP_test_E, size = TP_test_E + FN_test_E,
                             prob = hypervolume_prob)
      } else {
        warning(paste("No valid raster cells for CBP calculation in", hv_name))
        CBP_test_E <- NA
      }
    }

    tryCatch({
      volume_env <- get_volume(hypervolume_model)
    }, error = function(e) {
      warning(paste("Error calculating environmental volume for", hv_name, ":", e$message))
      volume_env <- NA
    })
  }

  ### VOLUME METRICS

  tryCatch({
    total_cells <- sum(!is.na(values(projected_raster)))
    suitable_cells <- sum(values(projected_raster) == 1, na.rm = TRUE)

    if (total_cells > 0) {
      volume_geo <- suitable_cells / total_cells
    } else {
      warning(paste("No valid raster cells for geographic volume calculation in", hv_name))
      volume_geo <- NA
    }
  }, error = function(e) {
    warning(paste("Error calculating geographic volume for", hv_name, ":", e$message))
    volume_geo <- NA
  })

  ### RETURN METRICS

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
