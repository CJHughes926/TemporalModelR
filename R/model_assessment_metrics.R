#' Calculate Model Assessment Metrics in Geographic and Environmental Space
#'
#' Modeling function that computes presence-only evaluation metrics for
#' hypervolume models in both geographic space (G-space) and environmental
#' space (E-space). Calculates sensitivity, continuous binomial probability, and
#' other metrics.
#'
#' @param hypervolume_model Hypervolume object from
#'   \code{\link[hypervolume]{hypervolume_gaussian}} or
#'   \code{\link[hypervolume]{hypervolume_svm}}.
#' @param projected_raster SpatRaster. Binary projection raster (0/1) from
#'   \code{\link[hypervolume]{hypervolume_project}}.
#' @param test_points_current_year sf, SpatVector, or data frame. Test points
#'   for the current time period with model variables and coordinates.
#' @param test_points_all_years sf, SpatVector, or data frame. Test points
#'   across all time periods with model variables.
#' @param variable_patterns Named character vector mapping variable names to
#'   raster filename patterns. Used to derive model variable names.
#'
#' @return List containing assessment metrics:
#' \itemize{
#'   \item G_volume: Geographic space volume (proportion of suitable area)
#'   \item TP_test_G: True positives in geographic space
#'   \item FN_test_G: False negatives in geographic space
#'   \item sensitivity_test_G: Sensitivity in geographic space
#'   \item omission_test_G: Omission rate in geographic space
#'   \item CBP_test_G: Continuous binomial probability in geographic space
#'   \item n_test_points_G: Number of test points in geographic space
#'   \item E_volume: Environmental space volume (hypervolume size)
#'   \item TP_test_E: True positives in environmental space
#'   \item FN_test_E: False negatives in environmental space
#'   \item sensitivity_test_E: Sensitivity in environmental space
#'   \item omission_test_E: Omission rate in environmental space
#'   \item CBP_test_E: Continuous binomial probability in environmental space
#'   \item n_test_points_E: Number of test points in environmental space
#' }
#'
#' @details
#' Metrics in geographic space (G-space) evaluate model predictions on the
#' landscape by extracting raster values at test point locations. Metrics in
#' environmental space (E-space) evaluate predictions using hypervolume
#' inclusion tests in niche dimensions.
#'
#' This function is called internally by
#' \code{\link{generate_spatiotemporal_predictions}} to produce comprehensive
#' evaluation tables across all models and time periods.
#'
#' @seealso
#' Modeling: \code{\link{build_hypervolume_models}},
#' \code{\link{generate_spatiotemporal_predictions}}
#'
#' Postprocessing: \code{\link{plot_model_assessment}}
#'
#' External: \code{\link[hypervolume]{hypervolume_inclusion_test}},
#' \code{\link[hypervolume]{get_volume}},
#' \code{\link[hypervolume]{hypervolume_project}}
#'
#' @examples
#' \dontrun{
#' metrics <- model_assessment_metrics(
#'   hypervolume_model = hv,
#'   projected_raster = pred_raster,
#'   test_points_current_year = test_data_current,
#'   test_points_all_years = test_data_all,
#'   variable_patterns = c("bio1" = "bio1_YEAR", "bio12" = "bio12_YEAR")
#' )
#' }
#'
#' @export
#' @importFrom hypervolume hypervolume_inclusion_test get_volume
#' @importFrom terra extract values crs ncell vect
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom stats dbinom
model_assessment_metrics <- function(hypervolume_model,
                                     projected_raster,
                                     test_points_current_year,
                                     test_points_all_years,
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

  if (missing(test_points_current_year)) {
    stop("ERROR: test_points_current_year is required")
  }

  if (missing(test_points_all_years)) {
    stop("ERROR: test_points_all_years is required")
  }

  ### Get hypervolume name for messages
  hv_name <- if (!is.null(hypervolume_model@Name) && hypervolume_model@Name != "untitled") {
    hypervolume_model@Name
  } else {
    "Unknown"
  }

  ### Geographic space metrics (G)
  if (inherits(test_points_current_year, "sf")) {
    test_points_proj <- suppressWarnings(terra::vect(test_points_current_year))
    terra::crs(test_points_proj) <- terra::crs(projected_raster)
  } else if (inherits(test_points_current_year, "SpatVector")) {
    test_points_proj <- test_points_current_year
    terra::crs(test_points_proj) <- terra::crs(projected_raster)
  } else {
    stop("test_points_current_year must be an sf or SpatVector object")
  }

  n_test_points_G <- nrow(test_points_proj)
  if (n_test_points_G == 0) {
    warning(paste0("WARNING: ", hv_name, " has 0 year-specific test points - G-space metrics will not be calculated"))
    TP_test_G <- 0; FN_test_G <- 0; sensitivity_test_G <- NA; omission_test_G <- NA; CBP_test_G <- NA
  } else {
    if (n_test_points_G < 10) {
      warning(paste0("WARNING: ", hv_name, " has only ", n_test_points_G,
                     " year-specific test points in geographic space - results may be unreliable"))
    } else {
      print(paste("  ", hv_name, "using", n_test_points_G, "year-specific test points in geographic space"))
    }

    hv_projected_values <- terra::extract(projected_raster, test_points_proj)[,2]
    hv_projected_values[is.na(hv_projected_values)] <- 0

    TP_test_G <- sum(hv_projected_values == 1, na.rm = TRUE)
    FN_test_G <- sum(hv_projected_values == 0, na.rm = TRUE)

    if ((TP_test_G + FN_test_G) == 0) {
      warning(paste0("WARNING: No valid test point extractions for ", hv_name, " in geographic space"))
      sensitivity_test_G <- NA; omission_test_G <- NA; CBP_test_G <- NA
    } else {
      sensitivity_test_G <- TP_test_G / (TP_test_G + FN_test_G)
      omission_test_G <- 1 - sensitivity_test_G

      total_area <- terra::ncell(projected_raster) - sum(is.na(values(projected_raster)))
      total_suitable_area <- sum(values(projected_raster) == 1, na.rm = TRUE)

      CBP_test_G <- if (total_area > 0) {
        dbinom(TP_test_G, size = TP_test_G + FN_test_G, prob = total_suitable_area / total_area)
      } else {
        warning(paste0("WARNING: No valid raster cells for CBP calculation in ", hv_name))
        NA
      }
    }
  }

  ### Environmental space metrics (E)
  if (inherits(test_points_all_years, "sf")) {
    test_data <- st_drop_geometry(test_points_all_years)
  } else if (inherits(test_points_all_years, "SpatVector")) {
    test_data <- as.data.frame(test_points_all_years)
  } else if (is.data.frame(test_points_all_years)) {
    test_data <- test_points_all_years
  } else {
    stop(paste0("ERROR: test_points_all_years must be sf, SpatVector, or data.frame"))
  }

  missing_vars <- clean_model_vars[!clean_model_vars %in% names(test_data)]
  if (length(missing_vars) > 0) {
    stop(paste0("ERROR: Missing variables in test_points_all_years: ",
                paste(missing_vars, collapse = ", "),
                " Available: ", paste(names(test_data), collapse = ", ")))
  }

  test_env <- dplyr::select(test_data, all_of(clean_model_vars))
  n_before <- nrow(test_env)
  test_env <- test_env[complete.cases(test_env), , drop = FALSE]
  n_after <- nrow(test_env)

  if (n_before > n_after) {
    warning(paste0("WARNING: Removed ", n_before - n_after, " test points with missing environmental data for ", hv_name))
  }

  if (nrow(test_env) == 0) {
    warning(paste0("WARNING: No complete test points for ", hv_name, " in environmental space"))
    TP_test_E <- 0; FN_test_E <- 0; sensitivity_test_E <- NA; omission_test_E <- NA; CBP_test_E <- NA; volume_env <- NA
  } else {
    if (nrow(test_env) < 10) {
      warning(paste0("WARNING: ", hv_name, " has only ", nrow(test_env),
                     " test points in environmental space - results may be unreliable"))
    }

    hv_inclusion <- hypervolume_inclusion_test(hypervolume_model, test_env)
    TP_test_E <- sum(hv_inclusion == TRUE)
    FN_test_E <- sum(hv_inclusion == FALSE)

    if ((TP_test_E + FN_test_E) == 0) {
      warning(paste0("WARNING: No valid inclusion tests for ", hv_name, " in environmental space"))
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
      volume_env <- get_volume(hypervolume_model)
    }, error = function(e) {
      warning(paste0("WARNING: Error calculating environmental volume for ", hv_name, ": ", e$message))
      volume_env <- NA
    })
  }

  ### Volume metrics (G)
  tryCatch({
    total_cells <- terra::ncell(projected_raster) - sum(is.na(values(projected_raster)))
    suitable_cells <- sum(values(projected_raster) == 1, na.rm = TRUE)
    volume_geo <- if (total_cells > 0) suitable_cells / total_cells else NA
  }, error = function(e) {
    warning(paste0("WARNING: Error calculating geographic volume for ", hv_name, ": ", e$message))
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
