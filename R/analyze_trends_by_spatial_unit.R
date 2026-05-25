utils::globalVariables(c("x", "y", "radius"))

#' Summarize Temporal Patterns and Trends by Spatial Unit
#'
#' Aggregates temporal pattern classifications and change metrics
#' across user-defined spatial units (e.g., states, counties, watersheds).
#' Returns summary tables and optionally generates simple visualizations.
#'
#' @usage
#' analyze_trends_by_spatial_unit(shapefile_path, name_field, binary_stack = NULL,
#'                                pattern_raster = NULL, time_decrease_raster = NULL,
#'                                time_increase_raster = NULL, time_steps = NULL,
#'                                output_dir = NULL, overwrite = FALSE,
#'                                create_plot = TRUE, pie_scale = 0.15,
#'                                verbose = TRUE)
#'
#' @param shapefile_path Character, sf object, or sfc object. Path to a
#'   shapefile or directory containing one, an sf object, or an sfc geometry.
#'   Shapefiles spatial units will be used as the units for the data summary.
#' @param name_field Character. Attribute field to use as spatial unit labels.
#' @param binary_stack \code{SpatRaster} or character. Optional stack of binary
#'   prediction rasters across time. Required for per-unit habitat summaries
#'   and time series plots.
#' @param pattern_raster \code{SpatRaster} or character. Optional pattern
#'   classification raster from \code{\link[TemporalModelR]{analyze_temporal_patterns}}.
#'   Required for pattern composition summaries and scatterpie map.
#' @param time_decrease_raster \code{SpatRaster} or character. Optional raster
#'   of first decrease time step from \code{\link[TemporalModelR]{analyze_temporal_patterns}}.
#' @param time_increase_raster \code{SpatRaster} or character. Optional raster
#'   of first increase time step from \code{\link[TemporalModelR]{analyze_temporal_patterns}}.
#' @param time_steps Vector. Time labels for layers in \code{binary_stack}.
#'   Required when \code{binary_stack} or change rasters are provided.
#' @param output_dir Character or \code{NULL}. Directory for CSV outputs and
#'   plots. If \code{NULL}, results are returned in memory only and no files
#'   are written. Default is \code{NULL}.
#' @param overwrite Logical. If \code{TRUE}, recomputes and overwrites existing
#'   CSVs. Default is \code{FALSE}.
#' @param create_plot Logical. If \code{TRUE} (default), generates and saves
#'   plots.
#' @param pie_scale Numeric in (0, 1]. Largest pie's radius as a fraction of
#'   the smaller map dimension. Default is \code{0.15}, meaning the biggest
#'   pie spans roughly 30\% of the smaller map dimension (radius = 15\%, so
#'   diameter = 30\%). Smaller pies scale so their area is proportional to
#'   the unit's pixel count. The fraction is converted internally to
#'   coordinate units, so this argument works identically across CRSes
#'   (degrees, meters, etc.) without manual rescaling.
#' @param verbose Logical. If \code{TRUE} (default), prints progress
#'   messages during processing. Includes details on raster extraction and
#'   per-spatial-unit summaries.
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item \code{overall_summary}: Pattern composition per spatial unit
#'     (present when \code{pattern_raster} is supplied).
#'   \item \code{timestep_summary}: Suitable pixel counts per unit per time step
#'     (present when \code{binary_stack} is supplied).
#'   \item \code{change_by_timestep}: Gain and loss pixel counts per unit per time
#'     step (present when both change rasters are supplied).
#'   \item \code{plots}: Named list of recorded plot objects (present when
#'     \code{create_plot = TRUE}).
#' }
#'
#' @details
#' Summarizes results from modeling and post-processing at the scale of specific
#' spatial blocks to allow for a nuanced look at spatiotemporal patterns.
#'
#' @seealso
#' Post-processing: \code{\link{summarize_raster_outputs}},
#'   \code{\link[TemporalModelR]{analyze_temporal_patterns}}
#'
#' @examples
#' \donttest{
#'   con_file <- system.file("extdata/binary/consensus_stack.tif",
#'                           package = "TemporalModelR")
#'
#'   binary_stack <- terra::rast(con_file)
#'
#'   study_crs <- sf::st_crs(binary_stack)
#'
#'   zones_sf <- rbind(
#'     sf::st_sf(ZONE = "West",
#'               geometry = sf::st_sfc(sf::st_polygon(list(
#'                 matrix(c(0, 0, 1500, 1500, 0,
#'                          0, 1500, 1500, 0, 0), ncol = 2)
#'               )), crs = study_crs)),
#'     sf::st_sf(ZONE = "East",
#'               geometry = sf::st_sfc(sf::st_polygon(list(
#'                 matrix(c(1500, 1500, 3000, 3000, 1500,
#'                          0,    1500, 1500, 0,    0),    ncol = 2)
#'               )), crs = study_crs))
#'   )
#'
#'   time_steps <- expand.grid(
#'     year             = 1:15,
#'     season           = "Spring",
#'     stringsAsFactors = FALSE
#'   )
#'
#'   analyze_trends_by_spatial_unit(
#'     shapefile_path = zones_sf,
#'     name_field     = "ZONE",
#'     binary_stack   = binary_stack,
#'     time_steps     = time_steps,
#'     create_plot    = FALSE,
#'     verbose        = FALSE
#'   )
#' }

#' @export
#' @importFrom sf st_read st_transform st_coordinates st_point_on_surface st_sf
#' @importFrom exactextractr exact_extract
#' @importFrom terra rast nlyr values crs
#' @importFrom graphics plot barplot lines points abline axis legend mtext par
#' @importFrom grDevices dev.cur dev.off png rainbow recordPlot
#' @importFrom utils head read.csv setTxtProgressBar txtProgressBar write.csv
analyze_trends_by_spatial_unit <- function(shapefile_path,
                                           name_field,
                                           binary_stack         = NULL,
                                           pattern_raster       = NULL,
                                           time_decrease_raster = NULL,
                                           time_increase_raster = NULL,
                                           time_steps           = NULL,
                                           output_dir           = NULL,
                                           overwrite            = FALSE,
                                           create_plot          = TRUE,
                                           pie_scale            = 0.15,
                                           verbose              = TRUE) {

  if (missing(shapefile_path)) {
    stop(paste0("ERROR: 'shapefile_path' is required. ",
                "Please provide a file path, sf object, or sfc geometry."))
  }
  if (missing(name_field)) {
    stop(paste0("ERROR: 'name_field' is required. Please specify the attribute ",
                "column to use as spatial unit labels."))
  }

  for (pkg in c("sf", "exactextractr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("ERROR: The '", pkg, "' package is required. ",
                  "Install with: install.packages('", pkg, "')"))
    }
  }

  has_binary  <- !is.null(binary_stack)
  has_pattern <- !is.null(pattern_raster)
  has_change  <- !is.null(time_decrease_raster) && !is.null(time_increase_raster)

  if (!has_binary && !has_pattern && !has_change) {
    stop("ERROR: At least one raster input must be provided.")
  }
  if ((has_binary || has_change) && is.null(time_steps)) {
    stop("ERROR: 'time_steps' must be provided when using binary_stack or change rasters.")
  }

  if (!is.null(time_steps)) {
    ts_norm           <- .normalize_time_steps(time_steps, "analyze_trends_by_spatial_unit")
    time_steps        <- ts_norm$time_steps
    secondary_filters <- ts_norm$secondary_filters
  } else {
    secondary_filters <- list()
  }
  if (xor(!is.null(time_decrease_raster), !is.null(time_increase_raster))) {
    stop("ERROR: Both 'time_decrease_raster' and 'time_increase_raster' must be provided together.")
  }

  if (is.character(shapefile_path) && dir.exists(shapefile_path)) {
    shp_files <- list.files(shapefile_path, pattern = "\\.shp$", full.names = TRUE)
    if (length(shp_files) == 0) stop(paste0("ERROR: No .shp file found in: ", shapefile_path))
    if (length(shp_files) > 1) stop(paste0("ERROR: Multiple .shp files in: ", shapefile_path))
    spatial_units <- sf::st_read(shp_files[1], quiet = TRUE)
  } else {
    spatial_units <- .load_shapefile_input(shapefile_path, "shapefile_path")
  }

  if (has_pattern) pattern_raster       <- .load_raster_input(pattern_raster,       "pattern_raster")
  if (has_change) {
    time_decrease_raster <- .load_raster_input(time_decrease_raster, "time_decrease_raster")
    time_increase_raster <- .load_raster_input(time_increase_raster, "time_increase_raster")
  }
  if (has_binary)  binary_stack <- .load_raster_input(binary_stack, "binary_stack")

  if (has_binary && length(secondary_filters) > 0) {
    lyr_names <- names(binary_stack)
    keep_idx  <- seq_len(terra::nlyr(binary_stack))
    for (sc in names(secondary_filters)) {
      val      <- as.character(secondary_filters[[sc]])
      keep_idx <- keep_idx[grepl(val, lyr_names[keep_idx], ignore.case = TRUE)]
    }
    if (length(keep_idx) == 0) {
      stop(paste0(
        "ERROR: No binary raster layers matched the secondary filter(s): ",
        paste(mapply(function(k, v) paste0(k, "=", v),
                     names(secondary_filters), secondary_filters), collapse = ", "),
        ". Layer names: ",
        paste(head(lyr_names, 5), collapse = ", "),
        if (length(lyr_names) > 5) paste0(" ... (", length(lyr_names), " total)") else ""
      ))
    }
    if (length(keep_idx) != length(time_steps)) {
      stop(paste0(
        "ERROR: After filtering by secondary column(s), found ", length(keep_idx),
        " matching raster layer(s) but time_steps has ", length(time_steps), " value(s). ",
        "Ensure time_steps values match the primary time steps in your predictions."
      ))
    }
    binary_stack <- binary_stack[[keep_idx]]
    if (verbose) message(paste0(
      "Filtered binary_stack to ", terra::nlyr(binary_stack), " layer(s) matching: ",
      paste(mapply(function(k, v) paste0(k, "=", v),
                   names(secondary_filters), secondary_filters), collapse = ", ")
    ))
  } else if (has_binary && !is.null(time_steps)) {

    if (terra::nlyr(binary_stack) > length(time_steps)) {
      match_counts <- vapply(as.character(time_steps), function(ts) {
        sum(grepl(paste0("(^|[^0-9])", ts, "([^0-9]|$)"), lyr_names))
      }, integer(1))
      ambiguous <- time_steps[match_counts > 1]
      if (length(ambiguous) > 0) {
        stop(paste0(
          "ERROR: The binary stack has ", terra::nlyr(binary_stack),
          " layers but time_steps only has ", length(time_steps), " value(s). ",
          "Time step value(s) ", paste(ambiguous, collapse = ", "),
          " each match more than one raster layer, suggesting the stack contains ",
          "predictions for multiple secondary time values (e.g. seasons). ",
          "Supply time_steps as a multi-column data frame specifying the secondary ",
          "column value to use. For example: ",
          "time_steps <- expand.grid(time_step = 1:15, season = \"Spring\", stringsAsFactors = FALSE)"
        ))
      }
    }
  }

  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  plots_dir <- if (!is.null(output_dir)) file.path(output_dir, "plots") else NULL
  if (create_plot && !is.null(plots_dir) && !dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!name_field %in% names(spatial_units)) {
    stop(paste0("ERROR: '", name_field, "' not found. Available: ",
                paste(names(spatial_units), collapse = ", ")))
  }

  if (verbose) message(paste("Loaded", nrow(spatial_units), "spatial units."))

  ref_raster  <- if (has_pattern) pattern_raster else if (has_binary) binary_stack else time_decrease_raster
  raster_crs  <- terra::crs(ref_raster, proj = TRUE)
  if (is.na(raster_crs) || is.null(raster_crs)) stop("ERROR: Reference raster has no CRS.")

  if (verbose) message("Transforming spatial units to raster CRS...")
  spatial_units_proj <- sf::st_transform(spatial_units, raster_crs)

  all_units <- trimws(as.character(spatial_units_proj[[name_field]]))
  n_units   <- length(all_units)

  overall_summary_file <- if (!is.null(output_dir)) file.path(output_dir, "pattern_summary.csv")     else NULL
  timestep_summary_file  <- if (!is.null(output_dir)) file.path(output_dir, "timestep_habitat.csv")       else NULL
  change_by_timestep_file  <- if (!is.null(output_dir)) file.path(output_dir, "change_by_timestep.csv")       else NULL

  overall_summary <- NULL
  timestep_summary  <- NULL
  change_by_timestep  <- NULL
  plots_list      <- list()

  if (has_pattern) {
    if (verbose) message("Extracting pattern classifications...")

    load_existing_pattern <- !is.null(overall_summary_file) &&
      file.exists(overall_summary_file) && !overwrite

    if (!load_existing_pattern) {

      pat_ext <- suppressWarnings(suppressMessages(
        exactextractr::exact_extract(pattern_raster, spatial_units_proj,
                                     progress = FALSE)
      ))

      pat_counts <- lapply(pat_ext, function(x) {
        v <- x$value
        vapply(seq_len(7), function(k) sum(v == k, na.rm = TRUE), integer(1))
      })

      pm <- do.call(rbind, pat_counts)

      overall_summary <- data.frame(
        Spatial_Unit   = all_units,
        Always_Absent  = pm[, 1],
        Always_Present = pm[, 2],
        No_Pattern     = pm[, 3],
        Increasing     = pm[, 4],
        Decreasing     = pm[, 5],
        Fluctuating    = pm[, 6],
        Failed         = pm[, 7],
        stringsAsFactors = FALSE
      )
      overall_summary[is.na(overall_summary)] <- 0
      overall_summary$Total_Pixels <- rowSums(overall_summary[, 2:8])

      tp <- overall_summary$Total_Pixels
      overall_summary$Pct_Always_Absent  <- round(100 * overall_summary$Always_Absent  / tp, 2)
      overall_summary$Pct_Always_Present <- round(100 * overall_summary$Always_Present / tp, 2)
      overall_summary$Pct_No_Pattern     <- round(100 * overall_summary$No_Pattern     / tp, 2)
      overall_summary$Pct_Increasing     <- round(100 * overall_summary$Increasing     / tp, 2)
      overall_summary$Pct_Decreasing     <- round(100 * overall_summary$Decreasing     / tp, 2)
      overall_summary$Pct_Fluctuating    <- round(100 * overall_summary$Fluctuating    / tp, 2)

      denom_no_absent  <- tp - overall_summary$Always_Absent
      denom_no_present <- tp - overall_summary$Always_Present
      overall_summary$Prop_Increasing       <- round(100 * overall_summary$Increasing     / denom_no_absent,  2)
      overall_summary$Prop_Stable_Suitable  <- round(100 * overall_summary$Always_Present / denom_no_absent,  2)
      overall_summary$Prop_Decreasing       <- round(100 * overall_summary$Decreasing     / denom_no_present, 2)
      overall_summary$Prop_Stable_Unsuitable <- round(100 * overall_summary$Always_Absent / denom_no_present, 2)
      overall_summary$Prop_Increasing[denom_no_absent  == 0] <- NA
      overall_summary$Prop_Stable_Suitable[denom_no_absent  == 0] <- NA
      overall_summary$Prop_Decreasing[denom_no_present == 0] <- NA
      overall_summary$Prop_Stable_Unsuitable[denom_no_present == 0] <- NA

      if (!is.null(overall_summary_file)) {
        utils::write.csv(overall_summary, overall_summary_file, row.names = FALSE)
        if (verbose) message(paste("Saved:", basename(overall_summary_file)))
      }
      gc(verbose = FALSE)

    } else {
      overall_summary <- utils::read.csv(overall_summary_file,
                                         stringsAsFactors = FALSE)
      overall_summary$Spatial_Unit <- trimws(as.character(overall_summary$Spatial_Unit))
      if (nrow(overall_summary) == 0) {
        warning("Loaded pattern_summary.csv is empty - rerun with overwrite = TRUE.")
        overall_summary <- NULL
      } else {
        count_cols <- c("Always_Absent", "Always_Present", "No_Pattern",
                        "Increasing", "Decreasing", "Fluctuating", "Failed")
        missing_cols <- count_cols[!count_cols %in% names(overall_summary)]
        if (length(missing_cols) > 0) {
          warning(paste0("pattern_summary.csv is missing columns: ",
                         paste(missing_cols, collapse = ", "),
                         " - rerun with overwrite = TRUE."))
          overall_summary <- NULL
        } else {
          if (!"Total_Pixels" %in% names(overall_summary))
            overall_summary$Total_Pixels <- rowSums(
              overall_summary[, count_cols, drop = FALSE], na.rm = TRUE
            )
          if (verbose) message(paste("Loaded existing:", basename(overall_summary_file),
                                     paste0("(", nrow(overall_summary), " units)")))
        }
      }
    }
  }

  if (has_binary) {
    if (verbose) message("Extracting per-unit habitat counts per time step...")

    load_existing_timestep <- !is.null(timestep_summary_file) &&
      file.exists(timestep_summary_file) && !overwrite

    if (!load_existing_timestep) {

      pb <- utils::txtProgressBar(min = 0, max = length(time_steps), style = 3, width = 50)
      timestep_results <- vector("list", length(time_steps))

      for (i in seq_along(time_steps)) {
        ts    <- time_steps[i]
        layer <- binary_stack[[i]]
        hab   <- suppressWarnings(suppressMessages(
          exactextractr::exact_extract(
            layer, spatial_units_proj,
            function(values, coverage_fractions)
              sum((values == 1) * coverage_fractions, na.rm = TRUE),
            progress = FALSE)
        ))
        timestep_results[[i]] <- data.frame(
          Spatial_Unit    = all_units,
          Time_Step       = ts,
          Pixels_Suitable = hab,
          stringsAsFactors = FALSE
        )
        utils::setTxtProgressBar(pb, i)
        if (i %% 5 == 0) gc(verbose = FALSE)
      }
      close(pb)
      if (verbose) message("")

      timestep_summary <- do.call(rbind, timestep_results)
      if (!is.null(timestep_summary_file)) {
        utils::write.csv(timestep_summary, timestep_summary_file, row.names = FALSE)
        if (verbose) message(paste("Saved:", basename(timestep_summary_file)))
      }
      gc(verbose = FALSE)

    } else {
      timestep_summary <- utils::read.csv(timestep_summary_file,
                                          stringsAsFactors = FALSE)
      timestep_summary$Spatial_Unit <- trimws(as.character(timestep_summary$Spatial_Unit))
      ### Support both old (Year) and new (Time_Step) column names
      if ("Year" %in% names(timestep_summary) && !"Time_Step" %in% names(timestep_summary)) {
        names(timestep_summary)[names(timestep_summary) == "Year"] <- "Time_Step"
      }
      timestep_summary$Time_Step       <- as.numeric(timestep_summary$Time_Step)
      timestep_summary$Pixels_Suitable <- as.numeric(timestep_summary$Pixels_Suitable)
      if (nrow(timestep_summary) == 0) {
        warning("Loaded timestep_habitat.csv is empty  -  rerun with overwrite = TRUE.")
        timestep_summary <- NULL
      } else {
        if (verbose) message(paste("Loaded existing:", basename(timestep_summary_file)))
      }
    }
  }

  if (has_change) {
    if (verbose) message("Extracting change events per time step...")

    load_existing_change <- !is.null(change_by_timestep_file) &&
      file.exists(change_by_timestep_file) && !overwrite

    if (!load_existing_change) {

      pb <- utils::txtProgressBar(min = 0, max = length(time_steps), style = 3, width = 50)
      dec_results <- vector("list", length(time_steps))
      inc_results <- vector("list", length(time_steps))

      for (i in seq_along(time_steps)) {
        ts  <- time_steps[i]
        dec <- suppressWarnings(suppressMessages(
          exactextractr::exact_extract(
            time_decrease_raster, spatial_units_proj,
            function(values, coverage_fractions) sum(values == ts, na.rm = TRUE),
            progress = FALSE)
        ))
        inc <- suppressWarnings(suppressMessages(
          exactextractr::exact_extract(
            time_increase_raster, spatial_units_proj,
            function(values, coverage_fractions) sum(values == ts, na.rm = TRUE),
            progress = FALSE)
        ))
        dec_results[[i]] <- data.frame(Spatial_Unit = all_units, Time_Step = ts,
                                       Decrease_Pixels = dec, stringsAsFactors = FALSE)
        inc_results[[i]] <- data.frame(Spatial_Unit = all_units, Time_Step = ts,
                                       Increase_Pixels = inc, stringsAsFactors = FALSE)
        utils::setTxtProgressBar(pb, i)
        if (i %% 5 == 0) gc(verbose = FALSE)
      }
      close(pb)
      if (verbose) message("")

      change_by_timestep <- merge(
        do.call(rbind, dec_results),
        do.call(rbind, inc_results),
        by = c("Spatial_Unit", "Time_Step"), all.x = TRUE
      )
      if (!is.null(change_by_timestep_file)) {
        utils::write.csv(change_by_timestep, change_by_timestep_file, row.names = FALSE)
        if (verbose) message(paste("Saved:", basename(change_by_timestep_file)))
      }
      gc(verbose = FALSE)

    } else {
      change_by_timestep <- utils::read.csv(change_by_timestep_file,
                                            stringsAsFactors = FALSE)
      change_by_timestep$Spatial_Unit <- trimws(as.character(change_by_timestep$Spatial_Unit))
      if ("Year" %in% names(change_by_timestep) && !"Time_Step" %in% names(change_by_timestep)) {
        names(change_by_timestep)[names(change_by_timestep) == "Year"] <- "Time_Step"
      }
      change_by_timestep$Time_Step       <- as.numeric(change_by_timestep$Time_Step)
      change_by_timestep$Decrease_Pixels <- as.numeric(change_by_timestep$Decrease_Pixels)
      change_by_timestep$Increase_Pixels <- as.numeric(change_by_timestep$Increase_Pixels)
      if (nrow(change_by_timestep) == 0) {
        warning("Loaded change_by_timestep.csv is empty  -  rerun with overwrite = TRUE.")
        change_by_timestep <- NULL
      } else {
        if (verbose) message(paste("Loaded existing:", basename(change_by_timestep_file)))
      }
    }
  }

  has_pattern <- !is.null(overall_summary)
  has_binary  <- !is.null(timestep_summary)
  has_change  <- !is.null(change_by_timestep)

  if (verbose) message(paste0("Data available - pattern: ", has_pattern,
                              " | binary: ", has_binary,
                              " | change: ", has_change))

  if (!has_pattern && !has_binary && !has_change) {
    stop(paste0(
      "ERROR: All cached CSVs in '", output_dir, "' are empty. ",
      "Rerun with overwrite = TRUE to recompute from rasters."
    ))
  }

  if (!is.null(overall_summary) && verbose) {
    message("\nPattern composition by spatial unit available in result$overall_summary.")
  }
  if (!is.null(timestep_summary) && verbose) {
    message("\nPer-unit habitat summary available in result$timestep_summary.")
  }
  if (!is.null(change_by_timestep) && verbose) {
    message("\nChange by time step available in result$change_by_timestep.")
  }

  if (!create_plot) {
    if (verbose) message("\nSpatial unit analysis complete (create_plot = FALSE, no plots generated).")
    result <- list()
    if (!is.null(overall_summary)) result$overall_summary <- overall_summary
    if (!is.null(timestep_summary))  result$timestep_summary  <- timestep_summary
    if (!is.null(change_by_timestep))  result$change_by_timestep  <- change_by_timestep
    return(invisible(result))
  }

  if (verbose) message("Generating visualizations...")

  c25 <- c(
    "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
    "green1", "yellow4", "yellow3", "darkorange4", "brown"
  )
  unit_colors <- if (n_units <= 25) c25[seq_len(n_units)] else
    grDevices::rainbow(n_units, s = 0.8, v = 0.8)
  names(unit_colors) <- all_units

  pat_cols_hex <- c(
    Always_Absent  = "#730000",
    Always_Present = "#267300",
    No_Pattern     = "#B2B2B2",
    Increasing     = "#A3FF73",
    Decreasing     = "#FF7F7F",
    Fluctuating    = "#A900E6",
    Failed         = "#000000"
  )

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  if (has_pattern) {
    has_gg <- requireNamespace("ggplot2",    quietly = TRUE)
    has_sp <- requireNamespace("scatterpie", quietly = TRUE)

    if (!has_gg || !has_sp) {
      missing_pkgs <- c(if (!has_gg) "ggplot2", if (!has_sp) "scatterpie")
      if (verbose) message(paste0(
        "NOTE: Scatterpie map skipped  -  install missing package(s): ",
        paste(missing_pkgs, collapse = ", "),
        "  ->  install.packages(c(",
        paste(paste0('"', missing_pkgs, '"'), collapse = ", "), "))"
      ))
    } else {
      if (verbose) message("Creating scatterpie map...")

      centroids       <- suppressWarnings(
        sf::st_coordinates(sf::st_point_on_surface(spatial_units_proj))
      )
      pie_data        <- overall_summary
      pie_data$x      <- centroids[match(overall_summary$Spatial_Unit, all_units), 1]
      pie_data$y      <- centroids[match(overall_summary$Spatial_Unit, all_units), 2]
      pie_cols        <- c("Always_Absent", "Always_Present", "No_Pattern",
                           "Increasing", "Decreasing", "Fluctuating", "Failed")
      for (col in pie_cols) if (!col %in% names(pie_data)) pie_data[[col]] <- 0

      ### pie_scale is a fraction of the smaller map dimension - convert to
      ### coordinate units here. This keeps the user-facing argument CRS-
      ### independent: pie_scale = 0.15 looks the same on a degrees map as on
      ### a meters map.
      if (!is.numeric(pie_scale) || length(pie_scale) != 1 ||
          pie_scale <= 0 || pie_scale > 1) {
        stop("ERROR: 'pie_scale' must be a single numeric value in (0, 1]. ",
             "It is interpreted as the largest pie's radius expressed as a ",
             "fraction of the smaller map dimension.")
      }

      bb_map        <- sf::st_bbox(spatial_units_proj)
      map_width     <- as.numeric(bb_map["xmax"] - bb_map["xmin"])
      map_height    <- as.numeric(bb_map["ymax"] - bb_map["ymin"])
      map_min_dim   <- min(map_width, map_height)
      r_max         <- pie_scale * map_min_dim

      ### Radius proportional to sqrt(pixel count) so pie *area* is
      ### proportional to count. Normalised by the largest unit so the
      ### biggest pie has radius = r_max.
      max_unit_px     <- max(overall_summary$Total_Pixels, na.rm = TRUE)
      pie_data$radius <- sqrt(overall_summary$Total_Pixels / max_unit_px) * r_max

      p_map <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = spatial_units_proj, fill = NA,
                         color = "black", linewidth = 0.3) +
        scatterpie::geom_scatterpie(
          ggplot2::aes(x = x, y = y, r = radius),
          data = pie_data, cols = pie_cols, color = NA
        ) +
        ggplot2::scale_fill_manual(
          values = pat_cols_hex,
          breaks = pie_cols,
          labels = c("Never Suitable", "Always Suitable", "No Pattern",
                     "Increasing", "Decreasing", "Fluctuating", "Failed")
        ) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(colour = NA))) +
        ggplot2::coord_sf() +
        ggplot2::labs(title = "Pattern Composition by Spatial Unit", fill = "Pattern",
                      x = NULL, y = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
          plot.title      = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
          legend.position = "right",
          axis.title      = ggplot2::element_blank(),
          axis.text       = ggplot2::element_blank(),
          axis.ticks      = ggplot2::element_blank(),
          axis.line       = ggplot2::element_blank()
        )

      if (!is.null(plots_dir)) {
        map_png <- file.path(plots_dir, "pattern_map.png")
        tryCatch({
          ggplot2::ggsave(map_png, p_map, width = 12, height = 8, dpi = 300)
          if (verbose) message(paste("Saved:", basename(map_png)))
        }, error = function(e) warning(paste("Could not save pattern map:", e$message)))
      }

      tryCatch(
        print(p_map),
        error = function(e) warning(paste("Could not render pattern map inline:", e$message))
      )
      plots_list$pattern_map <- p_map
    }
  }

  if (has_binary) {
    if (verbose) message("Creating time series plot...")

    ts_vals  <- sort(unique(as.numeric(timestep_summary$Time_Step)))
    ts_y_all <- as.numeric(timestep_summary$Pixels_Suitable)
    ts_y_max <- if (any(is.finite(ts_y_all))) max(ts_y_all, na.rm = TRUE) else 1
    ts_x_rng <- if (length(ts_vals) > 1) range(ts_vals) else c(ts_vals[1] - 1, ts_vals[1] + 1)
    x_ticks  <- pretty(ts_x_rng, n = 8)

    if (!is.null(plots_dir)) {
      ts_png <- file.path(plots_dir, "time_series.png")
      tryCatch({
        grDevices::png(ts_png, width = 1400, height = 700, res = 150)
        .plot_trend_timeseries(timestep_summary, all_units, unit_colors,
                               ts_x_rng, ts_y_max, x_ticks)
        grDevices::dev.off()
        if (verbose) message(paste("Saved:", basename(ts_png)))
      }, error = function(e) {
        warning(paste("Could not save time series:", e$message))
        if (grDevices::dev.cur() > 1) grDevices::dev.off()
      })
    }
    tryCatch({
      .plot_trend_timeseries(timestep_summary, all_units, unit_colors,
                             ts_x_rng, ts_y_max, x_ticks)
      plots_list$time_series <- grDevices::recordPlot()
    }, error = function(e) {
      warning(paste("Could not render time_series plot:", e$message))
    })
  }

  if (has_change) {

    change_has_data <- any(change_by_timestep$Increase_Pixels > 0, na.rm = TRUE) ||
      any(change_by_timestep$Decrease_Pixels > 0, na.rm = TRUE)

    if (!change_has_data) {
      dec_vals <- tryCatch(
        sort(unique(as.integer(terra::values(time_decrease_raster, mat = FALSE)),
                    na.rm = TRUE)),
        error = function(e) NULL
      )
      warning(paste0(
        "change_by_timestep has no non-zero pixel counts  -  skipping change plots.\n",
        "  time_steps range: ", min(time_steps), "-", max(time_steps), "\n",
        if (!is.null(dec_vals))
          paste0("  time_decrease_raster unique values (sample): ",
                 paste(head(dec_vals[!is.na(dec_vals)], 10), collapse = ", "))
        else
          "  Could not read time_decrease_raster values.",
        "\n  If ranges differ, check that time_steps matches the raster time step values."
      ))
    } else {

      yrs_chg <- sort(unique(change_by_timestep$Time_Step))

      if (verbose) message("Creating annual gains/losses plot...")

      gain_by_yr <- tapply(change_by_timestep$Increase_Pixels, change_by_timestep$Time_Step,
                           sum, na.rm = TRUE)
      loss_by_yr <- tapply(change_by_timestep$Decrease_Pixels, change_by_timestep$Time_Step,
                           sum, na.rm = TRUE)
      gain_v <- as.numeric(gain_by_yr[as.character(yrs_chg)])
      loss_v <- as.numeric(loss_by_yr[as.character(yrs_chg)])
      gain_v[!is.finite(gain_v)] <- 0
      loss_v[!is.finite(loss_v)] <- 0
      ac_ymax <- max(c(gain_v, loss_v), na.rm = TRUE)
      if (!is.finite(ac_ymax) || ac_ymax == 0) ac_ymax <- 1
      ac_ylim <- c(-ac_ymax * 1.15, ac_ymax * 1.15)

      if (!is.null(plots_dir)) {
        ac_png <- file.path(plots_dir, "annual_change.png")
        tryCatch({
          grDevices::png(ac_png, width = 1200, height = 700, res = 150)
          .plot_trend_change_per_step(gain_v, loss_v, yrs_chg, ac_ylim)
          grDevices::dev.off()
          if (verbose) message(paste("Saved:", basename(ac_png)))
        }, error = function(e) {
          warning(paste("Could not save annual change plot:", e$message))
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
        })
      }
      tryCatch({
        .plot_trend_change_per_step(gain_v, loss_v, yrs_chg, ac_ylim)
        plots_list$annual_change <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste("Could not render annual_change plot:", e$message))
      })

      if (verbose) message("Creating total change by unit plot...")

      unit_gain_tot <- tapply(change_by_timestep$Increase_Pixels,
                              change_by_timestep$Spatial_Unit, sum, na.rm = TRUE)
      unit_loss_tot <- tapply(change_by_timestep$Decrease_Pixels,
                              change_by_timestep$Spatial_Unit, sum, na.rm = TRUE)

      shared_units  <- intersect(names(unit_gain_tot), names(unit_loss_tot))
      shared_units  <- shared_units[order(
        as.numeric(unit_gain_tot[shared_units]) -
          as.numeric(unit_loss_tot[shared_units])
      )]
      g_tot <- as.numeric(unit_gain_tot[shared_units])
      l_tot <- as.numeric(unit_loss_tot[shared_units])
      g_tot[!is.finite(g_tot)] <- 0
      l_tot[!is.finite(l_tot)] <- 0
      tcu_ymax <- max(c(g_tot, l_tot), na.rm = TRUE)
      if (!is.finite(tcu_ymax) || tcu_ymax == 0) tcu_ymax <- 1

      if (!is.null(plots_dir)) {
        tcu_png <- file.path(plots_dir, "total_change_by_unit.png")
        tryCatch({
          grDevices::png(tcu_png, width = 1000, height = 800, res = 150)
          .plot_trend_total_unit(g_tot, l_tot, shared_units, time_steps, tcu_ymax)
          grDevices::dev.off()
          if (verbose) message(paste("Saved:", basename(tcu_png)))
        }, error = function(e) {
          warning(paste("Could not save total change plot:", e$message))
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
        })
      }
      tryCatch({
        .plot_trend_total_unit(g_tot, l_tot, shared_units, time_steps, tcu_ymax)
        plots_list$total_change_by_unit <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste("Could not render total_change_by_unit plot:", e$message))
      })

      if (verbose) message("Creating faceted time-step change plot...")

      n_cols <- min(3, n_units)
      n_rows <- ceiling(n_units / n_cols)

      if (!is.null(plots_dir)) {
        faceted_png <- file.path(plots_dir, "faceted_change.png")
        tryCatch({
          grDevices::png(faceted_png, width = 1400, height = 1000, res = 150)
          .plot_trend_facet_units(change_by_timestep, all_units, n_rows, n_cols)
          grDevices::dev.off()
          if (verbose) message(paste("Saved:", basename(faceted_png)))
        }, error = function(e) {
          warning(paste("Could not save faceted change plot:", e$message))
          if (grDevices::dev.cur() > 1) grDevices::dev.off()
        })
      }
      tryCatch({
        .plot_trend_facet_units(change_by_timestep, all_units, n_rows, n_cols)
        plots_list$faceted_change <- grDevices::recordPlot()
      }, error = function(e) {
        warning(paste("Could not render faceted_change plot:", e$message))
      })

    }
  }

  if (verbose) message("\nSpatial unit analysis complete.")
  if (!is.null(time_steps))
    if (verbose) message(paste0("Period: ", min(time_steps), "-", max(time_steps),
                                " | Spatial units: ", n_units))

  result <- list()
  if (!is.null(overall_summary)) result$overall_summary <- overall_summary
  if (!is.null(timestep_summary))  result$timestep_summary  <- timestep_summary
  if (!is.null(change_by_timestep))  result$change_by_timestep  <- change_by_timestep
  if (length(plots_list) > 0)    result$plots           <- plots_list

  invisible(result)
}
