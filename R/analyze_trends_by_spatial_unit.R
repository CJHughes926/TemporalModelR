#' Summarize Temporal Patterns and Trends by Spatial Unit
#'
#' Postprocessing function that aggregates temporal pattern classifications and
#' habitat change metrics across user-defined spatial units (e.g., states,
#' countries, watersheds). Generates summary tables and visualizations showing
#' how habitat suitability patterns and trends vary spatially.
#'
#' @param shapefile_path Character, sf object, or sfc object. Can be: (1) path
#'   to a shapefile, (2) path to a directory containing a single .shp file,
#'   (3) an sf object, or (4) an sfc geometry object. Spatial units are used
#'   for aggregation (e.g., administrative boundaries, watersheds).
#' @param name_field Character. Name of the attribute field in the shapefile to
#'   use as spatial unit identifiers.
#' @param binary_stack SpatRaster or character (optional). Stack of binary
#'   prediction rasters across time, or path to raster file. Required for
#'   yearly habitat summaries and time series plots.
#' @param pattern_raster SpatRaster or character (optional). Pattern
#'   classification raster from \code{\link{analyze_temporal_patterns}}, or
#'   path to raster file. Required for pattern composition summaries and
#'   scatterpie maps.
#' @param year_decrease_raster SpatRaster or character (optional). Raster
#'   showing time period of first decrease from
#'   \code{\link{analyze_temporal_patterns}}, or path to raster file. Required
#'   (along with year_increase_raster) for change event summaries.
#' @param year_increase_raster SpatRaster or character (optional). Raster
#'   showing time period of first increase from
#'   \code{\link{analyze_temporal_patterns}}, or path to raster file. Required
#'   (along with year_decrease_raster) for change event summaries.
#' @param time_steps Vector. Time period labels corresponding to layers in
#'   binary_stack. Required if binary_stack or change rasters are provided.
#' @param output_dir Character. Output directory for summary tables and plots.
#'   Default is "output/spatial_analysis".
#' @param overwrite Logical. If TRUE, overwrites existing output files. If
#'   FALSE, loads existing files when available. Default is FALSE.
#' @param pie_scale Numeric. Scaling factor for pie chart sizes in the map.
#'   Values > 1 make pies larger, values < 1 make them smaller. Default is 1.
#'
#' @return A list containing (elements present depend on inputs provided):
#' \itemize{
#'   \item overall_summary: Data frame with pattern composition statistics for
#'     each spatial unit, including raw counts, percentages of total pixels,
#'     and proportional metrics. Proportional metrics use adjusted denominators:
#'     Prop_Increasing and Prop_Stable_Suitable exclude always-absent pixels;
#'     Prop_Decreasing and Prop_Stable_Unsuitable exclude always-present pixels.
#'     (requires pattern_raster)
#'   \item yearly_summary: Data frame with habitat pixel counts by spatial unit
#'     and time period (requires binary_stack)
#'   \item change_by_year: Data frame with gain/loss pixel counts by spatial
#'     unit and time period (requires year_decrease_raster and
#'     year_increase_raster)
#'   \item plots: List of ggplot objects (contents depend on inputs provided)
#' }
#'
#' @details
#' At least one raster input (binary_stack, pattern_raster, or both change
#' rasters) must be provided. The function will generate outputs only for the
#' inputs that are supplied:
#' \itemize{
#'   \item pattern_raster: Generates pattern composition summary and scatterpie
#'     map
#'   \item binary_stack: Generates yearly habitat summary and time series plot
#'   \item year_decrease_raster + year_increase_raster: Generates change event
#'     summaries and all change-related plots
#' }
#'
#' For each spatial unit, the function extracts and summarizes (based on
#' available inputs):
#' \itemize{
#'   \item Pattern composition: Counts and percentages of pixels classified as
#'     Never Suitable, Always Suitable, No Pattern, Increasing, Decreasing, or
#'     Fluctuating
#'   \item Temporal habitat availability: Number of suitable pixels in each time
#'     period
#'   \item Change events: Number of pixels experiencing first increase or
#'     decrease in each time period
#' }
#'
#' @seealso
#' Postprocessing: \code{\link{summarize_raster_outputs}},
#' \code{\link{analyze_temporal_patterns}}
#'
#' @examples
#' \dontrun{
#' # Full analysis with all inputs
#' spatial_results <- analyze_trends_by_spatial_unit(
#'   shapefile_path = "admin_boundaries.shp",
#'   name_field = "STATE_NAME",
#'   binary_stack = "consensus_predictions/Binary_Rasters",
#'   pattern_raster = "temporal_patterns/Pattern_Classification.tif",
#'   year_decrease_raster = "temporal_patterns/Year_First_Decrease.tif",
#'   year_increase_raster = "temporal_patterns/Year_First_Increase.tif",
#'   time_steps = 2000:2020,
#'   output_dir = "spatial_analysis/"
#' )
#'
#' # Only binary stack for yearly summaries
#' spatial_results <- analyze_trends_by_spatial_unit(
#'   shapefile_path = "admin_boundaries.shp",
#'   name_field = "STATE_NAME",
#'   binary_stack = binary_rasters,
#'   time_steps = 2000:2020
#' )
#'
#' # Only pattern raster for composition analysis
#' spatial_results <- analyze_trends_by_spatial_unit(
#'   shapefile_path = spatial_units_sf,
#'   name_field = "STATE_NAME",
#'   pattern_raster = pattern_rast,
#'   pie_scale = .75
#' )
#' }
#'
#' @export
#' @importFrom sf st_read st_transform st_drop_geometry st_coordinates
#'   st_point_on_surface
#' @importFrom terra rast crs vect
#' @importFrom exactextractr exact_extract
#' @importFrom dplyr filter mutate select left_join group_by summarise bind_rows
#'   sym
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_sf geom_line geom_point geom_col geom_hline
#'   aes scale_color_manual scale_fill_manual scale_x_continuous coord_sf
#'   coord_flip facet_wrap labs theme_classic theme element_text element_rect
#'   ggsave
#' @importFrom scatterpie geom_scatterpie
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom grDevices rainbow
analyze_trends_by_spatial_unit <- function(shapefile_path,
                                           name_field,
                                           binary_stack = NULL,
                                           pattern_raster = NULL,
                                           year_decrease_raster = NULL,
                                           year_increase_raster = NULL,
                                           time_steps = NULL,
                                           output_dir = "output/spatial_analysis",
                                           overwrite = FALSE,
                                           pie_scale = 1) {

  require(sf)
  require(terra)
  require(exactextractr)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(scatterpie)
  require(utils)
  require(grDevices)


  ### Input validation

  has_binary <- !is.null(binary_stack)
  has_pattern <- !is.null(pattern_raster)
  has_change <- !is.null(year_decrease_raster) && !is.null(year_increase_raster)

  if (!has_binary && !has_pattern && !has_change) {
    stop("ERROR: At least one raster input must be provided (binary_stack, pattern_raster, or both year_decrease_raster and year_increase_raster)")
  }

  if ((has_binary || has_change) && is.null(time_steps)) {
    stop("ERROR: time_steps must be provided when using binary_stack or change rasters")
  }

  if (xor(!is.null(year_decrease_raster), !is.null(year_increase_raster))) {
    stop("ERROR: Both year_decrease_raster and year_increase_raster must be provided together, or neither")
  }

  if (inherits(shapefile_path, c("sf", "sfc"))) {
    spatial_units <- shapefile_path
    if (inherits(spatial_units, "sfc")) {
      spatial_units <- st_sf(geometry = spatial_units)
    }
  } else if (is.character(shapefile_path)) {
    if (dir.exists(shapefile_path)) {
      shp_files <- list.files(shapefile_path, pattern = "\\.shp$", full.names = TRUE)
      if (length(shp_files) == 0) {
        stop(paste0("ERROR: No .shp file found in directory: ", shapefile_path))
      }
      if (length(shp_files) > 1) {
        stop(paste0("ERROR: Multiple .shp files found in directory: ", shapefile_path,
                    ". Please specify a single file."))
      }
      print("Loading spatial units...")
      spatial_units <- st_read(shp_files[1], quiet = TRUE)
    } else if (file.exists(shapefile_path)) {
      print("Loading spatial units...")
      spatial_units <- st_read(shapefile_path, quiet = TRUE)
    } else {
      stop(paste0("ERROR: Shapefile not found: ", shapefile_path))
    }
  } else {
    stop("ERROR: shapefile_path must be either a file path (character), directory path, sf object, or sfc object")
  }

  ### Handle raster inputs as files or objects

  if (has_pattern) {
    if (is.character(pattern_raster)) {
      if (!file.exists(pattern_raster)) {
        stop(paste0("ERROR: Pattern raster file not found: ", pattern_raster))
      }
      pattern_raster <- rast(pattern_raster)
    } else if (!inherits(pattern_raster, "SpatRaster")) {
      stop("ERROR: pattern_raster must be a file path or SpatRaster object")
    }
  }

  if (has_change) {
    if (is.character(year_decrease_raster)) {
      if (!file.exists(year_decrease_raster)) {
        stop(paste0("ERROR: Year decrease raster file not found: ", year_decrease_raster))
      }
      year_decrease_raster <- rast(year_decrease_raster)
    } else if (!inherits(year_decrease_raster, "SpatRaster")) {
      stop("ERROR: year_decrease_raster must be a file path or SpatRaster object")
    }

    if (is.character(year_increase_raster)) {
      if (!file.exists(year_increase_raster)) {
        stop(paste0("ERROR: Year increase raster file not found: ", year_increase_raster))
      }
      year_increase_raster <- rast(year_increase_raster)
    } else if (!inherits(year_increase_raster, "SpatRaster")) {
      stop("ERROR: year_increase_raster must be a file path or SpatRaster object")
    }
  }

  if (has_binary) {
    if (is.character(binary_stack)) {
      if (!file.exists(binary_stack)) {
        stop(paste0("ERROR: Binary stack file not found: ", binary_stack))
      }
      binary_stack <- rast(binary_stack)
    } else if (!inherits(binary_stack, "SpatRaster")) {
      stop("ERROR: binary_stack must be a file path or SpatRaster object")
    }
  }

  ### Setup output directory

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  ### Prepare spatial units

  if (!name_field %in% names(spatial_units)) {
    stop(paste0("ERROR: Name field '", name_field, "' not found in shapefile. Available fields: ",
                paste(names(spatial_units), collapse = ", ")))
  }

  print(paste("Loaded", nrow(spatial_units), "spatial units"))

  ### Transform to raster CRS - use first available raster

  ref_raster <- if (has_pattern) pattern_raster else if (has_binary) binary_stack else year_decrease_raster
  raster_crs <- crs(ref_raster, proj = TRUE)

  if (is.na(raster_crs) || is.null(raster_crs)) {
    stop("ERROR: Reference raster has no CRS defined")
  }

  print("Transforming spatial units to raster CRS...")
  spatial_units_proj <- st_transform(spatial_units, raster_crs)

  ### Define output files

  overall_summary_file <- file.path(output_dir, "pattern_summary.csv")
  yearly_summary_file <- file.path(output_dir, "yearly_habitat.csv")
  change_by_year_file <- file.path(output_dir, "change_by_year.csv")

  ### Initialize return objects

  overall_summary <- NULL
  yearly_summary <- NULL
  change_by_year <- NULL
  plots_list <- list()

  ### Extract overall pattern classifications (if pattern_raster provided)

  if (has_pattern) {
    print("Extracting pattern classifications...")

    if (!file.exists(overall_summary_file) || overwrite) {

      pattern_extract <- suppressWarnings(suppressMessages(
        exact_extract(pattern_raster, spatial_units_proj, progress = FALSE)
      ))

      pattern_counts <- lapply(pattern_extract, function(x) {
        values <- x$value
        counts <- c(
          sum(values == 1, na.rm = TRUE),
          sum(values == 2, na.rm = TRUE),
          sum(values == 3, na.rm = TRUE),
          sum(values == 4, na.rm = TRUE),
          sum(values == 5, na.rm = TRUE),
          sum(values == 6, na.rm = TRUE),
          sum(values == 7, na.rm = TRUE)
        )
        return(counts)
      })

      pattern_matrix <- do.call(rbind, pattern_counts)

      overall_summary <- data.frame(
        Spatial_Unit = spatial_units_proj[[name_field]],
        Always_Absent = pattern_matrix[, 1],
        Always_Present = pattern_matrix[, 2],
        No_Pattern = pattern_matrix[, 3],
        Increasing = pattern_matrix[, 4],
        Decreasing = pattern_matrix[, 5],
        Fluctuating = pattern_matrix[, 6],
        Failed = pattern_matrix[, 7],
        stringsAsFactors = FALSE
      )

      overall_summary[is.na(overall_summary)] <- 0

      overall_summary$Total_Pixels <- rowSums(overall_summary[, -1])

      overall_summary$Pct_Always_Absent <- round(100 * overall_summary$Always_Absent / overall_summary$Total_Pixels, 2)
      overall_summary$Pct_Always_Present <- round(100 * overall_summary$Always_Present / overall_summary$Total_Pixels, 2)
      overall_summary$Pct_No_Pattern <- round(100 * overall_summary$No_Pattern / overall_summary$Total_Pixels, 2)
      overall_summary$Pct_Increasing <- round(100 * overall_summary$Increasing / overall_summary$Total_Pixels, 2)
      overall_summary$Pct_Decreasing <- round(100 * overall_summary$Decreasing / overall_summary$Total_Pixels, 2)
      overall_summary$Pct_Fluctuating <- round(100 * overall_summary$Fluctuating / overall_summary$Total_Pixels, 2)

      ### Proportional columns with adjusted denominators
      # For gains and stable suitable: exclude always absent (can't gain if never present)
      denom_excl_absent <- overall_summary$Total_Pixels - overall_summary$Always_Absent
      overall_summary$Prop_Increasing <- round(100 * overall_summary$Increasing / denom_excl_absent, 2)
      overall_summary$Prop_Stable_Suitable <- round(100 * overall_summary$Always_Present / denom_excl_absent, 2)

      # For losses and stable unsuitable: exclude always present (can't lose if always present)
      denom_excl_present <- overall_summary$Total_Pixels - overall_summary$Always_Present
      overall_summary$Prop_Decreasing <- round(100 * overall_summary$Decreasing / denom_excl_present, 2)
      overall_summary$Prop_Stable_Unsuitable <- round(100 * overall_summary$Always_Absent / denom_excl_present, 2)

      # Handle division by zero
      overall_summary$Prop_Increasing[denom_excl_absent == 0] <- NA
      overall_summary$Prop_Stable_Suitable[denom_excl_absent == 0] <- NA
      overall_summary$Prop_Decreasing[denom_excl_present == 0] <- NA
      overall_summary$Prop_Stable_Unsuitable[denom_excl_present == 0] <- NA

      write.csv(overall_summary, overall_summary_file, row.names = FALSE)
      print(paste("Saved:", basename(overall_summary_file)))

      gc(verbose = FALSE)

    } else {
      overall_summary <- read.csv(overall_summary_file)
      print(paste("Loaded existing:", basename(overall_summary_file)))
    }
  }

  ### Extract yearly habitat (if binary_stack provided)

  if (has_binary) {
    print("Extracting yearly habitat...")

    if (!file.exists(yearly_summary_file) || overwrite) {

      yearly_results <- list()

      pb <- txtProgressBar(min = 0, max = length(time_steps), style = 3, width = 50)

      for (i in 1:length(time_steps)) {
        year <- time_steps[i]

        year_layer <- binary_stack[[i]]

        habitat_pixels <- suppressWarnings(suppressMessages(
          exact_extract(year_layer, spatial_units_proj,
                        function(values, coverage_fractions) {
                          sum((values == 1) * coverage_fractions, na.rm = TRUE)
                        }, progress = FALSE)
        ))

        yearly_results[[i]] <- data.frame(
          Spatial_Unit = spatial_units_proj[[name_field]],
          Year = year,
          Pixels_Suitable = habitat_pixels,
          stringsAsFactors = FALSE
        )

        setTxtProgressBar(pb, i)

        if (i %% 5 == 0) gc(verbose = FALSE)
      }

      close(pb)
      print("")

      yearly_summary <- bind_rows(yearly_results)
      write.csv(yearly_summary, yearly_summary_file, row.names = FALSE)
      print(paste("Saved:", basename(yearly_summary_file)))

      gc(verbose = FALSE)

    } else {
      yearly_summary <- read.csv(yearly_summary_file)
      print(paste("Loaded existing:", basename(yearly_summary_file)))
    }
  }

  ### Extract change events (if both change rasters provided)

  if (has_change) {
    print("Extracting change events...")

    if (!file.exists(change_by_year_file) || overwrite) {

      decrease_results <- list()
      increase_results <- list()

      pb <- txtProgressBar(min = 0, max = length(time_steps), style = 3, width = 50)

      for (i in 1:length(time_steps)) {
        year <- time_steps[i]

        decrease_count <- suppressWarnings(suppressMessages(
          exact_extract(year_decrease_raster, spatial_units_proj,
                        function(values, coverage_fractions) {
                          sum(values == year, na.rm = TRUE)
                        }, progress = FALSE)
        ))

        increase_count <- suppressWarnings(suppressMessages(
          exact_extract(year_increase_raster, spatial_units_proj,
                        function(values, coverage_fractions) {
                          sum(values == year, na.rm = TRUE)
                        }, progress = FALSE)
        ))

        decrease_results[[i]] <- data.frame(
          Spatial_Unit = spatial_units_proj[[name_field]],
          Year = year,
          Decrease_Pixels = decrease_count,
          stringsAsFactors = FALSE
        )

        increase_results[[i]] <- data.frame(
          Spatial_Unit = spatial_units_proj[[name_field]],
          Year = year,
          Increase_Pixels = increase_count,
          stringsAsFactors = FALSE
        )

        setTxtProgressBar(pb, i)

        if (i %% 5 == 0) gc(verbose = FALSE)
      }

      close(pb)
      print("")

      decrease_df <- bind_rows(decrease_results)
      increase_df <- bind_rows(increase_results)
      change_by_year <- left_join(decrease_df, increase_df, by = c("Spatial_Unit", "Year"))

      write.csv(change_by_year, change_by_year_file, row.names = FALSE)
      print(paste("Saved:", basename(change_by_year_file)))

      gc(verbose = FALSE)

    } else {
      change_by_year <- read.csv(change_by_year_file)
      print(paste("Loaded existing:", basename(change_by_year_file)))
    }
  }

  ### Generate visualizations

  print("Generating visualizations...")

  all_units <- spatial_units_proj[[name_field]]
  n_units <- length(all_units)

  if (n_units <= 25) {
    c25 <- c(
      "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
      "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
      "maroon", "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
      "green1", "yellow4", "yellow3", "darkorange4", "brown"
    )
    unit_colors <- c25[1:n_units]
  } else {
    unit_colors <- rainbow(n_units, s = 0.8, v = 0.8)
  }

  names(unit_colors) <- all_units

  pattern_colors <- c(
    "Always_Absent" = "#730000",
    "Always_Present" = "#267300",
    "No_Pattern" = "#B2B2B2",
    "Increasing" = "#A3FF73",
    "Decreasing" = "#FF7F7F",
    "Fluctuating" = "#A900E6",
    "Failed" = "#000000"
  )

  ### Scatterpie map (if pattern_raster provided)

  if (has_pattern) {
    print("Creating scatterpie map...")

    pie_data <- overall_summary %>%
      left_join(
        spatial_units_proj %>%
          mutate(coords = suppressWarnings(st_coordinates(st_point_on_surface(geometry)))) %>%
          st_drop_geometry() %>%
          dplyr::select(!!sym(name_field), coords),
        by = c("Spatial_Unit" = name_field)
      ) %>%
      mutate(
        x = coords[,1],
        y = coords[,2]
      ) %>%
      dplyr::select(-coords)

    pie_cols <- c("Always_Absent", "Always_Present", "No_Pattern", "Increasing", "Decreasing", "Fluctuating", "Failed")

    for (col in pie_cols) {
      if (!col %in% names(pie_data)) {
        pie_data[[col]] <- 0
      }
    }

    total_pixels_all_units <- sum(overall_summary$Total_Pixels, na.rm = TRUE)
    pie_data$radius <- sqrt(overall_summary$Total_Pixels / total_pixels_all_units) * 0.3 * pie_scale

    p_map <- ggplot() +
      geom_sf(data = spatial_units_proj, fill = NA, color = "black", size = 0.3) +
      geom_scatterpie(aes(x = x, y = y, r = radius), data = pie_data, cols = pie_cols, color = NA) +
      scale_fill_manual(
        values = pattern_colors,
        breaks = c("Always_Absent", "Always_Present", "No_Pattern", "Increasing", "Decreasing", "Fluctuating", "Failed"),
        labels = c("Never Suitable", "Always Suitable", "No Pattern", "Increasing", "Decreasing", "Fluctuating", "Failed")
      ) +
      guides(fill = guide_legend(override.aes = list(colour = NA))) +
      coord_sf() +
      labs(title = "Pattern Composition by Spatial Unit", fill = "Pattern") +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "right"
      )

    ggsave(file.path(plots_dir, "pattern_map.png"), p_map, width = 12, height = 8, dpi = 300)
    plots_list$pattern_map <- p_map
  }

  ### Time series line plot (if binary_stack provided)

  if (has_binary) {
    print("Creating time series plot...")

    p_line <- ggplot(yearly_summary, aes(x = Year, y = Pixels_Suitable, color = Spatial_Unit)) +
      geom_line(size = 1.2, alpha = 0.8) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = unit_colors) +
      scale_x_continuous(breaks = seq(min(time_steps), max(time_steps), by = 5)) +
      labs(
        title = "Habitat Pixels Over Time",
        x = "Year",
        y = "Number of Pixels",
        color = "Spatial Unit"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(file.path(plots_dir, "time_series.png"), p_line, width = 14, height = 7, dpi = 300)
    plots_list$time_series <- p_line
  }

  ### Change-related plots (if both change rasters provided)

  if (has_change) {

    ### Overall gains/losses by year

    print("Creating annual gains/losses plot...")

    overall_change_year <- change_by_year %>%
      group_by(Year) %>%
      summarise(Total_Gain = sum(Increase_Pixels, na.rm = TRUE),
                Total_Loss = sum(Decrease_Pixels, na.rm = TRUE)) %>%
      pivot_longer(cols = c(Total_Gain, Total_Loss), names_to = "Change_Type", values_to = "Pixels") %>%
      mutate(Change_Type = ifelse(Change_Type == "Total_Gain", "Gain", "Loss"),
             Pixels = ifelse(Change_Type == "Loss", -Pixels, Pixels))

    p_overall_bar <- ggplot(overall_change_year, aes(x = Year, y = Pixels, fill = Change_Type)) +
      geom_col(width = 0.8) +
      geom_hline(yintercept = 0, color = "black", size = 0.5) +
      scale_fill_manual(values = c("Gain" = "#267300", "Loss" = "#E31A1C")) +
      scale_x_continuous(breaks = seq(min(time_steps), max(time_steps), by = 2)) +
      labs(
        title = "Annual Habitat Gains and Losses",
        x = "Year",
        y = "Pixels",
        fill = "Change Type"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    ggsave(file.path(plots_dir, "annual_change.png"), p_overall_bar, width = 12, height = 7, dpi = 300)
    plots_list$annual_change <- p_overall_bar

    ### Gains/losses by spatial unit stacked

    print("Creating stacked gains/losses plot...")

    change_by_unit <- change_by_year %>%
      pivot_longer(cols = c(Increase_Pixels, Decrease_Pixels),
                   names_to = "Change_Type", values_to = "Pixels") %>%
      mutate(Pixels = ifelse(Change_Type == "Decrease_Pixels", -Pixels, Pixels),
             Change_Type = gsub("_Pixels", "", Change_Type))

    p_unit_bar <- ggplot(change_by_unit, aes(x = factor(Year), y = Pixels, fill = Spatial_Unit)) +
      geom_col(position = "stack", color = NA) +
      geom_hline(yintercept = 0, color = "black", size = 0.8) +
      scale_fill_manual(values = unit_colors) +
      labs(
        title = "Annual Gains and Losses by Spatial Unit",
        x = "Year",
        y = "Pixels",
        fill = "Spatial Unit"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "right"
      )

    ggsave(file.path(plots_dir, "change_by_unit.png"), p_unit_bar, width = 14, height = 7, dpi = 300)
    plots_list$change_by_unit <- p_unit_bar

    ### Total gains/losses by unit

    print("Creating total change by unit plot...")

    overall_change_unit <- change_by_year %>%
      group_by(Spatial_Unit) %>%
      summarise(Total_Gain = sum(Increase_Pixels, na.rm = TRUE),
                Total_Loss = sum(Decrease_Pixels, na.rm = TRUE)) %>%
      pivot_longer(cols = c(Total_Gain, Total_Loss), names_to = "Change_Type", values_to = "Pixels") %>%
      mutate(Change_Type = ifelse(Change_Type == "Total_Gain", "Gain", "Loss"),
             Pixels = ifelse(Change_Type == "Loss", -Pixels, Pixels))

    p_overall_unit_bar <- ggplot(overall_change_unit, aes(x = reorder(Spatial_Unit, abs(Pixels)), y = Pixels, fill = Change_Type)) +
      geom_col(width = 0.8) +
      geom_hline(yintercept = 0, color = "black", size = 0.5) +
      scale_fill_manual(values = c("Gain" = "#267300", "Loss" = "#E31A1C")) +
      labs(
        title = paste("Total Gains and Losses by Unit", paste0("(", min(time_steps), "-", max(time_steps), ")")),
        x = "Unit",
        y = "Total Pixels",
        fill = "Change Type"
      ) +
      coord_flip() +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
      )

    ggsave(file.path(plots_dir, "total_change_by_unit.png"), p_overall_unit_bar, width = 10, height = 8, dpi = 300)
    plots_list$total_change_by_unit <- p_overall_unit_bar

    ### Faceted yearly gains/losses by unit

    print("Creating faceted yearly change plot...")

    p_faceted <- ggplot(change_by_unit, aes(x = Year, y = Pixels, fill = Change_Type)) +
      geom_col(width = 0.8) +
      facet_wrap(~ Spatial_Unit, scales = "free_y", ncol = 3) +
      geom_hline(yintercept = 0, color = "gray50", size = 0.3) +
      scale_fill_manual(values = c("Increase" = "#267300", "Decrease" = "#E31A1C")) +
      scale_x_continuous(breaks = seq(min(time_steps), max(time_steps), by = 10)) +
      labs(
        title = "Yearly Gains and Losses by Unit",
        x = "Year",
        y = "Pixels",
        fill = "Change Type"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "gray95", color = "gray80")
      )

    ggsave(file.path(plots_dir, "faceted_change.png"), p_faceted, width = 14, height = 10, dpi = 300)
    plots_list$faceted_change <- p_faceted
  }

  ### Summary

  print("Analysis complete")
  if (!is.null(time_steps)) {
    print(paste("Period:", min(time_steps), "-", max(time_steps)))
  }
  print(paste("Spatial units:", n_units))
  if (has_pattern) {
    print(paste("Total pixels analyzed:", format(sum(overall_summary$Total_Pixels), big.mark = ",")))
  }

  inputs_used <- c()
  if (has_pattern) inputs_used <- c(inputs_used, "pattern_raster")
  if (has_binary) inputs_used <- c(inputs_used, "binary_stack")
  if (has_change) inputs_used <- c(inputs_used, "change_rasters")
  print(paste("Inputs used:", paste(inputs_used, collapse = ", ")))

  result <- list()
  if (has_pattern) result$overall_summary <- overall_summary
  if (has_binary) result$yearly_summary <- yearly_summary
  if (has_change) result$change_by_year <- change_by_year
  if (length(plots_list) > 0) result$plots <- plots_list

  return(result)
}
