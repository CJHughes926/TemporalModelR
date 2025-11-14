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
#' @param binary_stack SpatRaster or character. Stack of binary prediction
#'   rasters across time, or path to raster file. Typically from
#'   \code{\link{summarize_raster_outputs}}.
#' @param pattern_raster SpatRaster or character. Pattern classification raster
#'   from \code{\link{analyze_temporal_patterns}}, or path to raster file.
#' @param year_decrease_raster SpatRaster or character. Raster showing time
#'   period of first decrease from \code{\link{analyze_temporal_patterns}}, or
#'   path to raster file.
#' @param year_increase_raster SpatRaster or character. Raster showing time
#'   period of first increase from \code{\link{analyze_temporal_patterns}}, or
#'   path to raster file.
#' @param time_steps Vector. Time period labels corresponding to layers in
#'   binary_stack.
#' @param output_dir Character. Output directory for summary tables and plots.
#'   Default is "output/spatial_analysis".
#' @param overwrite Logical. If TRUE, overwrites existing output files. If
#'   FALSE, loads existing files when available. Default is FALSE.
#' @param pie_scale Numeric. Scaling factor for pie chart sizes in the map.
#'   Values > 1 make pies larger, values < 1 make them smaller. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item overall_summary: Data frame with pattern composition statistics for
#'     each spatial unit
#'   \item yearly_summary: Data frame with habitat pixel counts by spatial unit
#'     and time period
#'   \item change_by_year: Data frame with gain/loss pixel counts by spatial
#'     unit and time period
#'   \item plots: List of ggplot objects including pattern map, time series,
#'     annual change, change by unit, total change, and faceted change plots
#' }
#'
#' @details
#' For each spatial unit, the function extracts and summarizes:
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
#' Generated visualizations include:
#' \enumerate{
#'   \item Scatterpie map showing pattern composition for each spatial unit
#'   \item Time series line plot of habitat availability over time
#'   \item Bar plot of total gains and losses by time period
#'   \item Stacked bar plot of gains and losses by spatial unit
#'   \item Bar plot of total gains and losses by spatial unit
#'   \item Faceted plots showing temporal trends within each spatial unit
#' }
#'
#' Spatial units with names containing "city" (case-insensitive) are
#' automatically filtered out during analysis.
#'
#' @seealso
#' Postprocessing: \code{\link{summarize_raster_outputs}},
#' \code{\link{analyze_temporal_patterns}}
#'
#' @examples
#' \dontrun{
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
#' spatial_units_sf <- st_read("admin_boundaries.shp")
#' spatial_results <- analyze_trends_by_spatial_unit(
#'   shapefile_path = spatial_units_sf,
#'   name_field = "STATE_NAME",
#'   binary_stack = binary_rasters,
#'   pattern_raster = pattern_rast,
#'   year_decrease_raster = decrease_rast,
#'   year_increase_raster = increase_rast,
#'   time_steps = 2000:2020,
#'   pie_scale = 1.5
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
                                           binary_stack,
                                           pattern_raster,
                                           year_decrease_raster,
                                           year_increase_raster,
                                           time_steps,
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

  if (is.character(pattern_raster)) {
    if (!file.exists(pattern_raster)) {
      stop(paste0("ERROR: Pattern raster file not found: ", pattern_raster))
    }
    pattern_raster <- rast(pattern_raster)
  } else if (!inherits(pattern_raster, "SpatRaster")) {
    stop("ERROR: pattern_raster must be a file path or SpatRaster object")
  }

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

  if (is.character(binary_stack)) {
    if (!file.exists(binary_stack)) {
      stop(paste0("ERROR: Binary stack file not found: ", binary_stack))
    }
    binary_stack <- rast(binary_stack)
  } else if (!inherits(binary_stack, "SpatRaster")) {
    stop("ERROR: binary_stack must be a file path or SpatRaster object")
  }

  ### Setup output directory

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  ### Prepare spatial units

  spatial_units <- spatial_units %>%
    filter(!grepl("city", !!sym(name_field), ignore.case = TRUE))

  if (!name_field %in% names(spatial_units)) {
    stop(paste0("ERROR: Name field '", name_field, "' not found in shapefile. Available fields: ",
                paste(names(spatial_units), collapse = ", ")))
  }

  print(paste("Loaded", nrow(spatial_units), "spatial units"))

  ### Transform to raster CRS

  raster_crs <- crs(pattern_raster, proj = TRUE)

  if (is.na(raster_crs) || is.null(raster_crs)) {
    stop("ERROR: Pattern raster has no CRS defined")
  }

  print("Transforming spatial units to raster CRS...")
  spatial_units_proj <- st_transform(spatial_units, raster_crs)

  ### Define output files

  overall_summary_file <- file.path(output_dir, "pattern_summary.csv")
  yearly_summary_file <- file.path(output_dir, "yearly_habitat.csv")
  change_by_year_file <- file.path(output_dir, "change_by_year.csv")

  ### Extract overall pattern classifications

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

    write.csv(overall_summary, overall_summary_file, row.names = FALSE)
    print(paste("Saved:", basename(overall_summary_file)))

    gc(verbose = FALSE)

  } else {
    overall_summary <- read.csv(overall_summary_file)
    print(paste("Loaded existing:", basename(overall_summary_file)))
  }

  ### Extract yearly habitat and change events

  print("Extracting yearly habitat and change events...")

  if ((!file.exists(yearly_summary_file) || !file.exists(change_by_year_file)) || overwrite) {

    yearly_results <- list()
    decrease_results <- list()
    increase_results <- list()

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

      yearly_results[[i]] <- data.frame(
        Spatial_Unit = spatial_units_proj[[name_field]],
        Year = year,
        Pixels_Suitable = habitat_pixels,
        stringsAsFactors = FALSE
      )

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

    yearly_summary <- bind_rows(yearly_results)
    write.csv(yearly_summary, yearly_summary_file, row.names = FALSE)
    print(paste("Saved:", basename(yearly_summary_file)))

    decrease_df <- bind_rows(decrease_results)
    increase_df <- bind_rows(increase_results)
    change_by_year <- left_join(decrease_df, increase_df, by = c("Spatial_Unit", "Year"))

    write.csv(change_by_year, change_by_year_file, row.names = FALSE)
    print(paste("Saved:", basename(change_by_year_file)))

    gc(verbose = FALSE)

  } else {
    yearly_summary <- read.csv(yearly_summary_file)
    print(paste("Loaded existing:", basename(yearly_summary_file)))

    change_by_year <- read.csv(change_by_year_file)
    print(paste("Loaded existing:", basename(change_by_year_file)))
  }

  ### Generate visualizations

  print("Generating visualizations...")

  all_units <- unique(change_by_year$Spatial_Unit)
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

  ### Scatterpie map

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

  total_pixels_all_units <- sum(overall_summary$Total_Pixels, na.rm = TRUE)
  pie_data$radius <- sqrt(overall_summary$Total_Pixels / total_pixels_all_units) * 0.3 * pie_scale

  p_map <- ggplot() +
    geom_sf(data = spatial_units_proj, fill = NA, color = "black", size = 0.3) +
    geom_scatterpie(aes(x = x, y = y, r = radius), data = pie_data, cols = pie_cols, color = NA) +
    scale_fill_manual(
      values = pattern_colors,
      labels = c("Never Suitable", "Always Suitable", "No Pattern", "Increasing", "Decreasing", "Fluctuating", "Failed")
    ) +
    coord_sf() +
    labs(title = "Pattern Composition by Spatial Unit", fill = "Pattern") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      legend.position = "right"
    )

  ggsave(file.path(plots_dir, "pattern_map.png"), p_map, width = 12, height = 8, dpi = 300)

  ### Time series line plot

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

  ### Summary

  print("Analysis complete")
  print(paste("Period:", min(time_steps), "-", max(time_steps)))
  print(paste("Spatial units:", n_units))
  print(paste("Total pixels analyzed:", format(sum(overall_summary$Total_Pixels), big.mark = ",")))

  return(list(
    overall_summary = overall_summary,
    yearly_summary = yearly_summary,
    change_by_year = change_by_year,
    plots = list(
      pattern_map = p_map,
      time_series = p_line,
      annual_change = p_overall_bar,
      change_by_unit = p_unit_bar,
      total_change_by_unit = p_overall_unit_bar,
      faceted_change = p_faceted
    )
  ))
}
