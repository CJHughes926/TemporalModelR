# Summarize Temporal Patterns and Trends by Spatial Unit

Aggregates temporal pattern classifications and change metrics across
user-defined spatial units (e.g., states, counties, watersheds). Returns
summary tables and optionally generates simple visualizations.

## Usage

``` r
analyze_trends_by_spatial_unit(shapefile_path, name_field, binary_stack = NULL,
                               pattern_raster = NULL, time_decrease_raster = NULL,
                               time_increase_raster = NULL, time_steps = NULL,
                               output_dir = NULL, overwrite = FALSE,
                               create_plot = TRUE, pie_scale = 0.15,
                               verbose = TRUE)
```

## Arguments

- shapefile_path:

  Character, sf object, or sfc object. Path to a shapefile or directory
  containing one, an sf object, or an sfc geometry. Shapefiles spatial
  units will be used as the units for the data summary.

- name_field:

  Character. Attribute field to use as spatial unit labels.

- binary_stack:

  `SpatRaster` or character. Optional stack of binary prediction rasters
  across time. Required for per-unit habitat summaries and time series
  plots.

- pattern_raster:

  `SpatRaster` or character. Optional pattern classification raster from
  [`analyze_temporal_patterns`](analyze_temporal_patterns.md). Required
  for pattern composition summaries and scatterpie map.

- time_decrease_raster:

  `SpatRaster` or character. Optional raster of first decrease time step
  from [`analyze_temporal_patterns`](analyze_temporal_patterns.md).

- time_increase_raster:

  `SpatRaster` or character. Optional raster of first increase time step
  from [`analyze_temporal_patterns`](analyze_temporal_patterns.md).

- time_steps:

  Vector. Time labels for layers in `binary_stack`. Required when
  `binary_stack` or change rasters are provided.

- output_dir:

  Character or `NULL`. Directory for CSV outputs and plots. If `NULL`,
  results are returned in memory only and no files are written. Default
  is `NULL`.

- overwrite:

  Logical. If `TRUE`, recomputes and overwrites existing CSVs. Default
  is `FALSE`.

- create_plot:

  Logical. If `TRUE` (default), generates and saves plots.

- pie_scale:

  Numeric in (0, 1\]. Largest pie's radius as a fraction of the smaller
  map dimension. Default is `0.15`, meaning the biggest pie spans
  roughly 30% of the smaller map dimension (radius = 15%, so diameter =
  30%). Smaller pies scale so their area is proportional to the unit's
  pixel count. The fraction is converted internally to coordinate units,
  so this argument works identically across CRSes (degrees, meters,
  etc.) without manual rescaling.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes details on raster extraction and per-spatial-unit
  summaries.

## Value

Invisibly returns a list containing:

- `overall_summary`: Pattern composition per spatial unit (present when
  `pattern_raster` is supplied).

- `timestep_summary`: Suitable pixel counts per unit per time step
  (present when `binary_stack` is supplied).

- `change_by_timestep`: Gain and loss pixel counts per unit per time
  step (present when both change rasters are supplied).

- `plots`: Named list of recorded plot objects (present when
  `create_plot = TRUE`).

## Details

Summarizes results from modeling and post-processing at the scale of
specific spatial blocks to allow for a nuanced look at spatiotemporal
patterns.

## See also

Post-processing:
[`summarize_raster_outputs`](summarize_raster_outputs.md),
[`analyze_temporal_patterns`](analyze_temporal_patterns.md)

## Examples

``` r
con_file <- system.file("extdata/binary/consensus_stack.tif",
                        package = "TemporalModelR")

binary_stack <- terra::rast(con_file)

study_crs <- sf::st_crs(binary_stack)

zones_sf <- rbind(
  sf::st_sf(ZONE = "West",
            geometry = sf::st_sfc(sf::st_polygon(list(
              matrix(c(0, 0, 1500, 1500, 0,
                       0, 1500, 1500, 0, 0), ncol = 2)
            )), crs = study_crs)),
  sf::st_sf(ZONE = "East",
            geometry = sf::st_sfc(sf::st_polygon(list(
              matrix(c(1500, 1500, 3000, 3000, 1500,
                       0,    1500, 1500, 0,    0),    ncol = 2)
            )), crs = study_crs))
)

time_steps <- expand.grid(
  year             = 1:15,
  season           = "Spring",
  stringsAsFactors = FALSE
)

analyze_trends_by_spatial_unit(
  shapefile_path = zones_sf,
  name_field     = "ZONE",
  binary_stack   = binary_stack,
  time_steps     = time_steps,
  create_plot    = FALSE,
  verbose        = FALSE
)
#>   |                                                          |                                                  |   0%  |                                                          |===                                               |   7%  |                                                          |=======                                           |  13%  |                                                          |==========                                        |  20%  |                                                          |=============                                     |  27%  |                                                          |=================                                 |  33%  |                                                          |====================                              |  40%  |                                                          |=======================                           |  47%  |                                                          |===========================                       |  53%  |                                                          |==============================                    |  60%  |                                                          |=================================                 |  67%  |                                                          |=====================================             |  73%  |                                                          |========================================          |  80%  |                                                          |===========================================       |  87%  |                                                          |===============================================   |  93%  |                                                          |==================================================| 100%
```
