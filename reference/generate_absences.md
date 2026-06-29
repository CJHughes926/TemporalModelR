# Generate Temporally Explicit Pseudoabsence Points

Generates pseudoabsence or background points for each fold produced by
[`spatiotemporal_partition`](spatiotemporal_partition.md), distributed
across time steps proportionally to the number of presence points in
each time step within each fold. Three generation methods are supported:
random sampling within the study area, buffer-constrained sampling
around presence points, and environmentally biased sampling that targets
areas outside the known environmental tolerance of the species.
Additionally, this function can be used to process user defined absences
to work with downflow operations. This may be any list of points for
which a user wants to count as 'absences' in downflow operations,
including negative occupancy surveys, or alternatively surveys for
similar species which may act as pseudoabsences for the presence of the
species of interest.

## Usage

``` r
generate_absences(partition_result, reference_shapefile_path, raster_dir,
                  variable_patterns, method = "random", ratio = 1,
                  buffer_distance = NULL, env_percentile = 0.05,
                  time_cols = NULL, pseudoabsence_times = NULL,
                  min_points_per_timestep = 1, user_absence_data = NULL,
                  xcol = NULL, ycol = NULL, points_crs = NULL,
                  create_plot = TRUE, plot_by_fold = FALSE,
                  plot_palette = "Dark 2", output_file = NULL,
                  verbose = TRUE)
```

## Arguments

- partition_result:

  List or character. Output from
  [`spatiotemporal_partition`](spatiotemporal_partition.md) or path to
  an `.rds` file containing that output.

- reference_shapefile_path:

  Character or sf object. Path to a polygon file or an `sf` polygon
  object defining the study area.

- raster_dir:

  Character. Directory containing environmental raster files (`.tif`),
  typically the output of [`raster_align`](raster_align.md) or
  [`scale_rasters`](scale_rasters.md). File names must follow the
  patterns supplied in `variable_patterns`, with any time placeholder
  substituted for the corresponding value from `time_cols`. Required for
  all methods.

- variable_patterns:

  Named character vector mapping clean variable names to raster filename
  patterns. For time-varying variables include the time placeholder in
  the pattern (e.g. `"forest_cover" = "forest_cover_YEAR"`); for static
  variables omit it (e.g. `"elevation" = "elevation"`). Time
  placeholders must match entries in `time_cols`.

- method:

  Character. Pseudoabsence generation method. One of `"random"`,
  `"buffer"`, `"environmental"`, or `"user_data"`. Default is
  `"random"`. When `"user_data"` is specified, `user_absence_data` is
  used as the source of absence locations instead of generating them.
  `ratio`, `buffer_distance`, `env_percentile`, and
  `pseudoabsence_times` are ignored. `raster_dir` and
  `variable_patterns` are required for all methods, including
  `"user_data"`, for temporally-matched environmental extraction at the
  absence points.

- ratio:

  Numeric. Number of pseudoabsence points to generate per presence
  point. Default is `1`. Values of 2, 10, 50, etc. are accepted. Points
  are always distributed proportionally across time steps within each
  fold. Set to `0` to disable proportional allocation and use a fixed
  number of points per time step instead, in which case
  `min_points_per_timestep` must be greater than 0. `ratio` and
  `min_points_per_timestep` cannot both be 0.

- buffer_distance:

  Numeric. Distance in the units of the CRS (typically meters for
  projected CRS) within which pseudoabsence points are sampled. Required
  when `method = "buffer"`. When `method = "environmental"`, supplying a
  value automatically applies a spatial buffer constraint before
  environmental profiling, following the three-step approach of Senay et
  al. (2013). If `NULL` for the environmental method, no spatial
  constraint is applied. Default is `NULL`.

- env_percentile:

  Numeric between 0 and 1. Quantile threshold used to define the
  boundary of the known environmental tolerance when
  `method = "environmental"`. Environmental cells within this quantile
  range across all variables are excluded from pseudoabsence sampling.
  Default is `0.05` (5th to 95th percentile envelope).

- time_cols:

  Character or character vector. Name of the column(s) containing the
  time step values. Must match `time_cols` used in
  [`spatiotemporal_partition`](spatiotemporal_partition.md) and the time
  placeholders used in `variable_patterns`. Default is `NULL`.

- pseudoabsence_times:

  Vector. Optional vector of specific time step values (for the first
  time column) at which to generate pseudoabsences. When `NULL`
  (default), all time steps present in the occurrence data are used.

- min_points_per_timestep:

  Integer. Minimum number of pseudoabsence points to generate per time
  step per fold. Default is `1`. When `ratio = 0`, this value sets the
  exact (fixed) number of points generated per time step per fold,
  independent of the number of presence points. `ratio` and
  `min_points_per_timestep` cannot both be 0.

- user_absence_data:

  Character, sf object, sfc object, Spatial object, or data frame. Path
  to occurrence data (`.csv`, `.shp`, `.geojson`, `.gpkg`) or a spatial
  object to be processed as absence data for downstream operations.
  Required when `method = "user_data"`. Should be preprocessed with
  [`spatiotemporal_rarefaction`](spatiotemporal_rarefaction.md) in the
  same formatting as presence data if not already thinned.

- xcol:

  Character. Name of the x-coordinate column in `user_absence_data`.
  Required when when `method = "environmental"` and `user_absence_data`
  is a CSV file or data frame.

- ycol:

  Character. Name of the y-coordinate column in `user_absence_data`.
  Required when when `method = "environmental"` and `user_absence_data`
  is a CSV file or data frame.

- points_crs:

  Character or CRS object. CRS of the `user_absence_data`. Required when
  when `method = "environmental"` and `user_absence_data` is a CSV file
  or data frame.

- create_plot:

  Logical. If `TRUE` (default), generates diagnostic plots showing the
  spatial and temporal distribution of generated pseudoabsence points
  alongside presence points.

- plot_by_fold:

  Logical. If `TRUE`, generates one map per fold. If `FALSE` (default),
  generates a single combined map.

- plot_palette:

  Character. Name of an HCL or RColorBrewer palette used to color folds
  in diagnostic plots. Accepts any HCL palette name (see
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html)) or, if
  RColorBrewer is installed, any Brewer palette name. Default is
  `"Dark 2"`.

- output_file:

  Character. Optional path to save the result as an `.rds` file. The
  parent directory will be created if it does not exist. Default is
  `NULL`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes per-fold and per-time-step pseudoabsence counts.

## Value

Invisibly returns a list containing:

- `pseudoabsences`: An sf object of all generated pseudoabsence points
  with columns `fold`, `temporal_block`, `presence` (always 0), the time
  column(s) if provided, and extracted environmental variable values
  matched to each point's time step.

- `plots`: A named list of recorded plot objects when
  `create_plot = TRUE`. Contains `temporal_distribution` and either
  `spatial_combined` or one `spatial_fold_N` entry per fold. Plots can
  be replayed with
  [`grDevices::replayPlot()`](https://rdrr.io/r/grDevices/recordplot.html).

- `summary`: A data frame summarising points generated per fold with
  columns `fold`, `n_presences`, `n_pseudoabsences`, and
  `ratio_achieved`.

## Details

Generates sets of background data based on user-specified methodology
that can be used as pseudoabsence data for the purposes of training
presence/absence models.

The four generation methods differ in how the absence locations are
obtained:

- **Random**: Points are sampled uniformly at random from the full study
  area, excluding a negligible buffer around presence locations to
  prevent exact overlap.

- **Buffer**: Points are sampled within buffers of radius
  `buffer_distance` drawn around all fold presences, clipped to the
  reference shapefile boundary.

- **Environmental**: Raster cells whose values fall outside the species
  tolerance envelope in at least one variable are identified as
  candidates. K-means clustering then selects a spatially representative
  subset. If `buffer_distance` is supplied the environmental filtering
  is applied only within that buffered region, implementing the full
  three-step approach of Senay et al. (2013).

- **User data**: Absence locations are taken directly from
  `user_absence_data` to be used when user has a predefined set of
  absense points. Points are assigned to folds by spatial join against
  the partition fold boundaries, with unmatched points routed to
  temporal folds by time value. Environmental values are then extracted
  at the supplied locations using the same time-matched logic as the
  generated methods.

## References

Senay SD, Worner SP, Ikeda T (2013) Novel Three-Step Pseudo-Absence
Selection Technique for Improved Species Distribution Modeling. PLoS ONE
8(8): e71218.

## See also

Preprocessing:
[`spatiotemporal_partition`](spatiotemporal_partition.md),
[`temporally_explicit_extraction`](temporally_explicit_extraction.md)

Modeling: [`build_temporal_hv`](build_temporal_hv.md),
[`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_rf`](build_temporal_rf.md)

## Examples

``` r
data(tmr_partition_small, package = "TemporalModelR")

scl_dir   <- system.file("extdata/rasters_scaled",
                         package = "TemporalModelR")

ref_file  <- system.file("extdata/rasters_raw/elevation.tif",
                         package = "TemporalModelR")

study_crs <- sf::st_crs(terra::rast(ref_file))

study_area_sf <- sf::st_as_sf(sf::st_as_sfc(
  sf::st_bbox(c(xmin = 0, xmax = 3000, ymin = 0, ymax = 1500),
              crs = study_crs)
))

generate_absences(
  partition_result         = tmr_partition_small,
  reference_shapefile_path = study_area_sf,
  raster_dir               = scl_dir,
  variable_patterns        = c(
    "elevation"    = "elevation",
    "forest_cover" = "forest_cover_YEAR"
  ),
  method                   = "random",
  ratio                    = 1,
  time_cols                = c("year"),
  create_plot              = FALSE,
  verbose                  = FALSE
)
```
