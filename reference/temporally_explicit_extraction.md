# Extract Time-Aligned Environmental Values at Species Occurrences

Preprocessing function that extracts raster values to species occurrence
records based on temporal components. Matches environmental layers to
occurrence timestamps and optionally computes scaling parameters for
standardization.

## Usage

``` r
temporally_explicit_extraction(points_sp, raster_dir, variable_patterns,
                               time_cols, xcol = NULL, ycol = NULL,
                               points_crs = NULL, output_dir,
                               output_prefix = "temp_explicit_df",
                               save_raw = TRUE, save_scaled = TRUE,
                               save_scaling_params = TRUE,
                               verbose = TRUE)
```

## Arguments

- points_sp:

  sf object, SpatialPointsDataFrame, file path to
  .csv/.shp/.geojson/.gpkg, or data frame with coordinate columns.

- raster_dir:

  Character. Directory containing environmental raster files (`.tif`),
  typically the output of [`raster_align`](raster_align.md). File names
  must follow the patterns supplied in `variable_patterns`, with any
  time placeholder substituted for the corresponding value from
  `time_cols`.

- variable_patterns:

  Named character vector mapping clean variable names to raster filename
  patterns. For time-varying variables include the time placeholder in
  the pattern (e.g. `"forest_cover" = "forest_cover_YEAR"`); for static
  variables omit it (e.g. `"elevation" = "elevation"`). Time
  placeholders must match entries in `time_cols`.

- time_cols:

  Character vector of time column names present in the point data (e.g.,
  c("YEAR"), c("YEAR", "MONTH")).

- xcol:

  Character. Name of the x-coordinate column. Required when `points_sp`
  is a CSV file or data frame.

- ycol:

  Character. Name of the y-coordinate column. Required when `points_sp`
  is a CSV file or data frame.

- points_crs:

  Character or CRS object. CRS of the input points. Required when
  `points_sp` is a CSV file or data frame.

- output_dir:

  Character. Directory to write output files.

- output_prefix:

  Character. Prefix for output filenames. Default is "temp_explicit_df".

- save_raw:

  Logical. If `TRUE` (default), writes raw extracted values CSV. If
  `FALSE`, skips raw values output.

- save_scaled:

  Logical. If `TRUE` (default), writes z-scaled values CSV. If `FALSE`,
  skips scaled values output.

- save_scaling_params:

  Logical. If `TRUE` (default), writes CSV of per-variable means and
  standard deviations. If `FALSE`, skips scaling parameters output.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes file loading, extraction progress, and file-save
  confirmation.

## Value

Invisibly returns a list containing:

- `raw_values`: Data frame of raw extracted values at each occurrence
  record (when `save_raw = TRUE`; `NULL` otherwise).

- `scaled_values`: Data frame of z-scaled extracted values (when
  `save_scaled = TRUE`; `NULL` otherwise).

- `scaling_params`: Data frame of per-variable means and standard
  deviations used for scaling (when `save_scaling_params = TRUE`; `NULL`
  otherwise). Pass this to [`scale_rasters`](scale_rasters.md).

- `files_created`: Named list of file paths written, with elements
  `raw`, `scaled`, and `scaling_params` (each `NULL` when the
  corresponding save flag is `FALSE`).

## Details

Extracts raster values to species occurrence records based on matched
temporal components. Matches environmental layers to occurrence
timestamps and optionally computes scaling parameters for
standardization.

Output CSV files are written to output_dir containing raw values, scaled
values, and scaling parameters.

Scaling parameters (mean and standard deviation) are optionally computed
across all occurrence records for each variable. These parameters should
be used with [`scale_rasters`](scale_rasters.md) to standardize
prediction layers.

## See also

Preprocessing:
[`spatiotemporal_rarefaction`](spatiotemporal_rarefaction.md),
[`scale_rasters`](scale_rasters.md),
[`spatiotemporal_partition`](spatiotemporal_partition.md)

## Examples

``` r
pts_file <- system.file(
  "extdata/points/synthetic_occurrence_points.csv",
  package = "TemporalModelR"
)

aln_dir  <- system.file("extdata/rasters_aligned",
                        package = "TemporalModelR")

ref_file <- system.file("extdata/rasters_raw/elevation.tif",
                        package = "TemporalModelR")

out_dir  <- file.path(tempdir(), "extracted")

temporally_explicit_extraction(
  points_sp           = pts_file,
  raster_dir          = aln_dir,
  variable_patterns   = c(
    "elevation"    = "elevation",
    "forest_cover" = "forest_cover_YEAR",
    "prseas"       = "prseas_YEAR_SEASON"
  ),
  time_cols           = c("year", "season"),
  xcol                = "x",
  ycol                = "y",
  points_crs          = terra::crs(terra::rast(ref_file)),
  output_dir          = out_dir,
  save_raw            = TRUE,
  save_scaled         = FALSE,
  save_scaling_params = TRUE,
  verbose             = FALSE
)
```
