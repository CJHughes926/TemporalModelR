# Spatiotemporal Rarefaction of Species Occurrence Data

Preprocessing function that rarefies species occurrence data to one
point per raster pixel, optionally accounting for temporal components.
Reduces sampling bias and spatial autocorrelation in occurrence
datasets.

## Usage

``` r
spatiotemporal_rarefaction(points_sp, output_dir, reference_raster,
                           time_cols = NULL, xcol = NULL, ycol = NULL,
                           points_crs = NULL, output_prefix = "Pts_Database",
                           verbose = TRUE)
```

## Arguments

- points_sp:

  Input point data. Accepts an sf object, data frame,
  SpatialPointsDataFrame, or file path to a `.csv`, `.shp`, `.geojson`,
  or `.gpkg` file.

- output_dir:

  Character. Directory where output CSV files will be saved.

- reference_raster:

  Character, SpatRaster, or RasterLayer. Raster used to define pixel
  boundaries for rarefaction. Accepts a file path, RasterLayer, or
  SpatRaster.

- time_cols:

  Character vector. Column names in `points_sp` defining temporal
  grouping for spatiotemporal rarefaction. When `NULL` (default), only
  spatial rarefaction is performed.

- xcol:

  Character. Name of the x-coordinate column. Required when `points_sp`
  is a CSV file or data frame.

- ycol:

  Character. Name of the y-coordinate column. Required when `points_sp`
  is a CSV file or data frame.

- points_crs:

  Character or CRS object. CRS of the input points. Required when
  `points_sp` is a CSV file or data frame.

- output_prefix:

  Character. Prefix for output file names. Default is `"Pts_Database"`.
  Output files are named `<output_prefix>_OnePerPix.csv` and, when
  `time_cols` are provided, `<output_prefix>_OnePerPixPerTimeStep.csv`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes file loading, missing-data removal counts,
  rarefaction summaries, and a final comparison of spatial vs
  spatiotemporal point counts.

## Value

Invisibly returns a list containing:

- `input_points`: Integer. Number of input points after CRS alignment
  and removal of rows with missing time column values.

- `spatial_points`: Integer. Number of points retained after
  spatial-only rarefaction (one per pixel).

- `spatiotemporal_points`: Integer. Number of points retained after
  spatiotemporal rarefaction (one per pixel per time combination). `NA`
  when `time_cols` is not provided.

- `time_cols_used`: Character vector of time columns used. `NULL` when
  `time_cols` is not provided.

- `spatial_table`: Data frame of spatially rarefied points with columns
  `pixel_id`, `X`, `Y`, and any `time_cols`.

- `spatiotemporal_table`: Data frame of spatiotemporally rarefied points
  with columns `pixel_id`, `X`, `Y`, and `time_cols`. `NULL` when
  `time_cols` is not provided.

- `files_created`: Named list of file paths written. Always contains
  `$spatial`; additionally contains `$spatiotemporal` when `time_cols`
  are provided.

## Details

The function assigns each point to a raster pixel using the resolution
and extent of `reference_raster`, then performs:

- Spatial rarefaction: retains one point per pixel, written to
  `<output_prefix>_OnePerPix.csv`.

- Spatiotemporal rarefaction (when `time_cols` are provided): retains
  one point per pixel per unique combination of time column values,
  written to `<output_prefix>_OnePerPixPerTimeStep.csv`.

Output CSV files are suitable as direct input to
[`temporally_explicit_extraction`](temporally_explicit_extraction.md).

## See also

Preprocessing:
[`temporally_explicit_extraction`](temporally_explicit_extraction.md),
[`spatiotemporal_partition`](spatiotemporal_partition.md)

## Examples

``` r
pts_file <- system.file(
  "extdata/points/synthetic_occurrence_points.csv",
  package = "TemporalModelR"
)

ref_file <- system.file("extdata/rasters_raw/elevation.tif",
                        package = "TemporalModelR")

out_dir  <- file.path(tempdir(), "rarefied")

spatiotemporal_rarefaction(
  points_sp        = pts_file,
  output_dir       = out_dir,
  reference_raster = ref_file,
  time_cols        = c("year", "season"),
  xcol             = "x",
  ycol             = "y",
  points_crs       = terra::crs(terra::rast(ref_file)),
  verbose          = FALSE
)
```
