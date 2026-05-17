# Spatiotemporal Cross-Validation Partitioning

Preprocesses species occurrence data by partitioning it into spatially
and temporally structured folds for cross-validation. Supports creation
of spatial-only folds, temporal-only folds, and random folds.

## Usage

``` r
spatiotemporal_partition(
  reference_shapefile_path,
  points_file_path,
  time_cols = NULL,
  xcol = NULL,
  ycol = NULL,
  points_crs = NULL,
  n_spatial_folds = 0,
  n_temporal_folds = 0,
  n_balanced_folds = 0,
  n_random_folds = 0,
  single_fold = FALSE,
  max_imbalance = 0.05,
  max_attempts = 10,
  create_plot = TRUE,
  plot_palette = "Dark 2",
  output_file = NULL,
  verbose = TRUE
)
```

## Arguments

- reference_shapefile_path:

  Character or sf object. Path to a polygon file or an `sf` polygon
  object defining the study area.

- points_file_path:

  Character, sf object, sfc object, Spatial object, or data frame. Path
  to occurrence data (`.csv`, `.shp`, `.geojson`, `.gpkg`) or a spatial
  object.

- time_cols:

  Character. Name of a single column containing temporal values (e.g.
  year). Used to define temporal blocks. Required when using temporal
  folds. Must be a single column name; does not support more than one
  time column unlike other functions in this package. Compound time
  representations (e.g. year + season) should be encoded into a single
  ordered numeric column before partitioning, or only one (e.g. year)
  should be used.

- xcol:

  Character. Name of the x-coordinate column. Required when
  `points_file_path` is a CSV file or data frame.

- ycol:

  Character. Name of the y-coordinate column. Required when
  `points_file_path` is a CSV file or data frame.

- points_crs:

  Character or CRS object. CRS of the input points. Required when
  `points_file_path` is a CSV file or data frame.

- n_spatial_folds:

  Integer. Number of spatially explicit folds. Ignored when using random
  folds. Default is `0`.

- n_temporal_folds:

  Integer. Number of temporally explicit folds. When used alone (with
  `n_spatial_folds = 0`), creates temporal-only folds where each fold
  spans the full study area but covers a distinct slice of the time
  series. When combined with `n_spatial_folds`, creates a spatiotemporal
  design. Ignored when using random folds. Default is `0`.

- n_balanced_folds:

  Integer. Reserved for future use. Default is `0` (disabled).

- n_random_folds:

  Integer. Number of random folds with no spatial or temporal structure.
  Overrides all other fold parameters. Default is `0`.

- single_fold:

  Logical. If `TRUE`, bypasses all partitioning and assigns all points
  to a single fold (fold 1). In this mode all points are used for both
  training and testing, producing a single model trained on the full
  dataset. All downstream functions accept the result identically to a
  standard multi-fold partition. Overrides all fold count parameters.
  Default is `FALSE`.

- max_imbalance:

  Numeric. Maximum allowed fold size imbalance as a proportion between 0
  and 1. Default is `0.05`.

- max_attempts:

  Integer. Maximum number of partitioning attempts for spatiotemporal
  and balanced modes. Each attempt re-runs the spatial block
  construction; the attempt with the lowest imbalance is returned.
  Ignored for random and spatial-only modes. Default is `10`.

- create_plot:

  Logical. If `TRUE` (default), generates diagnostic plots showing fold
  distributions.

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
  processing. Includes the partition mode, fold structure, per-fold
  point counts, and file-save confirmation.

## Value

Invisibly returns a list containing:

- `folds`: Data frame of fold assignments with a `fold` column
  identifying each point's cross-validation fold.

- `points_sf`: sf object of occurrence points with assigned folds.

- `voronoi_folds`: sf object of Voronoi polygons representing the
  spatial fold boundaries. `NULL` for random folds, temporal-only folds,
  and single-fold mode.

- `summary`: Data frame of partitioning summary statistics.

- `plots`: Named list of recorded plot objects when
  `create_plot = TRUE`. Empty list in single-fold mode.

## Details

Works better with smaller numbers of folds and may have difficulties
creating even folds for large numbers of groups or where sample sizes
are very small.

The function partitions data into folds using one of five modes:

- **Single fold**: All points are assigned to fold 1 and used for both
  training and testing. This produces a single model trained on the full
  dataset with no held-out validation. Useful when sample sizes are too
  small for cross-validation, or as a final production model step after
  cross-validation has already established model quality. Set
  `single_fold = TRUE`. All downstream functions accept the result
  identically to standard multi-fold output.

- **Random**: Points are assigned to folds by random shuffling with no
  spatial or temporal structure. Each fold is a simple random sample of
  the full dataset, intended as a naive baseline that makes no attempt
  to reduce spatial or temporal autocorrelation between training and
  test sets. Use `n_random_folds`.

- **Spatial-only**: The study area is divided into \\k\\ contiguous
  spatial regions using a recursive k-d tree bisection algorithm. At
  each step the point set is split along its longest spatial axis,
  recursively halving until the target number of folds is reached. A
  centroid reassignment pass then refines boundaries to improve balance.
  Each region becomes one fold, so training always occurs on data from
  geographically distinct areas relative to the test fold. No temporal
  separation is imposed, meaning that points from any time period may
  appear in any fold. Use `n_spatial_folds` alone.

- **Temporal-only**: Each fold covers the full spatial extent of the
  study area but is restricted to a distinct, non-overlapping slice of
  the time series. The global time series is divided into
  `n_temporal_folds` equal intervals using quantile-based breaks, and
  all points within each interval form one fold. This design tests model
  transferability across time while retaining full spatial coverage in
  every fold. Use `n_temporal_folds` alone (with `n_spatial_folds = 0`).
  Requires `time_cols`.

- **Spatiotemporal**: Folds are assigned using the same recursive k-d
  tree bisection as spatial-only mode, operating on the full point set
  to produce spatially contiguous groups. The resulting groups are then
  split into a spatial pool (`n_spatial_folds` folds drawn from
  geographically distinct regions) and a temporal pool
  (`n_temporal_folds` folds each restricted to a distinct slice of the
  time series but spanning the full study area). Together the two pools
  assess both geographic and temporal transferability in a single
  cross-validation design. Use `n_spatial_folds` and `n_temporal_folds`
  together. Requires `time_cols`.

Fold assignment uses a recursive k-d tree bisection algorithm that
splits points along their longest spatial axis at each step, followed by
a centroid reassignment pass to improve boundary regularity and
point-count balance. Voronoi tessellation on fold centroids is used only
for visualisation of the resulting spatial boundaries. For temporal
mode, temporal blocks are defined by dividing the global time series
into equal intervals using quantile-based breaks. For spatiotemporal
mode, the typical spatial assignment is done, but with one larger
spatial block made with enough points to represent all of the temporal
folds, then the temporal blocking is applied to those points.

Partitioned datasets are suitable for cross-validation in modeling
workflows, ensuring spatial and/or temporal independence between folds.

## See also

Preprocessing:
[`spatiotemporal_rarefaction`](spatiotemporal_rarefaction.md),
[`temporally_explicit_extraction`](temporally_explicit_extraction.md),
[`generate_absences`](generate_absences.md)

Modeling: [`build_temporal_hv`](build_temporal_hv.md),
[`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_rf`](build_temporal_rf.md)

## Examples

``` r
# \donttest{
  pts_file <- system.file(
    "extdata/points/extracted_seasonal_Scaled_Values.csv",
    package = "TemporalModelR"
  )

  ref_file <- system.file("extdata/rasters_raw/elevation.tif",
                          package = "TemporalModelR")

  study_crs <- sf::st_crs(terra::rast(ref_file))

  study_area_sf <- sf::st_as_sf(sf::st_as_sfc(
    sf::st_bbox(c(xmin = 0, xmax = 3000, ymin = 0, ymax = 1500),
                crs = study_crs)
  ))

  spatiotemporal_partition(
    reference_shapefile_path = study_area_sf,
    points_file_path         = pts_file,
    xcol                     = "x",
    ycol                     = "y",
    points_crs               = study_crs,
    time_cols                = "year",
    n_spatial_folds          = 2,
    n_temporal_folds         = 2,
    create_plot              = FALSE,
    verbose                  = FALSE
  )
#> Warning: Could not achieve target balance within 10 attempts. Final imbalance: 14.67%. Returning best result achieved. Try increasing max_imbalance or adjusting the fold configuration.
# }
```
