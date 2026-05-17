# Generate Spatiotemporal Predictions from Temporal Models

Generates temporally explicit habitat suitability predictions by
projecting fitted models onto environmental raster stacks matching each
time period. Accepts output from any of the four TemporalModelR modeling
functions. Produces per-time-step assessment metrics showing how
predictions and test point coverage vary over time.

## Usage

``` r
generate_spatiotemporal_predictions(
  partition_result,
  model_result,
  pseudoabsence_result = NULL,
  raster_dir,
  variable_patterns,
  time_cols,
  time_steps,
  output_dir = file.path(tempdir(), "Predictions"),
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- partition_result:

  List or character. Output from
  [`spatiotemporal_partition`](spatiotemporal_partition.md) or path to
  an `.rds` file containing that output.

- model_result:

  List or character. Output from any of
  [`build_temporal_hv`](build_temporal_hv.md),
  [`build_temporal_glm`](build_temporal_glm.md),
  [`build_temporal_gam`](build_temporal_gam.md), or
  [`build_temporal_rf`](build_temporal_rf.md), or a path to an `.rds`
  file. Model type is detected automatically from the `model_type`
  field.

- pseudoabsence_result:

  List, character, or `NULL`. Optional. Output from
  [`generate_absences`](generate_absences.md) or path to an `.rds` file.
  When supplied for presence/absence models (GLM, GAM, RF), the held-out
  pseudoabsence test points for each fold are filtered to the current
  time step and used alongside presence test points to compute
  per-timestep TN, FP, Specificity, and TSS. These columns are added to
  the timestep metrics table when pseudoabsences are available. Ignored
  for hypervolume models. Default is `NULL`.

- raster_dir:

  Character. Directory containing environmental raster files (`.tif`),
  typically the output of [`raster_align`](raster_align.md) or
  [`scale_rasters`](scale_rasters.md). File names must follow the
  patterns supplied in `variable_patterns`, with any time placeholder
  substituted for the corresponding value from `time_cols`.

- variable_patterns:

  Named character vector mapping clean variable names to raster filename
  patterns. For time-varying variables include the time placeholder in
  the pattern (e.g. `"forest_cover" = "forest_cover_YEAR"`); for static
  variables omit it (e.g. `"elevation" = "elevation"`). Time
  placeholders must match entries in `time_cols`.

- time_cols:

  Character. Name of the column(s) containing year or time step values
  in the occurrence data. Must match `time_cols` used in
  [`spatiotemporal_partition`](spatiotemporal_partition.md) and the time
  placeholders used in `variable_patterns`.

- time_steps:

  Vector, data frame, or matrix of time periods for which to generate
  predictions.

- output_dir:

  Character. Directory to write prediction rasters and the assessment
  metrics CSV. Default is `file.path(tempdir(), "Predictions")`.

- overwrite:

  Logical. If `TRUE`, overwrites existing output files. If `FALSE`
  (default), existing files are skipped.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes per-time-step prediction progress.

## Value

Invisibly returns a list containing:

- `timestep_metrics`: data frame of per-fold per-time-step metrics.
  Columns always present: `Fold`, the time column(s), `Pct_Suitable`,
  `N_Pres`, `TP`, `FN`, `Sensitivity`, `CBP`. When
  `pseudoabsence_result` is supplied for parametric models,
  additionally: `N_Abs`, `TN`, `FP`, `Specificity`, `TSS`. Saved as
  `Timestep_Assessment_Metrics.csv`.

- `overall_summary`: data frame of pooled metrics per fold. Columns:
  `Fold`, `N_Timesteps`, `Mean_Pct_Suitable`, `Total_TP`, `Total_FN`,
  `Overall_Sensitivity`, `Overall_CBP`. When pseudoabsences available:
  additionally `Total_TN`, `Total_FP`, `Overall_Specificity`,
  `Overall_TSS`.

- `prediction_files`: character vector of paths to saved prediction
  rasters.

- `model_type`: character string recording the model type used.

## Details

G-space predictions are produced as rasters for each time-step and fold
of an input model.

Per-timestep metrics are computed per fold per time step by extracting
raster predictions at the held-out presence points that fall within that
time step. `Pct_Suitable` records the proportion of the study area
predicted suitable. `Sensitivity` shows how consistently test points are
captured. `CBP` (cumulative binomial probability) tests whether the
observed number of correctly predicted test points is better than
expected by random placement: under the null, each point has
`Pct_Suitable` probability of falling in suitable area, so
`dbinom(TP, N_Test_Pts, Pct_Suitable)` gives the probability of
observing exactly `TP` correct by chance. Small values indicate
predictions are better than random.

## See also

Preprocessing: [`scale_rasters`](scale_rasters.md),
[`spatiotemporal_partition`](spatiotemporal_partition.md)

Modeling: [`build_temporal_hv`](build_temporal_hv.md),
[`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_rf`](build_temporal_rf.md)

Postprocessing:
[`summarize_raster_outputs`](summarize_raster_outputs.md),
[`plot_model_assessment`](plot_model_assessment.md)

## Examples

``` r
# \donttest{
  data(tmr_partition, package = "TemporalModelR")

  data(tmr_glm,       package = "TemporalModelR")

  data(tmr_absences,  package = "TemporalModelR")

  scl_dir    <- system.file("extdata/rasters_scaled",
                            package = "TemporalModelR")

  time_steps <- expand.grid(
    year             = 1:15,
    season           = "Spring",
    stringsAsFactors = FALSE
  )

  generate_spatiotemporal_predictions(
    partition_result     = tmr_partition,
    model_result         = tmr_glm,
    pseudoabsence_result = tmr_absences,
    raster_dir           = scl_dir,
    variable_patterns    = c(
      "elevation"    = "elevation",
      "forest_cover" = "forest_cover_YEAR",
      "prseas"       = "prseas_YEAR_SEASON"
    ),
    time_cols            = c("year", "season"),
    time_steps           = time_steps,
    output_dir           = tempdir(),
    overwrite            = TRUE,
    verbose              = FALSE
  )
# }
```
