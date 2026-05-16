# Visualize Model Assessment Metrics Across Time

Generates diagnostic plots from the per-timestep assessment table
produced by
[`generate_spatiotemporal_predictions`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.md),
optionally overlaying overall reference values from the model result
object.

## Usage

``` r
plot_model_assessment(
  predictions,
  time_column,
  secondary_time_mode = "combine",
  model_result = NULL,
  cbp_threshold = 0.05,
  plot_palette = "Dark 2",
  verbose = TRUE
)
```

## Arguments

- predictions:

  List returned by
  [`generate_spatiotemporal_predictions`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.md),
  or a named list with at least a `timestep_metrics` element. The
  `timestep_metrics` element may also be a path to a
  `Timestep_Assessment_Metrics.csv` file produced by
  [`generate_spatiotemporal_predictions`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.md).

- time_column:

  Character. Name of the primary time column in `timestep_metrics` to
  use as the x axis (e.g. `"year"`). When predictions span multiple time
  columns (e.g. `"year"` and `"season"`), provide all relevant column
  names as a character vector and control how secondary columns are
  handled via `secondary_time_mode`.

- secondary_time_mode:

  Character. How to handle secondary time columns when `time_column` has
  length \> 1. One of:

  - `"combine"` (default): secondary time values are appended to the
    primary value to form a single ordered x-axis label (e.g.
    `1_Spring`, `1_Summer`, `2_Spring`, ...).

  - `"facet"`: a separate plot is produced for each unique combination
    of secondary time values, with the primary time column as the x axis
    on every panel.

- model_result:

  List or character. Optional. Output from a `build_temporal_*()`
  function or path to its `.rds` file. When supplied, overall
  sensitivity and specificity from `model_result$fold_test_metrics` are
  added as per-fold reference lines. Default is `NULL`.

- cbp_threshold:

  Numeric. Significance threshold for CBP. Default is `0.05`.

- plot_palette:

  Character. Name of an HCL or RColorBrewer palette used to color folds
  in diagnostic plots. Accepts any HCL palette name (see
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html)) or, if
  RColorBrewer is installed, any Brewer palette name. Default is
  `"Dark 2"`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing.

## Value

Invisibly returns a named list containing:

- `pct_suitable`: Recorded plot of proportion of study area predicted
  suitable per time step.

- `sensitivity`: Recorded plot of per-timestep sensitivity.

- `specificity`: Recorded plot of per-timestep specificity (only present
  when pseudoabsence data were used).

- `cbp`: Recorded plot of cumulative binomial probability per time step
  on a log scale.

- `tp_fn`: Recorded plot of true positives and false negatives per time
  step.

- `tn_fp`: Recorded plot of true negatives and false positives per time
  step (only present when pseudoabsence data were used).

- `timestep_summary`: Data frame of per-time-step cross-fold mean and SD
  for each metric.

- `overall_summary`: Data frame from `predictions$overall_summary`, when
  present.

## Details

Plots per-fold and per-timestep diagnostic plots for data produced by
[`generate_spatiotemporal_predictions`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.md).
These quick visuals can be used by users to assess model performance and
significance and decide if the model's performance warrants further
interpretation of the results through postprocessing analyses.

## See also

Preprocessing:
[`spatiotemporal_partition`](https://cjhughes926.github.io/TemporalModelR/reference/spatiotemporal_partition.md),
[`generate_absences`](https://cjhughes926.github.io/TemporalModelR/reference/generate_absences.md)

Modeling:
[`build_temporal_glm`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_glm.md),
[`build_temporal_gam`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_gam.md),
[`build_temporal_rf`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_rf.md),
[`build_temporal_hv`](https://cjhughes926.github.io/TemporalModelR/reference/build_temporal_hv.md),

Postprocessing:
[`generate_spatiotemporal_predictions`](https://cjhughes926.github.io/TemporalModelR/reference/generate_spatiotemporal_predictions.md),
[`summarize_raster_outputs`](https://cjhughes926.github.io/TemporalModelR/reference/summarize_raster_outputs.md)

## Examples

``` r
# \donttest{
  data(tmr_predictions, package = "TemporalModelR")

  plot_model_assessment(
    predictions         = tmr_predictions,
    time_column         = c("year", "season"),
    secondary_time_mode = "combine",
    verbose             = FALSE
  )






# }
```
