# Pre-built spatiotemporal predictions (seasonal workflow)

Output from
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)
projecting `tmr_glm` across 15 years for the Spring season only (one
prediction layer per year-season combination, 15 total). Loaded with
`data(tmr_predictions)`.

## Usage

``` r
tmr_predictions
```

## Format

A list as returned by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
containing `$timestep_metrics`, `$overall_summary`, `$prediction_files`,
and `$model_type`.

## Details

Note that `$prediction_files` stores the absolute paths used at build
time. To work with the per-timestep rasters on a user machine,
regenerate them via
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)
or use the bundled prediction set in `inst/extdata/predictions/`.
