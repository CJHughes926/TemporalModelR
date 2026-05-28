# Summarize Prediction Rasters into Consensus Outputs

Post-processing function that synthesizes per-time-step fold-vote
rasters (output from
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md))
into binary consensus predictions and a temporal frequency summary. Each
input raster contains integer vote counts per pixel the number of
cross-validation folds that classified that pixel as suitable. The
`consensus` threshold controls how many folds must agree for a pixel to
be classified as suitable in the output binary rasters.

## Usage

``` r
summarize_raster_outputs(predictions_dir, output_dir = NULL,
                         consensus = 1, file_pattern = "Prediction_.*\\.tif$",
                         overwrite = FALSE, verbose = TRUE)
```

## Arguments

- predictions_dir:

  Character. Directory containing prediction raster files, typically the
  `output_dir` used in
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- output_dir:

  Character. Directory for output files. Defaults to `predictions_dir`
  if `NULL`.

- consensus:

  Integer. Minimum number of folds that must agree on suitability for a
  pixel to be classified as suitable in the binary output. For example,
  `consensus = 1` marks any pixel suitable if at least one fold predicts
  it suitable; `consensus = 4` requires at least four folds to agree.
  Must be between 1 and the maximum vote count in the rasters (number of
  folds). Default is `1`.

- file_pattern:

  Character. Regular expression to match prediction raster files.
  Default is `"Prediction_.*\.tif$"`.

- overwrite:

  Logical. If `TRUE`, overwrites existing output files. If `FALSE`
  (default), existing files are skipped.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing.

## Value

Invisibly returns a list containing:

- `consensus_stack`: `SpatRaster` stack of per-time-step binary
  consensus rasters (one layer per input file).

- `frequency_raster`: `SpatRaster` showing the proportion of time steps
  during which each pixel met the consensus threshold. Values range from
  0 (never suitable) to 1 (always suitable).

- `consensus`: the consensus threshold used.

- `n_timesteps`: number of time steps processed.

- `consensus_dir`: path to the directory containing per-time-step binary
  consensus rasters.

- `frequency_file`: path to the written frequency raster file.

## Details

Input rasters are fold-vote-count rasters produced by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
where each pixel value is the number of cross-validation fold models
that predicted suitability at that location for that time step. The
`consensus` threshold is applied as:
`binary = as.integer(vote_count >= consensus)`.

The frequency raster (proportion of time steps suitable under the chosen
consensus) serves as input to
[`analyze_temporal_patterns`](analyze_temporal_patterns.md) for
identifying long-term trends in suitability.

## See also

Upstream:
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)

Downstream: [`analyze_temporal_patterns`](analyze_temporal_patterns.md)

## Examples

``` r
# \donttest{
  pred_dir <- system.file("extdata/predictions",
                          package = "TemporalModelR")

  summarize_raster_outputs(
    predictions_dir = pred_dir,
    output_dir      = tempdir(),
    consensus       = 3,
    overwrite       = TRUE,
    verbose         = FALSE
  )
# }
```
