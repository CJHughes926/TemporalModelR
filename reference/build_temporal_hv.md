# Build Temporal Hypervolume Models Across Cross-Validation Folds

Modeling function that constructs hypervolume models for each
cross-validation fold using either Gaussian kernel density estimation or
one-class SVM. Each hypervolume reserves one fold as testing data and
uses the remaining folds as training data. The returned object follows
the same structure as [`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md), and
[`build_temporal_rf`](build_temporal_rf.md), and is accepted directly by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

## Usage

``` r
build_temporal_hv(partition_result, model_vars, method,
                  hypervolume_params = list(),
                  output_dir = file.path(tempdir(), "Hypervolume_Models"),
                  create_plot = TRUE, overwrite = FALSE,
                  plot_palette = "Dark 2", verbose = TRUE)
```

## Arguments

- partition_result:

  List or character. Output from
  [`spatiotemporal_partition`](spatiotemporal_partition.md) or path to
  an `.rds` file containing that output.

- model_vars:

  Character vector. Names of predictor columns to use in hypervolume
  construction. All variables must be present as columns in the
  occurrence data produced by
  [`temporally_explicit_extraction`](temporally_explicit_extraction.md).

- method:

  Character. Hypervolume method. One of `"gaussian"` for Gaussian kernel
  density estimation or `"svm"` for one-class support vector machine.
  Required.

- hypervolume_params:

  Named list. Additional parameters passed to
  [`hypervolume_gaussian`](https://rdrr.io/pkg/hypervolume/man/hypervolume_gaussian.html)
  or
  [`hypervolume_svm`](https://rdrr.io/pkg/hypervolume/man/hypervolume_svm.html).
  For Gaussian models, valid keys are `kde.bandwidth`,
  `quantile.requested`, `quantile.requested.type`, `chunk.size`,
  `verbose`, and `samples.per.point`. For SVM models, valid keys are
  `svm.nu`, `svm.gamma`, `chunk.size`, `verbose`, and
  `samples.per.point`. Default is an empty list, which uses the built-in
  defaults.

- output_dir:

  Character. Directory to write output files including saved hypervolume
  objects and plots. Default is
  `file.path(tempdir(), "Hypervolume_Models")`.

- create_plot:

  Logical. If `TRUE`, generates pairplot visualisations of each fold's
  hypervolume and a combined comparison plot. Default is `TRUE`.

- overwrite:

  Logical. If `TRUE`, overwrites existing saved hypervolume files. If
  `FALSE`, loads existing files when available. Default is `FALSE`.

- plot_palette:

  Character. Name of an HCL or RColorBrewer palette used to color folds
  in diagnostic plots. Accepts any HCL palette name (see
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html)) or, if
  RColorBrewer is installed, any Brewer palette name. Default is
  `"Dark 2"`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes per-fold hypervolume construction progress,
  volume summaries, and overlap statistics.

## Value

A list with class `"TemporalHypervolume"` containing:

- `models`: Named list of fitted `Hypervolume` objects, one per fold,
  named `fold1`, `fold2`, etc. This naming convention matches the other
  `build_temporal_*()` functions and is required by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- `volumes`: Named numeric vector of hypervolume sizes (units depend on
  the number of dimensions and bandwidth).

- `overlaps`: Named list of pairwise percent volume overlaps between all
  fold combinations.

- `method`: Character string recording the method used.

- `model_vars`: Character vector of predictor names used.

- `fold_training_data`: Named list of training data frames used to fit
  each fold model, retained for consistency with the other model types
  and for downstream use.

- `fold_test_metrics`: Data frame of E-space inclusion metrics computed
  on the held-out test points for each fold. Columns: `fold`, `n_test`,
  `volume`, `tp`, `fn`, `sensitivity`, `Sensitivity`. Printed as a
  summary table at the end of model building.

- `output_dir`: Path to the output directory.

- `model_type`: Character string `"hypervolume"`, used by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- `plots`: List of recorded plot objects (if `create_plot = TRUE`).
  Plots can be replayed with
  [`grDevices::replayPlot()`](https://rdrr.io/r/grDevices/recordplot.html).

## Details

For N folds, constructs N hypervolumes where each fold reserves one
group of points for testing and uses the remaining N-1 groups for
training. Pairwise overlap statistics quantify similarity between
hypervolumes across folds, with low overlap indicating that the
environmental space sampled differs substantially across folds.

Hypervolume objects are saved as a combined `.rds` file in `output_dir`.
If the file already exists and `overwrite = FALSE`, the saved file is
loaded rather than re-fitting.

The returned object is accepted directly by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
which uses the `model_type` field to dispatch hypervolume-specific
projection logic via
[`hypervolume_project`](https://rdrr.io/pkg/hypervolume/man/hypervolume_project.html).

## See also

Preprocessing: [`spatiotemporal_partition`](spatiotemporal_partition.md)

Modeling: [`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_rf`](build_temporal_rf.md),
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)

External:
[`hypervolume_gaussian`](https://rdrr.io/pkg/hypervolume/man/hypervolume_gaussian.html),
[`hypervolume_svm`](https://rdrr.io/pkg/hypervolume/man/hypervolume_svm.html)

## Examples

``` r
data(tmr_partition_small, package = "TemporalModelR")

build_temporal_hv(
  partition_result = tmr_partition_small,
  model_vars       = c("elevation", "forest_cover"),
  method           = "svm",
  output_dir       = tempdir(),
  create_plot      = FALSE,
  verbose          = FALSE
)
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> Retaining 4928/4928 hypervolume random points for comparison with 19 test points.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> Retaining 4928/4928 hypervolume random points for comparison with 19 test points.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> Retaining 6136/6136 hypervolume random points for comparison with 22 test points.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
#> Retaining 4763/4763 hypervolume random points for comparison with 17 test points.
#> 
#> Building tree... 
#> done.
#> Ball query... 
#> 
#> done.
```
