# Build Temporal GAM Models Across Cross-Validation Folds

Modeling function that constructs binomial generalized additive models
(GAMs) for each cross-validation fold using presence and pseudoabsence
data. Each model reserves one fold as testing data and uses the
remaining folds as training data. The user supplies the model formula
directly using standard mgcv formula syntax, including smooth terms such
as `s()`, `te()`, and `ti()`. Supports automatic or manual probability
thresholding for converting continuous predictions to binary suitability
classifications necessary for downstream analyses. The returned object
follows the same structure as
[`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_hv`](build_temporal_hv.md), and
[`build_temporal_rf`](build_temporal_rf.md), and is accepted directly by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

## Usage

``` r
build_temporal_gam(partition_result, pseudoabsence_result, model_formula,
                   link = "logit", gam_params = list(method = "REML"),
                   threshold_method  = "tss",
                   output_dir = file.path(tempdir(), "GAM_Models"),
                   create_plot = TRUE, plot_palette = "Dark 2",
                   overwrite = FALSE, time_cols = NULL, verbose = TRUE)
```

## Arguments

- partition_result:

  List or character. Output from
  [`spatiotemporal_partition`](spatiotemporal_partition.md) or path to
  an `.rds` file containing that output.

- pseudoabsence_result:

  List or character. Output from
  [`generate_absences`](generate_absences.md) or path to an `.rds` file
  containing that output.

- model_formula:

  Formula or character. The right-hand side of the model formula
  supplied as either a formula object or a character string. The
  response variable (`presence`) is always added automatically on the
  left-hand side, so only the right-hand side needs to be provided. Both
  of the following are accepted and equivalent:

  - `~ s(Var1) + s(Var2) + Var3`

  - `"~ s(Var1) + s(Var2) + Var3"`

  Standard mgcv formula syntax applies. Smooth terms are specified with
  `s()` for univariate smooths, `te()` for tensor product smooths of two
  or more variables, and `ti()` for tensor product interaction terms.
  Parametric terms can be included alongside smooth terms using `+`. The
  basis type and dimension can be controlled via arguments to `s()`,
  e.g. `s(Var1, k = 5, bs = "tp")`. All predictor names referenced in
  the formula must be present as columns in both the presence and
  pseudoabsence data.

- link:

  Character. The link function for the binomial GAM. One of `"logit"`
  (default), `"probit"`, `"cloglog"`, or `"cauchit"`. See
  [`binomial`](https://rdrr.io/r/stats/family.html) for details on each
  link function.

- gam_params:

  Named list. Additional arguments passed to
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html), such as `method` for
  the smoothing parameter estimation method (e.g. `"REML"`) or `select`
  for additional shrinkage. Default is `list(method = "REML")`.

- threshold_method:

  Character or numeric. Method used to convert continuous predicted
  probabilities to binary suitability. Accepted values:

  - `"prevalence"`: Sets threshold equal to the prevalence (proportion
    of presences) in the training data for that fold.

  - `"tss"`: Selects the threshold that maximizes the True Skill
    Statistic (sensitivity + specificity - 1) on the training data.
    Default.

  - A numeric value between 0 and 1 (e.g. `0.4`): Uses that value as a
    fixed threshold for all folds directly.

- output_dir:

  Character. Directory to write output files including saved model
  objects and plots. Default is `file.path(tempdir(), "GAM_Models")`.

- create_plot:

  Logical. If `TRUE`, generates per-fold response curve plots and a
  combined ROC curve summary. Default is `TRUE`.

- plot_palette:

  Character. Name of an HCL or RColorBrewer palette used to color folds
  in diagnostic plots. Accepts any HCL palette name (see
  [`hcl.pals`](https://rdrr.io/r/grDevices/palettes.html)) or, if
  RColorBrewer is installed, any Brewer palette name. Default is
  `"Dark 2"`.

- overwrite:

  Logical. If `TRUE`, overwrites existing saved model files. If `FALSE`,
  loads existing files when available. Default is `FALSE`.

- time_cols:

  Character. Name of the column(s) containing year or time step values
  in the occurrence data. Must match `time_cols` used in
  [`spatiotemporal_partition`](spatiotemporal_partition.md). Default is
  `NULL`.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes per-fold training summaries and file-saved
  messages. The completion summary and metrics table are always printed
  regardless of this setting.

## Value

A list with class `"TemporalGAM"` containing:

- `models`: Named list of fitted `gam` objects, one per fold.

- `thresholds`: Named numeric vector of probability thresholds used for
  binary classification, one per fold.

- `threshold_method`: Character string recording the thresholding method
  used.

- `model_formula`: The formula object as passed to the fitting function.

- `link`: Character string recording the link function used.

- `model_vars`: Character vector of predictor names extracted from the
  formula right-hand side.

- `fold_training_data`: Named list of training data frames used to fit
  each fold model, retained for downstream prediction.

- `fold_test_metrics`: Data frame of held-out test fold metrics per
  fold: `Threshold`, `AUC`, `TSS`, `Kappa`, `Sensitivity`, and
  `Specificity`. Also written to `Fold_Test_Metrics.csv` in
  `output_dir`.

- `output_dir`: Path to the output directory.

- `model_type`: Character string `"gam"`, used by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- `plots`: Named list of recorded plot objects when
  `create_plot = TRUE`. Plots can be replayed with
  [`grDevices::replayPlot()`](https://rdrr.io/r/grDevices/recordplot.html).

## Details

GAMs are fit using [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html) from
the mgcv package with `family = binomial(link = link)`. Smooth terms
default to thin plate regression splines (`bs = "tp"`) with the basis
dimension `k` chosen automatically by mgcv unless specified in the
formula. Smoothing parameters are estimated by REML by default.

The returned object is recognized by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
which uses the `model_type` field to use the correct prediction and
evaluation logic.

## See also

Preprocessing:
[`spatiotemporal_partition`](spatiotemporal_partition.md),
[`generate_absences`](generate_absences.md)

Modeling: [`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_rf`](build_temporal_rf.md),
[`build_temporal_hv`](build_temporal_hv.md),
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)

External: [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html),
[`s`](https://rdrr.io/pkg/mgcv/man/s.html)

## Examples

``` r
data(tmr_partition, package = "TemporalModelR")

data(tmr_absences,  package = "TemporalModelR")

build_temporal_gam(
  partition_result     = tmr_partition,
  pseudoabsence_result = tmr_absences,
  model_formula        = ~ s(elevation) + s(forest_cover) + s(prseas),
  threshold_method     = "tss",
  output_dir           = tempdir(),
  create_plot          = FALSE,
  time_cols            = c("year", "season"),
  verbose              = FALSE
)
```
