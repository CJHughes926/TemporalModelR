# Build Temporal GLM Models Across Cross-Validation Folds

Modeling function that constructs binomial generalized linear models
(GLMs) for each cross-validation fold using presence and pseudoabsence
data. Each model reserves one fold as testing data and uses the
remaining folds as training data. The user supplies the model formula
directly, giving full control over predictor terms, polynomials, and
interactions. The link function can be set to logit, probit,
complementary log-log, or cauchit. Supports automatic or manual
probability thresholding for converting continuous predictions to binary
suitability classifications necessary for downstream analyses. The
returned object follows the same structure as
[`build_temporal_hv`](build_temporal_hv.md),
[`build_temporal_gam`](build_temporal_gam.md), and
[`build_temporal_rf`](build_temporal_rf.md), and is accepted directly by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

## Usage

``` r
build_temporal_glm(partition_result, pseudoabsence_result, model_formula,
          link = "logit", threshold_method = "tss",
          output_dir = file.path(tempdir(), "GLM_Models"),
          create_plot = TRUE, plot_palette = "Dark 2", overwrite = FALSE,
          time_cols = NULL, verbose = TRUE)
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

  - `~ Var1 + Var2 + I(Var1^2)`

  - `"~ Var1 + Var2 + I(Var1^2)"`

  Standard R formula syntax applies: `+` for additive terms, `*` for
  main effects plus interaction, `:` for interaction only,
  [`I()`](https://rdrr.io/r/base/AsIs.html) for arithmetic
  transformations, [`poly()`](https://rdrr.io/r/stats/poly.html) for
  orthogonal polynomials, [`log()`](https://rdrr.io/r/base/Log.html),
  [`sqrt()`](https://rdrr.io/r/base/MathFun.html), and any other base R
  function that can appear in a formula. All predictor names referenced
  in the formula must be present as columns in both the presence and
  pseudoabsence data.

- link:

  Character. The link function for the binomial GLM. One of `"logit"`
  (default), `"probit"`, `"cloglog"`, or `"cauchit"`. See
  [`binomial`](https://rdrr.io/r/stats/family.html) for details on each
  link function.

- threshold_method:

  Character or numeric. Method used to convert continuous predicted
  probabilities to binary suitability. Accepted values:

  - `"prevalence"`: Sets threshold equal to the prevalence (proportion
    of presences) in the training data for that fold.

  - `"tss"`: Selects the threshold that maximises the True Skill
    Statistic (sensitivity + specificity - 1) on the training data.
    Default.

  - A numeric value between 0 and 1 (e.g. `0.4`): Uses that value as a
    fixed threshold for all folds directly.

- output_dir:

  Character. Directory to write output files including saved model
  objects and plots. Default is `file.path(tempdir(), "GLM_Models")`.

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

A list with class `"TemporalGLM"` containing:

- `models`: Named list of fitted `glm` objects, one per fold.

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

- `model_type`: Character string `"glm"`, used by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- `plots`: Named list of recorded plot objects when
  `create_plot = TRUE`. Plots can be replayed with
  [`grDevices::replayPlot()`](https://rdrr.io/r/grDevices/recordplot.html).

## Details

The `model_formula` argument accepts any standard R formula right-hand
side. The response (`presence`) is prepended automatically. All R
formula operators are valid, including
[`I()`](https://rdrr.io/r/base/AsIs.html),
[`poly()`](https://rdrr.io/r/stats/poly.html),
[`log()`](https://rdrr.io/r/base/Log.html),
[`sqrt()`](https://rdrr.io/r/base/MathFun.html), `:`, and `*`. Variable
names must match column names in the data exactly. Predictor names for
response curve plots are extracted via
[`all.vars()`](https://rdrr.io/r/base/allnames.html), which correctly
unwraps terms such as `I(Var1^2)` to the base variable `Var1`.

All models are fit as `stats::glm(..., family = binomial(link = link))`.
Predicted values are always probabilities on the 0-1 scale.

The returned object is recognised by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
which uses the `model_type` field to use the correct prediction and
evaluation logic.

## See also

Preprocessing:
[`spatiotemporal_partition`](spatiotemporal_partition.md),
[`generate_absences`](generate_absences.md)

Modeling: [`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_rf`](build_temporal_rf.md),
[`build_temporal_hv`](build_temporal_hv.md),
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)

## Examples

``` r
# \donttest{
  data(tmr_partition, package = "TemporalModelR")

  data(tmr_absences,  package = "TemporalModelR")

  build_temporal_glm(
    partition_result     = tmr_partition,
    pseudoabsence_result = tmr_absences,
    model_formula        = ~ elevation + forest_cover + prseas,
    threshold_method     = "tss",
    output_dir           = tempdir(),
    create_plot          = FALSE,
    time_cols            = c("year", "season"),
    verbose              = FALSE
  )
# }
```
