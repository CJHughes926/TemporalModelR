# Build Temporal Random Forest Models Across Cross-Validation Folds

Modeling function that constructs Random Forest classification models
for each cross-validation fold using presence and pseudoabsence data.
Each model reserves one fold as testing data and uses the remaining
folds as training data. The user specifies predictors as a character
vector. Predicted probabilities of presence are extracted from the
out-of-bag or in-bag vote fractions and thresholded to produce binary
suitability classifications. Variable importance is recorded for each
fold. The returned object follows the same structure as
[`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md), and
[`build_temporal_hv`](build_temporal_hv.md), and is accepted directly by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

## Usage

``` r
build_temporal_rf(
  partition_result,
  pseudoabsence_result,
  model_vars,
  rf_params = list(),
  threshold_method = "tss",
  output_dir = file.path(tempdir(), "RF_Models"),
  create_plot = TRUE,
  plot_palette = "Dark 2",
  overwrite = FALSE,
  time_cols = NULL,
  verbose = TRUE
)
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

- model_vars:

  Character vector. Names of predictor columns to include in the Random
  Forest. All variables must be present as columns in both the presence
  and pseudoabsence data.

- rf_params:

  Named list. Additional arguments passed to
  [`randomForest`](https://rdrr.io/pkg/randomForest/man/randomForest.html),
  such as `ntree` (number of trees, default 500), `mtry` (number of
  variables tried at each split, default
  `floor(sqrt(length(model_vars)))`), and `nodesize` (minimum node size,
  default 1 for classification). Default is an empty list, which uses
  randomForest defaults.

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
  objects and plots. Default is `file.path(tempdir(), "RF_Models")`.

- create_plot:

  Logical. If `TRUE`, generates a per-fold variable importance plot,
  partial dependence curves for each predictor, and a combined ROC curve
  summary. Default is `TRUE`.

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
  messages.

## Value

A list with class `"TemporalRF"` containing:

- `models`: Named list of fitted `randomForest` objects, one per fold.

- `thresholds`: Named numeric vector of probability thresholds used for
  binary classification, one per fold.

- `threshold_method`: Character string recording the thresholding method
  used.

- `model_vars`: Character vector of predictor names used.

- `variable_importance`: Named list of importance data frames, one per
  fold, with mean decrease in accuracy for each predictor.

- `fold_training_data`: Named list of training data frames used to fit
  each fold model, retained for downstream prediction.

- `fold_test_metrics`: Data frame of held-out test fold metrics per
  fold: `Threshold`, `AUC`, `TSS`, `Kappa`, `Sensitivity`, and
  `Specificity`. Also written to `Fold_Test_Metrics.csv` in
  `output_dir`.

- `output_dir`: Path to the output directory.

- `model_type`: Character string `"rf"`, used by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md).

- `plots`: Named list of recorded plot objects when
  `create_plot = TRUE`. Plots can be replayed with
  [`grDevices::replayPlot()`](https://rdrr.io/r/grDevices/recordplot.html).

## Details

Random Forests are fit using
[`randomForest`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
from the randomForest package. The response is treated as a factor
(`0`/`1`) so the model runs in classification mode, which produces class
vote fractions used as predicted probabilities. Importance is computed
with `importance = TRUE` and `type = 1` (mean decrease in accuracy).

Predicted probabilities are the vote fraction for class `1` from
`predict(..., type = "prob")[, "1"]`. These are used for threshold
selection and ROC curve construction.

Diagnostic plots include: a variable importance bar chart (mean decrease
in accuracy across folds), partial dependence curves for each predictor
showing the marginal effect of each variable while averaging over all
others (with rug marks for presences and pseudoabsences), and a combined
ROC curve panel.

The returned object is recognised by
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md),
which uses the `model_type` field to use the correct prediction and
evaluation logic.

## See also

Preprocessing:
[`spatiotemporal_partition`](spatiotemporal_partition.md),
[`generate_absences`](generate_absences.md)

Modeling: [`build_temporal_glm`](build_temporal_glm.md),
[`build_temporal_gam`](build_temporal_gam.md),
[`build_temporal_hv`](build_temporal_hv.md),
[`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)

External:
[`randomForest`](https://rdrr.io/pkg/randomForest/man/randomForest.html)

## Examples

``` r
# \donttest{
  data(tmr_partition, package = "TemporalModelR")

  data(tmr_absences,  package = "TemporalModelR")

  build_temporal_rf(
    partition_result     = tmr_partition,
    pseudoabsence_result = tmr_absences,
    model_vars           = c("elevation", "forest_cover", "prseas"),
    rf_params            = list(ntree = 100),
    threshold_method     = "tss",
    output_dir           = tempdir(),
    create_plot          = FALSE,
    time_cols            = c("year", "season"),
    verbose              = FALSE
  )
# }
```
