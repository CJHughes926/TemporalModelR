# Pre-built GLM result (seasonal workflow)

Output from [`build_temporal_glm`](build_temporal_glm.md) fit to
`tmr_partition` and `tmr_absences` with the formula
`~ forest_cover + prseas + elevation`, a logit link, and TSS-based
threshold selection. One model per fold. Loaded with `data(tmr_glm)`.

## Usage

``` r
tmr_glm
```

## Format

A list of class `"TemporalGLM"` as returned by
[`build_temporal_glm`](build_temporal_glm.md), containing `$models`,
`$thresholds`, `$threshold_method`, `$model_formula`, `$link`,
`$model_vars`, `$fold_training_data`, `$fold_test_metrics`,
`$output_dir`, `$model_type`, and `$plots`.
