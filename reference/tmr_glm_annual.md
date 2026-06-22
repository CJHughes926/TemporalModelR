# Pre-built GLM result (annual workflow)

Annual variant of `tmr_glm`. Fit with the formula
`~ forest_cover + pr_ann + elevation`, a logit link, and TSS-based
threshold selection. One model per fold. Loaded with
`data(tmr_glm_annual)`.

## Usage

``` r
tmr_glm_annual
```

## Format

A list of class `"TemporalGLM"`.
