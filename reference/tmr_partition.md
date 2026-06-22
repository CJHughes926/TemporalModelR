# Pre-built spatiotemporal partition (seasonal workflow)

Output from [`spatiotemporal_partition`](spatiotemporal_partition.md)
run on the synthetic occurrence dataset bundled in `inst/extdata/`.
Built with `n_spatial_folds = 2` and `n_temporal_folds = 2`, producing
four cross-validation folds (2 spatial + 2 temporal). Points span 15
years and the seasons Spring, Summer, and Autumn on a 15 x 30 cell study
area (3000 m x 1500 m, 100 m pixels) in a custom synthetic local CRS.
Predictors attached to each point are `forest_cover`, `prseas` (seasonal
precipitation), and `elevation`, all z-scored. Loaded with
`data(tmr_partition)`.

## Usage

``` r
tmr_partition
```

## Format

A list as returned by
[`spatiotemporal_partition`](spatiotemporal_partition.md), containing
`$folds`, `$points_sf`, `$voronoi_blocks`, `$voronoi_folds`, `$summary`,
and `$plots`.
