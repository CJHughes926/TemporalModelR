# Pre-built spatiotemporal partition, small version (seasonal workflow)

A subsampled version of [`tmr_partition`](tmr_partition.md) retaining
approximately half the points per fold, for use in examples and tests
where runtime is a concern. Built by sampling `ceiling(n / 2)` rows from
each fold of `tmr_partition$points_sf`. Loaded with
`data(tmr_partition_small)`.

## Usage

``` r
tmr_partition_small
```

## Format

A list with the same structure as [`tmr_partition`](tmr_partition.md),
containing `$folds`, `$points_sf`, `$voronoi_blocks`, `$voronoi_folds`,
`$summary`, and `$plots`.

## See also

[`tmr_partition`](tmr_partition.md)
