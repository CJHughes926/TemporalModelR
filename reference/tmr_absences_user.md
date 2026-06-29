# Pre-built user_data pseudoabsences (seasonal workflow)

Output of [`generate_absences`](generate_absences.md) run with
`method = "user_data"` against [`tmr_partition`](tmr_partition.md),
using `inst/extdata/points/synthetic_user_presences.csv` as the supplied
absence locations. The CSV contains occurrence-formatted points for a
'related speecies' to that data of interest to be used as a proxy for
survey effort for the species of interest, serving as an alternative
method for generating pseudoabsences.

## Usage

``` r
tmr_absences_user
```

## Format

A list as returned by [`generate_absences`](generate_absences.md),
containing `$pseudoabsences` (an sf object with attached predictor
columns), `$plots`, and `$summary`.
