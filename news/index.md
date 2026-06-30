# Changelog

## TemporalModelR 0.2.0

CRAN release: 2026-06-30

- Initial CRAN submission.

## TemporalModelR 0.3.0

- Fixed bug in plot_model_assessment which caused some plots to not
  correctly render in markdown documents

- Added functionality to generate_absences which now allowes for users
  to define their own absence data, with this function now able to
  format user defined absence data to work with all downstream
  operations

- Fixed bug with spatiotemporal_partition which was outputting the
  incorrect bounderies as part of the spatial partition output. Fixed
  this output to make it more compatible with partitioning user defined
  absences after the fact to corrospond to updates to generate_absences
