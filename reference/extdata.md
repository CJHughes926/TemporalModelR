# Bundled rasters, point files, and prediction outputs

Several non-R objects ship in `inst/extdata/` for use in function
examples and the package vignette. They cannot be portably serialized as
`.rda` files, so they are stored as GeoTIFF and CSV files and loaded via
[`system.file`](https://rdrr.io/r/base/system.file.html).

## Details

- `extdata/rasters_raw/`:

  Raw synthetic environmental rasters before alignment. Includes
  `elevation.tif`, `forest_cover_1.tif` through `forest_cover_15.tif`,
  `pr_ann_1.tif` through `pr_ann_15.tif`, and
  `prseas_<year>_<season>.tif` for each of 15 years and 4 seasons (60
  files).

- `extdata/rasters_aligned/`:

  Same rasters after [`raster_align`](raster_align.md) has reprojected
  and masked them to the reference grid.

- `extdata/rasters_scaled/`:

  Aligned rasters z-scored using the scaling parameters from the
  seasonal extraction. Used by examples and modeling functions that work
  with the seasonal predictor set.

- `extdata/rasters_scaled_annual/`:

  Aligned rasters z-scored using the scaling parameters from the annual
  extraction. Used by examples and modeling functions that work with the
  annual predictor set (`pr_ann` instead of `prseas`).

- `extdata/points/synthetic_occurrence_points.csv`:

  The raw synthetic presence dataset: 150 points with `x`, `y`, `year`,
  `season`, and `pres = 1`.

- `extdata/points/synthetic_occurrence_points.shp`:

  Same points as a shapefile.

- `extdata/points/synthetic_user_presences.csv`:

  A second-species presence dataset for demonstrating
  `method = "user_data"` in [`generate_absences`](generate_absences.md).
  Points were derived from buffer- constrained environmental
  pseudoabsences of the primary species and reformatted as presences
  (`pres = 1`). Same column structure as
  `synthetic_occurrence_points.csv` (`x`, `y`, `year`, `season`,
  `pres`).

- `extdata/points/extracted_seasonal_*.csv`:

  Outputs from
  [`temporally_explicit_extraction`](temporally_explicit_extraction.md)
  using the seasonal predictor set: `_Raw_Values.csv`,
  `_Scaled_Values.csv`, and `_Scaling_Parameters.csv`.

- `extdata/points/extracted_annual_*.csv`:

  Same outputs but using the annual predictor set.

- `extdata/predictions/`:

  Per-timestep fold-vote rasters from the seasonal workflow's
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)
  call. Fifteen files (`Prediction_<year>_Spring.tif`) suitable for
  direct input to
  [`summarize_raster_outputs`](summarize_raster_outputs.md).

- `extdata/binary/consensus_stack.tif`:

  Multi-layer GeoTIFF of binary suitable / unsuitable rasters produced
  by [`summarize_raster_outputs`](summarize_raster_outputs.md) with
  `consensus = 3`. Fifteen layers, one per year, ordered 1 through 15.

- `extdata/binary/frequency_raster.tif`:

  Companion single-layer raster giving the proportion of years each
  pixel was classified as suitable.

- `extdata/precomputed/`:

  Precomputed prediction outputs read directly by the modeling vignettes
  (V3a-V3d) so that
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)
  does not need to be rerun at vignette build time. One subdirectory per
  model type (`glm/`, `gam/`, `rf/`, `hv/`), each containing `preds.rds`
  and a `pred_tifs/` folder. `preds.rds` is the list returned by
  [`generate_spatiotemporal_predictions`](generate_spatiotemporal_predictions.md)
  for that model (`$timestep_metrics`, `$overall_summary`,
  `$prediction_files`, and `$model_type`), with `$prediction_files`
  reduced to bare file names rather than absolute build-time paths.
  `pred_tifs/` holds the 60 per-timestep fold-vote rasters
  (`Prediction_<year>_<season>.tif`, 15 years x 4 seasons) projected
  from that model. The vignettes load `preds.rds` and rebuild
  `$prediction_files` from `pred_tifs/` via
  [`system.file`](https://rdrr.io/r/base/system.file.html).

Example load patterns:


      ### Aligned raster directory
      aln_dir <- system.file("extdata/rasters_aligned",
                             package = "TemporalModelR")

      ### Consensus stack (multi-layer)
      binary_stack <- terra::rast(system.file(
        "extdata/binary/consensus_stack.tif", package = "TemporalModelR"
      ))

      ### Frequency raster
      frequency_rast <- terra::rast(system.file(
        "extdata/binary/frequency_raster.tif", package = "TemporalModelR"
      ))

      ### Per-timestep prediction directory
      pred_dir <- system.file("extdata/predictions",
                              package = "TemporalModelR")

      ### Precomputed GLM prediction object and its rasters
      glm_preds <- readRDS(system.file(
        "extdata/precomputed/glm/preds.rds", package = "TemporalModelR"
      ))
      glm_pred_files <- list.files(
        system.file("extdata/precomputed/glm/pred_tifs",
                    package = "TemporalModelR"),
        pattern = "\.tif$", full.names = TRUE
      )
