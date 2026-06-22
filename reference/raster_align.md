# Align and Standardize Raster Files to a Reference Raster

Preprocessing function that aligns a batch of raster files to a
specified reference raster by performing reprojection, resampling, and
masking operations. Ensures all rasters share identical CRS, resolution,
and spatial extent for later analyses.

## Usage

``` r
raster_align(input_dir, output_dir, reference_raster,
             output_suffix = "_Masked_Updated", pattern = ".*\\.tif$",
             resample_method = "bilinear",
             overwrite = FALSE, verbose = TRUE)
```

## Arguments

- input_dir:

  Character. Directory containing the input raster files.

- output_dir:

  Character. Directory where processed rasters will be saved.

- reference_raster:

  Character, SpatRaster, or RasterLayer. File path or raster object used
  as the alignment reference.

- output_suffix:

  Character. Suffix appended to output filenames. Default is
  `"_Masked_Updated"`. Output files are named
  `<original_name><output_suffix>.tif`.

- pattern:

  Character. Regular expression used to match raster files within
  `input_dir`. Default is `".*\.tif$"`.

- resample_method:

  Character. Resampling method passed to
  [`resample`](https://rspatial.github.io/terra/reference/resample.html).
  Default is `"bilinear"`, which is appropriate for continuous
  variables. Use `"near"` for categorical rasters (e.g. land cover
  classes) to avoid interpolating between class codes. Other accepted
  values include `"cubic"` and `"lanczos"`.

- overwrite:

  Logical. If `TRUE`, overwrites existing output files. If `FALSE`
  (default), existing files are skipped.

- verbose:

  Logical. If `TRUE` (default), prints progress messages during
  processing. Includes file counts, overwrite mode, per-file progress,
  and a completion summary.

## Value

Invisibly returns a list containing:

- `output_files`: Character vector of file paths written to
  `output_dir`.

- `n_processed`: Integer. Number of rasters successfully processed in
  this call (excludes skipped files when `overwrite = FALSE`).

## Details

For each raster in `input_dir` matching `pattern`, the function:

1.  Reprojects to the CRS of the reference raster

2.  Resamples to match the reference resolution using `resample_method`

3.  Masks values outside the reference raster's non-NA extent

4.  Saves the result as a GeoTIFF to `output_dir`

This preprocessing step ensures spatial alignment before applying other
downstream analyses such as
[`temporally_explicit_extraction`](temporally_explicit_extraction.md) or
[`scale_rasters`](scale_rasters.md).

## See also

Preprocessing: [`scale_rasters`](scale_rasters.md),
[`temporally_explicit_extraction`](temporally_explicit_extraction.md)

## Examples

``` r
raw_dir <- system.file("extdata/rasters_raw", package = "TemporalModelR")

ref     <- file.path(raw_dir, "elevation.tif")

out_dir <- file.path(tempdir(), "aligned")

raster_align(
  input_dir        = raw_dir,
  output_dir       = out_dir,
  reference_raster = ref,
  overwrite        = TRUE,
  verbose          = FALSE
)
```
