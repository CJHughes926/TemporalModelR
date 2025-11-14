
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TemporalModelR

<!-- badges: start -->
<!-- badges: end -->

TemporalModelR is an R package for building temporally explicit species
distribution models using hypervolume-based methods. The package
provides a complete workflow from data preprocessing through model
building, prediction, and temporal pattern analysis.

## Overview

TemporalModelR enables researchers to:

- **Preprocess spatial and temporal occurrence data** with
  spatiotemporal rarefication
- **Extract and scale environmental variables** matched to temporal
  occurrence records
- **Build hypervolume models** using Gaussian kernel density or
  one-class SVM methods
- **Generate spatiotemporal predictions** with comprehensive model
  evaluation metrics
- **Analyze temporal patterns** in habitat suitability using changepoint
  detection
- **Summarize trends by spatial units** for regional conservation
  assessments

## Installation

You can install the development version of TemporalModelR from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("CJHughes926/TemporalModelR")
```

## Workflow Overview

The TemporalModelR workflow consists of three main phases:

### 1. Preprocessing

- Align environmental rasters to a common projection and extent
- Rarefy occurrence data to reduce spatial and temporal bias
- Extract environmental values at occurrence locations
- Scale environmental rasters for model training

### 2. Modeling

- Partition occurrence data into spatiotemporal cross-validation folds
- Build hypervolume models for each fold
- Generate predictions across space and time
- Evaluate model performance in geographic and environmental space

### 3. Postprocessing

- Summarize predictions into consensus binary outputs
- Identify temporal patterns (increasing, decreasing, stable,
  fluctuating)
- Analyze trends by spatial units (e.g., counties, states, watersheds)

## Example Workflow: DC Metro Region Analysis

Here’s a complete workflow demonstrating the main functions using bird
occurrence data from the DC Metro region (1986-2023):

> **Note on Code Evaluation**: This example uses a mixed evaluation
> strategy: - **Preprocessing and modeling steps (Steps 1-9)**: Set to
> `eval=TRUE, message=FALSE, results='hide'` to avoid lengthy
> computation during README generation. These steps typically take hours
> to complete and require large spatial datasets. - **Postprocessing
> visualization steps (Steps 10-12)**: Set to
> `eval=TRUE, message=FALSE, results='hide'` to automatically generate
> plots from existing analysis results when knitting the README.
>
> To run the complete workflow from scratch, execute all chunks
> interactively in R or set `eval=TRUE, message=FALSE, results='hide'`
> for all chunks (be prepared for extended run times).

``` r
library(TemporalModelR)
library(terra)
library(sf)
library(raster)
```

### Phase 1: Preprocessing

#### Step 1: Prepare Reference Raster and Study Area

First, create a reference raster for your study area. This example uses
Loudoun and Montgomery counties:

``` r
counties_path <- "./DC_Metro/Shapefiles/DCMetro_Counties.shp"
reference_raster <- "./DC_Metro/reverse_water_mask_DCMetro.tif"

r <- raster(reference_raster)
counties <- st_read(counties_path)

loudoun <- counties[counties$NAMELSAD %in% c("Loudoun County", "Montgomery County"), ]

if (st_crs(loudoun) != crs(r)) {
  loudoun <- st_transform(loudoun, crs(r))
}

r_loudoun <- crop(r, loudoun)
r_loudoun <- mask(r_loudoun, loudoun)

r_loudoun_agg <- aggregate(r_loudoun, fact = 9, fun = mean, na.rm = TRUE)
```

#### Step 2: Align Environmental Rasters

Align all environmental rasters to match your reference raster’s
projection, extent, and resolution:

``` r
raster_align(
  input_dir = "G:/My Drive/VS_Rasters/DC_Metro/",
  output_dir = "G:/My Drive/VS_Rasters/DC_Metro/Masked_Projected_Variables_Simple/",
  reference_raster = r_loudoun_agg,
  overwrite = F
)
```

#### Step 3: Rarefy Occurrence Data

Reduce spatial and temporal bias in your occurrence data:

``` r
spatiotemporal_rarefication(
  points_sp = "./DC_Metro/simple_occupied_DC_Metro_AOU_5010.csv",
  xcol = "LONGDD",
  ycol = "LATDD",
  points_crs = 4326,
  output_dir = "./DC_Metro/PointFiles_simple/",
  reference_raster = r_loudoun_agg,
  time_cols = "Year"
)
```

This creates `Pts_Database_OnePerPixPerTimeStep.csv` with
spatiotemporally rarefied occurrences.

#### Step 4: Extract Environmental Values

Extract land cover and elevation values at occurrence locations,
calculating scaling parameters:

``` r
variable_patterns <- c(
  "Developed_Percentage2" = "Developed_Percentage2_YEAR",
  "Open_Percentage2" = "Open_Percentage2_YEAR",
  "Forest_Percentage2" = "Forest_Percentage2_YEAR",
  "elevation" = "elevation"
)

temporally_explicit_extraction(
  points_sp = "./DC_Metro/PointFiles_simple/Pts_Database_OnePerPixPerTimeStep.csv",
  xcol = "X",
  ycol = "Y",
  points_crs = 4326,
  raster_dir = "G:/My Drive/VS_Rasters/DC_Metro/Masked_Projected_Variables_Simple/",
  variable_patterns = variable_patterns,
  time_cols = "Year",
  output_dir = "./DC_Metro/PointFiles_simple/",
  output_prefix = "temp_explicit_df"
)
```

The function matches land cover rasters to occurrence years (e.g.,
`Developed_Percentage2_1990.tif` for 1990 occurrences), while static
variables like elevation are extracted once.

#### Step 5: Scale Environmental Rasters

Standardize all environmental rasters using the scaling parameters from
your occurrence data:

``` r
scale_rasters(
  input_dir = "G:/My Drive/VS_Rasters/DC_Metro/Masked_Projected_Variables_Simple/",
  output_dir = "G:/My Drive/VS_Rasters/DC_Metro/Scaled_simple/",
  scaling_params_file = "./DC_Metro/PointFiles_simple/temp_explicit_df_Scaling_Parameters.csv",
  variable_patterns = variable_patterns,
  time_cols = "Year",
  overwrite = TRUE
)
```

### Phase 2: Modeling

#### Step 6: Create Spatiotemporal Cross-Validation Folds

Partition occurrences into spatially and temporally independent folds:

``` r
partition_results <- spatiotemporal_partition(
  reference_shapefile_path = st_geometry(loudoun),
  points_file_path = "./DC_Metro/PointFiles_simple/temp_explicit_df_Scaled_Values.csv",
  time_col = "Year",
  xcol = "x",
  ycol = "y",
  points_crs = 4326,
  total_folds = 4,
  n_temporal = 2,
  n_spatial = 8,
  blocking_priority = "balanced",
  max_imbalance = 0.025,
  generate_plots = TRUE,
  output_file = "./DC_Metro/Partitioning_simple/partitioning_results_spatial.rds"
)
```

<img src="man/figures/README-partition_data-1.png" width="100%" />

#### Step 7: Build Hypervolume Models

Construct Gaussian hypervolumes for each cross-validation fold using
land cover variables:

``` r
hv_results <- build_hypervolume_models(
  partition_results = "./DC_Metro/Partitioning_simple/partitioning_results_spatial.rds",
  model_vars = c("Developed_Percentage2", "Open_Percentage2", "Forest_Percentage2"),
  method = "gaussian",
  output_dir = "./DC_Metro/VirtualSpecies_Results/Hypervolume_Gaussian_simple",
  hypervolume_params = list(
    quantile.requested = 0.95,
    quantile.requested.type = "probability"
  ),
  create_plot = TRUE,
  overwrite = TRUE
)
```

<img src="man/figures/README-build_hypervolumes-1.png" width="100%" /><img src="man/figures/README-build_hypervolumes-2.png" width="100%" /><img src="man/figures/README-build_hypervolumes-3.png" width="100%" /><img src="man/figures/README-build_hypervolumes-4.png" width="100%" /><img src="man/figures/README-build_hypervolumes-5.png" width="100%" />

#### Step 8: Generate Spatiotemporal Predictions

Project hypervolumes across all years (1986-2023) to create habitat
suitability predictions:

``` r
time_steps <- 1986:2023

variable_patterns <- c(
  "Developed_Percentage2" = "Developed_Percentage2_YEAR",
  "Open_Percentage2" = "Open_Percentage2_YEAR",
  "Forest_Percentage2" = "Forest_Percentage2_YEAR"
)

predictions <- generate_spatiotemporal_predictions(
  partition_results = "./DC_Metro/Partitioning_simple/partitioning_results_spatial.rds",
  hypervolume_results = "./DC_Metro/VirtualSpecies_Results/Hypervolume_Gaussian_simple/all_hypervolumes_gaussian.rds",
  time_col = "Year",
  time_steps = time_steps,
  variable_patterns = variable_patterns,
  raster_dir = "G:/My Drive/VS_Rasters/DC_Metro/Scaled_simple/",
  output_dir = "./predictions_from_package/",
  overwrite = FALSE
)
```

#### Step 9: Visualize Model Performance

Examine model evaluation metrics across the 38-year time series:

``` r
plot_model_assessment(
  data_file_path = "./predictions_from_package/Model_Assessment_Metrics.csv",
  time_column = "Year",
  separate_cbp = TRUE,
  cbp_threshold = 0.05
)
```

<img src="man/figures/README-plot_assessment-1.png" width="100%" /><img src="man/figures/README-plot_assessment-2.png" width="100%" /><img src="man/figures/README-plot_assessment-3.png" width="100%" /><img src="man/figures/README-plot_assessment-4.png" width="100%" />

### Phase 3: Postprocessing

#### Step 10: Create Consensus Predictions

Summarize predictions across all models to identify areas of agreement:

``` r
summary_results <- summarize_raster_outputs(
  predictions_dir = "./predictions_from_package/",
  output_dir = "./DC_Metro/VirtualSpecies_Results/Binary_Summaries_simple",
  overwrite = TRUE
)
```

#### Step 11: Identify Temporal Patterns

Apply changepoint detection to classify pixels into temporal trend
categories:

``` r
time_steps <- 1986:2023

pattern_results <- analyze_temporal_patterns(
  binary_stack = summary_results$binary_stack,
  summary_raster = summary_results$summary_raster,
  time_steps = time_steps,
  fastcpd_params = list(method = "BIC"),
  output_dir = "./output3",
  spatial_autocorrelation = TRUE,
  n_tiles_x = 2,
  n_tiles_y = 2,
  show_progress = TRUE,
  estimate_time = TRUE,
  overwrite = TRUE
)
```

<img src="man/figures/README-analyze_patterns-1.png" width="100%" />

The classification identifies seven pattern types:

- **Never Suitable** (dark red): Areas that remained unsuitable
  throughout
- **Decreasing** (green): Areas experiencing habitat loss
- **No Pattern** (gray): Areas with stable conditions
- **Increasing** (light green): Areas gaining suitable habitat
- **Always Suitable** (light red): Persistently suitable areas
- **Fluctuating** (purple): Areas with complex temporal dynamics
- **No Data** (yellow): Areas outside the study region or lacking data

#### Step 12: Summarize by County

Aggregate temporal patterns and trends for Loudoun and Montgomery
counties:

``` r
analyze_trends_by_spatial_unit(
  shapefile_path = st_as_sf(loudoun),
  name_field = "NAMELSAD",
  binary_stack = summary_results$binary_stack,
  pattern_raster = "./output3/pattern_raster_1986_2023.tif",
  year_decrease_raster = "./output3/year_first_decrease_1986_2023.tif",
  year_increase_raster = "./output3/year_first_increase_1986_2023.tif",
  output_dir = "./output3/spatial_analysis",
  time_steps = time_steps,
  pie_scale = 0.5
)
```

<img src="man/figures/README-spatial_trends-1.png" width="100%" /><img src="man/figures/README-spatial_trends-2.png" width="100%" /><img src="man/figures/README-spatial_trends-3.png" width="100%" /><img src="man/figures/README-spatial_trends-4.png" width="100%" /><img src="man/figures/README-spatial_trends-5.png" width="100%" /><img src="man/figures/README-spatial_trends-6.png" width="100%" />

The `analyze_trends_by_spatial_unit()` function generates multiple
output plots and saves them to the specified output directory:

- Pattern composition maps with pie charts
- Habitat availability time series by county
- Annual gains and losses
- Regional comparison plots

This analysis reveals county-level differences in habitat suitability
dynamics, showing how urbanization and land use change have affected
habitat availability from 1986 to 2023.

## Key Features

### Temporally Explicit Modeling

TemporalModelR is designed specifically for datasets where environmental
conditions change over time. The package:

- Matches environmental layers to occurrence timestamps
- Handles both static (e.g., elevation) and dynamic (e.g., land cover)
  variables
- Generates predictions for each time period in your study

### Hypervolume-Based Approach

Rather than traditional presence-background methods, TemporalModelR uses
hypervolumes to characterize species’ environmental niches:

- **Gaussian kernel density estimation** for smooth, probabilistic niche
  boundaries
- **One-class SVM** for complex, non-linear niche shapes
- Evaluation in both geographic space (G-space) and environmental space
  (E-space)

### Spatiotemporal Cross-Validation

The package implements sophisticated spatial and temporal partitioning:

- Creates spatially and temporally independent folds
- Generates Voronoi tessellations for spatial blocks
- Balances fold sizes across space and time
- Produces diagnostic plots for validation strategy assessment

### Comprehensive Temporal Analysis

Advanced changepoint detection identifies temporal trends:

- Classifies pixels as Never Suitable, Always Suitable, Increasing,
  Decreasing, No Pattern, or Fluctuating
- Accounts for spatial autocorrelation in pattern detection
- Identifies time periods when significant changes occurred
- Visualizes patterns across user-defined spatial units

## Function Reference

### Preprocessing Functions

- `raster_align()` - Align rasters to common
  projection/extent/resolution
- `spatiotemporal_rarefication()` - Rarefy occurrence data
  spatiotemporally
- `temporally_explicit_extraction()` - Extract environmental values at
  occurrences
- `scale_rasters()` - Standardize rasters using occurrence data
  statistics
- `spatiotemporal_partition()` - Create spatiotemporal cross-validation
  folds

### Modeling Functions

- `build_hypervolume_models()` - Build hypervolumes across CV folds
- `model_assessment_metrics()` - Calculate evaluation metrics (G-space
  and E-space)
- `generate_spatiotemporal_predictions()` - Project models across space
  and time

### Postprocessing Functions

- `summarize_raster_outputs()` - Create consensus predictions
- `analyze_temporal_patterns()` - Identify temporal trends via
  changepoint detection
- `plot_model_assessment()` - Visualize model performance over time
- `analyze_trends_by_spatial_unit()` - Aggregate trends by spatial units

## Data Requirements

### Occurrence Data

Your occurrence data should include:

- Geographic coordinates (latitude/longitude or projected coordinates)
- Temporal information (e.g., year, month)
- Species/observation identifiers

Example format:

    species,longitude,latitude,Year
    Species_A,-77.5,38.9,2015
    Species_A,-77.6,39.0,2016

### Environmental Data

Environmental rasters should:

- Cover your study area and time period
- Use consistent naming conventions with temporal placeholders
- Be in a format readable by `terra` (GeoTIFF recommended)

Example naming:

    Developed_Percentage2_1990.tif
    Developed_Percentage2_2000.tif
    Forest_Percentage2_1990.tif
    elevation.tif  (static variables have no time component)

### Spatial Units (Optional)

For regional summaries, provide:

- Polygon shapefile of administrative/ecological boundaries
- Attribute field with unique names for each unit

## Citation

If you use TemporalModelR in your research, please cite:

\[Citation information to be added\]

## Getting Help

- Report bugs or request features: [GitHub
  Issues](https://github.com/CJHughes926/TemporalModelR/issues)
- Ask questions: [GitHub
  Discussions](https://github.com/CJHughes926/TemporalModelR/discussions)

## License

\[License information to be added\]

## Acknowledgments

TemporalModelR builds on several excellent packages:

- [hypervolume](https://cran.r-project.org/package=hypervolume) for
  niche modeling
- [terra](https://cran.r-project.org/package=terra) for spatial data
  processing
- [sf](https://cran.r-project.org/package=sf) for vector spatial data
- [fastcpd](https://cran.r-project.org/package=fastcpd) for changepoint
  detection
- [exactextractr](https://cran.r-project.org/package=exactextractr) for
  spatial aggregation
