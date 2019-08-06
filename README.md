# enmSdm
<img align="right" src="enmSdm.png" height="190"/>
This package is a complement to the popular `dismo` package for R by Robert Hijmans. Its contains a suite of efficiency functions for preparing data, training and evaluating species distribution models, and comparing ecological niches.

You can install this package in R using these commands:

`install.packages('devtools') # if you haven't done this already`  
`library(devtools)`  
`install_github('adamlilith/omnibus')`  
`install_github('adamlilith/statisfactory')`  
`install_github('adamlilith/legendary')`  
`install_github('adamlilith/enmSdm')`  

NB: If for some reason these commands don't work, you can install the package(s) by downloading the latest zip/tar file from the `zipTarFiles` directory and installing the package(s) manually. If you do this, you will also have to install the `omnibus`,  `statisfactory`, and `legendary` packages, which are on GitHub also under my account (`adamlilith`).

## Data preparation ##
* `elimCellDups`: Eliminate duplicate points in each cell of a raster
* `geoFold`: Generate geographically distinct k-folds
* `geoThin` and `geoThinApprox`: Geographically thin points

## Model training ###
* `trainByCrossValid`: Wrapper for implementing any `trainXYZ` function across cross-validation folds
* `trainBrt`: Boosted regression trees (BRTs)
* `trainCrf`: Conditional regression trees (CRFs)
* `trainGam`: Generalized additive models (GAMs)
* `trainGlm`: Generalized linear models (GLMs)
* `trainGlmDredge`: Generalized linear models (GLMs)
* `trainLars`: Least-angle regression models (LARS)
* `trainMaxEnt` and `trainMaxNet`: Maxent models
* `trainNs`: Splines
* `trainRf`: Random forests (RFs)  

## Model evaluation ##
* `aucWeighted`: AUC (with/out site weights)
* `aucMultiWeighted`: Multivariate version of AUC (with/out site weight)
* `contBoyce`: Continuous Boyce Index (with/out site weights)
* `contBoyce2x`: "2X coverage" version of the Continuous Boyce Index (with/out site weights)
* `fpb`: Fpb (with/out site weights)
* `thresholdWeighted`: Thresholds to convert continuous predictions to binary predictions (with/out site weights)
* `thresholdStats`: Model performance statistics based on thresholded predictions (with/out site weights)
* `tssWeighted`: True Skill Statistic (TSS) (with/out site weights)
* `modelSize`: Number of response values in a model object

## Niche overlap ##
* `compareNiches`: Niche overlap metrics
* `compareResponse`: Compare niche model responses to a single variable
* `mop`: Calculate mobility-oriented parity, a measure of multivariate distance
* `nicheOverlap`: Calculate niche overlap as per Broennimann et al. Global Ecology and Biogeography 21:481-497
* `randPointsRespectingSelf`: Randomize geographic points while approximately respecting observed spatial autocorrelation structure between points
* `randPointsRespectingSelfOther2`: Randomize two sets of geographic points while approximately respecting observed spatial autocorrelation structure between and within sets
* `randPointsBatch`: Call `randPointsRespectingSelf` or `randPointsRespectingSelfOther2` multiple times
* `randPointsBatchExtract`: Extract environment from a set of rasters for sets of randomized points generated using `randPointsBatch`
* `randPointsBatchSampled`: Collate all sets of randomized points generated using `randPointsBatch`
* `randPointsBatchNicheOverlap`: Calculate niche overlap between sets of randomized points that were generated using `randPointsBatch`

## Spatial autocorrelation ##
* `localSpatialCorrForValues`: Calculate site-specific characteristic distance of local spatial autocorrelation for values associated with points or rasters
* `spatialCorrForPoints`: Calculate pairwise distance-based measure of global spatial autocorrelation between geographic points
* `spatialCorrForPointsSummary`: Characteristic cluster size of spatial points (distance of global autocorrelation)
* `spatialCorrForPointsPlot`: Plot observed and null distributions of pairwise distance-based measure of global spatial autocorrelation
* `spatialCorrForPointsWeight`: Assign weights to points based on pairwise distance-based measure of global spatial autocorrelation

## Functions for rasters
* `bioticVelocity`: Velocity of movement across a series of rasters
* `interpolateRasters`: Interpolate a stack of rasters
* `longLatRasters`: Generate rasters with values of longitude/latitude for cell values
* `sampleRast` and `sampleRastStrat`: Sample raster with/out replacement and possibly in a stratified manner

## Range area based on minimum convex polygons
* `mcpFromPolygons`: Minimum convex polygon from a set of polygons and points
* `areaFromPointsOrPoly`: Area of a spatial polygon or set of points

## Geographic utility functions ##
* `convertTropicosCoords`: Convert coordinates from the TROPICOS database
* `coordPrecision`: Calculate maximum possible coordinate precision
* `dmsToDecimal`: Convert degrees-minutes-seconds to decimal
* `getCRS`: Return a proj4string (coordinate reference system string) using a nickname
* `pointDist`: Geographic distance between set(s) of points
* `xToCoords`: Extract geographic coordinates from a data frame, matrix, or SpatialPoints* object
