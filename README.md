# enmSdm
This package is a complement to the popular **dismo** package for R by Robert Hijmans. Its contains a suite of efficiency functions for preparing data, training and evaluating *species distribution models*, and comparing ecological niches.

You can install this package in R using these commands:

`install.packages('devtools') # if you haven't done this already`  
`library(devtools)`  
`install_github('adamlilith/omnibus')`  
`install_github('adamlilith/statisfactory')`  
`install_github('adamlilith/enmSdm')`  

NB: If for some reason these commands don't work, you can install the package(s) by downloading the latest zip/tar file from the `zipTarFiles` directory and installing the package(s) manually.

## Data preparation ##
* `geoFold`: Generate geographically distinct k-folds
* `geoThin` and `geoThinApprox`: Geographically thin points

## Model training ###
* `trainBrt`: boosted regression trees (BRTs)
* `trainCrf`: conditional regression trees (CRFs)
* `trainGam`: generalized additive models (GAMs)
* `trainGlm`: generalized linear models (GLMs)
* `trainGlmDredge`: generalized linear models (GLMs)
* `trainLars`: least-angle regression models (LARS)
* `trainMaxEnt` and `trainMaxNet`: Maxent models
* `trainNs`: splines
* `trainRf`: random forests (RFs)  

## Model evaluation ##
* `aucWeighted`: AUC (with/out site weights)
* `aucMultiWeighted`: Multivariate version of AUC (with/out site weight)
* `contBoyce`: Continuous Boyce Index (with/out site weights)
* `Fpb`: Fpb (with/out site weights)

## Model information
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
* `spatialCorrForPoints`: Calculate pairwise distance-based measure of spatial autocorrelation between geographic points
* `spatialCorrForPointsSummary`: Characteristic cluster size of spatial points (distance of autocorrelation)
* `spatialCorrForPointsPlot`: Plot observed and null distributions of pairwise distance-based measure of spatial autocorrelation
* `spatialCorrForPointsWeight`: Assign weights to points based on pairwise distance-based measure of spatial autocorrelation

## Geographic utility functions ##
* `convertTropicosCoords`: Convert coordinates from the TROPICOS database
* `coordPrecision`: Calculate maximum possible coordinate precision
* `dmsToDecimal`: Convert degrees-minutes-seconds to decimal
* `elimCellDups`: Eliminate duplicate points in each cell of a raster
* `getCRS`: Return a proj4string (coordinate reference system string) using a nickname
* `longLatRasters`: Generate rasters with values of longitude/latitude for cell values
* `pointDist`: Geographic distance between set(s) of points
* `sampleRast` and `sampleRastStrat`: Sample raster with/out replacement and possibly in a stratified manner
* `xToCoords`: Extract geographic coordinates from a data frame, matrix, or SpatialPoints* object
