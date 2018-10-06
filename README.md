# enmSdm
This package is a complement to the popular **dismo** package by Robert Hijmans. Its contains a suite of efficiency functions for preparing data, training and evaluating *species distribution models*, and comparing ecological niches.

## Data preperation ##
* `geoFold`: Generate geographically distinct k-folds
* `geoThin` and `geoThinApprox`: Geographically thin points

## Model training ###
* `trainBrt`: boosted regression trees (BRTs)
* `trainCrf`: conditional regression trees (CRFs)
* `trainGam`: generalized additive models (GAMs)
* `trainGlm`: generalized linear models (GLMs)
* `trainLars`: least-angle regression models (LARS)
* `trainMaxEnt` and `trainMaxNet`: Maxent models
* `trainNs`: splines
* `trainRf`: random forests (RFs)  

## Model evaluation ##
* `aucWeighted`: AUC (with/out weights)
* `contBoyce`: Continuous Boyce Index (with/out weights)
* `Fpb`: Fpb

## Niche overlap ##
* `compareNiches`: Niche overlap metrics
* `compareResponse`: Compare responses to a single variable
* `randGeoBySelf`: Randomize points while respecting observed spatial autocorrelation structure between points
* `mop`: Calculate mobility-oriented parity, a measure of multivariate distance

## Geographic utility functions ##
* `convertTropicosCoords`: Convert coordinates from the TROPICOS database
* `coordPrecision`: Calculate maximum possible coordinate precision
* `dmsToDecimal`: Convert degrees-minutes-seconds to decimal
* `elimCellDups`: Eliminate duplicate points in each cell of a raster
* `getCRS`: Return a proj4string (coordinate reference system string) using a nikname
* `longLatRasters`: Generate rasters with values of longitude/latitude for cell values
* `pointDist`: Geographic distance between set(s) of points
* `sampleRast` and `sampleRastStrat`: Sample raster with/out replacement and possibly in a stratified manner
* `xToCoords`: Extract geographic coordinates from a data frame, matrix, or SpatialPoints* object
