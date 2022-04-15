#' enmSdm: Species distribution modeling and ecological niche modeling
#'
#' This package contains tools for modeling the distributions and niches of species or species-like entities. Its main features are a set of functions for training ENMs, evaluating niche overlap, and correcting for sampling bias.
#'
#' Create an issue on \href{https://github.com/adamlilith/enmSdm/issues}{GitHub}.
#'
#' @details
#' @section Data preparation:
#' 		\code{\link{elimCellDups}} Eliminate duplicate points in each cell of a raster \cr
#' 		\code{\link{geoFold}} Generate geographically distinct k-folds \cr
#' 		\code{\link{geoThin}} and \code{\link{geoThinApprox}} Geographically thin points \cr
#'
#' @section Model training:
#' 		\code{\link{trainByCrossValid}} Wrapper for implementing some \code{trainXYZ} function across cross-validation folds (see also \code{summaryByCrossValid}). \cr
#' 		\code{\link{trainBrt}} Boosted regression trees (BRTs) \cr
#' 		\code{\link{trainCrf}} Conditional regression trees (CRFs) \cr
#' 		\code{\link{trainGam}} Generalized additive models (GAMs) \cr
#' 		\code{\link{trainGlm}} Generalized linear models (GLMs) \cr
#' 		\code{\link{trainGlmDredge}} Generalized linear models (GLMs) \cr
#' 		\code{\link{trainMaxEnt}} and \code{\link{trainMaxNet}} Maxent models \cr
#' 		\code{\link{trainNs}} Natural splines (NSs) \cr
#' 		\code{\link{trainRf}} Random forests (RFs) \cr
#'
#' @section Model prediction:
#' 		\code{\link{predictEnmSdm}} Predict most model types using default settings \cr
#' 		\code{\link{predictMaxNet}} Predict MaxNet (MaxEnt) model \cr
#'
#' @section Model evaluation:
#' 		\code{\link{aucWeighted}} AUC (with/out site weights) \cr
#' 		\code{\link{aucMultiWeighted}} Multivariate version of AUC (with/out site weight) \cr
#' 		\code{\link{contBoyce}} Continuous Boyce Index (with/out site weights) \cr
#' 		\code{\link{contBoyce2x}} "2X coverage" version of the Continuous Boyce Index (with/out site weights) \cr
#' 		\code{\link{fpb}} Fpb (with/out site weights) \cr
#' 		\code{\link{thresholdWeighted}} Thresholds to convert continuous predictions to binary predictions (with/out site weights) \cr
#' 		\code{\link{thresholdStats}} Model performance statistics based on thresholded predictions (with/out site weights) \cr
#' 		\code{\link{tjursR2Weighted}} Tjur's R2 (with/out site weights) \cr
#' 		\code{\link{tssWeighted}} True Skill Statistic (TSS) (with/out site weights) \cr
#' 		\code{\link{modelSize}} Number of response values in a model object \cr
#'
#' @section Niche overlap:
#' 		\code{\link{compareNiches}} Niche overlap metrics \cr
#' 		\code{\link{compareResponse}} Compare niche model responses to a single variable \cr
#' 		\code{\link{mop}} Calculate mobility-oriented parity, a measure of multivariate distance as per Saupe et al. 2012. \cr
#' 		\code{\link{nicheOverlap}} Calculate niche overlap as per Broennimann et al. Global Ecology and Biogeography 21:481-497. \cr
#' 		\code{\link{randPointsRespectingSelf}} Randomize geographic points while approximately respecting observed spatial autocorrelation structure between points \cr
#' 		\code{\link{randPointsRespectingSelfOther2}} Randomize two sets of geographic points while approximately respecting observed spatial autocorrelation structure between and within sets \cr
#' 		\code{\link{randPointsBatch}} Call \code{\link{randPointsRespectingSelf}} or \code{\link{randPointsRespectingSelfOther2}} multiple times \cr
#' 		\code{\link{randPointsBatchExtract}} Extract environment from a set of rasters for sets of randomized points generated using `randPointsBatch` \cr
#' 		\code{\link{randPointsBatchSampled}} Collate all sets of randomized points generated using `randPointsBatch` \cr
#' 		\code{\link{randPointsBatchNicheOverlap}} Calculate niche overlap between sets of randomized points that were generated using `randPointsBatch` \cr
#' @section Spatial autocorrelation:
#' 		\code{\link{spatialCorrForPoints}} Calculate pairwise distance-based measure of global spatial autocorrelation between geographic points \cr
#' 		\code{\link{spatialCorrForPointsSummary}} Characteristic cluster size of spatial points (distance of global autocorrelation) \cr
#' 		\code{\link{spatialCorrForPointsPlot}} Plot observed and null distributions of pairwise distance-based measure of global spatial autocorrelation \cr
#' 		\code{\link{spatialCorrForPointsWeight}} Assign weights to points based on pairwise distance-based measure of global spatial autocorrelation \cr
#'
#' @section Functions for rasters:
#' 		\code{\link{bioticVelocity}} Velocity of movement across a series of rasters \cr
#'		\code{\link{getCores}} Get number of processor cores \cr
#' 		\code{\link{interpolateRasters}} Interpolate a stack of rasters \cr
#' 		\code{\link{longLatRasters}} Generate rasters with values of longitude/latitude for cell values \cr
#' 		\code{\link{rastWithSquareCells}} Create a raster with square cells \cr
#' 		\code{\link{sampleRast}} and \code{\link{sampleRastStrat}} Sample raster with/out replacement and possibly in a stratified manner \cr
#' 		\code{\link{squareRastCells}} Resample a raster so cells are "square" \cr
#'
#' @section Range area based on minimum convex polygons:
#' 		\code{\link{mcpFromPolygons}} Minimum convex polygon from a set of polygons \emph{and} points \cr
#' 		\code{\link{areaFromPointsOrPoly}} Area of a spatial polygon or set of points \cr
#'
#' @section Geographic utility functions:
#' 		\code{\link{convertTropicosCoords}} Convert coordinates from the TROPICOS database \cr
#' 		\code{\link{coordPrecision}} Calculate maximum possible coordinate precision \cr
#'		\code{\link{createPlotPoly}} Create a `SpatialPolygon` the same size as a plot region \cr
#' 		\code{\link{decimalToDms}} Convert decimal coordinate to degrees-minutes-seconds \cr
#' 		\code{\link{dmsToDecimal}} Convert degrees-minutes-seconds coordinate to decimal \cr
#' 		\code{\link{getCRS}} Return a proj4string (coordinate reference system string) using a nickname \cr
#' 		\code{\link{getEPSG}} Return a EPSG code0 (coordinate reference system code) using a nickname \cr
#' 		\code{\link{pointDist}} Geographic distance between set(s) of points \cr
#'		\code{\link{svToSpatial}} Convert \code{SpatVector} object to a \code{Spatial}* object. \cr
#' 		\code{\link{xToCoords}} Extract geographic coordinates from a data frame, matrix, or SpatialPoints* object \cr
#' @section Data:
#' 		\code{\link{lemurs}}: Lemur occurrences \cr
#' 		\code{\link{mad0}}: Madagascar spatial object \cr
#' @docType package
#' @author Adam B. Smith
#' @name enmSdm
#' @keywords internal
NULL
