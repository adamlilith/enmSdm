#' Velocity of shifts in densities across a series of rasters
#'
#' This function calculates several metrics of biotic velocity for a stack of rasters or an array representing a gridded landscape.
#' @param x Either a \code{RasterStack}, \code{RasterBrick}, or 3-dimensional array. Values should be either \code{NA} or >= 0.
#' \itemize{
#' 	\item If \code{x} is a \code{RasterStack} or \code{RasterBrick} then each layer is assumed to represent a time slice and the rasters \emph{must} be in an equal-area projection.
#' 	\item If \code{x} is an array then each "layer" in the third dimension is assumed to represent a map at a particular time slice in an equal-area projection. Note that if this is an array you should probably also specify the arguments \code{longitude} and \code{latitude}.
#' }
#' @param times Numeric vector \emph{or} with the same number of layers in \code{x} or \code{NULL} (default). This specifies the time period represented by each layer in \code{x} from oldest (top layer) to most recent (bottom layer). Times \emph{must} appear in sequential order. For example, if time periods are 24 kybp, 23 kybp, 22 kybp, use \code{c(-24, -23, -22)}, \emph{not} \code{c(24, 23, 22)}. If \code{NULL} (default), values are assigned starting at 1 and ending at the total number of rasters in \code{x}.
#' @param atTimes Numeric, values of \code{times} across which to calculate biotic velocity. You can use this to calculate biotic velocities across selected time periods (e.g., just the first and last time periods). Note that \code{atTimes} \emph{must} be the same as or a subset of \code{times}. The default is \code{NULL}, in which case velocity is calculated across all time slices (i.e., between period 1 and 2, 2 and 3, 3 and 4, etc.).
#' @param longitude Numeric matrix or \code{NULL} (default):
#' \itemize{
#'	\item If \code{x} is a \code{RasterStack} or \code{RasterBrick} then this is ignored (longitude is ascertained directly from the rasters, which \emph{must} be in equal-area projection for velocities to be valid).
#'	\item If \code{x} is an array and \code{longitude} is \code{NULL} (default), then longitude will be ascertained from column numbers in \code{x} and velocities will be in arbitrary spatial units (versus, for example, meters). Alternatively, this can be a two-dimensional matrix whose elements represent the longitude coordinates of the centers of cells of \code{x}. The matrix must have the same number of rows and columns as \code{x}. Coordinates must be from an equal-area projection for results to be valid.
#' }
#' @param latitude Numeric matrix or \code{NULL} (default):
#' \itemize{
#'	\item If \code{x} is a \code{RasterStack} or \code{RasterBrick} then this is ignored (latitude is obtained directly from the rasters, which \emph{must} be in equal-area projection for velocities to be valid).
#'	\item If \code{x} is an array and \code{latitude} is \code{NULL} (default), then latitude will be obtained from row numbers in \code{x} and velocities will be in arbitrary spatial units (versus, for example, meters). Alternatively, this can be a two-dimensional matrix whose elements represent the latitude coordinates of the centers of cells of \code{x}. The matrix must have the same number of rows and columns as \code{x}. Coordinates must be from an equal-area projection for results to be valid.
#' }
#' @param elevation Either \code{NULL} (default) or a raster or matrix representing elevation. If this is supplied, eleational change in range centroid is calculated (no need to specify this in \code{metrics}).
#' @param metrics Biotic velocity metrics to calculate (default is to calculate them all). All metrics ignore \code{NA} cells in \code{x}. Here "starting time period" represents one layer in \code{x} and "end time period" the next layer.
#' \itemize{
#' 	\item \code{centroid}: Movement of mass-weighted centroid.
#'  \item \code{nsCentroid} or \code{ewCentroid}: Movement in the north-south or east-west directions of the mass-weighted centroid. For north-south cardinality, positive values represent movement northward and negative southward. For east-west cardinality, positive values represent movement eastward and negative westward.
#'  \item \code{nCentroid}, \code{sCentroid}, \code{eCentroid}, and \code{wCentroid}: Movement of mass-weighted centroid of the portion of the raster north/south/east/west of the landscape-wide weighted centroid of the starting time period.
#'  \item \code{nsQuants} or \code{ewQuants}: Movement of the location of the \emph{n}th quantile of mass in the north-south or east-west directions. The quantiles can be specified in \code{quants}. For example, this could be the movement of the 5th, 50th, and 95th quantiles of population size going from south to north. The 0th quantile would measure the velocity of the southernmost or easternmost cell(s) with values >0, and the 100th quantile the northernmost or westernmost cell(s) with non-zero values.
#'  \item \code{mean}: Mean value across all cells (this is not really a measure of "velocity").
#'  \item \code{sum}: Total across all cells (this is not really a measure of "velocity").
#'  \item \code{quants}: \emph{N}th quantile values across all cells (this is not really a measure of "velocity"). Quantiles are given by \code{quants}.
#'  \item \code{prevalence}: Number of cells with values > 0 (this is not really a measure of "velocity").
#'  \item \code{similarity}: Several metrics of similarity between each time period. Some of these make sense for cases where values in \code{x} are in the range [0, 1], but not if some values are outside this range. The metrics are the simple mean difference, mean absolute difference, the Expected Fraction of Shared Presences or ESP (Godsoe 2014), the D statistic (Schoener 1968), and I (Warren et al. 2008). The last three metrics have been modified from their original versions by dividing the values by the total number of cells (see \emph{Details}) so they fall in the range [0, 1] provided that cell values are in the range [0, 1] or \code{NA}.
#'  \item \code{elevCentroid}: Velocity of the centroid of mass in elevation (up or down). Note that this metric is not calculated unless \code{elevation} is supplied.
#' 	\item \code{elevQuants}: Velocity of the emph{n}th quantile of mass in elevation (up or down). The quantiles to be evaluated are given by \code{quants}. The lowest elevation with mass >0 is the 0th quantile, and the highest elevation with mass >0 is the 100th. Note that this metric is not calculated unless \code{elevation} is supplied.
#' }
#' @param quants Numeric vector indicating the quantiles at which biotic velocity is calculated for the "\code{quant}" and "\code{Quants}" metrics. Default is \code{c(0.05, 0.10, 0.5, 0.9, 0.95)}.
#' @param onlyInSharedCells Logical, if \code{TRUE}, calculate biotic velocity using only those cells that are not \code{NA} in the start and end of each time period. This is useful for controlling for shifting land mass due to sea level rise, for example, when calculating biotic velocity for an ecosystem or a species. The default is \code{FALSE}.
#' @param cores Positive integer. Number of processor cores to use. Note that if the number of time steps at which velocity is calculated is small, using more cores may not always be faster.
#' @param warn Logical, if \code{TRUE} (default) then display function-specific warnings.
#' @param ... Other arguments (not used).
#' @return A data frame with biotic velocities and related values. Fields are as follows:
#' \itemize{
#' 	\item \code{timeFrom}: Start time of interval
#' 	\item \code{timeTo}: End time of interval
#' 	\item \code{timeSpan}: Duration of interval
#' 	\item \code{propSharedCellsNotNA}: Proportion of cells that were not \code{NA} in both starting and ending time periods
#' 	\item \code{timeFromPropNotNA} and \code{timeToPropNotNA}: Proportion of cells that were not \code{NA} in the starting and ending time periods
#' }
#' Depending on \code{metrics} that are specified, additional fields are as follows. All measurements of velocity are in distance units (typically meters) per time unit (which is the same as the units used for \code{times} and \code{atTimes}). For example, if the rasters are in an Albers equal-area projection and times are in years, then the output will be meters per year.
#' \itemize{
#' 	\item If \code{metrics} has \code{'centroid'}: Columns named \code{centroidVelocity}, \code{centroidLong}, \code{centroidLat} -- Velocity of weighted centroid, plus its longitude and latitude (in the "to" time period of each time step).
#' 	\item If \code{metrics} has \code{'nsCentroid'}: Columns named \code{nsCentroid}, \code{nsCentroidLat} -- Velocity of weighted centroid in north-south direction, plus its latitude (in the "to" time period of each time step).
#' 	\item If \code{metrics} has \code{'ewControid'}: \code{ewCentroid}, \code{ewCentroidLong} -- Velocity of weighted centroid in east-west direction, plus its longitude (in the "to" time period of each time step).
#' 	\item If \code{metrics} has \code{'nCentroid'}, \code{'sCentroid'}, \code{'eCentroid'}, and/or \code{'wCentroid'}: Columns named \code{nCentroidVelocity} and \code{nCentroidAbund}, \code{sCentroid} and \code{sCentroidAbund}, \code{eCentroid} and \code{eCentroidAbund}, and/or \code{wCentroid} and \code{wCentroidAbund}: Velocity of weighted centroid of all cells that fall north, south, east, or west of the landscape-wide centroid, plus a column indicating the weight (abundance) of all such populations.
#' 	\item if \code{metrics} contains any of \code{nsQuants} or \code{ewQuants}: Columns named \code{nsQuantVelocity_quantN} and \code{nsQuantLat_quantN}, or \code{ewQuantVelocity_quantN} and \code{ewQuantLat_quantN}: Velocity of the \emph{N}th quantile weight in the north-south or east-west directions, plus the latitude or longitude thereof (in the "to" time period of each time step). Quantiles are cumulated starting from the south or the east, so the 0.05th quantile, for example, is in the south or east of the range and the 0.95th in the north or west.
#' 	\item If \code{metrics} contains \code{mean}: A column named \code{mean} -- Mean weight in "timeTo" time step. In the same units as the values of the cells.
#' 	\item If \code{metrics} contains \code{quants}: A column named \code{quantile_quantN} -- The \emph{N}th quantile(s) of weight in the "timeTo" time step. In the same units as the values of the cells.
#' 	\item If \code{metrics} contains \code{prevalence}: A column named \code{prevalence} -- Proportion of non-\code{NA} cells with weight >0 in the "timeTo" time step. Unitless.
#' 	\item If \code{metrics} contains \code{similarity}: Columns named \code{simpleMeanDiff}, \code{meanAbsDiff}, \code{godsoeEsp}, \code{schoenersD}, and \code{warrensI} -- Measures of similarity (see \emph{Details}).
#' 	\item If \code{metrics} contains \code{elevCentroid}: Columns named \code{elevCentroidVelocity} and \code{elevCentroidElev} -- Velocity of the centroid in elevation (up or down) and the elevation in the "to" timestep. Velocity is signed (positive means up, negative means down).
#' 	\item If \code{metrics} contains \code{elevQuants}: Columns named \code{elevQuantVelocity_quantN} and \code{elevQuantVelocityElev_quantN} -- Velocity of the \emph{N}th quantile of mass in elevation (up or down) and the elevation of this quantile in the "to" timestep. Velocity is signed (positive means up, negative means down).
#' }
#' @details
#' \emph{Attention:}  
#'   
#' This function may yield erroneous velocities if the region of interest is near or spans a pole or the international date line. It converts rasters to arrays before doing calculations, so using very large rasters may yield slow performance and may not even work, depending on memory requirements. Results using the "Quant" and "quant" metrics may be somewhat counterintuitive if just one cell >0, or one row or column all with the same values with all other values equal to 0 or NA because defining quantiles in these situations is not intuitive. Results may also be counterintuitive if some cells have negative values because they can "push" a centroid away from what would seem to be the center of mass as assessed by visual examination of a map.  
#'  
#' \emph{Similarity metrics:}  
#'   
#' The similarity metrics are defined for two rasters or matrices \code{x1} and \code{x2} and the mean number of non-\code{NA} cells between them (\code{n}):
#' \itemize{
#' 	\item \code{simpleMeanDiff}: \code{sum(x2 - x1, na.rm=TRUE) / n}
#' 	\item \code{absMeanDiff}: \code{sum(abs(x2 - x1), na.rm=TRUE) / n}
#' 	\item \code{godsoeEsp}: \code{1 - sum(2 * (x1 * x2), na.rm=TRUE) / sum(x1 + x2, na.rm=TRUE)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar
#' 	\item \code{schoenerD}: \code{1 - (sum(abs(x1 - x2), na.rm=TRUE) / n)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar
#' 	\item \code{warrenI}: \code{1 - sqrt(sum((sqrt(x1) - sqrt(x2))^2, na.rm=TRUE) / n)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar
#' }
#'   
#' \emph{Note:}  
#'   
#' For the \code{nsQuant} and \code{ewQuant} metrics it is assumed that the latitude/longitude assigned to a cell is at its exact center. If a desired quantile does not fall exactly on the cell center, it is interpolated linearly between the rows/columns of cells that bracket the given quantile. For quantiles that fall south/eastward of the first row/column of cells, the cell border is assumed to be at 0.5 * cell length south/west of the cell center.
#' @references
#' Schoener, T. W. 1968. \emph{Anolis} lizards of Bimini: Resource partitioning in a complex fauna. \emph{Ecology} 49:704–726.  
#' Godsoe, W. 2014. Inferring the similarity of species distributions using Species’ Distribution Models. \emph{Ecography} 37:130-136.  
#' Warren, D.L., Glor, R.E., and Turelli, M. 2008. Environmental niche equivalency versus conservatism: Quantitative approaches to niche evolution. \emph{Evolution} 62:2868-2883.
#' @examples
#'
#' ### SIMPLE EXAMPLES ###
#' #######################
#'
#' library(raster)
#' 
#' ### movement in north-south directions
#' mat <- matrix(0, nrow=5, ncol=5)
#' mat1 <- mat2 <- mat
#' mat1[3, 3] <- 1
#' mat2[2, 3] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # note that nCentroidVelocity is NaN because just one cell is occupied
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs)) # north movement
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs)) # south movement
#' 
#' 
#' ### movement in east-west directions
#' mat <- matrix(0, nrow=5, ncol=5)
#' mat1 <- mat2 <- mat
#' mat1[3, 3] <- 1
#' mat2[3, 2] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # movement east
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' # movement west
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' 
#' ### movement of portions of range northward/southward
#' mat <- matrix(0, nrow=11, ncol=11)
#' mat1 <- mat2 <- mat
#' mat1[6, 5] <- 1 # bottom
#' mat1[5, 5] <- 1 # center
#' mat1[5, 4] <- 1 # west
#' mat1[5, 6] <- 1 # east
#' mat1[4, 5] <- 1 # north
#' 
#' mat2[6, 5] <- 1 # bottom
#' mat2[5, 5] <- 1 # center
#' mat2[5, 4] <- 1 # west
#' mat2[5, 6] <- 1 # east
#' mat2[3, 5] <- 1 # north
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # northern section moves north
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs))
#' # northern section moves south
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' ### quantile velocities: north/south movement
#' mat <- matrix(0, nrow=11, ncol=11)
#' mat1 <- mat2 <- mat
#' 
#' mat1[2:10, 6] <- 1
#' mat2[1:9, 6] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # shift north
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs))
#' # shift south
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' 
#' ### quantile velocities: east/west movement
#' mat <- matrix(0, nrow=11, ncol=11)
#' mat1 <- mat2 <- mat
#' 
#' mat1[6, 2:10] <- 1
#' mat2[6, 3:11] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # eastward shift
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs))
#' # westward shift
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' ### big block test
#' mat <- matrix(0, nrow=7, ncol=7)
#' mat1 <- mat2 <- mat
#' 
#' mat1[3:5, 3:5] <- 1
#' mat2[1:3, 1:3] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' # shift northwest
#' (bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs))
#' # shift southeast
#' (bioticVelocity(mats[ , , 2:1], times=1:2,
#' latitude=lats, longitude=longs))
#' 
#' 
#' ### big block test: common frame
#' mat <- matrix(0, nrow=7, ncol=7)
#' mat1 <- mat2 <- mat
#' 
#' mat1[3:5, 3:5] <- 1
#' mat1[1, ] <- NA
#' mat1[ , 1] <- NA
#' mat2[1:3, 1:3] <- 1
#' 
#' mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' matsRast <- stack(mat1rast, mat2rast)
#' plot(matsRast, col=c('gray', 'darkgreen'))
#' 
#' v1 <- bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs)
#' v2 <- bioticVelocity(mats, times=1:2,
#' latitude=lats, longitude=longs, onlyInSharedCells=TRUE)
#' (rbind(v1, v2))
#' 
#' ### elevational centroid velocity
#' 
#' # single cell moving
#' # goes up, holds, then down
#' mat <- matrix(0, nrow=7, ncol=7)
#' mat1 <- mat2 <- mat3 <- mat4 <- elevation <- mat
#' 
#' mat1[4, 4] <- 1
#' mat2[3, 4] <- 1
#' mat3[3, 4] <- 1
#' mat4[2, 4] <- 1
#' 
#' elevation[4, 4] <- 1
#' elevation[3, 4] <- 4
#' elevation[2, 4] <- 2
#' 
#' mats <- array(c(mat1, mat2, mat3, mat4), dim=c(nrow(mat1), ncol(mat1), 4))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' mat3rast <- raster(mat3)
#' mat4rast <- raster(mat4)
#' 
#' rownames(elevation) <- paste0('lat', nrow(mats):1)
#' colnames(elevation) <- paste0('long', 1:ncol(mats))
#' 
#' elevRast <- raster(elevation)
#' plot(stack(elevRast, mat1rast, mat2rast, mat3rast, mat4rast))
#' 
#' (elevVel <- bioticVelocity(mats, metrics='elevCentroid',
#' times=1:4, latitude=lats, longitude=longs, elevation=elevation))
#' 
#' ### elevational centroid quantiles
#' 
#' # quantiles mostly move up, hold, then down
#' mat <- matrix(0, nrow=4, ncol=4)
#' mat1 <- elevation <- mat
#' 
#' mat1[] <- runif(16)
#' mat2[] <- sort(mat1)
#' mat3 <- mat2
#' mat4 <- matrix(rev(mat3), nrow=4)
#' 
#' elevation[] <- 1:16
#' 
#' mats <- array(c(mat1, mat2, mat3, mat4), dim=c(nrow(mat1), ncol(mat1), 4))
#' rownames(mats) <- paste0('lat', nrow(mats):1)
#' colnames(mats) <- paste0('long', 1:ncol(mats))
#' lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
#' longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)
#' 
#' mat1rast <- raster(mat1)
#' mat2rast <- raster(mat2)
#' mat3rast <- raster(mat3)
#' mat4rast <- raster(mat4)
#' 
#' names(mat1rast) <- 'time1'
#' names(mat2rast) <- 'time2'
#' names(mat3rast) <- 'time3'
#' names(mat4rast) <- 'time4'
#' 
#' rownames(elevation) <- paste0('lat', nrow(mats):1)
#' colnames(elevation) <- paste0('long', 1:ncol(mats))
#' 
#' elevRast <- raster(elevation)
#' names(elevRast) <- 'elev'
#' 
#' plot(stack(elevRast, mat1rast, mat2rast, mat3rast, mat4rast))
#' 
#' (elevVel <- bioticVelocity(mats, metrics='elevQuants',
#' times=1:4, latitude=lats, longitude=longs, elevation=elevation))
#' 
#' \donttest{
#' ### multi-core
#' mats <- array(runif(100 * 100 * 1000), dim=c(100, 100, 500))
#' 
#' mats <- brick(mats)
#' projection(mats) <- getCRS('wgs84')
#' 
#' mc <- bioticVelocity(x=mats, metrics='centroid', cores=4)
#' }
#' 
#' @export

bioticVelocity <- function(
	x,
	times = NULL,
	atTimes = NULL,
	longitude = NULL,
	latitude = NULL,
	elevation = NULL,
	metrics = c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentroid', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants', 'mean', 'sum', 'quants', 'prevalence', 'similarity'),
	quants = c(0.05, 0.10, 0.5, 0.9, 0.95),
	onlyInSharedCells = FALSE,
	cores = 1,
	warn = TRUE,
	...
) {

	xClass <- class(x)

	### time of each period and times across which to calculate velocity
	####################################################################
	
		# total number of time periods
		totalTimes <- if ('array' %in% xClass) {
			dim(x)[3]
		} else if (xClass %in% c('RasterStack', 'RasterBrick')) {
			raster::nlayers(x)
		}

		# time of each period
		if (is.null(times)) times <- 1:totalTimes
	
		# times across which to calculate velocity
		if (is.null(atTimes)) atTimes <- times

		### indices of layers to use (discard others)
		atIndices <- which(times %in% atTimes)
	
	### catch errors
	################
		
		if (!all(order(times) == seq_along(times))) {
			stop('Times assigned to each period (argument "times") are not\nsequential (e.g., {1, 2, 3} versus {3, 2, 1}).')
		}
	
		if (!all(order(atTimes) == seq_along(atTimes))) {
			stop('Values in argument "atTimes" are not sequential (e.g.,\n{1, 2, 3} versus {3, 2, 1}).')
		}
	
		if (length(times) != totalTimes) stop('The length of "times" does not match the total number of time periods represented by "x".')
		if (!all(atTimes %in% times)) stop('All time slices specified in "atTimes" must also appear in "times".')
		
		if (is.null(elevation) & ('elevCentroid' %in% metrics | 'elevQuants' %in% metrics)) {
			
			if (warn) warning('If argument "metrics" includes "elevCentroid" and/or "elevQuants",\nthen argument "elevation" must be non-NULL. These metrics will not\nbe calculated.')
			
			if (any('elevCentroid' == metrics)) metrics <- metrics[-which(metrics == 'elevCentroid')]
			if (any('elevQuants' == metrics)) metrics <- metrics[-which(metrics == 'elevQuants')]
			
		}
		
	### convert input to array and get geographic information
	#########################################################

		# if rasters
		if (xClass %in% c('RasterStack', 'RasterBrick')) {

			if (is.null(longitude) | is.null(latitude)) {
		
				ll <- enmSdm::longLatRasters(x)
				longitude <- ll[['longitude']]
				latitude <- ll[['latitude']]
				
			}
			x <- raster::subset(x, atIndices)

		}
		
		# if input is an array and extent/CRS are specified
		# then get longitude and latitude and convert to array
		if ('array' %in% xClass) {

			x <- x[ , , atIndices]
			nRows <- dim(x)[1]
			nCols <- dim(x)[2]
			
			if (is.null(longitude)) {
				longitude <- matrix(1:nCols, nrow=nRows, ncol=nCols, byrow=TRUE)
				if (warn) warning('Argument "longitude" is not specified so using column number instead of longitude. Velocities will be in arbitrary spatial units.')
			}
			if (is.null(latitude)) {
				latitude <- matrix(nRows:1, nrow=nRows, ncol=nCols, byrow=FALSE)
				if (warn) warning('Argument "latitude" is not specified so using row number instead of latitude. Velocities will be in arbitrary spatial units.')
			}
			
			# convert array to raster stack... assuming cell sizes from longitude/latitude
			xmin <- longitude[1, 1] - 0.5 * (longitude[1, 2] - longitude[1, 1])
			xmax <- longitude[1, ncol(longitude)] + 0.5 * (longitude[1, 2] - longitude[1, 1])
			ymax <- latitude[1, 1] + 0.5 * (latitude[1, 1] - latitude[2, 1])
			ymin <- latitude[nrow(latitude), 1] - 0.5 * (latitude[1, 1] - latitude[2, 1])
		
			x <- raster::brick(x, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax)
			longitude <- raster::raster(longitude, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax)
			latitude <- raster::raster(latitude, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax)
			
			if (!is.null(elevation)) {
				if ('matrix' %in% class(elevation)) {
					elevation <- raster(elevation, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax)
				}
			}
				
		}
		
		# if (any(raster::minValue(x) < 0, na.rm=TRUE) & warn) warning('Negative values appear in "x". Output may be unreliable or undesirable.')
			
	### calculate weighted longitude and latitudes for starting time period
	#######################################################################

		# x1 <- x[[1]]
		# if (!is.null(elevation)) elev <- elevation

		# # correction for shared non-NA cells with next time period
		# if (onlyInSharedCells) {

			# x2 <- x[[2]]
			# x1mask <- x1 * 0 + 1
			# x2mask <- x2 * 0 + 1
			# x1x2mask <- x1mask * x2mask
			# x1 <- x1 * x1x2mask
			
			# if (!is.null(elevation)) elev <- elev * x1x2mask
			
		# }
			
		# # weighted longitude/latitude... used for centroid calculations for velocities
		# if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {

			# # weight longitude/latitude
			# x1weightedLongs <- longitude * x1
			# x1weightedLats <- latitude * x1

			# # centroid
			# x1sum <- cellStats(x1, 'sum')
			# x1centroidLong <- cellStats(x1weightedLongs, 'sum') / x1sum
			# x1centroidLat <- cellStats(x1weightedLats, 'sum') / x1sum
			
		# }
		
		### pre-calculations
		
		# for faster reference in the ."interpolateLongFromMatrix" and ".interpolateLatFromMatrix" functions
		if (any(c('nsQuants') %in% metrics)) {
			latVect <- rev(latitude[ , 1])
		}
		if (any(c('ewQuants') %in% metrics)) {
			longVect <- longitude[1, ]
		}
		
		if (!is.null(elevation) & ('elevCentroid' %in% metric | 'elevQuants' %in% metric)) elevVect <- raster::getValues(elev)
		
	### calculate velocities
	########################
	
		if (cores > 1) cores <- min(cores, parallel::detectCores())
		
		### multi-core
		##############

		# strategy: divide the key indices atTimes and indicesFrom into sets that can be done by one core
		# then feed bioticVelocity() a subset of rasters corresponding to these indices
		
		if (cores > 1) {
		
			# divvy up time periods among cores
			repAtTimes <- repIndicesFrom <- list()
			repSize <- floor(length(atTimes) / cores)
			numReps <- floor(length(atTimes) / repSize)
			for (i in 1:numReps) {
				
				extra <- ifelse(i > 1, 1, 0)
				repAtTimes[[i]] <- atTimes[(1 + (i - 1) * repSize - extra):(i * repSize)]
				repIndicesFrom[[i]] <- (1 + (i - 1) * repSize - extra):(i * repSize)
			
			}
			
			# add times rounded out
			if (tail(repAtTimes[[length(repAtTimes)]], 1) < tail(atTimes, 1)) {
				repAtTimes[[numReps + 1]] <- atTimes[(i * repSize)]:tail(atTimes, 1)
				repIndicesFrom[[numReps + 1]] <- (i * repSize):length(atTimes)
			}

			# multi-core
			if (cores > 1) {
				`%makeWork%` <- foreach::`%dopar%`
				cl <- parallel::makeCluster(cores)
				doParallel::registerDoParallel(cl)
			} else {
				`%makeWork%` <- foreach::`%do%`
			}
			
			mcOptions <- list(preschedule=TRUE, set.seed=FALSE, silent=FALSE)
			
			export <- c('bioticVelocity', '.euclid', '.cardinalDistance', '.interpolateLatFromMatrix', '.interpolateLongFromMatrix')
			
			out <- foreach::foreach(
				i=seq_along(repAtTimes),
				.options.multicore=mcOptions,
				.combine='rbind',
				.inorder=FALSE,
				.export=export,
				.packages = c('raster'),
				.verbose=FALSE) %makeWork%
				bioticVelocity(
					x = raster::subset(x, repIndicesFrom[[i]]),
					times = repAtTimes[[i]],
					atTimes = repAtTimes[[i]],
					longitude = longitude,
					latitude = latitude,
					elevation = elevation,
					metrics = metrics,
					quants = quants,
					onlyInSharedCells = onlyInSharedCells,
					cores = 1,
					warn = FALSE,
					...
				)
					
			parallel::stopCluster(cl)
		
			out <- out[order(out$timeFrom), ]
			
		### single-core
		###############
		
		} else {
		
			# output: data frame with one column per metric
			out <- data.frame()
			
			indicesFrom <- 1:(length(atTimes) - 1)

			### by each time period
			for (indexFrom in indicesFrom) {
			
				### get start time/end period layers and correct for shared non-NA cells
				x1 <- x[[indexFrom]]
				x2 <- x[[indexFrom + 1]]
				if (!is.null(elevation)) elev <- elevation

				### time
				timeFrom <- atTimes[indexFrom]
				timeTo <- atTimes[indexFrom + 1]
				timeSpan <- timeTo - timeFrom
				
				### remember
				thisOut <- data.frame(
					timeFrom = timeFrom,
					timeTo = timeTo,
					timeSpan = timeSpan
				)

				# correction for shared non-NA cells with next time period
				if (onlyInSharedCells) {
				
					x1mask <- x1 * 0 + 1
					x2mask <- x2 * 0 + 1
					x1x2mask <- x1mask * x2mask

					x1 <- x1 * x1x2mask
					x2 <- x2 * x1x2mask
							
					if (!is.null(elevation)) elev <- elev * x1x2mask

				}
				
				# statistics about cells
				size <- nrow(x1) * ncol(x1)
				x1ones <- x1 * 0 + 1
				x2ones <- x2 * 0 + 1
				propSharedCellsNotNA <- cellStats(x1ones * x2ones, 'sum') / size
				timeFromPropNotNA <- cellStats(x1ones, 'sum')  / size
				timeToPropNotNA <- cellStats(x2ones, 'sum')  / size

				thisOut <- cbind(
					thisOut,
					data.frame(
						propSharedCellsNotNA = propSharedCellsNotNA,
						timeFromPropNotNA = timeFromPropNotNA,
						timeToPropNotNA = timeToPropNotNA
					),
					row.names=NULL
				)
					
				### weighted longitude/latitude... used for centroid calculations for velocities
				if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {

					# weight longitude/latitude
					x1weightedLongs <- longitude * x1
					x1weightedLats <- latitude * x1

					x2weightedLongs <- longitude * x2
					x2weightedLats <- latitude * x2
					
					# centroid
					x1sum <- cellStats(x1, 'sum')
					x1centroidLong <- cellStats(x1weightedLongs, 'sum') / x1sum
					x1centroidLat <- cellStats(x1weightedLats, 'sum') / x1sum
					
					x2sum <- cellStats(x2, 'sum')
					x2centroidLong <- cellStats(x2weightedLongs, 'sum') / x2sum
					x2centroidLat <- cellStats(x2weightedLats, 'sum') / x2sum
					
					if (!is.null(elevation)) {

						x1weightedElev <- elev * x1
						x2weightedElev <- elev * x2
						
						x1elev <- x1weightedElev / x1sum
						x2elev <- x2weightedElev / x2sum
						
					} else {
						x1elev <- x2elev <- 0
					}

				}
				
				### weighted centroid metric
				if ('centroid' %in% metrics) {

					metric <- .euclid(c(x1centroidLong, x1centroidLat, x1elev), c(x2centroidLong, x2centroidLat, x2elev))
					metricRate <- metric / timeSpan
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							centroidVelocity = metricRate,
							centroidLong = x2centroidLong,
							centroidLat = x2centroidLat
						),
						row.names=NULL
					)
					
				}

				### velocity of occupied cells NORTH of start
				if ('nCentroid' %in% metrics) {

					cardOut <- .cardinalDistance(
						direction='n',
						longOrLat=latitude,
						x1=x1,
						x2=x2,
						refLong=x1centroidLong,
						refLat=x1centroidLat,
						x1weightedLongs=x1weightedLongs,
						x1weightedLats=x1weightedLats,
						x2weightedLongs=x2weightedLongs,
						x2weightedLats=x2weightedLats,
						x1weightedElev=x1weightedElev,
						x2weightedElev=x2weightedElev
					)

					metric <- cardOut$distance
					metricRate <- metric / timeSpan
					
					abundance <- cardOut$abundance
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							nCentroidVelocity = metricRate,
							nCentroidAbund = abundance
						),
						row.names=NULL
					)
					
				}
					
				### velocity of occupied cells SOUTH of start
				if ('sCentroid' %in% metrics) {

					cardOut <- .cardinalDistance(
						direction='s',
						longOrLat=latitude,
						x1=x1,
						x2=x2,
						refLong=x1centroidLong,
						refLat=x1centroidLat,
						x1weightedLongs=x1weightedLongs,
						x1weightedLats=x1weightedLats,
						x2weightedLongs=x2weightedLongs,
						x2weightedLats=x2weightedLats,
						x1weightedElev=x1weightedElev,
						x2weightedElev=x2weightedElev
					)

					metric <- cardOut$distance
					metricRate <- metric / timeSpan
					
					abundance <- cardOut$abundance
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							sCentroidVelocity = metricRate,
							sCentroidAbund = abundance
						),
						row.names=NULL
					)
					
				}
					
				### velocity of occupied cells EAST of start
				if ('eCentroid' %in% metrics) {

					cardOut <- .cardinalDistance(
						direction='e',
						longOrLat=longitude,
						x1=x1,
						x2=x2,
						refLong=x1centroidLong,
						refLat=x1centroidLat,
						x1weightedLongs=x1weightedLongs,
						x1weightedLats=x1weightedLats,
						x2weightedLongs=x2weightedLongs,
						x2weightedLats=x2weightedLats,
						x1weightedElev=x1weightedElev,
						x2weightedElev=x2weightedElev
					)

					metric <- cardOut$distance
					metricRate <- metric / timeSpan
					
					abundance <- cardOut$abundance
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							eCentroidVelocity = metricRate,
							eCentroidAbund = abundance
						),
						row.names=NULL
					)
					
				}
					
				### velocity of occupied cells WEST of start
				if ('wCentroid' %in% metrics) {

					cardOut <- .cardinalDistance(
						direction='w',
						longOrLat=longitude,
						x1=x1,
						x2=x2,
						refLong=x1centroidLong,
						refLat=x1centroidLat,
						x1weightedLongs=x1weightedLongs,
						x1weightedLats=x1weightedLats,
						x2weightedLongs=x2weightedLongs,
						x2weightedLats=x2weightedLats,
						x1weightedElev=x1weightedElev,
						x2weightedElev=x2weightedElev
					)

					metric <- cardOut$distance
					metricRate <- metric / timeSpan
					
					abundance <- cardOut$abundance
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							wCentroidVelocity = metricRate,
							wCentroidAbund = abundance
						),
						row.names=NULL
					)
					
				}
					
				### weighted north/south velocity
				if ('nsCentroid' %in% metrics) {
				
					metricSgn <- sign(x2centroidLat - x1centroidLat)
					metric <- .euclid(c(x2centroidLat, x2elev), c(x1centroidLat, x1elev))
					metricRate <- metricSgn * metric / timeSpan
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							nsCentroidVelocity = metricRate,
							nsCentroidLat = x2centroidLat
						),
						row.names=NULL
					)
					
				}
					
				### weighted east/west velocity
				if ('ewCentroid' %in% metrics) {
				
					metricSgn <- sign(x2centroidLong - x1centroidLong)
					metric <- .euclid(c(x2centroidLong, x2elev), c(x1centroidLong, x1elev))
					metricRate <- metricSgn * metric / timeSpan
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							ewCentroidVelocity = metricRate,
							ewCentroidLong = x2centroidLong
						),
						row.names=NULL
					)
					
				}
					
				### north/south quantile velocities (entire range)
				if ('nsQuants' %in% metrics) {

					# match location of quantiles
					# if quantile falls between rows then linearly extrapolate latitude
					for (countQuant in seq_along(quants)) {

						thisQuant <- quants[countQuant]
					
						# latitude of this quantile
						x1lat <- .interpolateLatFromMatrix(prob=thisQuant, x=x1, latVect=latVect)
						x2lat <- .interpolateLatFromMatrix(prob=thisQuant, x=x2, latVect=latVect)

						metric <- x2lat - x1lat
						metricRate <- metric / timeSpan
						
						# remember
						thisQuantOut <- data.frame(
							nsQuantVelocity = metricRate,
							nsQuantLat = x2lat
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
					
				}
					
				### east/west quantile velocities (entire range)
				if ('ewQuants' %in% metrics) {
				
					# match location of quantiles
					for (countQuant in seq_along(quants)) {

						thisQuant <- quants[countQuant]
						
						# longitudes of this quantile
						x1long <- .interpolateLongFromMatrix(prob=thisQuant, x=x1, longVect=longVect)
						x2long <- .interpolateLongFromMatrix(prob=thisQuant, x=x2, longVect=longVect)
						
						metric <- x2long - x1long
						metricRate <- metric / timeSpan
						
						# remember
						thisQuantOut <- data.frame(
							ewQuantVelocity = metricRate,
							ewQuantLong = x2long
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
					
				}

				# elevational centroid
				if ('elevCentroid' %in% metrics) {

					### sum of x1 and of x2 if doing elevational metrics
					if (!exists('x1sum', inherits=FALSE) | !exists('x2sum', inherits=FALSE)) {
						x1sum <- cellStats(x1, 'sum')
						x2sum <- cellStats(x2, 'sum')
					}

					elevCentroid_timeFrom <- cellStats(elev * x1, 'sum') / x1sum
					elevCentroid_timeTo <- cellStats(elev * x2, 'sum') / x2sum
					
					metric <- elevCentroid_timeTo - elevCentroid_timeFrom
					metricRate <- metric / timeSpan
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							elevCentroidVelocity = metricRate,
							elevCentroidElev = elevCentroid_timeTo
						),
						row.names=NULL
					)
					
				}

				
				# elevational quantiles
				if ('elevQuants' %in% metrics) {

					# match location of quantiles
					for (countQuant in seq_along(quants)) {

						thisQuant <- quants[countQuant]
						
						# weighted elevations
						weightedElevFrom <- elevVect * c(raster::as.matrix(x1))
						weightedElevTo <- elevVect * c(raster::as.matrix(x2))
						
						# quantile
						weightedElevQuantFrom <- quantile(weightedElevFrom, thisQuant, na.rm=TRUE)
						weightedElevQuantTo <- quantile(weightedElevTo, thisQuant, na.rm=TRUE)
						
						# sort weighted elevations and elevations by this same order (necessary for "bracket()" function)
						weightedElevFromOrder <- order(weightedElevFrom)
						weightedElevToOrder <- order(weightedElevTo)
						
						weightedElevFrom <- weightedElevFrom[weightedElevFromOrder]
						weightedElevTo <- weightedElevTo[weightedElevToOrder]
						
						elevVectFromOrdered <- elevVect[weightedElevFromOrder]
						elevVectToOrdered <- elevVect[weightedElevToOrder]
						
						# find index of elevation(s) closest to this quantile value
						elevQuantIndexFrom <- omnibus::bracket(weightedElevQuantFrom, by=weightedElevFrom, index=TRUE, warn=FALSE)
						elevQuantIndexTo <- omnibus::bracket(weightedElevQuantTo, by=weightedElevTo, index=TRUE, warn=FALSE)
						
						# get the elevations
						if (length(elevQuantIndexFrom) == 1) {
							elevFrom <- elevVectFromOrdered[elevQuantIndexFrom]
						} else {
							elevFrom <- thisQuant * elevVectFromOrdered[elevQuantIndexFrom[1]] + (1 - thisQuant) * elevVectFromOrdered[elevQuantIndexFrom[2]]
						}
						
						if (length(elevQuantIndexTo) == 1) {
							elevTo <- elevVectToOrdered[elevQuantIndexTo]
						} else {
							elevTo <- thisQuant * elevVectToOrdered[elevQuantIndexTo[1]] + (1 - thisQuant) * elevVectToOrdered[elevQuantIndexTo[2]]
						}
						
						metric <- elevTo - elevFrom
						metricRate <- metric / timeSpan
						
						# remember
						thisQuantOut <- data.frame(
							elevQuantVelocity = metricRate,
							elevQuantElev = elevTo
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
					
				}

				
				### mean abundance
				if ('mean' %in% metrics) {
				
					metric <- cellStats(x2, 'mean')
					thisOut <- cbind(
						thisOut,
						data.frame(
							mean = metric
						),
						row.names=NULL
					)
					
				}

				### total abundance
				if ('sum' %in% metrics) {
				
					metric <- cellStats(x2, 'sum')
					thisOut <- cbind(
						thisOut,
						data.frame(
							sum = metric
						),
						row.names=NULL
					)
					
				}

				### abundance quantiles
				if ('quants' %in% metrics) {
				
					metric <- raster::quantile(x2, quants, na.rm=TRUE)
					
					# remember
					for (countQuant in seq_along(metric)) {
					
						thisQuantOut <- data.frame(
							quantile = metric[countQuant]
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
					
				}
				
				### prevalence
				if ('prevalence' %in% metrics) {
				
					metric <- cellStats(x2 > 0, 'sum') / cellStats(x2ones, 'sum')
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							prevalence = metric
						),
						row.names=NULL
					)
					
				}

				### similarities
				if ('similarity' %in% metrics) {
					
					x1sizeNotNa <- cellStats(x1ones, 'sum')
					x2sizeNotNa <- cellStats(x1ones, 'sum')
					n <- mean(c(x1sizeNotNa, x2sizeNotNa))
					
					x1x2sum <- x1 + x2
					x1x2diff <- x2 - x1
					x1x2absDiff <- abs(x1x2diff)
					x1x2prod <- x1 * x2
					
					simpleMeanDiff <- cellStats(x1x2diff, 'sum') / n
					absMeanDiff <- cellStats(x1x2absDiff, 'sum') / n
					godsoeEsp <- 1 - cellStats(2 * x1x2prod, 'sum') / cellStats(x1x2sum, 'sum')
					schoenersD <- 1 - cellStats(x1x2absDiff, 'sum') / n
					warrenI <- 1 - sqrt(cellStats((sqrt(x1) - sqrt(x2))^2, 'sum') / n)

					thisOut <- cbind(
						thisOut,
						data.frame(
							simpleMeanDiff = simpleMeanDiff,
							absMeanDiff = absMeanDiff,
							godsoeEsp = godsoeEsp,
							schoenerD = schoenersD,
							warrenI = warrenI
						),
						row.names=NULL
					)
					
				}
				
				### remember
				out <- rbind(out, thisOut)
					
			} # next time period
			
		} # if single-core

	out
	
}

#' Euclidean distance between a pair of points
#'
#' Euclidean distance between a pair of points or two points. Note that the output is unsigned if \code{x2} and \code{y2} are provided, but signed if not.
#' @param a Numeric vector from 1 to 3 elements long
#' @param b Numeric vector from 1 to 3 elements long
#' @keywords internal
.euclid <- compiler::cmpfun(function(a, b) {

	if (length(a) != length(b)) stop('Length of "a" must be same as length of "b".')
	sqrt(sum((a - b)^2))
	
})

#' Movement of occupied cells in a given direction of a fixed point
#' 
#' This function calculates the weighted distance moved by a mass represented by set of cells which fall north, south, east, or west of a given location (i.e., typically the centroid of the starting population). Values >0 confer movement to the north, south, east, or west of this location.
#' @param direction Any of: \code{'n'} (north), \code{'s'} (south), \code{'e'} (east), or \code{'w'} (west).
#' @param longOrLat Numeric matrix, latitude or longitudes. If \code{direction} is \code{'n'} or \code{'s'} this must be latitudes. If \code{direction} is \code{'e'} or \code{'w'} this must be longitudes.
#' @param x1 Matrix of weights in time 1 (i.e., population size).
#' @param x2 Matrix of weights in time 2 (i.e., population size).
#' @param refLong Numeric, longitude of reference point from which to partition the weights into a northern, southern, eastern, or western portion.
#' @param refLat Numeric, latitude of reference point.
#' @param x1weightedLongs Matrix of longitudes weighted (i.e., by population size, given by \code{x1}).
#' @param x1weightedLats Matrix of latitudes weighted (i.e., by population size, given by \code{x1}).
#' @param x2weightedLongs Matrix of longitudes weighted (i.e., by population size, given by \code{x2}).
#' @param x2weightedLats Matrix of latitudes weighted (i.e., by population size, given by \code{x2}).
#' @param x1weightedElev Matrix of elevations weighted by x1 or \code{NULL}.
#' @param x2weightedElev Matrix of elevations weighted by x2 or \code{NULL}.
#' @return a list object with distance moved and abundance of all cells north/south/east/west of reference point.
#' @keywords internal
.cardinalDistance <- function(
	direction,
	longOrLat,
	x1,
	x2,
	refLong,
	refLat,
	x1weightedLongs,
	x1weightedLats,
	x2weightedLongs,
	x2weightedLats,
	x1weightedElev = NULL,
	x2weightedElev = NULL
) {

	# mask out cells north/south/east/west of or at starting centroid
	maskCells <- if (direction == 'n') {
		longOrLat >= refLat
	} else if (direction == 's') {
		longOrLat < refLat
	} else if (direction == 'e') {
		longOrLat >= refLong
	} else if (direction == 'w') {
		longOrLat < refLong
	}

	# censor x2
	x2censored <- x2 * maskCells
	abundance <- cellStats(x2censored, 'sum')
	
	# centroid of uncensored part of distribution
	if (abundance == 0) {
		distance <- 0
	} else {

		# censored, weighted longs/lats for x1
		x1censored <- x1 * maskCells
		x1weightedLongsCensored <- x1weightedLongs * maskCells
		x1weightedLatsCensored <- x1weightedLats * maskCells

		# censored, weighted longs/lats for x2
		x2weightedLongsCensored <- x2weightedLongs * maskCells
		x2weightedLatsCensored <- x2weightedLats * maskCells

		# weights for x1, x2
		x1sumCens <- cellStats(x1censored, 'sum')
		x2sumCens <- cellStats(x2censored, 'sum')

		# weighted long/lats for x1
		x1centroidLongCensored <- cellStats(x1weightedLongsCensored, 'sum') / x1sumCens
		x1centroidLatCensored <- cellStats(x1weightedLatsCensored, 'sum') / x1sumCens

		# weighted long/lats for x2
		x2centroidLongCensored <- cellStats(x2weightedLongsCensored, 'sum') / x2sumCens
		x2centroidLatCensored <- cellStats(x2weightedLatsCensored, 'sum') / x2sumCens
		
		# censored, weighted elevation for x1
		if (!is.null(x1weightedElev)) {

			x1weightedElevCensored <- cellStats(x1weightedElev * maskCells, 'sum')
			x1weightedElevCensored <- x1weightedElevCensored / x1sumCens
		
		} else {
			x1weightedElevCensored <- 0
		}
		
		# censored, weighted elevation for x2
		if (!is.null(x2weightedElev)) {

			x2weightedElevCensored <- cellStats(x2weightedElev * maskCells, 'sum')
			x2weightedElevCensored <- x2weightedElevCensored / x2sumCens
		
		} else {
			x2weightedElevCensored <- 0
		}
		
		distance <- .euclid(c(x2centroidLongCensored, x2centroidLatCensored, x1weightedElevCensored), c(x1centroidLongCensored, x1centroidLatCensored, x2weightedElevCensored))
	
	}

	list(distance=distance, abundance=abundance)
	
}

#' Latitude of of a quantile of the geographic abundance distribution-
#'
#' This function returns the latitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species abundance distribution. If a quantile would occur somewhere between two rows, the latitude is linearly interpolated between the latitudes of the two rows bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param x Matrix of abundances.
#' @param latVect Vector of latitudes, one per row in \code{x}... it is assumed these are in "reverse" order so correspond to the bottom-most row, second-bottom row, etc.
#' @keywords internal
.interpolateLatFromMatrix <- compiler::cmpfun(
	function(prob, x, latVect) {

	# standardized, cumulative sums of rows starting at bottom of matrix
	xRowSum <- raster::rowSums(x, na.rm=TRUE)
	xRowSum <- rev(xRowSum)

	if (!is.null(x1weightedElev)) {
		x1weightedElevVect <- raster::rowSums(x1weightedElevVect)
		x2weightedElevVect <- rev(x2weightedElevVect)
		x1elevVect <- x1weightedElevVect / xRowSum
		x1elevVect <- c(0, x1elevVect)
	}

	if (!is.null(x2weightedElev)) {
		x2weightedElevVect <- raster::rowSums(x2weightedElevVect)
		x2weightedElevVect <- rev(x2weightedElevVect)
		x2elevVect <- x2weightedElevVect / xRowSum
		x2elevVect <- c(0, x2elevVect)
	}

	xRowSum <- c(0, xRowSum)
	xRowCumSum <- cumsum(xRowSum)
	xRowCumSumStd <- xRowCumSum / max(xRowCumSum)
	
	rowIndex <- which(prob == xRowCumSumStd)
	
	# exact latitude
	if (length(rowIndex) != 0) {
		lat <- latVect[rowIndex]
	# interpolate latitude
	} else {
		
		cellLength <- mean(latVect[2:length(latVect)] - latVect[1:(length(latVect) - 1)])
		latVect <- c(latVect[1] - cellLength, 0.5 * cellLength + latVect)
		
		# if all values are equally distant from the desired quantile
		diffs1 <- abs(prob - xRowCumSumStd)
		row1 <- if (sd(diffs1) == 0) {
			round(median(seq_along(diffs1)))
		} else {
			if (prob > 0.5) {
				which.min(diffs1)
			} else {
				length(diffs1) - which.min(rev(diffs1)) + 1
			}
		}
		if (is.na(row1)) return(NA)
		lat1 <- latVect[row1]
		rowsCumSumStd_row1NA <- xRowCumSumStd
		rowsCumSumStd_row1NA[row1] <- NA
		
		# if all values are equally distant from the desired quantile
		diffs2 <- abs(prob - rowsCumSumStd_row1NA)
		row2 <- if (sd(diffs2, na.rm=TRUE) == 0) {
			round(median(seq_along(diffs2)[-is.na(diffs2)]))
		} else {
			if (prob > 0.5) {
				which.min(diffs2)
			} else {
				length(diffs2) - which.min(rev(diffs2)) + 1
			}
		}
		if (is.na(row2)) return(NA)
		lat2 <- latVect[row2]
	
		latVect <- c(lat1, lat2)
		latVectOrder <- order(latVect)
		latVect <- latVect[latVectOrder]
		xRowSum <- xRowSum[c(row1, row2)]
		xRowSum <- xRowSum[latVectOrder]
		# note: creates problems if abundance is 0 in both latitudes
		lat <- if (all(xRowSum == 0)) {
			NA
		} else {
			((1 - prob) * xRowSum[1] * latVect[1] + prob * xRowSum[2] * latVect[2]) / sum(xRowSum * c(1 - prob, prob))
		}
		
	}
	
	lat
	
})

#' Longitude of of a quantile of the geographic abundance distribution-
#'
#' This function returns the longitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species abundance distribution. If a quantile would occur somewhere between two columns, the longitude is linearly interpolated between the latitudes of the two columns bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param x Matrix of abundances.
#' @param longVect Vector of longitudes, one per column in \code{x}.
#' @keywords internal
.interpolateLongFromMatrix <- compiler::cmpfun(
	function(prob, x, longVect) {

	# standardized, cumulative sums of cols starting at right side of matrix
	xColSum <- raster::colSums(x, na.rm=TRUE)
	xColSum <- c(xColSum, 0)
	xColCumSum <- cumsum(xColSum)
	xColCumSumStd <- xColCumSum / max(xColCumSum)
	
	xColCumSumStd <- c(0, xColCumSumStd)
	
	colIndex <- which(prob == xColCumSumStd)
	# longVect <- longitude[1, ]
	
	# exact longitude
	if (length(colIndex) != 0) {
		long <- longVect[colIndex]
	# interpolate longitude
	} else {
		
		cellLength <- mean(longVect[2:length(longVect)] - longVect[1:(length(longVect) - 1)])
		longVect <- c(longVect[1] - 0.5 * cellLength, longVect + 0.5 * cellLength)
		
		# if all values are equally distant from the desired quantile
		diffs1 <- abs(prob - xColCumSumStd)
		col1 <- if (sd(diffs1) == 0) {
			round(median(seq_along(diffs1)))
		} else {
			if (prob < 0.5) {
				which.min(diffs1)
			} else {
				length(diffs1) - which.min(rev(diffs1)) + 1
			}
		}
		if (is.na(col1)) return(NA)
		long1 <- longVect[col1]
		colsCumSumStd_col1NA <- xColCumSumStd
		colsCumSumStd_col1NA[col1] <- NA
		
		# if all values are equally distant from the desired quantile
		diffs2 <- abs(prob - colsCumSumStd_col1NA)
		col2 <- if (sd(diffs2, na.rm=TRUE) == 0) {
			round(median(seq_along(diffs2)[-is.na(diffs2)]))
		} else {
			if (prob < 0.5) {
				which.min(diffs2)
			} else {
				length(diffs2) - which.min(rev(diffs2)) + 1
			}
		}
		if (is.na(col2)) return(NA)
		long2 <- longVect[col2]
	
		longVect <- c(long1, long2)
		longVectOrder <- order(longVect)
		longVect <- longVect[longVectOrder]
		xColSum <- xColSum[c(col1, col2)]
		xColSum <- xColSum[longVectOrder]
		# note: creates problems if abundance is 0 in both longitudes
		long <- if (all(xColSum == 0)) {
			NA
		} else {
			((1 - prob) * xColSum[1] * longVect[1] + prob * xColSum[2] * longVect[2]) / sum(xColSum * c(1 - prob, prob))
		}
		
	}
	
	long
	
})
