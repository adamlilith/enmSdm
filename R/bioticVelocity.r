#' Velocity of shifts in densities across a series of rasters
#'
#' This function calculates several metrics of biotic velocity for a stack of rasters or an array representing a gridded landscape.
#' @param x Either a \code{RasterStack}, \code{RasterBrick}, or 3-dimensional array. Values should be either \code{NA} or >= 0.
#' \itemize{
#' 	\item If \code{x} is a \code{RasterStack} or \code{RasterBrick} then each layer is assumed to represent a time slice and the rasters \emph{must} be in an equal-area projection.
#' 	\item If \code{x} is an array then each "layer" in the third dimension is assumed to represent a map at a particular time slice in an equal-area projection. Note that if this is an array you should probably also specify the arguments \code{longitude} and \code{latitude}.
#' }
#' @param times Numeric vector with the same number of layers in \code{x} or \code{NULL} (default). This specifies the time represented by each layer in \code{x} from oldest (top layer) to most recent (bottom layer). Times \emph{must} appear in sequential order. For example, if time periods are 24 kybp, 23 kybp, 22 kybp, use \code{c(-24, -23, -22)}, not \code{c(24, 23, 22)}. If \code{NULL} (default), values are assigned starting at 1 and ending at the total number of layers in \code{x}.
#' @param atTimes Numeric, values of \code{times} across which to calculate biotic velocity. You can use this to calculate biotic velocities across selected time periods (e.g., just the first and last time periods). Note that \code{atTimes} must be the same as or a subset of \code{times}. The default is \code{NULL}, in which case velocity is calculated across all time slices (i.e., between \code{times} 1 and 2, 2 and 3, 3 and 4, etc.).
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
#' @param elevation Either \code{NULL} (default) or a raster or matrix representing elevation. If this is supplied, changes in elevation are incorporated into all velocity and speed metrics. Additionally, you can also calculate the metrics \code{elevCentrioid} and \code{elevQuants}.
#' @param metrics Biotic velocity metrics to calculate (default is to calculate them all). All metrics ignore \code{NA} cells in \code{x}. Here "starting time period" represents one layer in \code{x} and "end time period" the next layer.
#' \itemize{
#'  \item \code{summary}: This calculates a series of measures, none of which are measures of velocity:
#'		\itemize{
#'			\item Mean: Mean value across all cells.
#'			\item Sum: Total across all cells.
#'			\item Quantiles: \emph{N}th quantile values across all cells. Quantiles are given by \code{quants}.
#'			\item Prevalence: Number of cells with values > 0.
#'		}
#' 	\item \code{centroid}: Speed of mass-weighted centroid.
#'  \item \code{nsCentroid} or \code{ewCentroid}: Velocity in the north-south or east-west directions of the mass-weighted centroid. For north-south cardinality, positive values represent movement northward and negative southward.
#'  \item \code{nCentroid}, \code{sCentroid}, \code{eCentroid}, and \code{wCentroid}: Speed of mass-weighted centroid of the portion of the raster north/south/east/west of the landscape-wide weighted centroid of the starting time period.
#'  \item \code{nsQuants} or \code{ewQuants}: Velocity of the location of the \emph{N}th quantile of mass in the north-south or east-west directions. The quantiles can be specified in \code{quants}. For example, this could be the movement of the 5th, 50th, and 95th quantiles of population size going from south to north. The 0th quantile would measure the velocity of the southernmost or easternmost cell(s) with values >0, and the 100th quantile the northernmost or westernmost cell(s) with non-zero values.
#'  \item \code{similarity}: Several metrics of similarity between each time period. Some of these make sense only for cases where values in \code{x} are in the range [0, 1], but not if some values are outside this range. See \code{\link{compareNiches}} for more details. The metrics are:
#'		\itemize{
#'			\item Simple mean difference
#'			\item Mean absolute difference
#'			\item Root-mean squared difference
#'			\item Expected Fraction of Shared Presences or ESP (Godsoe 2014)
#'			\item D statistic (Schoener 1968)
#'			\item I statistic (Warren et al. 2008)
#'			\item Pearson correlation
#'			\item Spearman rank correlation
#'		}
#'  \item \code{elevCentroid}: Velocity of the centroid of mass in elevation (up or down). Argument \code{elevation} must be supplied.
#' 	\item \code{elevQuants}: Velocity of the emph{n}th quantile of mass in elevation (up or down). The quantiles to be evaluated are given by \code{quants}. The lowest elevation with mass >0 is the 0th quantile, and the highest elevation with mass >0 is the 100th. Argument \code{elevation} must be supplied.
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
#' }
#' Depending on \code{metrics} that are specified, additional fields are as follows. All measurements of velocity are in distance units (typically meters) per time unit (which is the same as the units used for \code{times} and \code{atTimes}). For example, if the rasters are in an Albers equal-area projection and \code{times} are in years, then the output will be meters per year.
#' \itemize{
#' 	\item If \code{metrics} contains \code{summary}:
#' 	\itemize{
#'		\item A column named \code{propSharedCellsNotNA}: Proportion of cells that are not \code{NA} in both the "from" and "to" time step.
#'		\item Columns named \code{timeFromPropNotNA} and \code{timeToPropNotNA}: Proportion of cells in the "from" time and "to" steps that are not \code{NA}.
#' 		\item A column named \code{mean}: Mean weight in "timeTo" time step. In the same units as the values of the cells.
#' 		\item Columns named \code{quantile_quant}\emph{N}: The \emph{N}th quantile(s) of weight in the "timeTo" time step. In the same units as the values of the cells.
#' 		\item A column named \code{prevalence}: Proportion of non-\code{NA} cells with weight >0 in the "timeTo" time step relative to all non-\code{NA} cells. Unitless.
#' }
#' 	\item If \code{metrics} has \code{'centroid'}: Columns named \code{centroidVelocity}, \code{centroidLong}, \code{centroidLat} -- Speed of weighted centroid, plus its longitude and latitude (in the "to" time period of each time step). Values are always >= 0.
#' 	\item If \code{metrics} has \code{'nsCentroid'}: Columns named \code{nsCentroid} and \code{nsCentroidLat} -- Velocity of weighted centroid in north-south direction, plus its latitude (in the "to" time period of each time step). Positive values connote movement north, and negative values south.
#' 	\item If \code{metrics} has \code{'ewControid'}: \code{ewCentroid} and \code{ewCentroidLong} -- Velocity of weighted centroid in east-west direction, plus its longitude (in the "to" time period of each time step).  Positive values connote movement east, and negative values west.
#' 	\item If \code{metrics} has \code{'nCentroid'}, \code{'sCentroid'}, \code{'eCentroid'}, and/or \code{'wCentroid'}: Columns named \code{nCentroidVelocity} and \code{nCentroidAbund}, \code{sCentroid} and \code{sCentroidAbund}, \code{eCentroid} and \code{eCentroidAbund}, and/or \code{wCentroid} and \code{wCentroidAbund} -- Speed of weighted centroid of all cells that fall north, south, east, or west of the landscape-wide centroid, plus a column indicating the weight (abundance) of all such populations. Values are always >= 0.
#' 	\item If \code{metrics} contains any of \code{nsQuants} or \code{ewQuants}: Columns named \code{nsQuantVelocity_quant}\emph{N} and \code{nsQuantLat_quant}\emph{N}, or \code{ewQuantVelocity_quant}\emph{N} and \code{ewQuantLat_quant}\emph{N}: Velocity of the \emph{N}th quantile weight in the north-south or east-west directions, plus the latitude or longitude thereof (in the "to" time period of each time step). Quantiles are cumulated starting from the south or the west, so the 0.05th quantile, for example, is in the far south or west of the range and the 0.95th in the far north or east. Positive values connote movement north or east, and negative values movement south or west.
#' \item If \code{metrics} contains \code{similarity}, metrics of similarity are calculated for each pair of successive landscapes, defined below as \code{x1} and \code{x2}, with the number of shared non-\code{NA} cells between them being \code{n}:
#'	\itemize{
#'		\item A column named \code{simpleMeanDiff}: \code{sum(x2 - x1, na.rm=TRUE) / n}
#'		\item A column named \code{meanAbsDiff}: \code{sum(abs(x2 - x1), na.rm=TRUE) / n}
#'		\item A column named \code{rmsd} (root-mean square difference): \code{sqrt(sum((x2 - x1)^2, na.rm=TRUE)) / n}
#'		\item A column named \code{godsoeEsp}: \code{1 - sum(2 * (x1 * x2), na.rm=TRUE) / sum(x1 + x2, na.rm=TRUE)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar.
#'		\item A column named \code{schoenersD}: \code{1 - (sum(abs(x1 - x2), na.rm=TRUE) / n)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar.
#'		\item A column named \code{warrensI}: \code{1 - sqrt(sum((sqrt(x1) - sqrt(x2))^2, na.rm=TRUE) / n)}, values of 1 ==> maximally similar, 0 ==> maximally dissimilar.
#'		\item A column named \code{cor}: Pearson correlation between values of \code{x1} and \code{x2}.
#'		\item A column named \code{rankCor}: Spearman rank correlation between values of \code{x1} and \code{x2}.
#'	}
#' 	\item If \code{metrics} contains \code{elevCentroid}: Columns named \code{elevCentroidVelocity} and \code{elevCentroidElev} -- Velocity of the centroid in elevation (up or down) and the elevation in the "to" timestep. Positive values of velocity connote movement upward, and negative values downward.
#' 	\item If \code{metrics} contains \code{elevQuants}: Columns named \code{elevQuantVelocity_quant}\emph{N} and \code{elevQuantVelocityElev_quant}\emph{N} -- Velocity of the \emph{N}th quantile of mass in elevation (up or down) and the elevation of this quantile in the "to" timestep. Positive values of velocity connote movement upward, and negative values downward.
#' }
#' @details
#' \emph{Attention:}  
#'   
#' This function may yield erroneous velocities if the region of interest is near or spans a pole or the international date line. Results using the "Quant" and "quant" metrics may be somewhat counterintuitive if just one cell is >0, or one row or column has the same values with all other values equal to 0 or \code{NA} because defining quantiles in these situations is not intuitive. Results may also be counterintuitive if some cells have negative values because they can "push" a centroid away from what would seem to be the center of mass as assessed by visual examination of a map.  
#'   
#' \emph{Note:}  
#'   
#' For the \code{nsQuants} and \code{ewQuants} metrics it is assumed that the latitude/longitude assigned to a cell is at its exact center. If a desired quantile does not fall exactly on the cell center, it is interpolated linearly between the rows/columns of cells that bracket the given quantile. For quantiles that fall south/westward of the first row/column of cells, the cell border is assumed to be at 0.5 * cell length south/west of the cell center.
#' @examples
#' \dontrun{
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
#' ### elevation
#' # first case: all elevations the same, abundance moves north
#' n <- 5
#' mat <- matrix(NA, nrow=n, ncol=n)
#' mat1 <- mat
#' mat1[c(1, 5), ] <- 0
#' mat1[c(2, 4), ] <- 1
#' mat1[c(3), ] <- 3
#' 
#' mat2 <- mat
#' mat2[c(4, 5), ] <- 0
#' mat2[c(1, 3), ] <- 1
#' mat2[c(2), ] <- 3
#' 
#' lat <- matrix(rep(n:1, n), nrow=n)
#' long <- matrix(rep(1:n, each=n), nrow=n)
#' 
#' mats <- array(c(mat1, mat2), dim=c(n, n, 2))
#' colnames(mats) <- paste0('long', long[1, ])
#' rownames(mats) <- paste0('lat', lat[ , 1])
#' 
#' elev1 <- mat
#' elev1[] <- 1
#' 
#' metrics <- c('centroid', 'elevCentroid')
#' bv1 <- bioticVelocity(x=mats, longitude=long, latitude=lat, metrics=metrics,
#' elevation=elev1)
#' 
#' # second case: elevations higher in "northern" portion, abundance moves north
#' elev2 <- mat
#' elev2[] <- 1
#' elev2[1:2, ] <- 2
#' bv2 <- bioticVelocity(x=mats, longitude=long, latitude=lat, metrics=metrics,
#' elevation=elev2)
#'
#' rbind(bv1, bv2)
#' }
#'
#' \dontrun{
#' ### multi-core
#' mats <- array(runif(100 * 100 * 1000), dim=c(100, 100, 500))
#' mats <- brick(mats)
#' projection(mats) <- getCRS('wgs84')
#'
#' # adding elevation increases velocity because it increases average
#' # distance between cells
#' elev <- matrix(rnorm(100 * 100, 0, 0.1), nrow=100)
#' elev <- raster(elev, crs=getCRS('wgs84'))
#'
#' mc1 <- bioticVelocity(x=mats, metrics='centroid', cores=4)
#' mc2 <- bioticVelocity(x=mats, elevation=elev, metrics='centroid', cores=4)
#' plot(mc1$centroidVelocity, mc2$centroidVelocity,
#' xlab='No elevation', ylab='Elevation')
#' abline(0, 1)
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
	metrics = c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentroid', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants', 'similarity', 'summary'),
	quants = c(0.05, 0.10, 0.5, 0.9, 0.95),
	onlyInSharedCells = FALSE,
	cores = 1,
	warn = TRUE,
	...
) {

	### debugging
	if (FALSE) {
		times <- NULL
		atTimes <- NULL
		longitude <- NULL
		latitude <- NULL
		elevation <- NULL
		metrics <- c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentroid', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants', 'summary', 'similarity')
		quants <- c(0.05, 0.10, 0.5, 0.9, 0.95)
		onlyInSharedCells <- FALSE
		cores <- 1
		warn <- TRUE
	}
	
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
					elevation <- raster::raster(elevation, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax)
				}
			}
				
		}
		
		# if (any(raster::minValue(x) < 0, na.rm=TRUE) & warn) warning('Negative values appear in "x". Output may be unreliable or undesirable.')
			
	### pre-calculations
	####################

		# for faster reference in ."interpCoordFromQuantile"
		if (any(c('nsQuants', 'nCentroid', 'sCentroid') %in% metrics)) latVect <- rev(latitude[ , 1])
		if (any(c('nsQuants', 'eCentroid', 'wCentroid') %in% metrics)) longVect <- longitude[1, ]
		
		if (!is.null(elevation) & any('elevCentroid' %in% metrics)) elevVect <- raster::getValues(elevation)
		
	### calculate velocities
	########################
	
		cores <- if (cores > 1 & length(atTimes) > 2) {
			min(cores, parallel::detectCores())
		} else {
			1
		}

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
			
			# add times excluded because of rounding
			lastRepAtTime <- utils::tail(repAtTimes[[length(repAtTimes)]], 1)
			lastAtTime <- utils::tail(atTimes, 1)
			if (lastRepAtTime < lastAtTime) {
				indicesRemainingAtTimes <- which(atTimes == lastRepAtTime):length(atTimes)
				repAtTimes[[numReps + 1]] <- atTimes[indicesRemainingAtTimes]
				repIndicesFrom[[numReps + 1]] <- indicesRemainingAtTimes
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
			
			export <- c('bioticVelocity', '.euclid', '.cardinalDistance', '.interpCoordFromQuantile', 'compareNiches')
			
			# parallel::clusterEvalQ(cl, .libPaths(.libPaths()))
			
			out <- foreach::foreach(
				i=seq_along(repAtTimes),
				.options.multicore=mcOptions,
				.combine='rbind',
				.inorder=FALSE,
				.export=export,
				.packages = c('raster', 'doParallel'),
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
					warn = FALSE#,
					# ...
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
				
				### mean abundance
				if ('summary' %in% metrics) {

					# shared cells and NA cells
					size <- nrow(x1) * ncol(x1)
					x1ones <- x1 * 0 + 1
					x2ones <- x2 * 0 + 1
					x1x2ones <- x1ones * x2ones
					numSharedNonNaCells <- raster::cellStats(x1x2ones, 'sum')
					propSharedCellsNotNA <- numSharedNonNaCells / size
					timeFromPropNotNA <- raster::cellStats(x1ones, 'sum')  / size
					timeToPropNotNA <- raster::cellStats(x2ones, 'sum')  / size

					thisOut <- cbind(
						thisOut,
						data.frame(
							propSharedCellsNotNA = propSharedCellsNotNA,
							timeFromPropNotNA = timeFromPropNotNA,
							timeToPropNotNA = timeToPropNotNA
						),
						row.names=NULL
					)
					
					# mean
					metric <- raster::cellStats(x2, 'mean')
					thisOut <- cbind(
						thisOut,
						data.frame(
							mean = metric
						),
						row.names=NULL
					)
				
					# sum
					metric <- raster::cellStats(x2, 'sum')
					thisOut <- cbind(
						thisOut,
						data.frame(
							sum = metric
						),
						row.names=NULL
					)

					# quantiles
					metric <- raster::quantile(x2, quants, na.rm=TRUE)
					
					for (countQuant in seq_along(metric)) {
					
						thisQuantOut <- data.frame(
							quantile = metric[countQuant]
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}

					# prevalence
					metric <- raster::cellStats(x2 > 0, 'sum') / raster::cellStats(x2ones, 'sum')
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							prevalence = metric
						),
						row.names=NULL
					)
					
				}

				### pre-calculate weighted longitude/latitude... used for centroid calculations for velocities
				if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants') %in% metrics)) {

					# weight longitude/latitude
					x1weightedLongs <- longitude * x1
					x1weightedLats <- latitude * x1

					x2weightedLongs <- longitude * x2
					x2weightedLats <- latitude * x2
					
					# centroid
					x1sum <- raster::cellStats(x1, 'sum')
					x1centroidLong <- raster::cellStats(x1weightedLongs, 'sum') / x1sum
					x1centroidLat <- raster::cellStats(x1weightedLats, 'sum') / x1sum
					
					x2sum <- raster::cellStats(x2, 'sum')
					x2centroidLong <- raster::cellStats(x2weightedLongs, 'sum') / x2sum
					x2centroidLat <- raster::cellStats(x2weightedLats, 'sum') / x2sum
					
					if (!is.null(elevation)) {

						x1weightedElev <- elev * x1
						x2weightedElev <- elev * x2
						
						x1centroidElev <- raster::cellStats(x1weightedElev, 'sum') / x1sum
						x2centroidElev <- raster::cellStats(x2weightedElev, 'sum') / x2sum
						
					} else {
						x1weightedElev <- x2weightedElev <- NULL
						x1centroidElev <- x2centroidElev <- NA
					}

				}
				
				### weighted centroid metric
				if ('centroid' %in% metrics) {

					metric <- .euclid(c(x1centroidLong, x1centroidLat, x1centroidElev), c(x2centroidLong, x2centroidLat, x2centroidElev))
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
						coordVect=latVect,
						x1=x1,
						x2=x2,
						refCoord=x1centroidLat,
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
						coordVect=latVect,
						x1=x1,
						x2=x2,
						refCoord=x1centroidLat,
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
						coordVect=longVect,
						x1=x1,
						x2=x2,
						refCoord=x1centroidLong,
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
						coordVect=longVect,
						x1=x1,
						x2=x2,
						refCoord=x1centroidLong,
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

					metricSign <- sign(x2centroidLat - x1centroidLat)
					metric <- .euclid(c(x2centroidLat, x2centroidElev), c(x1centroidLat, x1centroidElev))
					metricRate <- metricSign * metric / timeSpan
					
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
				
					metricSign <- sign(x2centroidLong - x1centroidLong)
					metric <- .euclid(c(x2centroidLong, x2centroidElev), c(x1centroidLong, x1centroidElev))
					metricRate <- metricSign * metric / timeSpan
					
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

					# latitude of all quantiles
					x1loc <- .interpCoordFromQuantile(latOrLong='latitude', quants=quants, x=x1, coordVect=latVect, weightedElev=x1weightedElev, warn=warn)
					x2loc <- .interpCoordFromQuantile(latOrLong='latitude', quants=quants, x=x2, coordVect=latVect, weightedElev=x2weightedElev, warn=warn)
					
					# calculate velocity and remember
					for (countQuant in seq_along(quants)) {
					
						metricSign <- sign(x2loc$coord[countQuant] - x1loc$coord[countQuant])
						metric <- .euclid(c(x1loc$coord[countQuant], x1loc$elev[countQuant]), c(x2loc$coord[countQuant], x2loc$elev[countQuant]))
						metric <- metricSign * metric
						metricRate <- metric / timeSpan
						
						# remember
						thisQuantOut <- data.frame(
							nsQuantVelocity = metricRate,
							nsQuantLat = x2loc$coord[countQuant]
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
				
				}
					
				### east/west quantile velocities (entire range)
				if ('ewQuants' %in% metrics) {
				
					# latitude of all quantiles
					x1loc <- .interpCoordFromQuantile(latOrLong='longitude', quants=quants, x=x1, coordVect=longVect, weightedElev=x1weightedElev, warn=warn)
					x2loc <- .interpCoordFromQuantile(latOrLong='longitude', quants=quants, x=x2, coordVect=longVect, weightedElev=x2weightedElev, warn=warn)
					
					# calculate velocity and remember
					for (countQuant in seq_along(quants)) {
					
						metricSign <- sign(x2loc$coord[countQuant] - x1loc$coord[countQuant])
						metric <- .euclid(c(x1loc$coord[countQuant], x1loc$elev[countQuant]), c(x2loc$coord[countQuant], x2loc$elev[countQuant]))
						metric <- metricSign * metric
						metricRate <- metric / timeSpan
						
						# remember
						thisQuantOut <- data.frame(
							ewQuantVelocity = metricRate,
							ewQuantLong = x2loc$coord[countQuant]
						)
						
						quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
						names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
						
						thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
					
					}
					
				}

				# elevational centroid
				if ('elevCentroid' %in% metrics) {

					# abundance-weighted elevation
					if (!exists('x1centroidElev', inherits=FALSE) | !exists('x2centroidElev', inherits=FALSE)) {
					
						x1weightedElev <- elev * x1
						x2weightedElev <- elev * x2
						
						x1sum <- raster::cellStats(x1, 'sum')
						x2sum <- raster::cellStats(x2, 'sum')
						
						x1centroidElev <- raster::cellStats(x1weightedElev, 'sum') / x1sum
						x2centroidElev <- raster::cellStats(x2weightedElev, 'sum') / x2sum
						
					}

					metric <- x2centroidElev - x1centroidElev
					metricRate <- metric / timeSpan
					
					thisOut <- cbind(
						thisOut,
						data.frame(
							elevCentroidVelocity = metricRate,
							elevCentroidElev = x2centroidElev
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
						weightedElevFrom <- elevVect * as.vector(x1)
						weightedElevTo <- elevVect * as.vector(x2)
						
						# quantile
						weightedElevQuantFrom <- stats::quantile(weightedElevFrom, thisQuant, na.rm=TRUE)
						weightedElevQuantTo <- stats::quantile(weightedElevTo, thisQuant, na.rm=TRUE)
						
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

				### similarities
				if ('similarity' %in% metrics) {

					x1vect <- c(raster::as.matrix(x1))
					x2vect <- c(raster::as.matrix(x2))
					x1x2sum <- x1 + x2
					x1x2diff <- x2 - x1
					x1x2absDiff <- abs(x1x2diff)
					x1x2prod <- x1 * x2

					if (!exists('naNonNaCells', inherits=FALSE)) {

						x1ones <- x1 * 0 + 1
						x2ones <- x2 * 0 + 1
						x1x2ones <- x1ones * x2ones
						naNonNaCells <- c(raster::as.matrix(x1x2ones))
					}

					sims <- compareNiches(x1vect, x2vect, w=naNonNaCells, na.rm=TRUE)

					thisOut <- cbind(
						thisOut,
						data.frame(
							simpleMeanDiff = sims[['meanDiff']],
							meanAbsDiff = sims[['meanAbsDiff']],
							rmsd = sims[['rmsd']],
							godsoeEsp = sims[['esp']],
							schoenerD = sims[['d']],
							warrenI = sims[['i']],
							cor = sims[['cor']],
							rankCor = sims[['rankCor']]
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
	sqrt(sum((a - b)^2, na.rm=TRUE))
	
})

#' Movement of occupied cells in a given direction of a fixed point
#' 
#' This function calculates the weighted distance moved by a mass represented by set of cells which fall north, south, east, or west of a given location (i.e., typically the centroid of the starting population). Values >0 confer movement to the north, south, east, or west of this location.
#' @param direction Any of: \code{'n'} (north), \code{'s'} (south), \code{'e'} (east), or \code{'w'} (west).
#' @param longOrLat Numeric matrix, latitude or longitudes. If \code{direction} is \code{'n'} or \code{'s'} this must be latitudes. If \code{direction} is \code{'e'} or \code{'w'} this must be longitudes.
#' @param coordVect Vector of latitude or longitude of cell centers, depending on value of \code{longOrLat}. If latitude, these \emph{must} go from south to north. If \code{longitude}, these \emph{must} go from west to east.
#' @param x1 Matrix of weights in time 1 (i.e., population size).
#' @param x2 Matrix of weights in time 2 (i.e., population size).
#' @param refCoord Numeric, latitude or longitude (depending on \code{longOrLat}) of reference point from which to partition the weights into a northern, southern, eastern, or western portion.
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
	coordVect,
	x1,
	x2,
	refCoord,
	x1weightedLongs,
	x1weightedLats,
	x2weightedLongs,
	x2weightedLats,
	x1weightedElev = NULL,
	x2weightedElev = NULL
) {
	
	### calculate cell weightings
	
	# weight row/column containing reference lat/long by 0.5
	cellLength <- mean(coordVect[2:length(coordVect)] - coordVect[1:(length(coordVect) - 1)])
	cellSides <- c(coordVect[1] - 0.5 * cellLength, 0.5 * cellLength + coordVect)
	maskCellRows <- omnibus::bracket(refCoord, by=cellSides, inner=TRUE, index=TRUE)
	coord1 <- cellSides[maskCellRows[1]]
	coord2 <- cellSides[maskCellRows[2]]
	maskCellsCenter <- 0.5 * (longOrLat >= coord1) * (longOrLat <= coord2)
		
	# weight rows/columns east/west/north/south of reference by 1
	maskCells <- if (direction == 'n') {
		longOrLat >= coord2
	} else if (direction == 's') {
		longOrLat < coord1
	} else if (direction == 'e') {
		longOrLat >= coord2
	} else if (direction == 'w') {
		longOrLat < coord1
	}
	
	# combine masks
	maskCells <- maskCells + maskCellsCenter

	### censor cell values
	
	# censor x2
	x2censored <- x2 * maskCells
	abundance <- raster::cellStats(x2censored, 'sum')
	
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

		# censored weights for x1, x2
		x1sumCens <- raster::cellStats(x1censored, 'sum')
		x2sumCens <- raster::cellStats(x2censored, 'sum')

		# weighted, censored centroid for x1
		x1centroidLongCensored <- raster::cellStats(x1weightedLongsCensored, 'sum') / x1sumCens
		x1centroidLatCensored <- raster::cellStats(x1weightedLatsCensored, 'sum') / x1sumCens

		# weighted, censored centroid for x2
		x2centroidLongCensored <- raster::cellStats(x2weightedLongsCensored, 'sum') / x2sumCens
		x2centroidLatCensored <- raster::cellStats(x2weightedLatsCensored, 'sum') / x2sumCens
		
		# censored, weighted elevation for x1
		if (!is.null(x1weightedElev)) {

			x1weightedElevCensored <- x1weightedElev * maskCells
			x1elevCensored <- raster::cellStats(x1weightedElevCensored, 'sum') / x1sumCens
		
		} else {
			x1elevCensored <- NA
		}
		
		# censored, weighted elevation for x2
		if (!is.null(x2weightedElev)) {

			x2weightedElevCensored <- x2weightedElev * maskCells
			x2elevCensored <- raster::cellStats(x2weightedElevCensored, 'sum') / x2sumCens
		
		} else {
			x2elevCensored <- NA
		}
		
		distance <- .euclid(c(x2centroidLongCensored, x2centroidLatCensored, x1elevCensored), c(x1centroidLongCensored, x1centroidLatCensored, x2elevCensored))
	
	}

	list(distance=distance, abundance=abundance)
	
}

#' Latitude of quantile(s) of the geographic abundance distribution
#'
#' This function returns the latitude or longitude of quantile(s) of the geographic abundance distribution. The input is derived from a rasterized map of the species abundance distribution. If a quantile would occur somewhere between two cells, the latitude/longitude is linearly interpolated between the cells bracketing its value.
#' @param latOrLong Either 'latitude' or 'longitude'
#' @param quants Quantile value(s) (i.e., in the range [0, 1])
#' @param x Matrix of abundances.
#' @param coordVect Vector of latitudes, one per row in \code{x} (from south to north!!!) **OR** vector or longitudes, one per column in \code{x} (from west to east!!!).
#' @param weightedElev Raster of elevations weighted by x1 or x2 or \code{NULL}.
#' @param warn Logical. Show warnings.
#' @keywords internal
.interpCoordFromQuantile <- compiler::cmpfun(
	function(latOrLong, quants, x, coordVect, weightedElev=NULL, warn=TRUE) {

	# PRELIMINARY: standardized, cumulative sums of cells starting at bottom/left of matrix
	xSumVect <- raster::rowSums(x, na.rm=TRUE)

	# elevation weighted by x1 and x2
	if (!is.null(weightedElev)) {
		
		weightedElevVect <- raster::rowSums(weightedElev)
		weightedElevVect <- weightedElevVect / xSumVect
		if (latOrLong == 'latitude') weightedElevVect <- rev(weightedElevVect)
		weightedElevVect <- c(NA, weightedElevVect) # adding NA because first row of lat vector is for southern-most edge of cell
		
	}

	if (latOrLong == 'latitude') xSumVect <- rev(xSumVect)
	xSumVect <- c(0, xSumVect)
	xCumulSumVect <- cumsum(xSumVect)
	xCumulSumVect <- xCumulSumVect / max(xCumulSumVect)
	
	# coordinates of cell edges starting at southernmost edge of southernmost cell... exception: the last latitude will be the northernmost edge of the northernmost cell
	cellLength <- mean(coordVect[2:length(coordVect)] - coordVect[1:(length(coordVect) - 1)])
	cellSides <- c(coordVect[1] - 0.5 * cellLength, 0.5 * cellLength + coordVect)
	
	# holds output, one row per quantile value
	quantOut <- data.frame()
	
	# cycle through quantiles
	for (quant in quants) {

		# all abundances are 0
		if (all(xSumVect == 0) | all(is.na(xSumVect))) {
			coord <- NA
			elev <- NA
		} else {

			# find the row(s) that bracket this quantile
			cells <- omnibus::bracket(quant, by=xCumulSumVect, index=TRUE, inner=TRUE, warn=warn)

			# get bracketing latitudes and abundances
			cell1 <- cells[1]
			coord1 <- cellSides[cell1]
			xCumulSum1 <- xCumulSumVect[cell1]
			
			if (length(cells) > 1) {
				cell2 <- cells[2]
				coord2 <- cellSides[cell2]
				xCumulSum2 <- xCumulSumVect[cell2]
			}

			# interpolate latitude
			coord <- if (length(cells) == 1) { # exact match
				coord1
			} else if (length(cells) == 2) { # bracketed
				coord1 + ((quant - xCumulSum1) / (xCumulSum2 - xCumulSum1)) * (coord2 - coord1)
			}

			# interpolate elevation
			if (!is.null(weightedElev)) {
				
				if (length(cells) == 1) { # exact match
					elev <- weightedElevVect[cell1]
				} else if (length(cells) == 2) { # bracketed
					
					elev1 <- weightedElevVect[cell1]
					elev2 <- weightedElevVect[cell2]
					
					elev <- elev1 + ((quant - xCumulSum1) / (xCumulSum2 - xCumulSum1)) * (elev2 - elev1)
				}
				
			} else {
				elev <- NA
			}
			
		} # if all abundances are 0

		# remember output
		quantOut <- rbind(
			quantOut,
			data.frame(
				quant=quant,
				coord=coord,
				elev=elev
			)
		)
		
	} # next quantile
		
	quantOut
	
})
