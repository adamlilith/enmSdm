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
#' }
#' @param quants Numeric vector indicating the quantiles at which biotic velocity is calculated for the "\code{quant}" and "\code{Quants}" metrics. Default is \code{c(0.05, 0.10, 0.5, 0.9, 0.95)}.
#' @param onlyInSharedCells Logical, if \code{TRUE}, calculate biotic velocity using only those cells that are not \code{NA} in the start and end of each time period. This is useful for controlling for shifting land mass due to sea level rise, for example, when calculating biotic velocity for an ecosystem or a species. The default is \code{FALSE}.
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
#' \donttest{ 
#' ### COMPLEX EXAMPLES ###
#' ########################
#'
#' ### setup for all complex examples
#' ##################################
#'
#' library(raster)
#' 
#' # simulate species starting at center of map
#' # will be using a spatially Gaussian distribution
#' # that uses only longitude and latitude as predictors
#' gauss <- function(x1, x2, mu1=0, mu2=0, sigma1=1, sigma2=0, rho=0) {
#' 
#'		first <- ((x1 - mu1) / sigma1)^2
#'		prod <- ((2 * rho * (x1 - mu1) * (x2 - mu2)) / (sigma1 * sigma2))
#'		second <- ((x2 - mu2) / sigma2)^2
#'		denom <- 2 * (1 - rho^2)
#' 
#'		inside <- first - prod + second
#'		inside <- (-1 * inside) / denom
#' 
#'		expo <- exp(inside)
#'		expo
#' 
#' }
#' 
#' r <- raster()
#' r[] <- 1
#' mollweide <- enmSdm::getCRS('mollweide', TRUE)
#' r <- raster::projectRaster(r, crs=mollweide)
#' ll <- enmSdm::longLatRasters(r, m=TRUE)
#' 
#' # simulated population trajectory:
#' # time 0 to 1: no change
#' # time 1 to 2: size doubles, no shift
#' # time 2 to 3: size halves, no shift
#' # time 3 to 4: size same, moves north
#' # time 4 to 5: size same, moves south
#' # time 5 to 6: size same, moves east
#' # time 6 to 7: size same, moves west
#' # times 7 to 8: size same, moves northwest
#' # times 8 to 10: size same, moves northwest, double time
#' # times 10 to 11: size doubles, moves northwest
#' # times 11 to 12: same size, moves northwest
#' x1 <- ll[['longitude']]
#' x2 <- ll[['latitude']]
#' 
#' s1 <- 2000000
#' s2 <- 1000000
#' shift <- 1500000
#' 
#' species_t0 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t1 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t2 <- gauss(x1=x1, x2=x2, sigma1=2 * s1, sigma2=2 * s2)
#' species_t3 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t4 <- gauss(x1=x1, x2=x2, mu2=shift, sigma1=s1, sigma2=s2)
#' species_t5 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t6 <- gauss(x1=x1, x2=x2, mu1=shift, sigma1=s1, sigma2=s2)
#' species_t7 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t8 <- gauss(x1=x1, x2=x2, mu1=-1 * shift, mu2=1 * shift, sigma1=s1, sigma2=s2)
#' species_t10 <- gauss(x1=x1, x2=x2, mu1=-2 * shift, mu2=2 * shift, sigma1=s1, sigma2=s2)
#' species_t11 <- gauss(x1=x1, x2=x2, mu1=-3 * shift, mu2=3 * shift, sigma1=2 * s1, sigma2=s2)
#' species_t12 <- gauss(x1=x1, x2=x2, mu1=-4 * shift, mu2=4 * shift, sigma1=2 * s1, sigma2=s2,
#' 	rho=-0.5)
#' 	
#' rasts <- raster::stack(species_t0, species_t1, species_t2, species_t3,
#' species_t4, species_t5, species_t6, species_t7, species_t8,
#' 	species_t10, species_t11, species_t12)
#' times <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12)
#' names(rasts) <- paste0('t', times)
#' plot(rasts)
#' 
#' ### example with stationary land mass
#' #####################################
#' 
#' # across just start and end time periods
#' (vels <- bioticVelocity(rasts, times=times, atTimes=c(0, 12),
#' metrics='centroid'))
#' 
#' # across each time interval
#' (vels <- bioticVelocity(rasts, times=times, metrics='centroid'))
#' 
#' # plot movement of centroid
#' weightedLong <- species_t0 * ll[['longitude']]
#' weightedLat <- species_t0 * ll[['latitude']]
#' startLong <- cellStats(weightedLong, 'sum') / cellStats(species_t0, 'sum')
#' startLat <- cellStats(weightedLat, 'sum') / cellStats(species_t0, 'sum')
#' 
#' plot(species_t0)
#' 
#' for (i in 1:(nrow(vels) - 1)) {
#' 
#' x1pos <- vels$centroidLong[i]
#' y1pos <- vels$centroidLat[i]
#' x2pos <- vels$centroidLong[i + 1]
#' y2pos <- vels$centroidLat[i + 1]
#' 
#' 	move <- sqrt((x1pos - x2pos)^2 + (y1pos - y2pos)^2)
#' if (move > 10) arrows(x1pos, y1pos, x2pos, y2pos, angle=15, length=0.05)
#' 
#' }
#' 
#' # all metrics
#' (vels <- bioticVelocity(x=rasts, times=times))
#' 
#' ### example with shifting land mass
#' ###################################
#' 
#' species_t0 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' rasts <- stack(species_t0, species_t0, species_t0, species_t0, species_t0)
#' # simulated population trajectory:
#' # time 0 to 1: species no change, lose eastward land
#' # time 1 to 2: species no change, lose southward land
#' # time 2 to 3: species no change, gain southward land
#' # time 3 to 4: species no change, reset land to t0
#' 
#' # create masks
#' ext <- extent(species_t0)
#' 
#' long <- ext@xmin + 0.55 * (ext@xmax - ext@xmin)
#' longMask <- extent(long, ext@xmax, ext@ymin, ext@ymax)
#' longMask <- as(longMask, 'SpatialPolygons')
#' projection(longMask) <- projection(species_t0)
#' 
#' lat <- ext@ymin + 0.45 * (ext@ymax - ext@ymin)
#' latMask <- extent(ext@xmin, ext@xmax, lat, ext@ymax)
#' latMask <- as(latMask, 'SpatialPolygons')
#' projection(latMask) <- projection(species_t0)
#' 
#' species_t1s <- mask(species_t0, longMask, inverse=TRUE) # lose eastward
#' species_t2s <- mask(species_t1s, latMask) # lose southward
#' species_t3s <- species_t1s # gain southward
#' species_t4s <- species_t0 # gain eastward
#' 
#' rasts <- stack(species_t0, species_t1s, species_t2s, species_t3s,
#' species_t4s)
#' 
#' plot(rasts)
#' 
#' # is this what you want? could be considered an artifact, could
#' # be considered "true" biotic velocity... range shift is extrinsic
#' (vels <- bioticVelocity(x=rasts, times=0:4, metrics=c('centroid', 'mean')))
#' 
#' # calculate velocity only for cells shared by start and end rasters
#' # in each time step
#' (vels <- bioticVelocity(x=rasts, times=0:4, metrics=c('centroid', 'mean'),
#' onlyInSharedCells=TRUE))
#' }
#'
#' @export
bioticVelocity <- function(
	x,
	times = NULL,
	atTimes = NULL,
	longitude = NULL,
	latitude = NULL,
	metrics = c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentroid', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants', 'mean', 'sum', 'quants', 'prevalence', 'similarity'),
	quants = c(0.05, 0.10, 0.5, 0.9, 0.95),
	onlyInSharedCells = FALSE,
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
	
		if (!all(order(times) == seq_along(times)) & warn) {
			warning('Times assigned to each period are not sequential (e.g., {1, 2, 3}\nversus {3, 2, 1}). Velocities may have incorrect signs.')
		}
	
		# times across which to calculate velocity
		if (is.null(atTimes)) atTimes <- times
		atTimes <- sort(atTimes)

		### indices of layers to use (discard others)
		atIndices <- which(times %in% atTimes)
	
	### catch errors
	################
		
		if (length(times) != totalTimes) stop('The length of "times" does not match the total number of time periods represented by "x".')
		if (!all(atTimes %in% times)) stop ('All time slices specified in "atTimes" must also appear in "times".')

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
				
		}
		
		if (any(raster::minValue(x) < 0, na.rm=TRUE) & warn) warning('Negative values appear in "x". Output may be unreliable or undesirable.')
			
	### calculate weighted longitude and latitudes for starting time period
	#######################################################################

		x1 <- x[[1]]

		# correction for shared non-NA cells with next time period
		if (onlyInSharedCells) {

			x2 <- x[[2]]
			x1mask <- x1 * 0 + 1
			x2mask <- x2 * 0 + 1
			x1x2mask <- x1mask * x2mask
			x1 <- x1 * x1x2mask
			
		}
			
		# weighted longitude/latitude... used for centroid calculations for velocities
		if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {

			# weight longitude/latitude
			x1weightedLongs <- longitude * x1
			x1weightedLats <- latitude * x1

			# centroid
			x1sum <- cellStats(x1, 'sum')
			x1centroidLong <- cellStats(x1weightedLongs, 'sum') / x1sum
			x1centroidLat <- cellStats(x1weightedLats, 'sum') / x1sum
			
		}
		
		# for faster reference in the ."interpolateLongFromMatrix" and ".interpolateLatFromMatrix" functions
		if (any(c('nsQuants') %in% metrics)) {
			latVect <- rev(latitude[ , 1])
		}
		if (any(c('ewQuants') %in% metrics)) {
			longVect <- longitude[1, ]
		}
		
	### calculate velocities
	########################
	
		# output: data frame with one column per metric
		out <- data.frame()
		
		indicesFrom <- 1:(length(atTimes) - 1)

		### by each time period
		for (indexFrom in indicesFrom) {
		
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

			### get start time/end period layers and correct for shared non-NA cells
			x1 <- x[[indexFrom]]
			x2 <- x[[indexFrom + 1]]

			# correction for shared non-NA cells with next time period
			if (onlyInSharedCells) {
			
				x1mask <- x1 * 0 + 1
				x2mask <- x2 * 0 + 1
				x1x2mask <- x1mask * x2mask

				x1 <- x1 * x1x2mask
				x2 <- x2 * x1x2mask
						
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
				
			}

			### weighted centroid metric
			if ('centroid' %in% metrics) {

				metric <- .euclid(x1centroidLong, x1centroidLat, x2centroidLong, x2centroidLat)
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
					x1=x1,
					x2=x2,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
					x1weightedLongs=x1weightedLongs,
					x1weightedLats=x1weightedLats,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					longOrLat=latitude
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
					x1=x1,
					x2=x2,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
					x1weightedLongs=x1weightedLongs,
					x1weightedLats=x1weightedLats,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					longOrLat=latitude
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
					x1=x1,
					x2=x2,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
					x1weightedLongs=x1weightedLongs,
					x1weightedLats=x1weightedLats,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					longOrLat=longitude
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
					x1=x1,
					x2=x2,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
					x1weightedLongs=x1weightedLongs,
					x1weightedLats=x1weightedLats,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					longOrLat=longitude
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
			
				metric <- .euclid(x2centroidLat, x1centroidLat)
				metricRate <- metric / timeSpan ### metricRate <- -1 * metric / timeSpan
				
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
			
				metric <- .euclid(x2centroidLong, x1centroidLong)
				metricRate <- metric / timeSpan
				
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

					metric <- .euclid(x2lat, x1lat)
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
					
					metric <- .euclid(x2long, x1long)
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

	out
	
}

#' Euclidean distance between a pair of points
#'
#' Euclidean distance between a pair of points or two points. Note that the output is unsigned if \code{x2} and \code{y2} are provided, but signed if not.
#' @param x1 Numeric
#' @param y1 Numeric
#' @param x2 Numeric
#' @param y2 Numeric
#' @details If \code{x2} and \code{y2} are \code{NULL} then the output is simply \code{x1 - y2}.
#' @keywords internal
.euclid <- compiler::cmpfun(function(x1, y1, x2=NULL, y2=NULL) {

	if (is.null(x2) & is.null(y2)) {
		x1 - y1
	} else if (!is.null(x2) & !is.null(y2)) {
		sqrt((x1 - x2)^2 + (y1 - y2)^2)
	} else {
		stop('Please specify arguments "x1" and "x2" OR "x1", "x2", "y1", and "y2".')
	}
	
})

#' Movement of occupied cells in a given direction of a fixed point
#' 
#' This function calculates the weighted distance moved by a mass represented by set of cells which fall north, south, east, or west of a given location (i.e., typically the centroid of the starting population). Values >0 confer movement to the north, south, east, or west of this location. #' @param direction Any of: \code{'n'} (north), \code{'s'} (south), \code{'e'} (east), or \code{'w'} (west).
#' @param refLong Numeric, longitude of reference point from which to partition the weights into a northern, southern, eastern, or western portion.
#' @param refLat Numeric, latitude of reference point.
#' @param x1 Matrix of weights in time 1 (i.e., population size).
#' @param x2 Matrix of weights in time 2 (i.e., population size).
#' @param x1weightedLongs Matrix of longitudes weighted (i.e., by population size, given by \code{x1}).
#' @param x1weightedLats Matrix of latitudes weighted (i.e., by population size, given by \code{x1}).
#' @param x2weightedLongs Matrix of longitudes weighted (i.e., by population size, given by \code{x2}).
#' @param x2weightedLats Matrix of latitudes weighted (i.e., by population size, given by \code{x2}).
#' @param longOrLat Numeric matrix, latitude or longitudes. If \code{direction} is \code{'n'} or \code{'s'} this must be latitudes. If \code{direction} is \code{'e'} or \code{'w'} this must be longitudes.
#' @return a list object with distance moved and abundance of all cells north/south/east/west of reference point.
#' @keywords internal
.cardinalDistance <- function(
	direction,
	refLong,
	refLat,
	x1,
	x2,
	x1weightedLongs,
	x1weightedLats,
	x2weightedLongs,
	x2weightedLats,
	longOrLat
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
	
	x1censored <- x1 * maskCells
	x1weightedLongsCensored <- x1weightedLongs * maskCells
	x1weightedLatsCensored <- x1weightedLats * maskCells

	x2censored <- x2 * maskCells
	x2weightedLongsCensored <- x2weightedLongs * maskCells
	x2weightedLatsCensored <- x2weightedLats * maskCells

	abundance <- cellStats(x2censored, 'sum')
	
	# centroid of uncensored part of distribution
	if (abundance == 0) {
		distance <- 0
	} else {

		x1centroidLongCensored <- cellStats(x1weightedLongsCensored, 'sum') / cellStats(x1censored, 'sum')
		x1centroidLatCensored <- cellStats(x1weightedLatsCensored, 'sum') / cellStats(x1censored, 'sum')

		x2centroidLongCensored <- cellStats(x2weightedLongsCensored, 'sum') / cellStats(x2censored, 'sum')
		x2centroidLatCensored <- cellStats(x2weightedLatsCensored, 'sum') / cellStats(x2censored, 'sum')
		
		distance <- .euclid(x1=x2centroidLongCensored, y1=x2centroidLatCensored, x2=x1centroidLongCensored, y2=x1centroidLatCensored)
	
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
	xRowSum <- c(xRowSum, 0)
	xRowSum <- rev(xRowSum)
	xRowCumSum <- cumsum(xRowSum)
	xRowCumSumStd <- xRowCumSum / max(xRowCumSum)
	
	rowIndex <- which(prob == xRowCumSumStd)
	# latVect <- rev(latitude[ , 1])
	
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
