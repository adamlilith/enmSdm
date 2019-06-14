#' Shifts in densities across rasters
#'
#' This function calculates several metrics of biotic velocity for a stack of rasters, an array representing a gridded landscape, or a "pops" (list) object.
#' @param x Either a \code{RasterStack}, a 3-dimensional array, or list representing colonization histories (a "pops" object). Regardless of the class of \code{x}, it is assumed values representing the entity for which to calculate velocities are either \code{NA} or >= 0 (i.e., no negative population size).
#' \itemize{
#' 	\item If \code{x} is a RasterStack then each raster is assumed to represent a time slice and the rasters \emph{must} be in an equal-area projection.
#' 	\item If \code{x} is an array then each matrix in the third dimension is assumed to represent a map at a particular time slice in an equal-area projection. Note that if this is an array you can also specify the arguments \code{longitude} and \code{latitude}.
#'  \item If \code{x} is a list (a "pops" object), then it is assumed to have at minimum the field \code{x$Nvecs}. To ensure velocities are in spatial units (usually meters--as opposed to arbitrary units), it must have the fields \code{x$pophist$longitude} and \code{x$pophist$latitude} (if not, then it must have It must also have \code{x$pophist$row} and \code{x$pophist$col}). To ensure velocities are in units of time (vs arbitrary units) it also needs \code{x$pophist$time}.
#' }
#' @param times Numeric vector \emph{or} \code{NULL} (default). This specifies the time period represented by each layer in \code{x}. Times \emph{must} appear in sequential order. For example, if time periods are 24 kybp, 23 kybp, 22 kypb, use \code{c(-24, -23, -22)}, \emph{not} \code{c(24, 23, 22)}. If \code{NULL} (default), default values are assigned:
#' \itemize{
#'	\item If \code{x} is a \code{RasterStack} and \code{times} is \code{NULL}, then \code{times} will be assigned values starting at 1 and ending at the total number of rasters in \code{x}. Alternatively, you can supply a numeric vector with the same number of values as layers in \code{x}.
#'	\item If \code{x} is an array and \code{times} is \code{NULL}, then \code{times} will be assigned values starting at 1 and ending at \code{dim(x)[3]}. Alternatively, you can supply a numeric vector with the same number of values as layers in the third dimension of \code{x}.
#'	\item If \code{x} is a list (a "pops" object) and contains a field named \code{x$pophist$time} and \code{times} is \code{NULL}, then \code{.$pophist$time} will be used to assign times to each period. Alternatively, you can supply a numeric vector with the same number of values as columns in \code{x$Nvecs}.
#' }
#' @param atTimes Numeric, values of \code{times} across which to calculate biotic velocity. You can use this to calculate biotic velocities across selected time periods (e.g., just the first and last time periods). Note that \code{across} \emph{must} be the same as or a subset of \code{times}. The default is \code{NULL}, in which case biotic velocity is calculated across all time slices (i.e., between period 1 and 2, 2 and 3, 3 and 4, etc.).
#' @param longitude Numeric matrix or \code{NULL} (default).
#' \itemize{
#'	\item If \code{x} is \code{RasterStack} then this is ignored (longitude is ascertained directly from the rasters, which \emph{must} be in equal-area projection for velocities to be valid!).
#'	\item If \code{x} is an array and \code{longitude} is \code{NULL} (default), then longitude will be ascertained from column numbers in \code{x} and velocities will be in arbitrary spatial units (versus, for example, meters). Alternatively, this can be a two-dimensional matrix whose elements represent the longitude coordinates of the centers of cells of \code{x}. The matrix must have the same number of rows and columns as \code{x}. Note that the coordinates must be from an equal-area projection for results to be valid!
#' 	\item If \code{x} is a list (a "pops" object), the function uses values in \code{x$pophist$longitude}. If \code{x$pophist$longitude} does not exist, then the \code{x$pophist$col} is used as longitude (with a warning). Alternatively, you can supply a two-dimensional matrix with longitude values for each cell. Note that the coordinates must be from an equal-area projection for results to be valid!
#' }
#' @param latitude Numeric matrix or \code{NULL} (default).
#' \itemize{
#'	\item If \code{x} is \code{RasterStack} then this is ignored (latitude is ascertained directly from the rasters, which \emph{must} be in equal-area projection for velocities to be valid!).
#'	\item If \code{x} is an array and \code{latitude} is \code{NULL} (default), then latitude will be ascertained from row numbers in \code{x} and velocities will be in arbitrary spatial units (versus, for example, meters). Alternatively, this can be a two-dimensional matrix whose elements represent the latitude coordinates of the centers of cells of \code{x}. The matrix must have the same number of rows and columns as \code{x}. Note that the coordinates must be from an equal-area projection for results to be valid!
#' 	\item If \code{x} is a list (a "pops" object), the function uses values in \code{x$pophist$latitude}. If \code{x$pophist$latitude} does not exist, then the \code{x$pophist$row} is used as latitude (with a warning). Alternatively, you can supply a two-dimensional matrix with latitude values for each cell. Note that the coordinates must be from an equal-area projection for results to be valid!
#' }
#' @param metrics Biotic velocity metrics to calculate (default is to calculate them all). The units for all metrics will be the spatial units of the map (usually meters) divided by the temporal units (same units used for \code{times}---most often years). All metrics ignore \code{NA} cells in \code{x}.
#' \itemize{
#' 	\item \code{centroid}: Movement of mass-weighted centroid.
#'  \item \code{nsCentroid} or \code{ewCentroid}: Movement of mass-weighted centroid in the north-south or east-west directions. For north-south cardinality, positive values represent movement northward and negative southward. For east-west cardinality, positive values represent movement eastward and negative westward.
#'  \item \code{nCentroid}, \code{sCentroid}, \code{eCentroid}, and \code{wCentroid}: Movement of mass-weighted centroid in of the portion of the range that is north, south, east, or west of the mass-weighted centroid of the previous time period.
#'  \item \code{nsQuants} or \code{ewQuants}: Movement of the location of the \emph{x}th quantile of mass in the north-south or east-west directions. For example, this could be the movement of the 5th, 50th, and 95th quantiles of population size from south to north. The quantiles can be specified in \code{quants}.
#'  \item \code{mean}: Mean value across all cells (this is not really a measure of "velocity").
#'  \item \code{quants}: \emph{x}th quantile values across all cells (this is not really a measure of "velocity"). Quantiles are given by \code{quants}.
#'  \item \code{prevalence}: Number of cells with values > 0 (this is not really a measure of "velocity").
#' }
#' @param quants Numeric vector indicating the quantiles at which biotic velocity is calculated for the "quant" and "Quants" metrics. Default is \code{c(0.05, 0.10, 0.5, 0.9, 0.95)}.
#' @param warn Logical, if \code{TRUE} then display function-specific warnings.
#' @param ... Other arguments (not used).
#' @return A data frame with biotic velocities. Fields are as follows:
#' \itemize{
#' 	\item \code{timeFrom}: Start time of interval
#' 	\item \code{timeTo}: End time of interval
#' 	\item \code{timeSpan}: Duration of interval
#' }
#' Depending on \code{metrics} that are specified, additional fields are as follows:
#' \itemize{
#' 	\item \code{centroidVelocity}, \code{centroidLong}, \code{centroidLat}: Velocity of weighted centroid and its longitude and latitude.
#' 	\item \code{nsCentroid}, \code{nsCentroidLat}: Velocity of weighted centroid in north-south direction and its latitude.
#' 	\item \code{ewCentroid}, \code{ewCentroidLong}: Velocity of weighted centroid in east-west direction and its longitude.
#' 	\item \code{nCentroid} and \code{nCentroidAbund}, \code{sCentroid} and \code{sCentroidAbund}, \code{eCentroid} and \code{eCentroidAbund}, or \code{wCentroid} and \code{wCentroidAbund}: Velocity of weighted centroid of all cells with weight >0 that fall north, south, east, or west of centroid of initial population (population in first time step), plus weight of all such populations.
#' 	\item \code{nsCentroidVelocity_quantX} and \code{nsCentroidLat_quantX}, or \code{ewCentroidVelocity_quantX} and \code{ewCentroidLat_quantX}: Velocity of the \emph{x}th quantile weight in the north-south or east-west directions, plus the latitude or longitude thereof.
#' 	\item \code{mean}: Mean weight in "timeTo" time step.
#' 	\item \code{quantile_quantX}: The \emph{X}th quantile(s) of weight in the "timeTo" time step.
#' 	\item \code{prevalence}: Proportion of non-\code{NA} cells with weight >0 in the "timeTo" time step.
#' }
#' @details This function may yield erroneous velocities if the region of interest is near a pole and will yield erroneous results if the region spans the international date line. It converts rasters to arrays before doing calculations, so using very large rasters may yield slow performance and may not even work, depending on memory requirements.
#' @examples
#' # simulate species starting at center of map
#' # will be using a spatially Gaussian distribution
#' # that uses only longitude and latitude as predictors
#' gauss <- function(x1, x2, mu1=0, mu2=0, sigma1=1, sigma2=0, rho=0) {
#' 
#' 	first <- ((x1 - mu1) / sigma1)^2
#' 	prod <- ((2 * rho * (x1 - mu1) * (x2 - mu2)) / (sigma1 * sigma2))
#' 	second <- ((x2 - mu2) / sigma2)^2
#' 	denom <- 2 * (1 - rho^2)
#' 
#' 	inside <- first - prod + second
#' 	inside <- (-1 * inside) / denom
#' 	
#' 	expo <- exp(inside)
#' 	
#' 	expo
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
#' 	species_t4, species_t5, species_t6, species_t7, species_t8,
#' 	species_t10, species_t11, species_t12)
#' times <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12)
#' names(rasts) <- paste0('t', times)
#' plot(rasts)
#' 
#' # across just start and end time periods
#' (vels <- bioticVelocity(rasts, times=times, atTimes=c(0, 12),
#' 	metrics='centroid'))
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
#' 	x1 <- vels$centroidLong[i]
#' 	y1 <- vels$centroidLat[i]
#' 	x2 <- vels$centroidLong[i + 1]
#' 	y2 <- vels$centroidLat[i + 1]
#' 	
#' 	move <- sqrt((x1 - x2)^2 + (y1 - y2)^2)
#' 	if (move > 10) arrows(x1, y1, x2, y2, angle=15, length=0.05)
#' 
#' }
#' 
#' # all metrics
#' (vels <- bioticVelocity(x=rasts, times=times))
#' @export
bioticVelocity <- function(
	x,
	times = NULL,
	atTimes = NULL,
	longitude = NULL,
	latitude = NULL,
	metrics = c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentroid', 'eCentroid', 'wCentroid', 'nsQuants', 'ewQuants', 'mean', 'quants', 'prevalence'),
	quants = c(0.05, 0.10, 0.5, 0.9, 0.95),
	warn = TRUE,
	...
) {

	xClass <- class(x)

	### time of each period and times across which to calculate velocity
	####################################################################
	
		# total number of time periods
		totalTimes <- if ('array' %in% xClass) {
			dim(x)[3]
		} else if ('RasterStack' %in% xClass) {
			raster::nlayers(x)
		} else if ('list' %in% xClass) {
			ncol(x$Nvecs)
		}

		# time of each period
		if (is.null(times)) {
		
			# get from "pops" object
			if ('list' %in% xClass) {
				if (!is.null(x$pophist$time)) {
					times <- x$pophist$time
				} else {
					if (warn) warning('The population history object does not contain a field named ".$pophist$time". Velocities will be in arbitrary time units.')
					times <- 1:totalTimes
				}
			} else {
				times <- 1:totalTimes
			}
			
		}
	
		if (!all(order(times) == 1:totalTimes) & warn) {
			warning('Times assigned to each period are not sequential (e.g., {1, 2, 3} versus {3, 2, 1}). Velocities can have incorrect sign.')
		}
	
		# times across which to calculate velocity
		if (is.null(atTimes)) atTimes <- times
		atTimes <- sort(atTimes)

		### indices of layers to use (discard others)
		atIndices <- which(times %in% atTimes)
	
	### catch errors
	################
		
		if (length(times) != totalTimes) stop('The length of "times" does not match the total number of time periods represented by "x".')
		if (!all(atTimes %in% times)) stop ('All time slices specified in "atTimes" must also be in "times".')

	### convert input to array and get geographic information
	#########################################################

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
				
		}
		
		# convert from raster stack
		if ('RasterStack' %in% xClass) {

			ll <<- enmSdm::longLatRasters(x)
			longitude <- raster::as.matrix(ll[['longitude']])
			latitude <- raster::as.matrix(ll[['latitude']])
			x <- raster::subset(x, atIndices)
			x <- raster::as.array(x)

		}

		# convert from list ("pops" object)
		if ('list' %in% xClass) {
			
			nRows <- max(x$pophist$row)
			nCols <- max(x$pophist$col)
			xNew <- array(NA, dim=c(nRows, nCols, totalTimes))
			for (i in 1:totalTimes) xNew[, , i] <- matrix(x[['Nvecs']][ , i], nrow=nRows, ncol=nCols, byrow=TRUE)
			
			if (is.null(x$pophist$longitude)) {
				longitude <- matrix(x$pophist$col, nrow=nRows, ncol=nCols, byrow=TRUE)
				if (warn) warning('The population history object does not contain a field named ".$pophist$longitude". Velocities will be in arbitrary spatial units.')
			} else if (is.null(longitude)) {
				longitude <- matrix(x$pophist$longitude, nrow=nRows, ncol=nCols, byrow=TRUE)
			}
			
			if (is.null(x$pophist$latitude)) {
				# note: reversing row numbers for latitude so high values occur in the "north" and low in the "south" so they match trends in latitude
				latitude <- matrix(rev(x$pophist$row), nrow=nRows, ncol=nCols, byrow=TRUE)
				if (warn) warning('The population history object does not contain a field named ".$pophist$latitude". Velocities will be in arbitrary spatial units.')
			} else if (is.null(latitude)) {
				latitude <- matrix(x$pophist$latitude, nrow=nRows, ncol=nCols, byrow=TRUE)
			}
			
			# subset to time periods of interest
			x <- xNew[ , , atIndices]

		}

		if (any(x < 0, na.rm=TRUE) & warn) warning('Negative values appear in "x". Output may be unreliable or undesirable.')

	### calculate weighted longitude and latitudes for starting time period
	#######################################################################

		x1 <- x[ , , 1]

		# weighted longitude/latitude... used for centroid calculations for velocities
		if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {
			
			# mask longitude/latitude by NA cases in weights
			x1maskedLong <- longitude * !is.na(x1)
			x1maskedLat <- latitude * !is.na(x1)
			
			# weight longitude/latitude
			x1weightedLongs <- x1maskedLong * x1
			x1weightedLats <- x1maskedLat * x1

			x1centroidLong <- sum(x1weightedLongs, na.rm=TRUE) / sum(x1, na.rm=TRUE)
			x1centroidLat <- sum(x1weightedLats, na.rm=TRUE) / sum(x1, na.rm=TRUE)
			
			if (any(c('nCentroid', 'sCentroid', 'eCentroid', 'wCentroid') %in% metrics)) {
			
				x1sum <- sum(x1, na.rm=TRUE)
				x2centroidLong <- sum(x1weightedLongs, na.rm=TRUE) / x1sum
				x2centroidLat <- sum(x1weightedLats, na.rm=TRUE) / x1sum
			
			}
			
		}

	### calculate velocities
	########################
	
		# output: data frame with on column per metric
		out <- data.frame()
		
		indicesFrom <- 1:(length(atTimes) - 1)

		### by each time period
		for (indexFrom in indicesFrom) {
		
			### time
			timeFrom <- atTimes[indexFrom]
			timeTo <- atTimes[indexFrom + 1]
			timeSpan <- timeTo - timeFrom
			
			### remember
			timeSpan <- timeTo - timeFrom
			
			thisOut <- data.frame(
				timeFrom = timeFrom,
				timeTo = timeTo,
				timeSpan = timeSpan
			)
			
			### "to" layer
			x2 <- x[ , , indexFrom + 1]

			# weighted longitude/latitude... used for centroid calculations for velocities
			if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {

				# mask longitude/latitude by NA cases in weights
				x2maskedLong <- longitude * !is.na(x2)
				x2maskedLat <- latitude * !is.na(x2)
				
				# weight longitude/latitude
				x2weightedLongs <- x2maskedLong * x2
				x2weightedLats <- x2maskedLat * x2
						
				x2sum <- sum(x2, na.rm=TRUE)
				x2centroidLong <- sum(x2weightedLongs, na.rm=TRUE) / x2sum
				x2centroidLat <- sum(x2weightedLats, na.rm=TRUE) / x2sum
				
			}
			
			### weighted centroid
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

				cardOut <- .cardinalVelocity(
					direction='n',
					x2=x2,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
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

				cardOut <- .cardinalVelocity(
					direction='s',
					x2=x2,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
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

				cardOut <- .cardinalVelocity(
					direction='e',
					x2=x2,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
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

				cardOut <- .cardinalVelocity(
					direction='w',
					x2=x2,
					x2weightedLongs=x2weightedLongs,
					x2weightedLats=x2weightedLats,
					refLong=x1centroidLong,
					refLat=x1centroidLat,
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
				metricRate <- metric / timeSpan
				
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

				# standardized, cumulative sums of rows starting at bottom of matrix
				x1rowSum <- apply(x1, 1, sum, na.rm=TRUE)
				x2rowSum <- apply(x2, 1, sum, na.rm=TRUE)

				x1rowCumSum <- cumsum(rev(x1rowSum))
				x2rowCumSum <- cumsum(rev(x2rowSum))
				
				x1rowCumSumStd <- x1rowCumSum / max(x1rowCumSum)
				x2rowCumSumStd <- x2rowCumSum / max(x2rowCumSum)
							
				# match location of quantiles
				# if quantile falls between rows then linearly extrapolate latitude
				for (countQuant in seq_along(quants)) {

					thisQuant <- quants[countQuant]
				
					# latitude of this quantile
					x1lat <- .interpolateLatFromMatrix(prob=thisQuant, rowsCumSumStd=x1rowSum, latitude=latitude)
					x2lat <- .interpolateLatFromMatrix(prob=thisQuant, rowsCumSumStd=x2rowSum, latitude=latitude)

					metric <- .euclid(x2lat, x1lat)
					metricRate <- metric / timeSpan
					
					# remember
					thisQuantOut <- data.frame(
						nsCentroidVelocity = metricRate,
						nsCentroidLat = x2lat
					)
					
					quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
					names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
					
					thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
				
				}
				
			}
				
			### east/west quantile velocities (entire range)
			if ('ewQuants' %in% metrics) {
			
				# standardized, cumulative sums of cols starting at bottom of matrix
				x1colSum <- apply(x1, 2, sum, na.rm=TRUE)
				x2colSum <- apply(x2, 2, sum, na.rm=TRUE)

				x1colCumSum <- cumsum(x1colSum)
				x2colCumSum <- cumsum(x2colSum)
				
				x1colCumSumStd <- x1colCumSum / max(x1colCumSum)
				x2colCumSumStd <- x2colCumSum / max(x2colCumSum)
				
				# match location of quantiles
				for (countQuant in seq_along(quants)) {

					thisQuant <- quants[countQuant]
					
					# longitudes of this quantile
					x1long <- .interpolateLongFromMatrix(prob=thisQuant, colsCumSumStd=x1colSum, longitude=longitude)
					x2long <- .interpolateLongFromMatrix(prob=thisQuant, colsCumSumStd=x2colSum, longitude=longitude)
					
					metric <- .euclid(x2long, x1long)
					metricRate <- metric / timeSpan
					
					# remember
					thisQuantOut <- data.frame(
						ewCentroidVelocity = metricRate,
						ewCentroidLong = x2long
					)
					
					quantCharacter <- gsub(quants[countQuant], pattern='[.]', replacement='p')
					names(thisQuantOut) <- paste0(names(thisQuantOut), '_quant', quantCharacter)
					
					thisOut <- cbind(thisOut, thisQuantOut, row.names=NULL)
				
				}
				
			}
			
			### mean abundance
			if ('mean' %in% metrics) {
			
				metric <- mean(x2, na.rm=TRUE)
				thisOut <- cbind(
					thisOut,
					data.frame(
						mean = metric
					),
					row.names=NULL
				)
				
			}

			### abundance quantiles
			if ('quants' %in% metrics) {
			
				metric <- quantile(x2, quants, na.rm=TRUE)
				
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
			
				metric <- sum(x2 > 0, na.rm=TRUE) / sum(!is.na(x2))
				
				thisOut <- cbind(
					thisOut,
					data.frame(
						prevalence = metric
					),
					row.names=NULL
				)
				
			}
			
			### re-assign to obviate re-calculation
			x1 <- x2
			
			# weighted longitude/latitude... used for centroid calculations for velocities
			if (any(c('centroid', 'nsCentroid', 'ewCentroid', 'nCentroid', 'sCentriod', 'eCentroid', 'wCentroid') %in% metrics)) {

				x1weightedLongs <- x2weightedLongs
				x1weightedLats <- x2weightedLats
				
				x1centroidLong <- x2centroidLong
				x1centroidLat <- x2centroidLat

			}

			### remember
			out <- rbind(out, thisOut)
				
		} # next time period

	out
	
}
