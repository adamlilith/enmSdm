#' Calculate the precision of a geographic coordinate
#'
#' This function calculates the maximum precision of geographic coordinates in decimal-degrees format. Coordinate precision is the maximum distance between any point defined by increasing/decreasing the value of each longitude/latitude pair by one more digit past the decimal place then increasing/decreasing this value by 5. For example, if the longitude is 82.3 and latitude 45.75 then the precision is the distance between the point pairs (82.25, 45.745) and (82.35, 45.755).
#' @param x Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distVincentyEllipsoid}} is used to calculate distances. Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param ... Arguments to pass to \code{distVincentyEllipsoid}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code} or \code{crs = <object of class CRS (see sp package)>}. Otherwise WGS84 is assumed.
#' @return Numeric values (by default in units of meters).
#' @seealso \code{\link[geosphere]{distGeo}}
#' @examples
#' # coarse-precision cases
#' long <-	c(45, 1, 45.1, 5.1)
#' lat <-  c(45, 1, 45.1, 5.10)
#' ll <- cbind(long, lat)
#' precision_m <- coordPrecision(ll)
#' cbind(ll, precision_m)
#'
#' # fine-precision cases
#' long <-	rep(45, 8)
#' lat <-  c(45, 45.1, 45.11, 45.111, 45.1111, 45.11111, 45.111111, 45.1111111)
#' ll <- cbind(long, lat)
#' precision_m <- coordPrecision(ll)
#' cbind(ll, precision_m)
#'
#' # precision varies with latitude
#' long <- rep(45, 181)
#' lat <- seq(-90, 90)
#' ll <- cbind(long, lat)
#' precision_m <- coordPrecision(ll)
#' cbind(ll, precision_m)
#' plot(lat, precision_m / 1000, xlab='Latitude', ylab='Precision (km)')
#' 
#' # dateline/polar cases
#' long <-	c(0, 180, 45, 45)
#' lat <-  c(45, 45, 90, -90)
#' ll <- cbind(long, lat)
#' precision_m <- coordPrecision(ll)
#' cbind(ll, precision_m)
#'
#' \donttest{
#' # known error when longitude is negative and latitude is -90
#' long <- -45
#' lat <- -90
#' ll <- cbind(long, lat)
#' coordPrecision(ll)
#' }
#' @export
coordPrecision <- function(
	x,
	distFunct=NULL,
	...
) {

	if (is.null(distFunct)) distFunct <- geosphere::distVincentyEllipsoid

	ellipses <- list(...)
	if ('crs' %in% omnibus::ellipseNames(list)) crs <- ellipses$crs

	# convert SpatialPointsDataFrame to SpatialPoints
	if (class(x) == 'SpatialPointsDataFrame') x <- sp::SpatialPoints(coordinates(x), getCRS(raster::projection(x)))

	# convert matrix/data frame to SpatialPoints
	x <- if (class(x) %in% c('matrix', 'data.frame')) {

		if (exists('crs', inherits=FALSE)) {
			sp::SpatialPoints(x[ , 1:2, drop=FALSE], sp::CRS(crs))
		} else {
			sp::SpatialPoints(x[ , 1:2, drop=FALSE], getCRS('wgs84', TRUE))
		}

	}

	# correct CRS
	if (raster::projection(x) != getCRS('wgs84') & raster::projection(x) != getCRS('nad83')) {

		warning('Coordinates are not in WGS84 or NAD83. Projecting them to WGS84.')
		x <- sp::spTransform(x, getCRS('wgs84', TRUE))

	}

	out <- rep(NA, length(x))

	# calculate precision
	for (i in seq_along(x)) {

		## get coordinates and add deltas
		xx <- sp::coordinates(x[i])
		digits <- omnibus::countDecDigits(c(xx))

		deltaLong <- 5 * 10^-(1 + digits[1])
		deltaLat <- 5 * 10^-(1 + digits[2])

		plusPlus <- xx + cbind(deltaLong, deltaLat)
		minusMinus <- xx - cbind(deltaLong, deltaLat)

		### correct datelines/polar cases

		mult <- 1

		# crosses dateline but not at pole
		if (plusPlus[1] > 180 & (minusMinus[2] > -90 | plusPlus[2] < 90)) {

			minusMinus[1] <- minusMinus[1] - 10
			plusPlus[1] <- plusPlus[1] - 10

		# crosses dateline but not at pole
		} else if (minusMinus[1] < 0 & (minusMinus[2] > -90 | plusPlus[2] < 90)) {

			minusMinus[1] <- minusMinus[1] + 10
			plusPlus[1] <- plusPlus[1] + 10

		# crosses dataline and at a pole
		} else if (plusPlus[1] > 180 | minusMinus[1] < 0) {

			mult <- 2

			# north pole
			if (plusPlus[2] > 90) {

				plusPlus <- c(2, 90)
				minusMinus <- c(1, 89.5)

			# south pole
			} else if (minusMinus[2] < -90) {

				plusPlus <- c(2, -90)
				minusMinus <- c(1, -89.5)

			}

		# at north pole but does not cross dateline
		} else if (plusPlus[2] > 90) {

			mult <- 2

			plusPlus <- c(plusPlus[1], 90)
			minusMinus <- c(minusMinus[1], 89.5)

		# at south pole but does not cross dateline
		} else if (minusMinus[2] < -90) {

			mult <- 2

			plusPlus <- c(plusPlus[1], -90)
			minusMinus <- c(minusMinus[1], -89.5)

		}

		# calculate precision
		out[i] <- mult * distFunct(plusPlus, minusMinus, ...)

	}

	out

}
