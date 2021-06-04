#' Calculate distance-based weights for spatial points
#'
#' This function calculates a set of weights for point data based on inter-point distances. Weights are intended to reflect degree of independence between points as a function of their proximity. The weighting scheme is based on the following set of assumptions:
#' \itemize{
#'		\item	An observation that has no other observations within a given distance \code{dm} from it provides one "unit" of information on the process being observed.
#' 		\item	If two observations were repeated at the same exact site (same coordinates) under the same conditions, each observation would provide half a unit of information on the process because they are redundant. Similarly, if three redundant observations were made, each would provide one-third of a unit, and so on.
#'		\item   If two observations are made at a distance \code{d} apart that is less than \code{dm}, then each should receive a weight <1, depending on \code{d}. The more distant they are, the more independent the information they provide.
#' }
#' If we define the "neighborhood" of a focal point as the area circumscribed by a circle with radius \code{dm}, and the distance between the focal point and all other points in the neighborhood is stored in a vector \code{d}, the weight assigned to the focal point is: \cr
#' 
#' \deqn{latex}{\frac{1}{n + 1}}
#' 
#' @param x An object of class \code{SpatialPoints} or \code{SpatialPointsDataFrame}, or a two-column matrix or data frame representing coordinates and assumed to be unprojected (WGS84), or a symmetrical distance matrix representing distances between points.
#' @param dm Maximum distance in meters beyond which observations are considered independent.
#'
#' @return Numeric vector.
#' @examples
#'

data(lemurs)
ll <- c('longitude', 'latitude')
er <- lemurs[lemurs$species=='Eulemur rubriventer', ll]

w <- weightByDist(er, dm=10000)
hist(w)

data(mad0)
plot(mad0)
points(er, cex=10 * w) # point size ~ weight


# plot weights as a function of distance assuming
# just ONE point is in the neighborhood
dm <- 3000 # neighborhood size
focalCoords <- cbind(-99, 37)

scenarios <- 21
w <- dist_m <- rep(NA, scenarios) # stores weights and inter-point distances

for (i in 1:scenarios) {

	neighPointCoords <- focalCoords + cbind(0, (i - 1) * 0.0025)
	x <- rbind(focalCoords, neighPointCoords)
	w[i] <- weightByDist(x, dm=dm)[1]
	dist_m[i] <- geosphere::distm(x)[1, 2]

}

plot(dist_m, w, xlab='Distance (m)', ylab='Weight', ylim=c(0, 1))
abline(v=dm, lty='dotted')

# plot weights as a function of distance assuming
# TWO points in the neighborhood
dm <- 3000 # neighborhood size
focalCoords <- cbind(-99, 37)

scenarios <- 21
w <- dist_m <- rep(NA, scenarios) # stores weights and inter-point distances

for (i in 1:scenarios) {

	neighPointCoords1 <- neighPointCoords2 <-
		focalCoords + cbind(0, (i - 1) * 0.0025)
	x <- rbind(focalCoords, neighPointCoords1, neighPointCoords2)
	w[i] <- weightByDist(x, dm=dm)[1]
	dist_m[i] <- geosphere::distm(x)[1, 2]

}

plot(dist_m, w, xlab='Distance (m)', ylab='Weight', ylim=c(0, 1))
abline(v=dm, lty='dotted')




#' @export
weightByDist <- function(x, dm) {

	cl <- class(x)

	# convert x into distance matrix
	if (any(c('SpatialPoints', 'SpatialPointsDataFrame') %in% cl)) {
		x <- geosphere::distm(x)
	} else if (ncol(x) == 2) {
		x <- sp::SpatialPoints(x, getCRS('wgs84', TRUE))
		x <- geosphere::distm(x)
	}
	
	if (dim(x)[1L] > 1) diag(x) <- NA
	
	# weighting
	w <- rep(NA, nrow(x))

	for (i in 1:nrow(x)) {
	
		potentialNeighs <- x[i, ]
		d <- potentialNeighs[which(potentialNeighs <= dm)]
		n <- length(d)
		
		if (n == 0) { # no neighbors
			w[i] <- 1
		} else { # neighbors
			
			(1 / (n + 1)) * ((sum(1 - exp(-d)) / (n * (1 - exp(-dm)))) + 1)
			
			(1 / (n + 1)) * ((sum(1 - exp(-d / dm)) / (n * (1 - exp(-1)))) + 1)
		
		}

	}
	
	w

}
