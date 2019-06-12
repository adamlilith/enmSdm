#' Movement of occupied cells in a given direction of a fixed point
#' 
#' This function calculates the movement of a set of occupied cells that are north, south, east, or west of a given point (i.e., typically the centroid of the starting population). Values >0 confer movement to the north, south, east, or west of this portion of a species' distribution. Values equal to 0 mean no occupied cells were north, south, east, or west of the reference point.
#' @param direction Any of: \code{'n'} (north), \code{'s'} (south), \code{'e'} (east), or \code{'w'} (west).
#' @param x2 Matrix of weights (i.e., population size).
#' @param x2weightedLongs Matrix of longitudes weighted (i.e., by population size, given by \code{x2}).
#' @param x2weightedLats Matrix of latitudes weighted (i.e., by population size, given by \code{x2}).
#' @param refLong Numeric, longitude of reference point.
#' @param refLat Numeric, latitude of reference point.
#' @param longOrLat Numeric matrix, latitude or longitudes. If \code{direction} is \code{'n'} or \code{'s'} this must be latitudes. If \code{direction} is \code{'e'} or \code{'w'} this must be longitudes.
#' @return a list object with distance moved and abundance of all cells north/south/east/west of reference point.
#' @keywords internal
.cardinalVelocity <- function(
	direction,
	x2,
	x2weightedLongs,
	x2weightedLats,
	refLong,
	refLat,
	longOrLat
) {

	# mask out cells north/south/east/west of or at starting centroid
	maskCells <- if (direction == 'n') {
		longOrLat > refLat
	} else if (direction == 's') {
		longOrLat < refLat
	} else if (direction == 'e') {
		longOrLat > refLong
	} else if (direction == 'w') {
		longOrLat < refLong
	}
		
	x2censored <- x2 * maskCells
	x2weightedLongsCensored <- x2weightedLongs * maskCells
	x2weightedLatsCensored <- x2weightedLats * maskCells

	# centroid of uncensored part of distribution
	x2centroidLongCensored <- sum(x2weightedLongsCensored, na.rm=TRUE) / sum(x2censored, na.rm=TRUE)
	x2centroidLatCensored <- sum(x2weightedLatsCensored, na.rm=TRUE) / sum(x2censored, na.rm=TRUE)
	
	distance <- .euclid(x2centroidLongCensored, x2centroidLatCensored, refLong, refLat)
	abundance <- sum(x2censored, na.rm=TRUE)

	list(distance=distance, abundance=abundance)
	
}
