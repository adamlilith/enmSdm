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
		longOrLat > refLat
	} else if (direction == 's') {
		longOrLat < refLat
	} else if (direction == 'e') {
		longOrLat > refLong
	} else if (direction == 'w') {
		longOrLat < refLong
	}
	
	x1censored <- x1 * maskCells
	x1weightedLongsCensored <- x1weightedLongs * maskCells
	x1weightedLatsCensored <- x1weightedLats * maskCells

	x2censored <- x2 * maskCells
	x2weightedLongsCensored <- x2weightedLongs * maskCells
	x2weightedLatsCensored <- x2weightedLats * maskCells

	# centroid of uncensored part of distribution
	if (sum(x2censored, na.rm=TRUE) == 0) {
		distance <- 0
	} else {

		x1centroidLongCensored <- sum(x1weightedLongsCensored, na.rm=TRUE) / sum(x1censored, na.rm=TRUE)
		x1centroidLatCensored <- sum(x1weightedLatsCensored, na.rm=TRUE) / sum(x1censored, na.rm=TRUE)

		x2centroidLongCensored <- sum(x2weightedLongsCensored, na.rm=TRUE) / sum(x2censored, na.rm=TRUE)
		x2centroidLatCensored <- sum(x2weightedLatsCensored, na.rm=TRUE) / sum(x2censored, na.rm=TRUE)
		
		distance <- .euclid(x1=x2centroidLongCensored, y1=x2centroidLatCensored, x2=x1centroidLongCensored, y2=x1centroidLatCensored)
	
	}

	abundance <- sum(x2censored, na.rm=TRUE)
	list(distance=distance, abundance=abundance)
	
}
