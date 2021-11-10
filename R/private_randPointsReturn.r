#' Return output common to "randPointsRespecting~" functions
#'
#' Returns an object of same class as input unless user desired to strip non-geographic data from the output.
#' @param x Original input data (matrix, data frame, SpatialPoints, or SpatialPointsDataFrame).
#' @param randSites SpatialPoint object.
#' @param crs Character or object of class \code{CRS}.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned value will have just shuffled coordinates.
#' @param out \code{NULL} (should be \code{NULL} if \code{keepData} is \code{FALSE}), \emph{or} a copy of the original input data (the unshuffled \code{x}).
#' @keywords internal

.randPointsReturn <- function(
	x,
	randSites,
	crs,
	keepData,
	out = NULL
) {

	coords <- sp::coordinates(randSites)

	if (!keepData) {
		if (class(x) %in% c('SpatialPoints', 'SpatialPointDataFrame')) {
			out <- sp::SpatialPoints(coords, sp::CRS(crs))
		} else {
			out <- coords
			if (class(x) == 'data.frame') out <- as.data.frame(out)
		}
	} else {
		if (class(x) == 'SpatialPointsDataFrame') {
			out <- sp::SpatialPointsDataFrame(coords, data=as.data.frame(out), sp::CRS(crs))
		} else {
			out[ , 1:2] <- coords
		}
	}
	
	out
	
}
