#' Thin spatial points so that no more than one point falls within each cell of a raster
#'
#' This function thins spatial points such that no more than one point falls within each cell of a reference raster. If more than one point falls in a cell the first point in the input data is retained unless the user specifies a priority for keeping points.
#'
#' @param x Data frame, matrix, SpatialPoints*, or SpatVector object.
#' @param rast Raster* or SpatRaster object.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a SpatialPoints object then this is ignored.
#' @param priority Either \code{NULL}, in which case for every cell with more than one point the first point in \code{x} is chosen, or a numeric or character list indicating preference for some points over others when points occur in the same cell. There should be the same number of elements in \code{priority} as there are rows or points in \code{x}. Priority is assigned by the natural sort order of \code{priority}. For example, for 3 points in a cell for which \code{priority} is \code{c(2, 1, 3)}, the script will retain the second point and discard the rest. Similarly, if \code{priority} is \code{c('z', 'y', 'x')} then the third point will be chosen. Priorities assigned to points in other cells are ignored when thinning points in a particular cell.
#' @return Object of class \code{x}.
#' @examples
#'
#' x <- data.frame(long=c(-90.1, -90.1, -90.2, 20), lat=c(38, 38, 38, 38), point=letters[1:4])
#' rast <- raster::raster() # empty raster covering entire world with 1-degree resolution
#' elimCellDups(x, rast, longLat=c(1, 2))
#' elimCellDups(x, rast, longLat=c(1, 2), priority=c(3, 2, 1, 0))
#'
#' @export

elimCellDups <- function(
	x,
	rast,
	longLat = NULL,
	priority = NULL
) {

	# reformulate
	if (inherits(x, 'SpatVector')) {
		inWasVect <- TRUE
		x <- as(x, 'Spatial')
	} else {
		inWasVect <- FALSE
	}
	if (inherits(rast, 'SpatRaster')) rast <- as(rast, 'Raster')

	# get coordinates
	xy <- xToCoords(x, longLat, out='sp')

	# get cell numbers for each point and adjoin with data frame
	cellNum <- raster::cellFromXY(rast, xy)

	# remember original row names
	index <- 1:length(xy)

	# define priority
	if (is.null(priority)) priority <- 1:length(xy)

	# index of points to remove
	removeThese <- integer()

	# remove redundant points in each cell
	uniCells <- unique(cellNum)
	for (thisCell in uniCells) {

		# if more than one point per cell
		if (sum(cellNum == thisCell) > 1) {

			thisRow <- index[thisCell == cellNum]
			thisPriority <- priority[thisCell == cellNum]

			thisRow <- thisRow[order(thisPriority)]
			removeThese <- c(removeThese, thisRow[2:length(thisRow)])

		}

	}

	# remove redundant points
	if (length(removeThese) > 0) {
		x <- if (inherits(x, c('matrix', 'data.frame', 'SpatialPointsDataFrame'))) {
			x[-removeThese, ]
		} else {
			x[-removeThese]
		}
	}

	if (inWasVect) x <- terra::vect(x)

	x

}
