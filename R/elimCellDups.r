#' Thin spatial points so that no more than one point falls within each cell of a raster
#'
#' This function thins spatial points such that no more than one point falls within each cell of a reference raster. if more than one point falls in a cell the first point in the input data is retained unless the user specifies a priority for keeping points.
#' @param x Data frame, matrix, or SpatialPoints* object.
#' @param r Raster* object.
#' @param longLat Two-element character list **or** two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} **or** a two-element list of integers that correspond to longitide and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that repredsent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a SpatialPoints object then this is ignored.
#' @param priority Either \code{NULL}, in which case for every cell with more than one point the first point in \code{x} is chosen, or a numeric or character list indicating preference for some points over others when points occur in the same cell. There should be the same number of elements in \code{priority} as there are rows or points in \code{x}. Priority is assigned by the natureal sort order of \code{priority}. For example, for 3 points in a cell for which \code{priority} is \code{c(2, 1, 3)}, the script will retain the second point and discard the rest. Similarly, if \code{priority} is \code{c('x', 'y', 'z')} then the third point will be chosen. Priorities assigned to points in other cells are ignored when thinning points in a particular cell.
#' @return Object of class \code{x}.
#' @details It is assumed that the points and raster have the same coordinate reference system and projection (if projected). If the points \code{x} are a SpatialPoints* object then the function will abort with an error if the coordinate reference system (CRS) for the points and raster do not match. If the raster lacks a CRS then it will be assumed to have the same CRS as the points (and a warning will be displayed).
#' @examples
#' x <- data.frame(long=c(-90.1, -90.1, -90.2, 20), lat=c(38, 38, 38, 38), point=letters[1:4])
#' r <- raster::raster() # empty raster covering enture world with 1-degree resolution
#' elimCellDups(x, r, longLat=c(1, 2))
#' elimCellDups(x, r, longLat=c(1, 2), priority=c(3, 2, 1, 1))
#' @export

elimCellDups <- function(
	x,
	r,
	longLat = NULL,
	priority = NULL
) {

	# check CRS of points and raster
	if (class(x) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
		if (is.na(raster::projection(r))) {
			warning('Raster will be assumed to have same coordinate reference system as points.', .immediate=TRUE)
			raster::projection(r) <- raster::projection(x)
		} else if (raster::projection(r) != raster::projection(x)) {
			stop('Raster and points do not have the same coordinate reference system.')
		}
	}

	# get coordinates
	xy <- xToCoords(x, longLat, sp=FALSE)

	# get cell numbers for each point and adjoin with data frame
	cellNum <- raster::cellFromXY(r, xy)

	# remember original row names
	index <- 1:nrow(xy)

	# define priority
	if (is.null(priority)) priority <- 1:nrow(xy)

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
	x <- if (class(x) %in% c('data.frame', 'matrix', 'SpatialPointsDataFrame')) {
		x[-removeThese, ]
	} else {
		x[-removeThese]
	}

	x

}

