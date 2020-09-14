#' Create a raster with square cells
#' 
#' This function creates a raster from an object with an extent (i.e., another raster or a Spatial object) with square cells. The user can specify cell resolution (linear dimension) \emph{or} the approximate number of cells desired.
#' @param x An object with a spatial extent property (e.g., a raster, raster stack, raster brick, SpatialPoints*, SpatialLines*, or SpatialPolygons*).
#' @param numCells Positive integer, approximate number of cells desired. If this is specified, then \code{res} is ignored. If this number of cells cannot be fit into the desired extent exactly, then the actual number of cells will be larger.
#' @param res Positive numeric. Size of a cell in the units of the projection of \code{x} (typically meters). Ignored if \code{numCells} is not \code{NULL}.
#' @param vals Numeric, value to assign to cells. Note that if this is shorter than the number of cells in the output, then values will be recycled. If longer, then values will be truncated. The default is to generate random values in the range [0, 1].
#' @return Raster object. The raster will have an extent of the same size or larger than the extent of \code{x}.
#' @seealso \code{\link{squareRastCells}}
#' @examples
#' # project outline of Madagascar to equal-area:
#' data(mad0)
#' mad0Ea <- sp::spTransform(mad0, getCRS('madAlbers', TRUE))
#'
#' numCells <- 101
#' cellSize_meters <- 10E4
#' byNumCells <- rastWithSquareCells(mad0Ea, numCells=numCells)
#' byCellSize <- rastWithSquareCells(mad0Ea, res=cellSize_meters)
#' 
#' par(mfrow=c(1, 2))
#' main1 <- paste0('Cells: ', numCells, ' desired, ',
#' ncell(byNumCells), ' actual')
#' plot(byNumCells, main=main1)
#' plot(mad0Ea, add=TRUE)
#' main2 <- paste0('Cells ', cellSize_meters, ' m on a side')
#' plot(byCellSize, main=main2)
#' plot(mad0Ea, add=TRUE)
#' # Note that they look the same, but the one on the left has one less row
#' # than the one on the right.
#'
#' @export

rastWithSquareCells <- function(x, numCells = NULL, res = NULL, vals = runif(100), ...) {

	if (is.null(numCells) & is.null(res)) stop('Either "numCells" or "res" must be specified.')
	if (!is.null(numCells) & !is.null(res)) warning('Both "numCells" and "res" are specified. Ignoring argument "res".')
	
	crs <- raster::projection(x)
	ext <- raster::extent(x)
	
	longDist <- ext@xmax - ext@xmin
	latDist <- ext@ymax - ext@ymin

	# create raster with approximate number of desired cells
	if (!is.null(numCells)) {
	
		res <- sqrt((longDist * latDist) / numCells)

	}

	# number of rows/columns
	numRows <- ceiling(latDist / res)
	numCols <- ceiling(longDist / res)

	# expand extent to accommodate rows/columns
	padLat <- res * (numRows - (latDist / res)) / 2
	padLong <- res * (numCols - (longDist / res)) / 2
	
	ext <- c(
		ext@xmin - padLong,
		ext@xmax + padLong,
		ext@ymin - padLat,
		ext@ymax + padLat
	)
	
	ext <- raster::extent(ext)
	out <- raster(ext, nrows=numRows, ncols=numCols, crs=raster::projection(x))
	numCells <- raster::ncell(out)
	vals <- if (length(vals) > numCells) {
		vals[1:numCells]
	} else {
		rep(vals, length.out=numCells)
	}
	raster::values(out) <- vals
	out
	
}
