#' Interpolate a stack of rasters
#'
#' This function returns a stack of rasters interpolated from a stack of rasters. For example, the input might represent rasters of a process measured at times t, t + 1, and t + 4. The rasters at t + 2 and t + 3 could be interpolated based on the values in the other rasters.
#' @param rasts Raster stack.
#' @param interpFrom Numeric vector, one value per raster in \code{rasts}. Values represent "distance" along the set of rasters rasters (e.g., time).
#' @param interpTo Numeric vector, values of "distances" at which to interpolate the rasters.
#' @param type Either \code{'linear'} or \code{'spline'}. The \code{\link[stats]{approx}} function is used for linear interpolation and \code{\link[stats]{spline}} for spline-based interpolation.
#' @param ... Other arguments passed to \code{\link[stats]{approx}} or \code{\link[stats]{spline}}. \emph{Do not} include any of these arguments: \code{x}, \code{y}, or \code{xout}.
#' @return A raster stack, one per element in \code{interpTo}.
#' @seealso \code{\link[raster]{approxNA}}, \code{\link[stats]{approx}}, \code{\link[stats]{spline}}
#' @examples
#' interpFrom <- c(1, 3, 4, 8)
#' interpTo <- c(2, 5, 6, 7, 8, 9)
#' rx <- raster(nrows=10, ncols=10)
#' r1 <- setValues(rx, 1:100)
#' r3 <- setValues(rx, 100:1)
#' r4 <- setValues(rx, 100:1 - 30)
#' r8 <- setValues(rx, runif(100))
#' rasts <- stack(r1, r3, r4, r8)
#' names(rasts) <- paste0('rasts', interpFrom)
#' 
#' linear <- interpolateRasters(rasts, interpFrom, interpTo)
#' linear <- interpolateRasters(rasts, interpFrom, interpTo, rule=2)
#' splines <- interpolateRasters(rasts, interpFrom, interpTo, type='spline')
#' plot(rasts)
#' x11()
#' plot(linear)
#' x11()
#' plot(spline)
#' @export

interpolateRasters <- function(
	rasts,
	interpFrom,
	interpTo,
	type = 'linear',
	...
) {

	if (nlayers(rasts) < 2) stop('Argument "rasts" must have >1 raster layer.')
	if (length(interpFrom) != raster::nlayers(rasts)) stop('Argument "interpFrom" must have same length as number of rasters in argument "rasts".')
	
	# blank array for output
	width <- ncol(rasts)
	height <- nrow(rasts)
	numWhere <- length(interpTo)
	
	outArray <- array(NA, dim=c(width, height, numWhere))
	
	# by each cell
	for (countRow in 1:height) {
	
		yCoord <- raster::yFromRow(rasts, countRow)
	
		for (countCol in 1:width) {
		
			xCoord <- raster::xFromCol(rasts, countCol)
			xy <- cbind(xCoord, yCoord)
			
			y <- raster::extract(rasts, xy)

			cellInterpol <- if (type == 'linear') {
				stats::approx(x=interpFrom, y=y, xout=interpTo, ...)$y
			} else if (type == 'spline') {
				stats::spline(x=interpFrom, y=y, xout=interpTo, ...)$y
			}
			
			outArray[countRow, countCol, ] <- cellInterpol
		
		}
		
	}
	
	proj4 <- raster::projection(rasts)
	out <- raster::raster(outArray[ , , 1], template=rasts)
	
	for (countWhere in 2:length(interpTo)) {
	
		thisOut <- raster::raster(outArray[ , , countWhere], template=rasts)
		out <- raster::stack(out, thisOut)
	
	}

	names(out) <- paste0('interpFrom', interpTo)
	out
	
}
