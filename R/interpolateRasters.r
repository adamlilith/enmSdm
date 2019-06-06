#' Interpolate a stack of rasters
#'
#' This function returns a stack of rasters interpolated from a stack of rasters. For example, the input might represent rasters of a process measured at times t, t + 1, and t + 4. The rasters at t + 2 and t + 3 could be interpolated based on the values in the other rasters.
#' @param r Raster stack.
#' @param x Numeric vector, one value per raster in \code{r}. Values represent "distance" along the set of rasters rasters (e.g., time).
#' @param where Numeric vector, values of "distances" at which to interpolate the rasters.
#' @param method Either \code{'linear'} or \code{'spline'}. The \code{\link[stats]{approx}} function is used for linear interpolation and \code{\link[stats]{spline}} for spline-based interpolation.
#' @param ... Other arguments passed to \code{\link[stats]{approx}} or \code{\link[stats]{spline}}.
#' @return A raster stack, one per element in \code{where}.
#' @seealso \code{\link[raster]{approxNA}}, \code{\link[stats]{approx}}, \code{\link[stats]{spline}}
#' @examples
#' x <- c(1, 3, 4, 8)
#' 
#' rx <- raster(nrows=10, ncols=10)
#' r1 <- setValues(rx, 1:100)
#' r3 <- setValues(rx, 100:1)
#' r4 <- setValues(rx, 100:1 - 30)
#' r8 <- setValues(rx, runif(100))
#' r <- stack(r1, r3, r4, r8)
#' names(r) <- paste0('r', x)
#' 
#' where <- c(2, 5, 6, 7, 8, 9)
#' linear <- interpolateRasters(r, x, where)
#' splines <- interpolateRasters(r, x, where, method='spline')
#' plot(r)
#' x11()
#' plot(linear)
#' x11()
#' plot(spline)
#' @export

interpolateRasters <- function(
	r,
	x,
	where,
	method  = 'linear',
	...
) {

	if (nlayers(r) < 2) stop('Argument "r" must have >1 raster layer.')
	if (length(x) != raster::nlayers(r)) stop('Argument "x" must have same length as number of rasters in argument "r".')
	
	# blank array for output
	width <- ncol(r)
	height <- nrow(r)
	numWhere <- length(where)
	
	outArray <- array(NA, dim=c(width, height, numWhere))
	
	# by each cell
	for (countRow in 1:height) {
	
		yCoord <- raster::yFromRow(r, countRow)
	
		for (countCol in 1:width) {
		
			xCoord <- raster::xFromCol(r, countCol)
			xy <- cbind(xCoord, yCoord)
			
			y <- raster::extract(r, xy)

			cellInterpol <- if (method == 'linear') {
			
				# stats::approx(x=x, y=y, xout=where, ...)$y
				stats::approx(x=x, y=y, xout=where)$y
				
			} else if (method == 'spline') {
				
				# stats::spline(x=x, y=y, xout=where, ...)$y
				stats::spline(x=x, y=y, xout=where)$y
				
			}
			
			outArray[countRow, countCol, ] <- cellInterpol
		
		}
		
	}
	
	proj4 <- raster::projection(r)
	out <- raster::raster(outArray[ , , 1], template=r)
	
	for (countWhere in 2:length(where)) {
	
		thisOut <- raster::raster(outArray[ , , countWhere], template=r)
		out <- raster::stack(out, thisOut)
	
	}

	names(out) <- paste0('x', where)
	out
	
}
