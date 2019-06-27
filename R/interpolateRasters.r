#' Interpolate a stack of rasters
#'
#' This function returns a stack of rasters interpolated from a stack of rasters. For example, the input might represent rasters of a process measured at times t, t + 1, and t + 4. The rasters at t + 2 and t + 3 could be interpolated based on the values in the other rasters.
#' @param rasts Raster stack.
#' @param interpFrom Numeric vector, one value per raster in \code{rasts}. Values represent "distance" along the set of rasters rasters (e.g., time).
#' @param interpTo Numeric vector, values of "distances" at which to interpolate the rasters.
#' @param type Either \code{'linear'} or \code{'spline'}. The \code{\link[stats]{approx}} function is used for linear interpolation and \code{\link[stats]{spline}} for spline-based interpolation.
#' @param na.rm Logical, if \code{TRUE} (default) then ignore cases where all values in the same cells across rasters from which interpolations are made are \code{NA} (i.e., do not throw an error). If \code{FALSE}, then throw an error when this occurs.
#' @param ... Other arguments passed to \code{\link[stats]{approx}} or \code{\link[stats]{spline}}. \emph{Do not} include any of these arguments: \code{x}, \code{y}, or \code{xout}.
#' @return A raster stack, one per element in \code{interpTo}.
#' @details This function can be very memory-intensive for large rasters.  It may speed things up (and make them possible) to do interpolations piece by piece (e.g., instead of interpolating between times t0, t1, t2, t3, ..., interpolate between t0 and t1, then t1 and t2, etc. This may give results that differ from using the entire set, however. Note that using linear and splines will often yield very similar results except that in a small number of cases splines may produce very extreme interpolated values.
#' @seealso \code{\link[raster]{approxNA}}, \code{\link[stats]{approx}}, \code{\link[stats]{spline}}
#' @examples
#' interpFrom <- c(1, 3, 4, 8)
#' interpTo <- c(2, 5, 6, 7, 8, 9)
#' rx <- raster(nrows=10, ncols=10)
#' r1 <- setValues(rx, 1:100)
#' r3 <- setValues(rx, 100:1)
#' r4 <- setValues(rx, 100:1 - 30)
#' r8 <- setValues(rx, c(runif(95), rep(NA, 5)))
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
	na.rm = TRUE,
	...
) {

	### check for errors
	####################
		
		if (!any(c('linear', 'spline') %in% type)) stop('Argument "type" must be "linear" or "spline".')
		if (nlayers(rasts) < 2) stop('Argument "rasts" must have >1 raster layer.')
		if (length(interpFrom) != raster::nlayers(rasts)) stop('Argument "interpFrom" must have same length as number of rasters in argument "rasts".')
		
	### reserve blank array for output
	##################################

		rows <- ncol(rasts)
		cols <- nrow(rasts)
		numInterps <- length(interpTo)
		
		outArray <- array(NA, dim=c(rows, cols, numInterps))
		
		template <- rasts[[1]]
		rasts <- raster::as.matrix(rasts)

	### interpolate each cell
	#########################
		
		countCell <- 1
		for (countRow in 1:rows) {
		
			for (countCol in 1:cols) {
			
				y <- rasts[countCell, ]
				
				cellInterpol <- if (sum(!is.na(y)) < 2) {
					rep(NA, numInterps)
				} else if (type == 'linear') {
					stats::approx(x=interpFrom, y=y, xout=interpTo, ...)$y
					# stats::approx(x=interpFrom, y=y, xout=interpTo)$y
				} else if (type == 'spline') {
					stats::spline(x=interpFrom, y=y, xout=interpTo, ...)$y
					# stats::spline(x=interpFrom, y=y, xout=interpTo)$y
				}
				
				outArray[countRow, countCol, ] <- cellInterpol
				countCell <- countCell + 1
				
			} # next row of output
			
		} # next column of output

	### reconfigure back to raster format
	#####################################
	
		out <- raster::raster(outArray[ , , 1], template=template)
		
		for (countInterpTo in 2:length(interpTo)) {
		
			thisOut <- raster::raster(outArray[ , , countInterpTo], template=template)
			out <- raster::stack(out, thisOut)
		
		}

		interpToNames <- gsub(interpTo, pattern='-', replacement='Neg')
		interpToNames <- gsub(interpToNames, pattern='\\.', replacement='p')
		names(out) <- paste0('interpTo', interpToNames)
		out
	
}
