#' Longitude of of a quantile of the geographic abundance distribution
#'
#' This function returns the longitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species' abundance distribution. If a quantile would occur somewhere between two rows the longitude is linearly interpolated between the latitudes of the two rows bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param colsCumSumStd Vector representing summed rows of a matrix.
#' @param longitude Matrix of latitudes.
#' @keywords internal
.interpolateLongFromMatrix <- compiler::cmpfun(
	function(prob, colsCumSumStd, longitude) {

	colIndex <- which(prob == colsCumSumStd)
	
	# interpolate longitude
	if (length(colIndex) == 0) {
	
		col1 <- which.min(abs(prob - colsCumSumStd))
		colIndexNA <- colsCumSumStd
		colIndexNA[col1] <- NA
		col2 <- which.min(abs(prob - colIndexNA))
	
		longs <- longitude[1, c(col1, col2)]
		colsSummedCol1 <- colsCumSumStd[col1]
		colsSummedCol2 <- colsCumSumStd[col2]

		longs <- sort(longs)
		long <- longs[1] + prob * (longs[2] - longs[1])
	
	} else {
	
		long <- longitude[1, colIndex]
		
	}
	
	long
	
})
