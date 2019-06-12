#' Latitude of of a quantile of the geographic abundance distribution-
#'
#' This function returns the latitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species' abundance distribution. If a quantile would occur somewhere between two rows the latitude is linearly interpolated between the latitudes of the two rows bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param rowsCumSumStd Vector representing cumulative sum of rows of a matrix, standardized to a maximim value of 1.
#' @param latitude Matrix of latitudes.
#' @keywords internal
.interpolateLatFromMatrix <- compiler::cmpfun(
	function(prob, rowsCumSumStd, latitude) {

	rowIndex <- which(prob == rowsCumSumStd)

	# interpolate latitude
	if (length(rowIndex) == 0) {
	
		row1 <- which.min(abs(prob - rowsCumSumStd))
		rowIndexNA <- rowsCumSumStd
		rowIndexNA[row1] <- NA
		row2 <- which.min(abs(prob - rowIndexNA))
	
		lats <- latitude[c(row1, row2), 1]
		rowsSummedRow1 <- rowsCumSumStd[row1]
		rowsSummedRow2 <- rowsCumSumStd[row2]
		
		lats <- sort(lats)
		lat <- lats[1] + prob * (lats[2] - lats[1])
	
	} else {
	
		lat <- latitude[rowIndex, 1]
		
	}
	
	lat
	
})
