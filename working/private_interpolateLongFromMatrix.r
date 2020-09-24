#' Longitude of of a quantile of the geographic abundance distribution-
#'
#' This function returns the longitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species abundance distribution. If a quantile would occur somewhere between two columns, the longitude is linearly interpolated between the latitudes of the two columns bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param x Matrix of abundances.
#' @param longitude Matrix of latitudes.
#' @keywords internal
.interpolateLongFromMatrix <- compiler::cmpfun(
	function(prob, x, longitude) {

	# standardized, cumulative sums of rows starting at bottom of matrix
	xColSum <- colSums(x, na.rm=TRUE)
	xColCumSum <- cumsum(xColSum)
	xColCumSumStd <- xColCumSum / max(xColCumSum)

	colIndex <- which(prob == xColCumSumStd)

	# exact longitude
	if (length(colIndex) != 0) {
	
		long <- longitude[colIndex, 1]
		
	# interpolate longitude
	} else {
	
		longVect <- longitude[1, ]
		
		# if all values are equally distant from the desired quantile
		diffs1 <- abs(prob - xColCumSumStd)
		col1 <- if (sd(diffs1) == 0) {
			round(median(seq_along(diffs1)))
		} else {
			if (prob > 0.5) {
				which.min(diffs1)
			} else {
				length(diffs1) - which.min(rev(diffs1)) + 1
			}
		}
		if (is.na(col1)) return(NA)
		long1 <- longVect[col1]
		colsCumSumStd_row1NA <- xColCumSumStd
		colsCumSumStd_row1NA[col1] <- NA
		
		# if all values are equally distant from the desired quantile
		diffs2 <- abs(prob - colsCumSumStd_row1NA)
		col2 <- if (sd(diffs2, na.rm=TRUE) == 0) {
			round(median(seq_along(diffs2)[-is.na(diffs2)]))
		} else {
			if (prob > 0.5) {
				which.min(diffs2)
			} else {
				length(diffs2) - which.min(rev(diffs2)) + 1
			}
		}
		if (is.na(col2)) return(NA)
		long2 <- longVect[col2]
	
		longs <- sort(c(long1, long2))
		# note: could calculate location weighted by abundance, but this creates problems if abundance is 0 in both latitudes or if it is negative in one cell
		long <- longs[1] + prob * (longs[2] - longs[1])
		
	}
	
	long
	
})
