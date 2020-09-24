#' Latitude of of a quantile of the geographic abundance distribution-
#'
#' This function returns the latitude of a quantile of the geographic abundance distribution. The input is derived from a rasterized map of the species abundance distribution. If a quantile would occur somewhere between two rows, the latitude is linearly interpolated between the latitudes of the two rows bracketing its value. It will probably return \code{NA} if the quantile falls outside the range of the values.
#' @param prob Quantile value (i.e., in the range [0, 1])
#' @param x Matrix of abundances.
#' @param latitude Matrix of latitudes.
#' @keywords internal
.interpolateLatFromMatrix <- compiler::cmpfun(
	function(prob, x, latitude) {

	# standardized, cumulative sums of rows starting at bottom of matrix
	xRowSum <- rowSums(x, na.rm=TRUE)
	xRowCumSum <- cumsum(rev(xRowSum))
	xRowCumSumStd <- xRowCumSum / max(xRowCumSum)
	
	rowIndex <- which(prob == xRowCumSumStd)

	# exact latitude
	if (length(rowIndex) != 0) {
	
		lat <- rev(latitude[ , 1])[rowIndex]
		
	# interpolate latitude
	} else {
	
		latVect <- rev(latitude[ , 1])
		
		# if all values are equally distant from the desired quantile
		diffs1 <- abs(prob - xRowCumSumStd)
		row1 <- if (sd(diffs1) == 0) {
			round(median(seq_along(diffs1)))
		} else {
			if (prob > 0.5) {
				which.min(diffs1)
			} else {
				length(diffs1) - which.min(rev(diffs1)) + 1
			}
		}
		if (is.na(row1)) return(NA)
		lat1 <- latVect[row1]
		rowsCumSumStd_row1NA <- xRowCumSumStd
		rowsCumSumStd_row1NA[row1] <- NA
		
		# if all values are equally distant from the desired quantile
		diffs2 <- abs(prob - rowsCumSumStd_row1NA)
		row2 <- if (sd(diffs2, na.rm=TRUE) == 0) {
			round(median(seq_along(diffs2)[-is.na(diffs2)]))
		} else {
			if (prob > 0.5) {
				which.min(diffs2)
			} else {
				length(diffs2) - which.min(rev(diffs2)) + 1
			}
		}
		if (is.na(row2)) return(NA)
		lat2 <- latVect[row2]
	
		lats <- sort(c(lat1, lat2))
		# note: could calculate location weighted by abundance, but this creates problems if abundance is 0 in both latitudes or if it is negative in one cell
		lat <- lats[1] + prob * (lats[2] - lats[1])
		
	}
	
	lat
	
})
