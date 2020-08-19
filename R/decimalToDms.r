#' Convert geographic coordinates in decimal format to degrees-minutes-second
#'
#' This function converts geographic coordinates in decimal format to degrees-minutes-seconds (DD-MM-SS) format.
#' @param x Numeric or vector of numeric values, longitude or latitude in decimal format.
#' @return A numeric matrix with three columns: degrees, seconds, and seconds. Note that the hemisphere (i.e., indicated by the sign of x) is not returned since it could be either north/south or east/west.
#' @examples
#' decimalToDms(38.56123) # latitude of St. Louis, Missouri, USA
#' decimalToDms(90.06521) # longitude of St. Louis, Missouri, USA
#' @export
decimalToDms <- function(x) {

	n <- length(x)
	
	out <- matrix(NA, ncol=3, nrow=n)
	colnames(out) <- c('dd', 'mm', 'ss')
	
	for (i in 1:n) {
		
		xx <- abs(x[i])
		dd <- trunc(xx)
		mm <- (xx - dd) * 60
		ss <- (mm - trunc(mm)) * 60
		
		mm <- trunc(mm)
		
		out[i, ] <- c(dd, mm, ss)
		
	}
		
	out
	
}
