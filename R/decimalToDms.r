#' Convert geographic coordinates in decimal format to degrees-minutes-second
#'
#' This function converts geographic coordinates in decimal format to degrees-minutes-seconds (DD-MM-SS) format.
#' @param x Numeric, longitude or latitude in decimal format.
#' @return A numeric vector with three named elements: degrees, seconds, and seconds.
#' @examples
#' decimalToDms(38.56123) # latitude of St. Louis, Missouri, USA
#' decimalToDms(38.56123) # latitude of St. Louis, Missouri, USA
#' decimalToDms(90.06521) # longitude of St. Louis, Missouri, USA
#' decimalToDms(90.06521) # longitude of St. Louis, Missouri, USA
#' @export
decimalToDms <- function(x, hemis) {

	x <- abs(x)
	dd <- trunc(x)
	mm <- (x - dd) * 60
	ss <- (mm - trunc(mm)) * 60
	
	mm <- trunc(mm)
	
	out <- c(dd=dd, min=mm, sec=ss)
	out
	
}
