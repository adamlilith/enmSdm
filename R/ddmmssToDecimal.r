#' Convert geographic coordinates in degrees-minutes-second to decimal format
#'
#' This function converts geographic coordinates in degrees-minutes-seconds (DD-MM-SS) format to decimal format.
#' @param dd Numeric. Degrees longitude or latitude. Can be a decimal value.
#' @param mm Numeric. Minutes longitude or latitude. Can be a decimal value.
#' @param ss Numeric. Second longitude or latitude. Can be a decimal value.
#' @param eastNorth Logical. Set to TRUE if the location is east of the Prime Meridian (longitude) or north of the equator (latitude).
#' @return Numeric.
#' @examples
#' ddmmssToDecimal(38, 37, 38, TRUE) # latitude of St. Louis, Missouri, USA
#' ddmmssToDecimal(90, 11, 52.1, FALSE) # longitude of St. Louis, Missouri, USA
#' @export

ddmmssToDecimal <- function(dd, mm, ss, eastNorth) {

	if (eastNorth) {
		dd + mm / 60 + ss / 3600
	} else {
		-dd - mm / 60 - ss / 3600
	}
	
}
