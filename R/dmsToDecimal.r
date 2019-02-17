#' Convert geographic coordinates in degrees-minutes-second to decimal format
#'
#' This function converts geographic coordinates in degrees-minutes-seconds (DD-MM-SS) format to decimal format.
#' @param dd Numeric. Degrees longitude or latitude. Can be a decimal value.
#' @param mm Numeric. Minutes longitude or latitude. Can be a decimal value.
#' @param ss Numeric. Second longitude or latitude. Can be a decimal value.
#' @param hemis Logical or character. If logical, set equal to \code{TRUE} if the location is east of the Prime Meridian (longitude) or north of the equator (latitude). If a character, this can be "N" (north), "S" (south), "E" (east), or "W" (west).
#' @return Numeric.
#' @examples
#' dmsToDecimal(38, 37, 38, TRUE) # latitude of St. Louis, Missouri, USA
#' dmsToDecimal(38, 37, 38, 'N') # latitude of St. Louis, Missouri, USA
#' dmsToDecimal(90, 11, 52.1, FALSE) # longitude of St. Louis, Missouri, USA
#' dmsToDecimal(90, 11, 52.1, 'W') # longitude of St. Louis, Missouri, USA
#' @export
dmsToDecimal <- function(dd, mm, ss, hemis) {

	sgn <- if ((class(hemis) == 'logical' && !hemis) | (hemis == 'S' | hemis == 'W')) {
		-1
	} else {
		1
	}
		
	out <- sgn * (dd + mm / 60 + ss / 3600)
	out
	
}
