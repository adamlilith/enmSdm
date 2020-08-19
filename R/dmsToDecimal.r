#' Convert geographic coordinates in degrees-minutes-second to decimal format
#'
#' This function converts geographic coordinates in degrees-minutes-seconds (DD-MM-SS) format to decimal format.
#' @param dd Numeric. Degrees longitude or latitude. Can be a decimal value.
#' @param mm Numeric. Minutes longitude or latitude. Can be a decimal value.
#' @param ss Numeric. Second longitude or latitude. Can be a decimal value.
#' @param hemis Character or \code{NULL} (default). "N" (north), "S" (south), "E" (east), or "W" (west). If left as \code{NULL}, then the value returned will always be positive, even if it is in the western or southern hemisphere.
#' @return Numeric.
#' @examples
#' dmsToDecimal(38, 37, 38) # latitude of St. Louis, Missouri, USA
#' dmsToDecimal(38, 37, 38, 'N') # latitude of St. Louis, Missouri, USA
#' dmsToDecimal(90, 11, 52.1) # longitude of St. Louis, Missouri, USA
#' dmsToDecimal(90, 11, 52.1, 'W') # longitude of St. Louis, Missouri, USA
#' @export
dmsToDecimal <- function(dd, mm, ss, hemis = NULL) {

	n <- length(dd)

	if (all(is.null(hemis))) {
		sgn <- rep(1, n)
		warning('Hemisphere not specified. Assuming northern/western hemisphere.')
	} else {
		sgn <- ifelse(toupper(hemis) == 'S' | toupper(hemis) == 'W', -1, ifelse(toupper(hemis) == 'E' | toupper(hemis) == 'N', 1, 1))
		if (!any(toupper(hemis) %in% c('N', 'S', 'E', 'W'))) warning('Hemisphere invalid for at least some coordinates. Assuming northern/western hemisphere for these coordinates.')
	}
	
	out <- sgn * (dd + mm / 60 + ss / 3600)
	out
	
}
