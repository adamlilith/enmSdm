#' Euclidean distance between a pair of points
#'
#' Euclidean distance between a pair of points or two points. 
#' @param x1 Numeric
#' @param y1 Numeric
#' @param x2 Numeric
#' @param y2 Numeric
#' @details If \code{x2} and \code{y2} are \code{NULL} then the output is simply \code{x1 - y2}.
#' @keywords internal
.euclid <- compiler::cmpfun(function(x1, y1, x2=NULL, y2=NULL) {

	if (is.null(x2) | is.null(y2)) {
		x1 - y1
	} else {
		sqrt((x1 - x2)^2 + (y1 - y2)^2)
	}
	
})
