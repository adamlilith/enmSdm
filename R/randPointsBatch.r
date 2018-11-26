#' Iterate "randPointsRespecting..." functions.
#'
#' This function is a wrapper function for any of the \code{randPointsRespecting~} functions. It is useful for calling the function multiple times using the same arguments. The output is a list with one element per call.
#' @param randFunct Character, any of: \code{'randPointsRespectingSelf'}, \code{'randPointsRespectingSelfOther1'}, or \code{'randPointsRespectingSelfOther2'}.
#' @param ... Arguments to pass to \code{randPointsRespectingSelf}, \code{randPointsRespectingSelfOther1}, or \code{randPointsRespectingSelfOther2}.
#' @param iterations Positive integer, number of times to call the function.
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned value will have just shuffled coordinates.
#' @param verbose Logical, if \code{TRUE} (default) show progress along iterations.
#' @param verboseEach Logical, if \code{FALSE} (default), then suppress progress indicators for each call of the \code{randPointsRespecting~} function.
#' @details Note that if you use the \code{randPointsRespectingSelfOther2} function and intend to summarize across iterations using subsequent functions (e.g., \code{\link[enmSdm]{randPointsSampled}}), it is highly advisable to ensure that the \code{x1} and \code{x2} arguments to that function have the same class (matrix, data frame, SpatialPoints, or SpatialPointsDataFrame). You should also either set the argument \code{keepData} to \code{FALSE} or ensure that the column names of \code{x1} and \code{x2} are exactly the same (which allows you to use \code{keepData = TRUE}).
#' @return A list with \code{iterations} elements.
#' @seealso \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther1}}, \code{\link[enmSdm]{randPointsRespectingSelfOther2}}, \code{\link[enmSdm]{randPointsSampled}}
#' @examples
#' # madagascar
#' library(dismo)
#' madElev <- getData('alt', country='MDG')
#' par(layout(matrix(c(1, 2), nrow=1)))
#' plot(madElev, main='Madagascar')
#' data(lemur)
#' points(lemur, pch=16)
#' rands <- randPointsRespectingSelf(lemur, mad, verbose=TRUE)
#' par(fig=1, new=FALSE)
#' points(rand, col='red')
#' @export

randPointsBatch <- function(
	randFunct,
	iterations = 100,
	...,
	rast,
	distFunct = NULL,
	keepData = FALSE,
	verbose=TRUE,
	verboseEach=FALSE
) {

	if (verboseEach) verbose <- TRUE

	if (verbose) omnibus::say('Calling ', iterations, ' of ', randFunct)

	out <- list()

	funct <- match.fun(randFunct)
	
	# for each iteration
	for (iter in 1:iterations) {
	
		if (verbose) {
			level <- if (verboseEach) { 1} else { NULL } 
			omnibus::say('Iteration ', iter, ' of ', iterations, level=level)
		}

		out[[iter]] <- funct(..., rast = rast, distFunct = distFunct, verbose = verboseEach)
		
	} # next iteration
	
	attr(out, 'randFunct') <- funct
	attr(out, 'randFunctName') <- randFunct
	attr(out, 'iterations') <- iterations
	
	out

}
