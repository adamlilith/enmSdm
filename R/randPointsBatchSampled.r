#' Returns collated set of sampled sites from a set of "randPoints~" calls
#'
#' This function is called using a list object typically generated using the \code{\link[enmSdm]{randPontsMaster}} function. It returns a \code{SpatialPoints} object that represents all of the randomized points from all of the iterations. These can be used as the "available environment" in a niche overlap test.
#' @return A \code{\link[sp]{SpatialPoints}} object.
#' @seealso \code{\link[enmSdm]{sampleRast}}, \code{\link[enmSdm]{randPointsMaster}}, \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther1}}, \code{\link[enmSdm]{randPointsRespectingSelfOther2}}
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

randPointsBatchSampled <- function(rands) {

	if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
	
		out <- rands[[1]]
		
		for (i in 2:rands$iterations) {
			out <- rbind(out, rands[[i]])
		}
	
	} else if (attr(rands, 'randFunctName') == 'randPointsRespectingSelfOther2') {
	
		out <- rands[[1]]$x1rand
		
		if (class(out) != class(rands[[1]]$x2rand)) stop('"x1" and "x2" must have the same class (SpatialPoints, SpatialPointsDataFrame, data frame, or matrix.')

		out <- rbind(out, rands[[1]]$x2rand)
		
		for (i in 2:length(rands)) {
			out <- rbind(out, rands[[i]]$x1rand)
			out <- rbind(out, rands[[i]]$x2rand)
		}
		
	}
	
	out

}
