#' Updates a list of random points with environmental data from rasters
#'
#' This function is called using a list object typically generated using the \code{\link[enmSdm]{randPointsMaster}} function. To each set of random points represented in that list it adds environmental data extracted from a raster stack or brick.
#' @param rands A list object typically generated using the \code{\link[enmSdm]{randPontsMaster}} function.
#' @param rast A raster, raster stack, or raster brick from which to extract data.
#' @param verbose Logical, if \code{TRUE} display progress.
#' @return A list.
#' @seealso \code{\link[enmSdm]{randPointsMaster}}, \code{\link[enmSdm]{randPointsSampled}}
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

randPointsBatchExtract <- function(rands, rast, verbose = FALSE) {

	if (verbose) {
		randsSize <- length(rands)
		prog <- utils::txtProgressBar(min=0, max=randsSize, width=32, style=3)
	}

	if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
	
		for (i in seq_along(rands)) {

			x <- rands[[i]]
			
			coords <- if (class(x) %in% c('data.frame', 'matrix')) {
				x[ , 1:2]
			} else {
				sp::coordinates(x)
			}
			
			env <- raster::extract(rast, coords)
			if (class(rands[[i]]) == 'SpatialPoints') rands[[i]] <- as(rands[[i]], 'SpatialPointsDataFrame')
			rands[[i]] <- cbind(rands[[i]], env)
		
			if (verbose) utils::setTxtProgressBar(prog, i)
		
		}
	
	} else if (attr(rands, 'randFunctName') == 'randPointsRespectingSelfOther2') {
	
		for (i in seq_along(rands)) {
		
			x1 <- rands[[i]]$x1rand
			x2 <- rands[[i]]$x2rand
			
			coords1 <- if (class(x1) %in% c('data.frame', 'matrix')) {
				x1[ , 1:2]
			} else {
				sp::coordinates(x1)
			}
			
			coords2 <- if (class(x2) %in% c('data.frame', 'matrix')) {
				x2[ , 1:2]
			} else {
				sp::coordinates(x2)
			}
			
			env1 <- raster::extract(rast, coords1)
			env2 <- raster::extract(rast, coords2)
			
			if (class(rands[[i]]$x1rand) == 'SpatialPoints') {
				rands[[i]]$x1rand <- as(rands[[i]]$x1rand, 'SpatialPointsDataFrame')
				rands[[i]]$x1rand@data <- as.data.frame(env1)
			} else {
				rands[[i]]$x1rand <- cbind(rands[[i]]$x1rand, env1)
			}
			
			if (class(rands[[i]]$x2rand) == 'SpatialPoints') {
				rands[[i]]$x2rand <- as(rands[[i]]$x2rand, 'SpatialPointsDataFrame')
				rands[[i]]$x2rand@data <- as.data.frame(env2)
			} else {
				rands[[i]]$x2rand <- cbind(rands[[i]]$x2rand, env2)
			}
			
			if (verbose) utils::setTxtProgressBar(prog, i)
			
		}
		
	}
	
	if (verbose) close(prog)
	
	rands

}
