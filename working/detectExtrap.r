#' Identify environments requiring extrapolation
#'
#' This function re-creates the procedure described in Mesgaran et al. 2014, wherein the authors describe a method to classify environments according to whether they require no extrapolation, fall within the univariate range of each variable used in a training set but have unique combinations of variables requiring extrapolation, or fall outside the univariate range of variables. The method relies on calculating Mahalanobis distance between the centroid of a "reference" set of environments and each site. An "extrapolation" score is then assigned to each site.
#' @param x Either a RasterStack or RasterBrick, or a data frame used to define a reference set of conditions. Each column/layer is assumed to represent a variable.
#' @param y Either a RasterStack or RasterBrick, or a data frame for which to assess "distance" from the environment in \code{x}.
#' @param ... Other arguments (not used).
#' @return If \code{y} is a raster object, the output will be a list object. The first element will be a raster stack, one per variable, plus a raster with the identify of the variable most different from the reference set, plus a final raster with the combined (smallest) values of the extrapolation statistic across variables. If \code{y} is a data frame, the output is a data frame with the score for each variable, the identify of the most dissimilar variable, and the combined (smallest) values of the extrapolation statistic.
#' @seealso \code{\link[dismo]{mess}}
#' @examples
#' @export

detectExtrap <- function(
	x,
	y,
	na.rm = TRUE
) {

	if (class(x) %in% c('RasterBrick', 'RasterStack')) {
		x <- as.data.frame(x)
	}
	
	yWasSpatial <- FALSE
	if (class(y) %in% c('RasterBrick', 'RasterStack')) {
		yWasSpatial <- TRUE
		yOrig <- y
		y <- as.data.frame(y)
	}

	out <- matrix(NA, ncol=ncol(y), nrow=nrow(y)
	colnames(out) <- colnames(x)

	## cases of simple extraction
	for (thisVar in 1:ncol(x)) {
	
		xx <- x[ , thisVar, drop=TRUE]
		xMedian <- median(xx, na.rm=TRUE)
		xMax <- max(xx, na.rm=TRUE)
		xMin <- min(xx, na.rm=TRUE)
		
		dist <- 
	
	}

}
