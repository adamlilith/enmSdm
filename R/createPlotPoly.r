#' Create SpatialPolygon same size as a plot
#'
#' This function creates a "rectangular" \code{SpatialPolygon} object with the same dimensions as a plot window. It is especially useful for cropping subsequent rasters or \code{SpatialPolygon} objects that would benefit from cropping. A plot must be made before calling this function.
#'
#' @param Either \code{NULL} or an object of class \code{CRS}, a coordinate reference string (proj4 string), or an object with a CRS. If any of these is provided, the \code{SpatialPolygon} object will have this CRS.
#' @return \code{SpatialPolygon}
#' @examples
#' data(mad0)
#' plot(mad0)
#' poly <- createPlotPoly(mad0)
#' plot(poly, border='blue')
#' plot(mad0, add=TRUE)
#'
#' @export

createPlotPoly <- function(x = NULL) {

	usr <- par('usr')
	ext <- raster::extent(usr)
	ext <- as(ext, 'SpatialPolygons')
	if (!is.null(x)) {
		if (inherits(x, 'CRS') | inherits(x, 'character')) {
			raster::projection(ext) <- x
		} else {
			raster::projection(ext) <- raster::projection(x)
		}
	}
	
	ext
	
}
