#' Select random points from a portion of a raster
#'
#' This function selects random points from a randomly-selected subregion of raster. The subregion is selected by randomly locating a point on the any non-\code{NA} cell in the raster, buffering it by a given width, rasterizing the buffered area, then selecting from any non-\code{NA} cells in this region.
#' @param rast Raster, RasterStack, or RasterBrick demarcating the area in which randomized sites are to be placed. If a RasterStack or a RasterBrick is used then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param n Positive integer, number of random points to select.
#' @param width Positive numeric, size of buffer radius (in map units--usually meters).
#' @param restrict Logical, if \code{TRUE} (default), then select random points from any non-\code{NA} cell in a restricted area. If \code{FALSE}, then select from any non-\code{NA} cell on the raster. If \code{FALSE}, then is exactly the same as using \code{\link[dismo]{randomPoints}} except that it automatically returns a SpatialPoints object for the random points.
#' @param circleCent A SpatialPoints or SpatialPointsDataFrame object. The default value is \code{NULL}, in which case a centroid for the circle is randomly chosen. If supplied, then this will be the center of the circle.
#' @param ... Other arguments (not used).
#' @return A list object with four elements: 1) A SpatialPoints object of randomly selected points; 2) a SpatialPoints object representing the  center of the circle chosen to represent the subregion; 3) a numeric value representing the buffer width; and 4) a numeric value representing the number of non-\code{NA} cells in the region from which the points were selected.
#' @keywords internal

.sampleRastSubregion <- compiler::cmpfun(function(
	rast,
	n,
	width,
	restrict = TRUE,
	circleCent = NULL,
	...
) {

	crs <- raster::projection(rast)

	# initial placement of random sites in a restricted area
	if (restrict) {
	
		# random circle center
		if (is.null(circleCent)) {
			circleCent <- enmSdm::sampleRast(rast, 1, prob=FALSE)
			circleCent <- SpatialPoints(circleCent, sp::CRS(crs))
		}

		# rasterize the restricted area
		focalArea <- raster::buffer(circleCent, width=width)
		rastCrop <- raster::crop(rast, focalArea)
		rast <- raster::rasterize(focalArea, rastCrop)
		
		# for some reason rasterizing can mask NA cells so they appear non-NA, so remove them
		rast <- rast * rastCrop
	
	}

	n <- if (raster::ncell(rast) > n) { ncell(rast) } else { n }
	randPoints <- dismo::randomPoints(rast, n)
	randPoints <- sp::SpatialPoints(randPoints, sp::CRS(crs))
	
	out <- list(randPoints = randPoints, n=n, circleCent = circleCent, width = width)

})
