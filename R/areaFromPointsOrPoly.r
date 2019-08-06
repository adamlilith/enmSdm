#' Area of a spatial polygon or set of points
#'
#' This function returns the area of a SpatialPolygon or SpatialPolygonsDataFrame object \emph{or} of the minimum convex polygon of a set of points. Input can be a set of coordinates or a polygon.
#' @param x Any of: SpatialPoints, SpatialPointsDataFrame, SpatialPolygons, or SpatialPolygonsDataFrame. Must be in an equal-area projection!
#' @return Numeric (area in km2).
#' @export
areaFromPointsOrPoly <- function(x) {

	# minimum convex polygon
	if (class(x) %in% c('SpatialPoints', 'SpatialPointsDataFrame')) {
		x <- adehabitatHR::mcp(x, percent=100, unout='km2')
	}

	# area of MCP
	area_m2 <- rgeos::gArea(x)
	area_km2 <- area_m2 / 1000^2

	area_km2
	
}
