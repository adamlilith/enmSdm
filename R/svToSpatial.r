#' Convert SpatVector to Spatial*
#' 
#' This function converts a \code{SpatVector} object from the \pkg{terra} package to a \code{Spatial} object of the appropriate class (\code{SpatialPoints}, \code{SpatialPointsDataFrame}, \code{SpatialPolygons}, or \code{SpatialPolygonsDataFrame}).
#'
#' @param x		\code{SpatVector} object.
#' @return Object of class \code{Spatial}.
#' @examples
#'
#' f <- system.file('ex/lux.shp', package='terra')
#' v <- terra::vect(f)
#' spat <- svToSpatial(v)
#' spat
#'
#' @export

svToSpatial <- function(x) {

	x <- sf::st_as_sf(x)
	x <- as(x, 'Spatial')
	x

}
