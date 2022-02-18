#' Extract geographic coordinates from data frame, matrix, or \code{SpatialPoints}* object
#'
#' This function extracts geographic coordinates from a data frame, matrix, or \code{SpatialPoints}* object
#' @param x A data frame, matrix, \code{SpatialPoints*}, or \code{SpatVector} object. If a data frame or matrix then the coordinate reference system is assumed to be unprojected (WGS84).
#' @param longLat Two-element character vector \emph{or} two-element integer vector. If \code{x} is a data frame then this should be a character list specifying the names of the fields in \code{x} \emph{or} a two-element vector of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element vector indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} or \code{SpatVector} object then this is ignored.
#' @param out Type of object to return. Options are:
#' \itemize{
#'		\item \code{'default'}: Same type as input
#' 		\item \code{'sp'}: \code{SpatialPoints} object
#' 		\item \code{'vect'}: \code{SpatVaector} object
#' }
#'
#' @return Object of class\code{SpatialPoints}, class\code{SpatVector}, a 2-column \code{matrix}, or a 2-column \code{data.frame}.
#' @seealso \code{\link[sp]{SpatialPoints}}, \code{\link[terra]{SpatVector}}
#' @export
xToCoords <- function(x, longLat = NULL, out = 'default') {

	if (inherits(x, c('data.frame', 'matrix'))) {
		if (is.null(longLat) & ncol(x) == 2) longLat <- 1:2
		x <- x[ , longLat, drop=FALSE]
		crs <- getCRS('wgs84', asCRS=(out == 'sp'))
	} else if (inherits(x, 'Spatial')) {
		crs <- sp::CRS(raster::projection(x))
		if (out == 'vect') crs <- terra:crs(crs)
		x <- sp::coordinates(x)
	} else if (inherits(x, 'SpatVector')) {
		crs <- terra::crs(x)
		if (out == 'sp') crs <- as(crs, 'CRS')
		x <- terra::geom(x)
	} else {
		stop('The class of the input object is not recognized.')
	}

	if (out=='sp') {
		x <- sp::SpatialPoints(x, proj4string=crs)
	} else if (out=='vect') {
		x <- terra::vect(as.matrix(x), crs=crs)
	}
	x

}

