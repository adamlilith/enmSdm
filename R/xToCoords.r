#' Extract geographic coordinates from data frame, matrix, or \code{SpatialPoints}* object
#'
#' This function extracts geographic coordinates from a data frame, matrix, or \code{SpatialPoints}* object
#' @param x A data frame, matrix, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object. If a data frame or matrix then the coordinate reference system is assumed to be unprojected (WGS84).
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this is ignored.
#' @param sp Logical. If \code{TRUE} then return object of class\code{SpatialPoints}. If \code{FALSE} then return object of class \code{matrix}.
#' @return Object of class\code{SpatialPoints} or a 2-column \code{matrix}.
#' @seealso \code{\link[sp]{SpatialPoints}}, \code{\link[sp]{SpatialPointsDataFrame}}
#' @export
xToCoords <- function(x, longLat = NULL, sp = TRUE) {

	if (class(x) %in% c('matrix', 'data.frame')) {
		if (is.null(longLat) & ncol(x) == 2) longLat <- 1:2
		x <- x[ , longLat, drop=FALSE]
		crs <- getCRS('wgs84', asCRS=TRUE)
	} else {
		crs <- sp::CRS(raster::projection(x))
		x <- sp::coordinates(x)
	}

	if (sp) x <- sp::SpatialPoints(x, proj4string=crs)
	x

}

