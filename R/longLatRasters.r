#' Generate rasters with cell values equal to cell longitude or latitude
#'
#' This function generates a raster stack with two rasters, one with cell values equal to the cell's longitude and the other with cell values equal to cthe cell's latitude.
#' @param x Raster* object. The output will have the same resolution, extent, and coordinate projection system as \code{x}.
#' @param m Logical or raster* object. If TRUE then the output rasters will have \code{NA}s placed in cells where \code{NA}s exist in \code{x} (or the first layer in \code{x} if \code{x} is a raster stack or brick. If \code{FALSE} then longitude and latitude will be calcualted for all cells. If \code{m} is a raster then values in the output will be \code{NA} wherever they are \code{NA} in \code{m}.
#' @param filePath String or \code{NULL}. If a string, then this is the path (not including file name) to which to write the raster stack with longitude/latitude rasters. If \code{NULL} then no file is written.
#' @param ... Arguments to pass to \code{writeRaster} (if \code{filePath} is not \code{NULL}).
#' @return Object of class \code{RasterStack}.
#' @examples
#' \donttest{
#' x <- raster::raster() # raster with 1 deg resolution and extent equal to entire world
#' ll <- longLatRasters(x)
#' plot(ll)
#' }
#' @export

longLatRasters <- function(
	x,
	m = NULL,
	filePath = NULL,
	...
) {

	x <- x[[1]]

	# get mask raster
	if (!is.null(m) && class(m) == 'logical') {
		if (m) {
			m <- x * 0 + 1
		} else if (!m) {
			m[] <- 1
		}
	} else {
		m <- x
		m[] <- 1
	}

	

	##########
	## MAIN ##
	##########

	# initiate lat/long rasters
	lat <- raster::raster(raster::extent(x), nrow=nrow(x), ncol=ncol(x))
	long <- raster::raster(raster::extent(x), nrow=nrow(x), ncol=ncol(x))

	if (!is.null(filePath)) {

		long <- raster::writeStart(x=long, filename=paste0(filePath, '/longitude'), ...)
		lat <- raster::writeStart(x=lat, filename=paste0(filePath, '/latitude'), ...)

	}

	# write rasters
	if (!is.null(filePath)) {

		# for each block, calculate latitude and longitude
		for (countRow in 1:nrow(x)) {

			# initiate list to store vectors of lat/long
			theseLong <- rep(NA, ncol(x))
			theseLat <- rep(NA, ncol(x))

			# assign latitudes and longitudes
			theseLong <- raster::xFromCol(object=x, col=1:ncol(x))
			theseLat <- rep(raster::yFromRow(object=x, row=countRow), ncol(x))

			# mask
			theseLong <- theseLong * m[countRow, ]
			theseLat <- theseLat * m[countRow, ]

			# remember output
			raster::writeValues(x=long, v=theseLong, start=countRow)
			raster::writeValues(x=lat, v=theseLat, start=countRow)

		} # for each block

		# stop writing
		raster::projection(lat) <- raster::projection(long) <- raster::projection(x)

		names(long) <- 'longitude'
		names(lat) <- 'latitude'

		long <- raster::writeStop(long)
		lat <- raster::writeStop(lat)

	# do in memory--no writing!
	} else {

		long[] <- rep(raster::xFromCol(object=x, col=1:ncol(x)), nrow(x))
		lat[] <- rep(raster::yFromRow(object=x, row=1:nrow(x)), each=ncol(x))

		long <- long * m
		lat <- lat * m

	}

	names(long) <- 'longitude'
	names(lat) <- 'latitude'

	raster::projection(long) <- raster::projection(x)
	raster::projection(lat) <- raster::projection(x)

	ll <- raster::stack(long, lat)
	ll

}
