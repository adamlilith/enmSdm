#' Memory-wise extraction of data from a raster or spatial object
#'
#' This function is the same as the \code{\link[terra]{extract}} function in the \pkg{terra} package, except that it is built to better handle cases when the points/lines/polygons where the extractions are are very large in memory. Essentially, the function is a wrapper for \code{extract} with a for loop that allows you to discard unwanted columns. Note that this function is not intended to make things fast... just feasible for cases where memory may be a problem. A few tips for memory management:
#' \itemize{
#'	\item	If you don't need specific fields in an object, remove them (e.g., \code{x$doNotNeed <- NULL}).
#'	\item	Likewise, removing rasters in a raster stack that you don't need may help.
#'	\item	On a Windows machine, you can give R more memory using \code{memory.limit(memory.limit() * 2^x)} where \code{x} is an integer. On my machine, which has 32 GB of RAM, I use \code{x = 29}.
#'	\item	Run R through the terminal, not through RStudio or the R GUI. This can reduce overhead because these programs take more memory.
#' }
#'
#' @param x \code{SpatRaster} or \code{SpatVector} from which data will be out.
#' @param y \code{SpatVector} with points/lines/polygons denoting locations where data will be out.
#' @param keep Names or indices of the fields in \code{x} to keep. If \code{NULL} (default), then all columns will be kept and there will be very little advantage to using this function over just using \code{extract}.
#' @param atATime Number of records in \code{y} to extract at one time. If the object \code{x} has just a few things you are extracting (e.g., a few layers for a raster set or a few fields from a polygon), you can probably use a larger number; maybe, say, 100000. If your raster set or polygon has a lot of fields, then a smaller number may be required.
#' @param verbose Report progress (\code{TRUE}) or not (\code{FALSE}).
#' @param ... Arguments to pass to \code{\link[terra]{extract}}.
#'
#' @return Data frame.
#' @examples
#'
#' # This example does *not* require this function because the raster
#' # and points are so small in memory.
#'
#' data(mad0)
#' data(lemurs)
#'
#' ll <- c('longitude', 'latitude')
#' lemurs <- sp::SpatialPoints(lemurs[ , ll], getCRS('wgs84', TRUE))
#' 
#' # using extractLarge()
#' ex1 <- extractLarge(mad0, lemurs, atATime=10)
#'
#' # using just extract()
#' ex2 <- extract(mad0, lemurs)
#'
#' @export

extractLarge <- function(x, y, keep = NULL, atATime = 10000, verbose = TRUE, ...) {

	if (inherits(x, 'Raster')) x <- terra::rast(x)
	if (inherits(x, 'Spatial')) {
		x <- sf::st_as_sf(x)
		x <- terra::vect(x)
	}

	if (inherits(y, 'data.frame')) y <- as.matrix(y)
	if (inherits(y, 'matrix')) y <- terra::vect(y, type='points')
	if (inherits(y, 'sf')) y <- terra::vect(y)
	if (inherits(y, 'Spatial')) {
		y <- sf::st_as_sf(y)
		y <- terra::vect(y)
	}

	# pre-allocate output
	thisExt <- terra::extract(x, y[1L], ...)
	if (!is.null(keep)) thisExt <- thisExt[ , keep, drop=FALSE]
	for (i in 1L:ncol(thisExt)) thisExt[ , i] <- NA
	out <- thisExt[rep(1L, n), ]
	row.names(out) <- NULL

	# sets
	n <- nrow(y)
	nPerSet <- ceiling(n / atATime)
	sets <- ceiling(n / nPerSet)
	
	# progress
	if (verbose) pb <- txtProgressBar(min=0, max=sets, width=30, style=3)

	# extract each set
	for (set in 1:sets) {
	
		if (verbose)
	
		start <- nPerSet * (set - 1) + 1
		end <- min(nPerSet * (set), n)
		yy <- y[start:end, ]
	
		thisExtracted <- terra::extract(x, yy, ...)
		if (!is.null(keep)) thisExtracted <- thisExtracted[ , keep, drop=FALSE]
		out[start:end, ] <- thisExtracted
		gc()
		
		if (verbose) setTxtProgressBar(pb, set)
	
	}
	
	if (verbose) close(pb)
		
	out

}
