#' Memory-wise extraction of data from a raster or spatial object
#'
#' This function is the same as the \code{\link[terra]{extract}} function in the \pkg{terra} package, except that it is built to better handle cases when the points/lines/polygons where the extractions are are very large in memory. Essentially, the function is a wrapper for \code{extract} with a for loop that allows you to discard unwanted columns. Note that this function is not intended to make things fast... just feasible for cases where memory may be a problem. A few tips for memory management:
#' \itemize{
#'	\item	If you do not need need specific fields from an object you will be extracting from, remove them (e.g., \code{x$doNotNeed <- NULL}).
#'	\item	Likewise, removing rasters in a raster stack that you don not need may help.
#'	\item	On a Windows machine, you can give R more memory using \code{memory.limit(memory.limit() * 2^x)} where \code{n} is an integer. On my machine, which has 32 GB of RAM, I use \code{n = 29}.
#'	\item	Run R through the terminal, not through RStudio or the R GUI. This can reduce overhead because these programs take more memory.
#'	\item	Close other programs you may not be using.
#'	\item	Start in a fresh R session and load only what you need to get the extraction done. Save the output before proceeding.
#'	\item	Go get a coffee, because it may take a while.
#' }
#'
#' @param	x	\code{sf}, \code{Raster}*, \code{SpatRaster}, \code{Spatial}*, or \code{SpatVector} object from which data will be extracted from.
#' @param	y \code{sf}, \code{Spatial}* or \code{SpatVector} with points/lines/polygons denoting locations where data will be extracted from.
#' @param	keep Names or indices of the fields or values in \code{x} to keep. If \code{NULL} (default), then all columns will be kept. If you are unsure of these, try just extracting a small number items from \code{y} and see what the names in the output are.
#' @param	atATime Number of records in \code{y} to extract at one time. If the object \code{x} has just a few things you are extracting (e.g., a few layers for a raster set or a few fields from a polygon), you can probably use a larger number; maybe, say, 100000. If your raster set or polygon has a lot of fields, then a smaller number may be required.
#' @param	verbose Report progress (\code{TRUE}) or not (\code{FALSE}).
#' @param	... Arguments to pass to \code{\link[terra]{extract}}.
#'
#' @return Data frame.
#' @examples
#'
#' # This example does *not* require this function because the raster
#' # and points are so small in memory. But it illustrates the process.
#'
#' data(mad0)
#' data(lemurs)
#'
#' ll <- c('longitude', 'latitude')
#' lemurs <- lemurs[ , ll]
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
	if (inherits(x, 'sf')) x <- terra::vect(x)
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

	# sets
	n <- nrow(y)
	sets <- ceiling(n / atATime)
	
	# pre-allocate output
	thisExt <- terra::extract(x, y[1L], ...)
	if (!is.null(keep)) thisExt <- thisExt[ , keep, drop=FALSE]
	for (i in 1L:ncol(thisExt)) thisExt[ , i] <- NA
	out <- thisExt[rep(1L, n), , drop=FALSE]
	row.names(out) <- NULL

	# progress
	if (verbose) pb <- txtProgressBar(min=0, max=sets, width=30, style=3)

	# extract each set
	for (set in 1:sets) {
	
		if (verbose)
	
		start <- atATime * (set - 1) + 1
		end <- min(atATime * set, n)
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
