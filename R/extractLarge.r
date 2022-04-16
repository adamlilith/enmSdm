#' Memory-wise extraction of data from a raster or spatial object
#'
#' This function is the same as the \code{\link[terra]{extract}} function in the \pkg{terra} package, except that it is built to better handle cases when the points/lines/polygons where the extractions are are very large in memory. Essentially, the function is a wrapper for \code{extract} with a for loop that allows you to discard unwanted columns. Note that this function is not intended to make things fast... just feasible
#'
#' @param x \code{SpatRaster} or \code{SpatVector} from which data will be out.
#' @param y \code{SpatVector} with points/lines/polygons denoting locations where data will be out.
#' @param keep Names or indices of the fields in \code{x} to keep. If \code{NULL} (default), then all columns will be kept and there will be very little advantage to using this function over just using \code{extract}.
#' @param atATime Number of records in \code{y} to extract at one time.
#' @param verbose Report progress (\code{TRUE}) or not (\code{FALSE}).
#' @param ... Arguments to pass to \code{\link[terra]{extract}}.
#'
#' @return Data frame.
#' @examples
#'
#' @export

extractLarge <- function(x, y, keep = NULL, atATime = 10000, verbose = TRUE, ...) {

	n <- nrow(y)
	nPerSet <- ceiling(n / atATime)
	sets <- ceiling(n / nPerSet)
	
	if (verbose) pb <- txtProgressBar(min=0, max=sets, width=50)
	
	out <- data.frame()
	for (set in 1:sets) {
	
		if (verbose)
	
		start <- nPerSet * (set - 1) + 1
		end <- min(nPerSet * (set), n)
		yy <- y[start:end, ]
	
		thisExtracted <- terra::extract(x, yy, ...)
		thisExtracted <- thisExtracted[ , keep, drop=FALSE]
		out <- rbind(out, thisExtracted, make.row.names=FALSE)
		gc()
		
		if (verbose) setTxtProgressBar(pb, set)
	
	}
	
	if (verbose) close(pb)
		
	out

}
