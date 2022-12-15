#' Calculate number of cores to use for multi-core functions
#'
#' This function calculates the number of cores to use for multi-core functions based on user preferences, size of the raster object, number of available cores, and whether or not the user wants to use multi-core functionality even if it is possible to do the entire raster processing in one block.
#' @param rast Raster, RasterStack, or RasterBrick
#' @param cores Integer > 0. Number of cores desired to be used. Default is 2.
#' @param forceMulti Logical. If \code{TRUE} then force number of cores to be >1 even if the raster processing can be done in one chunk. Default is \code{FALSE}.
#' @return Integer.
#' @export
getCores <- function(rast, cores = 2, forceMulti = FALSE) {

	cores <- min(cores, parallel::detectCores(logical = FALSE))
	if (!forceMulti) cores <- min(cores, raster::blockSize(rast, minblocks = 1)$n)
	cores

}
