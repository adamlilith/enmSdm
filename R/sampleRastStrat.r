#' Sample random points from a raster in manner stratified by raster values with/without replacement.
#'
#' This function is a robust version of \code{sampleStratified()} in the \code{raster} package where sampling can be performed with replacement in each stratum. Each unique value in a raster is identified as a stratum. Each stratum is sampled a given number of times.
#' @param x Raster* object.
#' @param n number of points to sample from *each* stratum
#' @param adjArea Logical. If TRUE then adjust sample probability of cells by cell area.
#' @param replace Logical. If TRUE then sample with replacement.
#' @return 2-column matrix with longitude and latitude of random points.
#' @seealso \code{\link[raster]{sampleStratified}}, \code{\link[randomPoints]{randomPoints}}, \code{\link{sampleRast}}
#' @export
sampleRastStrat <- function(x, n, adjArea = TRUE, replace = TRUE) {

	# calculate cell area
	if (adjArea) cellArea <- raster::area(x[[1]], na.rm=TRUE)

	# convert raster to vector with no NAs
	val <- as.vector(x[[1]])

	# stratifications
	strats <- sort(stats::na.omit(unique(val)))
	rm(val)

	# sites (cell numbers)
	cellNum <- integer()

	for (i in seq_along(strats)) {

		maskStrat <- raster::calc(x, fun=function(x) ifelse(x == strats[i], 1, NA))
		if (adjArea) maskStrat <- maskStrat * cellArea

		# get just cells from this stratum
		valStrat <- as.vector(maskStrat)
		cellNumStrat <- 1:raster::ncell(maskStrat)
		cellNumStrat <- cellNumStrat[!is.na(valStrat)]
		valStrat <- valStrat[!is.na(valStrat)]

		probs <- if (adjArea) { valStrat } else { rep(1, length(valStrat)) }

		thisCellNum <- cellNumStrat[sample(seq_along(cellNumStrat), size=n, replace=replace, prob=probs)]
		cellNum <- c(cellNum, thisCellNum)

	}

	xy <- raster::xyFromCell(x, cellNum)
	xy

}
