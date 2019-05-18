#' Sample random points from a raster with/out replacement
#'
#' This function returns coordinates randomly located on a raster where cells can be sampled more than once if desired (sampled with replacement) and where the probability of selection is proportinate to the cell value (plus maybe cell area--both if desired).
#' @param x Raster* object.
#' @param n Positive integer. Number of points to draw.
#' @param adjArea Logical. If TRUE then adjust probabilities so sampling accounts for cell area.
#' @param replace Logical. If TRUE then sample with replacement.
#' @param prob Logical. If TRUE then sample cells with probabilities proportional to cell values. If `adjArea` is also TRUE then probabilities are drawn proportional to the product of cell area * the value of the cell.
#' @return 2-column matrix with longitude and latitude of random points.
#' @seealso \code{\link[dismo]{randomPoints}}, \code{\link{sampleRastStrat}}
#' r <- raster::raster()
#' nc <- raster::ncell(r)
#' r[] <- 1:nc
#' rands1 <- sampleRast(r, 10000)
#' rands2 <- sampleRast(r, 10000, adjArea=FALSE)
#' rands3 <- sampleRast(r, 10000, prob=FALSE)
#' rands4 <- sampleRast(r, 10000, adjArea=FALSE, prob=FALSE)
#' par(mfrow=c(2, 2))
#' plot(r, main='adjArea = TRUE & prob = TRUE')
#' points(rands1, pch='.')
#' plot(r, main='adjArea = FALSE & prob = TRUE')
#' points(rands2, pch='.')
#' plot(r, main='adjArea = TRUE & prob = FALSE')
#' points(rands3, pch='.')
#' plot(r, main='adjArea = FALSE & prob = FALSE')
#' points(rands4, pch='.')
#' @export

## function to sample raster with/out replacement with probabilities proportional to raster values
sampleRast <- function(x, n, adjArea = TRUE, replace = TRUE, prob = TRUE) {

	val <- as.vector(x[[1]])

	# adjust probabilities for cell area and/or cell values
	if (adjArea) {

		areas <- raster::area(x[[1]], na.rm=TRUE)
		areas <- as.vector(areas)
		probs <- if (prob) {
			val * areas
		} else {
			areas
		}

	} else if (!adjArea & prob) {

		if (any(probs < 0)) warning('Some probabilities are < 0.', .immediate=TRUE)
		probs <- val

	} else if (!adjArea & !prob) {

		probs <- rep(1, length(val))

	}

	cellNum <- 1:raster::ncell(x)
	cellNum <- cellNum[!is.na(val)]

	probs <- probs[!is.na(val)]

	sites <- sample(cellNum, size=n, replace=replace, prob=probs)
	xy <- raster::xyFromCell(x, sites)

	xy

}
