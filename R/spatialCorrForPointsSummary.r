#' Summarize spatial autocorrelation for geographic points
#'
#' This function summarizes output from the function \code{\link[enmSdm]{spatialCorrForPoints}} which estimates a null distribution of pairwise distances between points for an observed set of points.
#' @param x Matrix generated by the function \code{\link[enmSdm]{spatialCorrForPoints}}.
#' @param perc Numeric value in the range [0, 100], indicates the upper quantile of randomized distance frequencies above which observed distance frequencies are considered significant. Default is 95.
#' @param dist Logical, if \code{TRUE} (default) then the midpoint of the first distance interval at which there is an insignificant difference is returned. If \code{FALSE} then the row number of the first such interval is returned.
#' @param verbose Logical, if \code{TRUE} display results.
#' @return Numeric, midpoint of the smallest distance interval at which the observed distance distribution is insignificantly different from the null distribution.
#' @details The idea behind this measure of spatial autocorrelation is that a set of geographic points is "independent" of one another if their pairwise distances are indistinguishable from pairwise distances of the same number of points randomly located across a landscape. Typically a set of points displays clustering (non-independence) across some distances but not all distances. Thus to identify the scale of clustering pairwise distances are tabulated into bins. (I suggest using overlapping bins, e.g., from 0 to 20000 m, 10000 to 30000 m, 20000 to 40000 m, etc. See the example below for how to do this).  
#' The function \code{spatialCorrForPoints} first calculates the observed distance distribution and tabulates the frequency of distances into bins. Then, it generates a set of randomly located points equal to the same number of points as in the observed set. It then calculates the randomized distance distribution and tabulates the distances. The randomization is repeated a large number of times (the default is 100). The observed frequency of distances can be compared to the set of random distances using \code{spatialCorrForPointsSummary} and \code{spatialCorrForPointsPlot}. The default values in those functions assume that clustering occurs if the observed pairwise distance is > the 95th quantile of the null frequency distribution for that bin (i.e., a 1-tailed test), but users can specify a different percentile to demarcate significance. In practice a series of distance bins often show clustering, but the one usually of interest is the first distance bin (the one closest to 0) that has a non-significant difference between observed and expected distances. This is the characteristic diameter of a cluster of points. Points closer than this distance can be considered non-independent of one another.  
#' Alternatively, one can specify a set of points using the \code{fixed} argument. In this case, the "observed" pairwise distance distribution is tabulated from the set of pairwise distances between the points specified by argument \code{pts} and \code{fixed}. The randomized distance distribution is calculated by randomly re-locating points in \code{pts} and calculating distances to \code{fixed}.  
#' The function \code{spatialCorrForPointsWeight} calculates weights for a set of points based on the characteristic scale of spatial autocorrelation.
#' @seealso \code{\link[enmSdm]{spatialCorrForPoints}}, \code{\link[enmSdm]{spatialCorrForPointsPlot}}, \code{\link[enmSdm]{spatialCorrForPointsWeight}}
#' @examples
#' \donttest{
#' # create raster of Madagascar
#' data(mad0)
#' rast <- raster::raster(mad0, res=c(1/12, 1/12))
#' rast[] <- 1
#' rast <- raster::crop(rast, mad0)
#' mad0rast <- raster::rasterize(mad0, rast)
#' rast <- rast * mad0rast
#' 
#' # lemur point data
#' data(lemurs)
#' fulvus <- lemurs[lemurs$species == 'Eulemur fulvus', c('longitude', 'latitude')]
#' 
#' # create overlapping bins for tabulating pairwise distances
#' ext <- extent(rast)
#' southwest <- c(ext@xmin, ext@ymin)
#' northeast <- c(ext@xmax, ext@ymax)
#' maxDist <- geosphere::distGeo(southwest, northeast)
#' 
#' binLength <- 60000 # in meters
#' maxDist <- binLength * ceiling(maxDist / binLength)
#' 
#' breaks <- data.frame(
#' 	lower=seq(0, maxDist - binLength, by=binLength / 2),
#' 	upper=seq(binLength, maxDist, by=binLength / 2)
#' )
#' 
#' # compare observed pairwise distance distribution to null distribution
#' # of pairwise values from randomly located points
#' obsAndNullDistrib <- spatialCorrForPoints(
#' 	pts = fulvus,
#' 	rast = rast,
#' 	breaks = breaks,
#' 	iters = 100,
#' 	verbose = TRUE
#' )
#' 
#' # summary and plot
#' sacDist <- spatialCorrForPointsSummary(obsAndNullDistrib)
#' main <- paste('Characteristic cluster size:', sacDist, 'meters')
#' spatialCorrForPointsPlot(obsAndNullDistrib, xlab='Distance (m)', main=main)
#'
#' # calculate weights
#' weight <- 4 * spatialCorrForPointsWeight(x=obsAndNullDistrib, pts=fulvus)
#' plot(mad0, main='Point size represents weight')
#' points(fulvus, pch=1, cex=weight)
#' }
#' @export

spatialCorrForPointsSummary <- function(
	x,
	perc = 95,
	dist = TRUE,
	verbose=TRUE
) {

	obs <- x[ , 'observedProportion']
	nullCols <- which(grepl(colnames(x), pattern='randProportion'))
	nulls <- x[ , nullCols]
	nullUpper <- apply(nulls, 1, quantile, probs=perc / 100)

	firstInsigIndex <- which.max(obs <= nullUpper)

	if (length(firstInsigIndex) == 0) {
		
		if (verbose) {
			omnibus::say('The observed distance distribution is significantly clustered')
			omnibus::say('across all distance intervals when using the upper ', perc, 'th')
			omnibus::say('quantile of null distances.')
		}
		firstInsigDist <- NA
	
	} else {

		lower <- x[firstInsigIndex, 'lower']
		middle <- x[firstInsigIndex, 'middle']
		upper <- x[firstInsigIndex, 'upper']
		
		if (verbose) {
			
			omnibus::say('The distance interval at which the observed distance')
			omnibus::say('distribution is first insignificantly different from')
			omnibus::say('the null distribution using the upper ', perc, 'th quantile')
			omnibus::say('occurs in the interval defined by:')
			omnibus::say('lower: ', lower) 
			omnibus::say('middle: ', middle) 
			omnibus::say('upper: ', upper)
		}
		
	}
	
	if (dist) {
		middle
	} else {
		firstInsigIndex
	}

}
