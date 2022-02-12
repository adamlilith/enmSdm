#' Velocity of climate relative to set starting locations
#'
#' This function calculates the velocity of climate from a given starting point, here defined as the distance traversed by a given climate over time divided by the duration of the time period. This is commonly called "gradient" climate velocity.  There are several ways to match climate at a point in a time period to climate in a different period:
#' \itemize{
#'		\item	\code{'exact'}: Finds the location of the climate value exactly the same as the starting climate value. If no cell has an exact match, the function attemps to find an exact match using interpolation. If none can be found, the function returns \code{NA}. If there are ties, the one closest in space is chosen.
#' 		\item	\code{'nearest'}: Finds the location of the climate value that is most similar to the starting value. If there are ties, the one closest in space is chosen.
#'		\item	\code{'delta'}: Finds the location of the climate value most similar to the starting value plus/minus a small tolerance value (provided by \code{delta}). If there are ties, the one closest in space is chosen.
#' }
#'
#' @param	loc Coordinates of each starting point at time 0. This must be a 2-column matrix or data frame with longitude in the first column and latitude in the second.
#' @param	clim A stack of two climate rasters, one at the starting time period and one at the ending time period. The raster "on top" must be the for "starting" and the raster under that for "ending".
#' @param 	direction Direction in which climate velociy is to be calculated. Currently, this function only supports calculating velocity in the cardinal directions (east-west or north-south). Valid values are \code{'ew'} (east-west) and/or \code{'ns'} (north-south).
#' @param	method Any of \code{'exact'}, \code{'nearest'}, or \code{'delta'}. This specifies the method used for finding a climatic match.
#' @param	delta Either \code{NULL} (default) or a positive numeric value. This is the tolerance value for the "delta" method.
#'
#' @return One numeric value for each point. These represent distances (i.e., not "velocities" per se). To convert them to velocity, divide by the time span between the two rasters. Values are positive if velocity is in the eastward direction for east-west velocity and negative in the westeward direction.  Values are positive if velocity is in the northward direction for north-south velocity and negative in the southward direction.

#' @export

climateVelocity <- function(
	loc,
	clim,
	direction = c('ew', 'ns'),
	method = 'exact',
	delta = NULL
) {

	for (thisDir in direction) {
		
		if (thisDir == 'ew') {
		}
	
	
	
	
	
	} # next direction
	
}
