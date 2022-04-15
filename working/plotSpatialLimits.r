#' Get x- and y-axis limits for plotting a raster or other spatial object
#'
#' This function returns the x- and y-axis limits for plotting a raster or spatial polygon. It is useful for plotting spatial objects with different spatial extents on the same plot. Using the \code{plot} interface, the first layer sets the extent. Using the \pkg{ggplot2} interface, the extent that includes all the objects is used. This can be problematic because if you want to focus the map on a specific object.

#' @param 	x Object of class \code{Raster}* (\pkg{raster} package, \code{Spatial}* (\pkg{sp} package) or \code{SpatVector} (\pkg{terra}).
#' @param	square If \code{TRUE}, then the function attempts to return x- and y-limits that would make the plot "square." In practice, this is not usualy possible because of teh curvature of the Earth, but it can help rectify very "tall" or "wide" plots.
#' @return	List
#' @examples

data(mad0)
data(lemurs)
ll <- c('longitude', 'latitude')
lemur <- lemurs[lemurs$species == 'Eulemur rufifrons', ll]

#' extent is all of Madagascar
plot(mad0)
points(lemur, col='blue')

#' extent is just that of lemurs
lims <- plotSpatialLimits(mad0)
plot(mad0, xlim=lims$xlim, ylim=lims$ylim)
points(lemur, col='blue')

#' @export

plotSpatialLimits <- function(x, square = FALSE) {

	# plot extent
	if (inherits(x, c('Raster', 'Spatial'))) {

		ext <- raster::extent(x)
		xlim <- c(ext@xmin, ext@xmax)
		ylim <- c(ext@ymin, ext@ymax)
	
	} else if (inherits(x, 'SpatVector')) {
		
		ext <- terra::ext(x)
		xlim <- ext@ptr$vector[1:2]
		ylim <- ext@ptr$vector[3:4]
		
	}
		
	# rescale to roughly square
	if (square) {

		xr <- diff(xlim)
		yr <- diff(ylim)
	
		if (yr > xr) {
		
			ratio <- yr / xr
			mid <- mean(xlim)
			xlim[1] <- xlim[1] - abs(xlim[1] - mid) * 0.5 * ratio
			xlim[2] <- xlim[2] + abs(xlim[2] - mid) * 0.5 * ratio
		
		} else {
		
			ratio <- xr / yr
			mid <- mean(ylim)
			ylim[1] <- ylim[1] - abs(ylim[1] - mid) * 0.5 * ratio
			ylim[2] <- ylim[2] + abs(ylim[2] - mid) * 0.5 * ratio
		
		}
		
	}
	
	list(xlim = xlim, ylim = ylim)
	
}
