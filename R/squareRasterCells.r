#' Resample a raster so cells are "square"
#' 
#' This function resamples a raster so that cells are "square" (have the same linear dimension on each side). It uses either the "length" or "width" as a baseline distance, then resamples the other dimension to be the same. Some things to note:
#' \itemize{
#' 	\item If a raster is in an unprojected coordinate system (e.g., if it has WGS84 or NAD83), then the new cells (as the old cells) will have the same length/width dimensions in degrees, but cells will not actually be "square" on the ground because degrees of longitude represent smaller absolute distances at higher latitudes.
#' 	\item The new cell dimensions will probably not be exactly the same due to the difficulty of fitting square frames onto a round Earth. However, they will differ at most by a fraction of a percent.
#'	\item Since re-sizing cells typically means that at the edges only a fractional cell can be fit, a template raster is created that has the same extent as the input raster is created. This is then extended by one cell in both directions in the dimension in which cells of the input raster are being expanded/contracted. The new raster is created by resampling to the extent and dimensions of this template raster. As a result, the new extent will be slightly larger than the old extent.
#'  \item Owing to the preceding reason, rasters that have edges near a pole and/or the international date line may not yield workable results if the new cells extend the raster "over" the pole or date line.
#' }
#' @param x Raster.
#' @param keepWidth Logical, if \code{TRUE}, then use the width (east-west direction) as the baseline distance and resample so that height (north-south distance) of a cell is the same. If \code{FALSE}, use height as the baseline.
#' @param ... Arguments to send to \code{\link[raster]{resample}}.
#' @return Raster object.
#' @examples
#' \donttest{
#' # get WORLDCLIM elevation data
#' elev <- raster::getData('worldclim', var='alt', res=10)
#' 
#' # crop to Madagascar
#' data(mad0)
#' elev <- crop(elev, mad0)
#'
#' # project to equal-area:
#' elev <- raster::projectRaster(elev, crs=getCRS('madAlbers'))
#'
#' # square the cells (in degrees... will not be equal-area!)
#' byWidth <- squareRasterCells(elev)
#' byHeight <- squareRasterCells(elev, keepWidth=FALSE)
#' print('Original cell dimensions:')
#' print(res(elev))
#' print('Cell dimensions, width fixed:')
#' print(res(byWidth))
#' print('Cell dimensions, height fixed:')
#' print(res(byHeight))
#' print('Relative difference in width:')
#' print(res(byWidth)[2] / res(byWidth)[1])
#' print('Relative difference in height:')
#' print(res(byHeight)[1] / res(byHeight)[2])
#' }
#' @export
squareRasterCells <- function(x, keepWidth = TRUE, ...) {

	crs <- raster::projection(x)
	
	width <- res(x)[1]
	height <- res(x)[2]
	
	ext <- raster::extent(x)
	ext <- c(ext@xmin, ext@xmax, ext@ymin, ext@ymax)
	
	northSouth <- abs(ext[4] - ext[3])
	eastWest <- abs(ext[2] - ext[1])
		
	if (keepWidth) {
		
		ncols <- ncol(x)
		nrows <- northSouth / width
		extra <- width * (nrows - trunc(nrows)) / 2
		ext <- c(
			ext[1],
			ext[2],
			ext[3] - extra,
			ext[4] + extra
		)
	
	} else if (!keepWidth) {
	
		nrows <- nrow(x)
		ncols <- eastWest / height
		extra <- height * (ncols - trunc(ncols)) / 2
		ext <- c(
			ext[1] - extra,
			ext[2] + extra,
			ext[3],
			ext[4]
		)
	
	}
		
	ext <- raster::extent(ext)
		
	template <- raster::raster(ext, nrows=nrows, ncols=ncols, crs=crs)
	rast <- raster::resample(x, template, ...)
	rast
	
}
