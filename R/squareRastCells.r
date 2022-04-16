#' Resample a raster so cells are "square"
#' 
#' This function resamples a raster so that cells are "square" (have the same linear dimension on each side). It uses either the "length" or "width" as a baseline distance, then resamples the other dimension to be the same. Some things to note:
#' \itemize{
#' 	\item Cells will not really be square on the ground because, after all, the Earth is spherical and every projection distorts shape.
#' 	\item If a raster is in an unprojected coordinate system (e.g., if it has WGS84 or NAD83), then the new cells (as the old cells) will have the same length/width dimensions in degrees, but cells will not actually be "square" on the ground because degrees of longitude represent smaller absolute distances at higher latitudes.
#'	\item Since re-sizing cells typically means that ony a fractional cell can be fit at the edges of the extent, a template raster is created that has the same extent as the pkg raster is created. This is then extended by one cell in both directions in the dimension in which cells of the pkg raster are being expanded/contracted. The new raster is created by resampling to the extent and dimensions of this template raster. As a result, the new extent will be slightly larger than the old extent.
#'  \item Owing to the preceding reason, rasters that have edges near a pole and/or the international date line may not yield workable results if the new cells extend the raster "over" the pole or date line.
#' }
#' @param x \code{Raster} or \code{SpatRaster}.
#' @param keepWidth Logical, if \code{TRUE}, then use the width (east-west direction) as the baseline distance and resample so that height (north-south distance) of a cell is the same. If \code{FALSE}, use height as the baseline.
#' @param ... Arguments to send to \code{\link[raster]{resample}} (\pkg{raster} package) or \code{\link[terra]{resample}} (\pkg{terra} package).
#' @return Raster object.
#' @seealso \code{\link{rastWithSquareCells}}, \code{\link[raster]{resample}} (\pkg{raster} package), \code{\link[terra]{resample}} (\pkg{terra} package), \code{\link[raster]{aggregate}} (\pkg{raster} package), \code{\link[terra]{aggregate}} (\pkg{terra} package)
#' @examples
#' \dontrun{
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
#' byWidth <- squareRastCells(elev)
#' byHeight <- squareRastCells(elev, keepWidth=FALSE)
#' print('Original cell dimensions:')
#' print(res(elev))
#' print('Cell dimensions, width fixed:')
#' print(res(byWidth))
#' print('Cell dimensions, height fixed:')
#' print(res(byHeight))
#' }
#' @export

squareRastCells <- function(x, keepWidth = TRUE, ...) {

	if (inherits(x, 'Raster')) {

		crs <- raster::projection(x)
		ext <- raster::extent(x)
		ext <- c(ext@xmin, ext@xmax, ext@ymin, ext@ymax)
		res <- raster::res(x)
		pkg <- 'raster'
		
	} else if (inherits(x, 'SpatRaster')) {
	
		crs <- terra::crs(x)
		ext <- terra::ext(x)
		ext <- ext@ptr$vector
		res <- terra::res(x)
		pkg <- 'terra'
	
	}
	
	if (keepWidth) {
		
		# pad extent so it can accommodate all necessary rows
		width <- res[1]
		northSouth <- ext[4] - ext[3]
		ncols <- ncol(x)

		nrows <- ceiling(northSouth / width)
		crammedRows <- northSouth / width
		pad <- (width * (nrows - crammedRows)) / 2
		ext <- c(
			ext[1],
			ext[2],
			ext[3] - pad,
			ext[4] + pad
		)
	
	} else if (!keepWidth) {
	
		# pad extent so it can accommodate all necessary columns
		height <- res[2]
		eastWest <- ext[2] - ext[1]
		nrows <- nrow(x)

		ncols <- ceiling(eastWest / height)
		crammedCols <- eastWest / height
		pad <- (height * (ncols - crammedCols)) / 2
		ext <- c(
			ext[1] - pad,
			ext[2] + pad,
			ext[3],
			ext[4]
		)
	
	}
		
	ext <- raster::extent(ext)
		
	template <- raster::raster(ext, nrows=nrows, ncols=ncols, crs=crs)
	if (input == 'terra') {
		rast <- terra::rast(rast)
		rast <- terra::resample(x, template, ...)
	} else {
		rast <- raster::resample(x, template, ...)
	}
	rast
	
}
