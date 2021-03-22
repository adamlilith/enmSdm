#' Return an EPSG code
#'
#' This function returns the EPSG code for a particular coordinate reference system and projection.
#' @param x Character. Name of EPSG code to return. Letter case is not important, but included below to aid memory (at least mine). Possible values include:
#'	\itemize{
#'	\item \code{albersNA} Albers equal-area conic for North America: \code{102008}
#'	\item \code{albersCONUS} Albers equal-area conic for coterminous United States: \code{5070}
#'	\item \code{chelsa} CHELSA climate dataset (actually WGS84): \code{4326}.
#'	\item \code{lccCONUS} Lambert conformal conic for conterminous United States: \code{102004}
#'	\item \code{lccEp} Lambert conformal conic for Europe: \code{3034}
#'	\item \code{lccNA} Lambert conformal conic for North America: \code{102009}
#'	\item \code{conusLccNA} Lambert conformal conic for conterminous United States: \code{102004}
#'	\item \code{lccNAsia} Lambert conformal conic for North Asia: \code{102027}
#'	\item \code{lccSA} Lambert conformal conic for South America: \code{102015}
#'	\item \code{lccSAsia} Lambert conformal conic for South Asia: \code{102030}
#'	\item \code{laeaN} Lambert azimuthal equal area for the North Pole: \code{102017}
#'	\item \code{laeaS} Lambert azimuthal equal area for the South Pole: \code{102020}
#'	\item \code{mollweide} Mollweide equal-area for the world: \code{54009}
#'	\item \code{nad27} NAD27 (unprojected): \code{4267}
#'	\item \code{nad83} NAD83 (unprojected): \code{4269}
#'	\item \code{plateCaree} Plate Caree using WGS84: \code{32662}
#'	\item \code{robinson} Robinson: \code{54030}
#'	\item \code{wgs72} CRS for WGS72 (unprojected): \code{4322}
#'	\item \code{wgs84} CRS for WGS84 (unprojected): EPSG: \code{4326}.
#'	\item \code{worldclim} WORLDCLIM climate data set (actually WGS84): \code{4326}
#' }
#' @return Integer.
#' @seealso \code{\link[enmSdm]{getCRS}}
#' @examples
#' getEPSG('wgs84')
#' getEPSG('albersNA')
#' getEPSG('mollweide')
#' @export

getEPSG <- function(x) {

	x <- tolower(x)

	out <- if (x == 'albersna') {
		102008L
	} else if (x == 'albersconus') {
		5070L # or 102003L
	} else if (x == 'chelsa') { # WGS84
		4326L
	} else if (x == 'lccnasia') {
		102027L
	} else if (x == 'lccsasia') {
		102030L
	} else if (x == 'lccsa') {
		102015L
	} else if (x == 'lccep') {
		3034L
	} else if (x == 'lccna') {
		102009L
	} else if (x == 'lccconus') {
		102004L
	} else if (x == 'laean') {
		102017L
	} else if (x == 'laeas') {
		102020L
	} else if (x == 'mollweide') {
		54009L
	} else if (x == 'nad27') {
		4267L
	} else if (x == 'nad83') {
		4269L
	} else if (x == 'platecaree') {
		32662L
	} else if (x == 'robinson') {
		54030L
	} else if (x == 'wgs72') {
		4322L
	} else if (x == 'wgs84') {
		4326L
	} else if (x == 'worldclim') {
		4326L
	} else {
		NA
	}

	out

}
