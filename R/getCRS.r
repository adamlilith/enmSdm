#' Return a "proj4" string for a datum and projection
#'
#' This function returns the "proj4" string for a particular datum and possily projection. It is intended to be easier to use than looking up particular proj4 strings (i.e., on the web).
#' @param x Character. Name of proj4 string to return. Possible values include:
#'	\itemize{
#'	\item \code{albersNA} CRS for Equal-area Albers for North America (\code{+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs})
#'	\item \code{climateNA} CRS for ClimateNA dataset (\code{+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0})
#'	\item \code{chelsa} CRS for CHELSA climate dataset (actually WGS84, or \code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}).
#'	\item \code{dayMet} CRS for DayMet dataset (\code{+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60})
#'	\item \code{mollweide} CRS for Mollweide (\code{+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}).
#'	\code{nad27} CRS for NAD27 (\code{+proj=longlat +datum=NAD27 +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat}).
#'	\item \code{nad83} CRS for NAD83 (\code{+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0}).
#'	\item \code{nad83harn} CRS for NAD83 HARN (\code{+proj=longlat +ellps=GRS80 +no_defs}).
#'	\item \code{plateCaree} CRS Plate Caree using WGS84 (\code{C:/ecology/Drive/GIS/Projection Templates/geographic_world_plateCaree.sbn}).
#'	\item \code{prism} CRS for PRISM dataset (\code{+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs}).
#'	\item \code{wgs72} CRS for WGS72 (\code{+proj=longlat +ellps=WGS72 +no_defs}).
#'	\item \code{wgs84} CRS for WGS84 (\code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}).
#'	\item \code{worldclim} CRS for the WORLDCLIM climate data set (actually WGS84, or \code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}).
#'	\item \code{madAlbers} CRS for Albers Equal-Area Conic projection of Madagascar (\code{+proj=aea +lat_1=-14 +lat_2=-24 +lat_0=-19 +lon_0=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,}).
#' \item UTM projections require the coordinate reference system (e.g., NAD27, NAD83, WGS84, etc.), the hemisphere (north or south) and the zone. The "code" is as:
#' \code{utm<coordinate system><hemisphere><zone>}. For example, \code{utmNad83north9} returns thr CRS for UTM zone 9 in the northern hemisphere using the NAD83 coordinate reference system. Likewise, \code{utmNad27south6} returns the CRS for UTM zone 6 in the southern hemisphere using the NAD27 coordinate system. So far the following have been implemented:
#' 	\itemize{
#' 		\item Coordinate Systems: \code{nad27}, \code{nad83}, \code{wgs84}
#' 		\item Hemispsheres: \code{north}
#' 	}
#' }
#' @param asCRS Logical. If TRUE then return object is of class \code{CRS}. If FALSE then returned object is of class \code{character}.
#' @return Object of class \code{CRS} or \code{character}. Default is FALSE.
#' @seealso \code{\link[sp]{CRS}}
#' @examples
#' getCRS('wgs84')
#' getCRS('prism')
#' getCRS('mollweide', TRUE)
#' @export

getCRS <- function(
	x,
	asCRS = FALSE
) {

	x <- tolower(x)

	out <- if (x == 'albersna') {
		'+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	} else if (x == 'chelsa') {
		'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
	} else if (x == 'climatena') {
		'+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
	} else if (x == 'daymet') {
		'+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60'
	} else if (x == 'madalbers') {
		'+proj=aea +lat_1=-14 +lat_2=-24 +lat_0=-19 +lon_0=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,'
	} else if (x == 'mollweide') {
		'+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
	} else if (x == 'nad27') {
		'+proj=longlat +datum=NAD27 +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat'
	} else if (x == 'nad83') {
		'+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'
	} else if (x == 'nad83harn') {
		'+proj=longlat +ellps=GRS80 +no_defs'
	} else if (x == 'platecaree') {
		'+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
	} else if (x == 'prism') {
		'+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs'
	} else if (x == 'wgs72') {
		'+proj=longlat +ellps=WGS72 +no_defs'
	} else if (x == 'wgs84') {
		'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
	} else if (x == 'worldlclim') {
		'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
	} else if (substr(x, 1, 13) == 'utmnad27north') {
		paste0('+proj=utm +zone=', substr(x, nchar(x) - 1, nchar(x)), ' +datum=NAD27 +units=m +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat')
	} else if (substr(x, 1, 13) == 'utmnad83north') {
		paste0('+proj=utm +zone=',  substr(x, nchar(x) - 1, nchar(x)), ' +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0')
	} else if (substr(x, 1, 13) == 'utmwgs84north') {
		paste0(' +proj=utm +zone=', substr(x, nchar(x) - 1, nchar(x)), ' +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
	} else {
		NA
	}

	if (asCRS) out <- sp::CRS(out)
	out

}
