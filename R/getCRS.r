#' Return a "proj4" string for a datum and projection
#'
#' This function returns the "proj4" string for a particular datum and possibly projection. It is intended to be easier to use than looking up particular proj4 strings (i.e., on the web).
#' @param x Character. Name of proj4 string to return. Possible values include:
#'	\itemize{
#'	\item \code{albersNA} Albers equal-area conic for North America: \code{+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs}.
#'	\item \code{madAlbers} Albers equal-area conic for Madagascar: \code{+proj=aea +lat_1=-14 +lat_2=-24 +lat_0=-19 +lon_0=47 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,}.
#'	\item \code{climateNA} ClimateNA dataset: \code{+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0}.
#'	\item \code{chelsa} CHELSA climate dataset (actually WGS84): or \code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}.
#'	\item \code{dayMet} DayMet dataset: \code{+proj=lcc +lon_0=-100 +lat_0=42.5 +x_0=0 +y_0=0 +a=6378137 +rf=298.257223563 +lat_1=25 +lat_2=60}.
#'	\item \code{lccCONUS} Lambert conformal conic for the contiguous United States: \code{+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs}.
#'	\item \code{lccEp} Lambert conformal conic for Europe: \code{+proj=lcc +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs}.
#'	\item \code{lccNA} Lambert conformal conic for North America: \code{+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs}.
#'	\item \code{lccNAsia} Lambert conformal conic for North Asia: \code{+proj=lcc +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}.
#'	\item \code{lccSA} Lambert conformal conic for South America: \code{+proj=lcc +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs}.
#'	\item \code{lccSAsia} Lambert conformal conic for South Asia: \code{+proj=lcc +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}.
#'	\item \code{laeaN} Lambert azimuthal equal area for the North Pole: \code{+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}.
#'	\item \code{laeaS} Lambert azimuthal equal area for the South Pole: \code{+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}.
#'	\item \code{mollweide} Mollweide equal-area for the world: \code{+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs}.
#'	\code{nad27} NAD27 (unprojected): \code{+proj=longlat +datum=NAD27 +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat}.
#'	\item \code{nad83} NAD83 (unprojected): \code{+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0}.
#'	\item \code{nad83harn} NAD83 HARN (unprojected): \code{+proj=longlat +ellps=GRS80 +no_defs}.
#'	\item \code{plateCaree} Plate Caree using WGS84: \code{C:/ecology/Drive/GIS/Projection Templates/geographic_world_plateCaree.sbn}.
#'	\item \code{prism} PRISM dataset (actually WGS84): \code{+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs}.
#'	\item \code{wgs72} CRS for WGS72 (unprojected): \code{+proj=longlat +ellps=WGS72 +no_defs}.
#'	\item \code{wgs84} CRS for WGS84 (unprojected): \code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}.
#'	\item \code{worldclim} WORLDCLIM climate data set (actually WGS84): \code{+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0}.
#'  \item UTM projections require the datum (e.g., NAD27, NAD83, WGS84, etc.), the hemisphere (north or south) and the zone. The "code" is as:
#' \code{utm<coordinate system><hemisphere><zone>}. For example, \code{utmNad83north9} returns the CRS for UTM zone 9 in the northern hemisphere using the NAD83 coordinate reference system. Likewise, \code{utmNad27south6} returns the CRS for UTM zone 6 in the southern hemisphere using the NAD27 coordinate system. So far the following have been implemented:
#' 	\itemize{
#' 		\item Coordinate Systems: \code{nad27}, \code{nad83}, \code{wgs84}
#' 		\item Hemispheres: \code{north}
#' 	}
#' }
#' @param asCRS Logical. If \code{TRUE} then return object is of class \code{CRS}. If \code{FALSE} (default) then returned object is of class \code{character}.
#' @return Object of class \code{CRS} or \code{character}.
#' @seealso \code{\link[birdsEye]{makeCRS}}
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
	} else if (x == 'lccnasia') {
		'+proj=lcc +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
	} else if (x == 'lccsasia') {
		'+proj=lcc +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
	} else if (x == 'lccsa') {
		'+proj=lcc +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs'
	} else if (x == 'lccep') {
		'+proj=lcc +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs'
	} else if (x == 'lccna') {
		'+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	} else if (x == 'lccconus') {
		'+proj=lcc +lat_1=33 +lat_2=45 +lat_0=39 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'
	} else if (x == 'laean') {
		'+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
	} else if (x == 'laeas') {
		'+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
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
