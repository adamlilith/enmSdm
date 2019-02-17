#' Convert coordinates from TROPICOS records to decimal degrees
#'
#' The TROPICOS (http://www.tropicos.org/) plant specimen database of the Missouri Botanical Garden returns records in degrees-minutes-seconds format with symbols for degrees, minutes, second, and hemisphere.  This function converts this format to signed decimal degrees format.  By default, values surrounded by square brackets (\code{[]}) are returned as \code{NA} because they are geolocated to the center of a political unit (usually county or equivalent).
#' @param x A data frame with data downloaded from TROPICOS. If \code{x} is supplied then the script will guess which fields pertain to longitude and latitude. The advantage of supplying \code{x} (versus \code{long} and \code{lat} is that the data in \code{x} is returned plus fields for decimal longitude and latitude and whether or not the record is approximate. The default value is \code{NULL}.
#' @param long Character list, longitudes in TROPICOS format. If this is supplied then do not supply \code{x}.
#' @param lat Character list, latitudes in TROPICOS format. If this is supplied then do not supply \code{x}.
#' @return A data five-column frame. The first two columns are the original TROPICOS coordinates, the next two the converted coordinates, and the last a flag to indicate if the TROPICOS coordinates were approximate (\code{TRUE}) or "exact" (\code{FALSE}) or missing (\code{NA}).
#' @seealso \code{\link[enmSdm]{dmsToDecimal}}
#' @examples
#' long <- c("078°14'05\"\"W", "091°49'39\"\"W", NA, "[091°53'38\"\"W]")
#' lat <- c("39°40'41\"\"N", "36°34'14\"\"S", NA, "[36°46'14\"\"N]")
#' convertTropicosCoords(long=long, lat=lat)
#' @export
convertTropicosCoords <- function(
	x=NULL,
	long=NULL,
	lat=NULL
) {

	if (!is.null(x) & (!is.null(long) | !is.null(lat))) stop('You must supply either "x" OR "long" plus "lat" to function convertTropicosCoords().')

	# define output structure
	if (!is.null(x)) {
		out <- x
		out$approximate <- out$latitude <- out$longitude <- NA
		tropicosLong <- x$Longitude
		tropicosLat <- x$Latitude
	} else {
		out <- data.frame(tropicosLong=long, tropicosLat=lat, longitude=NA, latitude=NA, approximate=NA)
		tropicosLong <- long
		tropicosLat <- lat
	}

	
	# convert each coordinate set
	for (i in 1:nrow(out)) {

		thisLong <- as.character(tropicosLong[i])
		thisLat <- as.character(tropicosLat[i])

		# if NA
		if (is.na(thisLong) || is.na(thisLat) || thisLong == '' || thisLat == '') {

			out$longitude[i] <- NA
			out$latitude[i] <- NA
			
		# convert
		} else {

			approx <- (substr(thisLong, 1, 1) == '[')

			if (approx) {

				thisLong <- substr(thisLong, 2, nchar(thisLong) - 1)
				thisLat <- substr(thisLat, 2, nchar(thisLat) - 1)
				out$approximate[i] <- TRUE

			} else {

				out$approximate[i] <- FALSE

			}

			thisLongDms <- c(substr(thisLong, 1, 3), substr(thisLong, 5, 6), substr(thisLong, 8, 9))
			thisLatDms <- c(substr(thisLat, 1, 2), substr(thisLat, 4, 5), substr(thisLat, 7, 8))
			
			hemisLong <- substr(thisLong, nchar(thisLong), nchar(thisLong))
			hemisLat <- substr(thisLat, nchar(thisLat), nchar(thisLat))
			
			thisLong <- dmsToDecimal(dd=as.numeric(thisLongDms[1]), mm=as.numeric(thisLongDms[2]), ss=as.numeric(thisLongDms[3]), hemis=hemisLong)
			thisLat <- dmsToDecimal(dd=as.numeric(thisLatDms[1]), mm=as.numeric(thisLatDms[2]), ss=as.numeric(thisLatDms[3]), hemis=hemisLat)

			out$longitude[i] <- thisLong
			out$latitude[i] <- thisLat

		}

	}

	out

}
