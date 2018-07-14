#' Convert coordinates from TROPICOS records to decimal degrees
#'
#' The TROPICOS (http://www.tropicos.org/) plant specimen database of the Missouri Botanical Garden returns records in degrees-minutes-seconds format with symbols for degrees, minutes, second, and hemisphere.  This function converts this format to signed decimal degrees format.  By default, values surrounded by square brackets (\code{[]}) are returned as \code{NA} because they are geolocated to the center of a political unit (usually county or equivalent).
#' @param long Character list, longitudes in TROPICOS format.
#' @param lat Character list, latitudes in TROPICOS format.
#' @param returnApprox Logical, if \code{TRUE} then approximate coordinates (surrounded by square brackets in TROPICOS format) are converted. If \code{FALSE} (default), then \code{NA} is returned.
#' @return A data five-column frame. The first two columns are the original TROPICOS coordinates, the next two the converted coordinates, and the last a flag to indicate if the TROPICOS coordinates were approximate (\code{TRUE}) or "exact" (\code{FALSE}).
#' @seealso \code{\link[enmSdm]{dmsToDecimal}}
#' @examples
#' long <- c("078°14'05\"\"W", "091°49'39\"\"W", NA, "[091°53'38\"\"W]")
#' lat <- c("39°40'41\"\"N", "36°34'14\"\"S", NA, "[36°46'14\"\"N]")
#' convertTropicosCoords(long, lat)
#' convertTropicosCoords(long, lat, TRUE)
#' @export
convertTropicosCoords <- function(
	long,
	lat,
	returnApprox = FALSE
) {

	out <- data.frame(tropicosLong=long, tropicosLat=lat, longitude=NA, latitude=NA, approximate=NA)

	# convert each coordinate set
	for (i in 1:nrow(out)) {

		thisLong <- as.character(out$tropicosLong[i])
		thisLat <- as.character(out$tropicosLat[i])

		# if NA
		if (is.na(thisLong) | is.na(thisLat)) {

			out$longitude[i] <- NA
			out$latitude[i] <- NA

		# if not converting approximate coordinates
		} else if (substr(thisLong, 1, 1) == '[' & !returnApprox) {

			out$longitude[i] <- NA
			out$latitude[i] <- NA
			out$approximate[i] <- TRUE

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

	# mask approximate coordinates
	if (!returnApprox && any(out$approximate)) {

		nas <- cbind(NA, NA)
		nas <- nas[rep(1, sum(out$approximate, na.rm=TRUE)), ]
		out[out$approximate & !is.na(out$approximate), c('longitude', 'latitude')] <- nas
	
	}
	
	out

}
