#' Convert coordinates from TROPICOS records to decimal degrees.
#'
#' The [TROPICOS](http://www.tropicos.org/) plant specimen database of the Missouri Botanical Garden returns records in degrees-minutes-seconds format with symbols for degrees, minutes, second, and hemisphere.  This function converts this format to signed decimal degrees format.  By default, values surrounded by square brackets (\code{[]}) are returned as \code{NA} because they are geolocated to the center of a political unit (usually county or equivalent).
#' @param long Character list, longitudes in TROPICOS format.
#' @param lat Character list, latitudes in TROPICOS format.
#' @param returnApprox Logical, if \code{TRUE} then approximate coordinates (surrounded by square brackets in TROPICOS format) are converted. If \code{FALSE} (default), then \code{NA} is returned.
#' @return A data five-column frame. The first two columns are the original TROPICOS coordinates, the next two the converted coordinates, and the last a flag to indicate if the TROPICOS coordinates were approximate (\code{TRUE}) or "exact" (\code{FALSE}).
#' @seealso \code{\link[enmSdm]{ddmmssToDecimal}}
#' @examples
#' long <- c("078°14'05\"\"W", "091°49'39\"\"W", NA, "[091°53'38\"\"W]")
#' lat <- c("39°40'41\"\"N", "36°34'14\"\"S", NA, "[36°46'14\"\"N]")
#' tropicosCoordConvert(long, lat)
#' tropicosCoordConvert(long, lat, TRUE)
#' @export
tropicosCoordConvert <- function(
	long,
	lat,
	returnApprox = FALSE
) {

	out <- data.frame(tropicosLong=long, tropicosLat=lat, longitude=NA, latitude=NA, approximate=NA)

	# convert each coordinate set
	for (i in 1:nrow(out)) {

		thisLong <- out$tropicosLong[i]
		thisLat <- out$tropicosLat[i]

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

			thisLong <- strsplit(thisLong, split='°')[[1]]
			thisLat <- strsplit(thisLat, split='°')[[1]]

			thisLong <- unlist(strsplit(thisLong, split='\''))
			thisLat <- unlist(strsplit(thisLat, split='\''))

			thisLong <- unlist(strsplit(thisLong, split='\"'))
			thisLat <- unlist(strsplit(thisLat, split='\"'))

			ew <- if (thisLong[5] == 'W') {
				-1
			} else {
				1
			}

			ns <- if (thisLat[5] == 'S') {
				-1
			} else {
				1
			}

			thisLong <- as.numeric(thisLong[1]) +
				as.numeric(thisLong[2]) / 60 +
				as.numeric(thisLong[3]) / 3600

			thisLat <- as.numeric(thisLat[1]) +
				as.numeric(thisLat[2]) / 60 +
				as.numeric(thisLat[3]) / 3600

			thisLong <- ew * thisLong
			thisLat <- ns * thisLat

			out$longitude[i] <- thisLong
			out$latitude[i] <- thisLat

		}

	}

	out

}
