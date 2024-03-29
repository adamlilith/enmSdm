% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertTropicosCoords.r
\name{convertTropicosCoords}
\alias{convertTropicosCoords}
\title{Convert coordinates from TROPICOS records to decimal degrees}
\usage{
convertTropicosCoords(x = NULL, long = NULL, lat = NULL)
}
\arguments{
\item{x}{A data frame with data downloaded from TROPICOS. If \code{x} is supplied then the script will guess which fields pertain to longitude and latitude. The advantage of supplying \code{x} (versus \code{long} and \code{lat} is that the data in \code{x} is returned plus fields for decimal longitude and latitude and whether or not the record is approximate. The default value is \code{NULL}.}

\item{long}{Character list, longitudes in TROPICOS format. If this is supplied then do not supply \code{x}.}

\item{lat}{Character list, latitudes in TROPICOS format. If this is supplied then do not supply \code{x}.}
}
\value{
A data five-column frame. The first two columns are the original TROPICOS coordinates, the next two the converted coordinates, and the last a flag to indicate if the TROPICOS coordinates were approximate (\code{TRUE}) or "exact" (\code{FALSE}) or missing (\code{NA}).
}
\description{
The TROPICOS (http://www.tropicos.org/) plant specimen database of the Missouri Botanical Garden returns records in degrees-minutes-seconds format with symbols for degrees, minutes, second, and hemisphere.  This function converts this format to signed decimal degrees format.  By default, values surrounded by square brackets (\code{[]}) are returned as \code{NA} because they are geolocated to the center of a political unit (usually county or equivalent).
}
\examples{
long <- c("078°14'05\"\"W", "091°49'39\"\"W", NA, "[091°53'38\"\"W]")
lat <- c("39°40'41\"\"N", "36°34'14\"\"S", NA, "[36°46'14\"\"N]")
convertTropicosCoords(long=long, lat=lat)
}
\seealso{
\code{\link[enmSdm]{dmsToDecimal}}
}
