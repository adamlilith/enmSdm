% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geoThinApprox.r
\name{geoThinApprox}
\alias{geoThinApprox}
\title{Thin geographic points stochastically}
\usage{
geoThinApprox(
  x,
  minDist,
  longLat = NULL,
  distFunct = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{Data frame, matrix, or SpatialPoints* object. See \code{Details} for further information on the coordinate reference system.}

\item{minDist}{Numeric. Minimum distance (usually in m) thinned points must be from their nearest neighbor.}

\item{longLat}{Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this argument is ignored.}

\item{distFunct}{Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and references therein). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.}

\item{verbose}{Logical. If \code{TRUE} then display progress.}

\item{...}{Extra arguments to pass to \code{distFunct}.}
}
\value{
Object of class \code{x}.
}
\description{
This function thins a set of geographic points so of the remainder, none are closer than a given distance. The function is nearly the same as the function \code{thin.algorithm} in the \pkg{spThin} package, except that it accepts a data frame, matrix, SpatialPoints, or SpatialPointsDataFrame as a main argument and the user can specify the distance function to be used and that it returns an object of the same kind.
}
\details{
If \code{x} is a data frame or a matrix then it will be assumed to be unprojected (WGS84 coordinate reference system) and \code{minDist} should be in units of meters.
}
\examples{
x <- data.frame(long=c(-90.1, -90.1, -90.15, -90.17, -90.2, -89),
   lat=c(38, 38, 38, 38, 38, 38), point=letters[1:6])
set.seed(123)
geoThinApprox(x, 10000, longLat=c(1, 2)) # run #1
geoThinApprox(x, 10000, longLat=c(1, 2)) # run #2
geoThinApprox(x, 10000, longLat=c(1, 2)) # run #3
geoThinApprox(x, 10, longLat=c(1, 2))

# example using SpatialPointsDataFrame
data(lemur)
data(mad0)
par(mfrow=c(1, 3))
for (count in 1:3) {
    plot(mad0, main='Madagascar')
    points(lemur, col='red')
    thinned <- geoThinApprox(lemur, 50000)
    points(thinned, pch=16)
    legend('topright', legend=c('retained', 'discarded'),
    col=c('black', 'red'), pch=c(16, 1))
}
}
\seealso{
\code{\link[geosphere]{distGeo}}, \code{\link{geoThin}}
}
