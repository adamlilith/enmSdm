% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcpFromPolys.r
\name{mcpFromPolys}
\alias{mcpFromPolys}
\title{Minimum convex polygon from a set of spatial polygons and points}
\usage{
mcpFromPolys(polys, pts = NULL)
}
\arguments{
\item{polys}{Object representing spatial polygons of (for example) counties in which a species is known to reside. \emph{This must be in an equal-area projection!}  This object can be either a \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} (\pkg{sp} package), \code{SpatVector} (\pkg{terra} package)), or \code{sf} \code{POLYGON} or \code{MULTIPOLYGON} (\pkg{sf} package) object.}

\item{pts}{Either \code{NULL} (default) or a set of spatial points. If provided, this must be a "points" object of a class from the same package used for \code{polys}. For example, if you supply \code{polys} a \code{SpatialPolygons} object from the \pkg{sp} package, \code{pts} must be a \code{SpatialPoints} or \code{SpatialPointsDataFrame}, also from the \pkg{sp} package. \emph{These must also be in an equal-area projection!}}
}
\value{
\code{SpatialPolygons}, \code{SpatVector}, or \code{sf POLYGON} representing a minimum convex polygon.
}
\description{
This function returns a minimum convex polygon constructed from a set of spatial polygons (and possibly points).
}
\details{
This function constructs a minimum convex polygon (MCP) from a set of spatial polygons. If points are provided in \code{pts}, then the MCP is constructed from the set of \code{pts} plus the set of points that are lie on the border of each polygon where it is closest to the centroid of \code{pts}. If \code{pts} is not supplied, then the MCP is constructed from the point on each polygon that lies on the border of the polygon that is closest to the centroid of the centroids of the polygons. The centroid of the centroid is used so as to weigh each polygon equally.
}
\examples{

# This is a contrived example based on red-bellied lemurs in Madagascar
# represented by points data and (pretend) Faritras-level occurrences.

### example using Spatial* inputs (sp package)
##############################################

# Tananarive (Paris) / Laborde Grid - EPSG:29701
madProj <- '+init=epsg:29701'
madProj <- sp::CRS(madProj)

data(mad1)
mad1 <- sp::spTransform(mad1, madProj)

data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
wgs84 <- enmSdm::getCRS('wgs84', TRUE)
redBelly <- sp::SpatialPoints(redBelly[ , ll], proj4string=wgs84)
redBelly <- sp::spTransform(redBelly, madProj)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
polys <- mad1[mad1$NAME_2 \%in\% faritras, ]

mcpPolys <- mcpFromPolys(polys)
mcpPolysPoints <- mcpFromPolys(polys, redBelly)

# extent of occurrence in m2
rgeos::gArea(mcpPolys)
rgeos::gArea(mcpPolysPoints)

plot(mad1)
plot(polys, col='gray80', add=TRUE)
plot(mcpPolysPoints, col=scales::alpha('green', 0.4), add=TRUE)
plot(mcpPolys, col=scales::alpha('purple', 0.4), add=TRUE)
plot(redBelly, pch=16, add=TRUE)
legend('bottomright', 
legend=c('Presences', '"Occupied" Faritras',
'MCP w/ polygons', 'MCP w/ polygons & points'),
fill=c(NA, 'gray', scales::alpha('purple', 0.4),
scales::alpha('green', 0.4)),
pch=c(16, NA, NA, NA),
border=c(NA, 'black', 'black', 'black'))

### example using sf* inputs (sf package)
#########################################

# Tananarive (Paris) / Laborde Grid - EPSG:29701
madProj <- sf::st_crs(29701)

data(mad1)
mad1 <- sf::st_as_sf(mad1)
mad1 <- sf::st_transform(mad1, madProj)

data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- sf::st_as_sf(redBelly[ , ll], crs=4326, coords=ll)
redBelly <- sf::st_transform(redBelly, madProj)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
polys <- mad1[mad1$NAME_2 \%in\% faritras, ]

mcpPolys <- mcpFromPolys(polys)
mcpPolysPoints <- mcpFromPolys(polys, redBelly)

# extent of occurrence in m2... Areas are slightly different from using "Spatial"
# because of different projection.
sf::st_area(mcpPolys)
sf::st_area(mcpPolysPoints)

plot(sf::st_geometry(mad1))
plot(sf::st_geometry(polys), col='gray80', add=TRUE)
plot(mcpPolysPoints, col=scales::alpha('green', 0.4), add=TRUE)
plot(mcpPolys, col=scales::alpha('purple', 0.4), add=TRUE)
plot(redBelly, pch=16, add=TRUE)
legend('bottomright', 
legend=c('Presences', '"Occupied" Faritras',
'MCP w/ polygons', 'MCP w/ polygons & points'),
fill=c(NA, 'gray', scales::alpha('purple', 0.4),
scales::alpha('green', 0.4)),
pch=c(16, NA, NA, NA),
border=c(NA, 'black', 'black', 'black'))

### example using SpatVect inputs (terra package)
#################################################

# Tananarive (Paris) / Laborde Grid - EPSG:29701
wgs84 <- '+init=epsg:4326'
madProj <- '+init=epsg:29701'

data(mad1)
mad1 <- terra::vect(mad1)
mad1 <- terra::project(mad1, madProj)

data(lemurs)
redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
ll <- c('longitude', 'latitude')
redBelly <- terra::vect(redBelly[ , ll], geom=ll, crs=wgs84)
redBelly <- terra::project(redBelly, madProj)

faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe',
'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana',
'Analamanga', 'Itasy')
polys <- mad1[mad1$NAME_2 \%in\% faritras, ]

mcpPolys <- mcpFromPolys(polys)
mcpPolysPoints <- mcpFromPolys(polys, pts=redBelly)

# extent of occurrence in m2
terra::expanse(mcpPolys)
terra::expanse(mcpPolysPoints)

plot(mad1)
plot(polys, col='gray80', add=TRUE)
plot(mcpPolysPoints, col=scales::alpha('green', 0.4), add=TRUE)
plot(mcpPolys, col=scales::alpha('purple', 0.4), add=TRUE)
plot(redBelly, pch=16, add=TRUE)
legend('bottomright', 
legend=c('Presences', '"Occupied" Faritras',
'MCP w/ polygons', 'MCP w/ polygons & points'),
fill=c(NA, 'gray', scales::alpha('purple', 0.4),
scales::alpha('green', 0.4)),
pch=c(16, NA, NA, NA),
border=c(NA, 'black', 'black', 'black'))

### NOTE
# Using SpatVector input (terra package) yields EOOs that are slightly
# larger than using Spatial* (sp) or sf (sf) objects (by about 0.03-0.07\%
# in this example). The difference arises because terra::expanse yields a
# different value from rgeos::gArea and sf::st_area.

}
\references{
Smith, A.B., Murphy, S., Henderson, D., and Erickson, K.D. Including imprecisely georeferenced specimens improves accuracy of species distribution models and estimates of niche breadth. \doi{10.1101/2021.06.10.447988}
}
