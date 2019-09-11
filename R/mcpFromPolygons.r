#' Minimum convex polygon from a set of polygons and points
#'
#' This function returns a minimum convex polygon constructed from a set of spatial polygons (and possibly points). See \emph{Details} for more information.
#' @param polys SpatialPolygons or SpatialPolygonsDataFrame object, representing (for example) counties in which a species is known to reside. These must be in an equal-area projection!
#' @param pts Either \code{NULL} or a \code{SpatialPoints} or \code{SpatialPointsDataFrame} object in an equal-area projection. These must be in an equal-area projection! See \emph{Details}.
#' @return SpatialPolygons object representing a minimum convex polygon.
#' @details This function constructs a minimum convex polygon (MCP) from a set of spatial polygons. The MCP is constructed from the point on each polygon that lies on the border of the polygon that is closest to the centroid of the point \code{pts}, if they are provided, or if not the centroid of the set of polygons if only the polygons are provided.
#' @examples
#' # red-bellied lemur in Madagascar
#' # represented by points data and (pretend) Frarita-level occurrences
#' mad <- raster::getData(name='GADM', country='MDG', level=2)
#' madEaProj <- sp::CRS('+init=epsg:32738')
#' mad <- sp::spTransform(mad, madEaProj)
#'
#' data(lemurs)
#' redBelly <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' ll <- c('longitude', 'latitude')
#' wgs84 <- enmSdm::getCRS('wgs84', TRUE)
#' redBelly <- sp::SpatialPoints(redBelly[ , ll], proj4string=wgs84)
#' redBelly <- sp::spTransform(redBelly, madEaProj)
#'
#' faritras <- c('Vakinankaratra', 'Haute matsiatra', 'Ihorombe', 'Vatovavy Fitovinany', 'Alaotra-Mangoro', 'Analanjirofo', 'Atsinanana', 'Analamanga', 'Itasy')
#' polys <- mad[mad$NAME_2 %in% faritras, ]
#'
#' mcpPolys <- mcpFromPolygons(polys)
#' mcpPolysPoints <- mcpFromPolygons(polys, redBelly)
#'
#' # range size in km2
#' areaFromPointsOrPoly(redBelly)
#' areaFromPointsOrPoly(mcpPolys)
#' areaFromPointsOrPoly(mcpPolysPoints)
#'
#' plot(mad)
#' plot(polys, col='gray80', add=TRUE)
#' plot(mcpPolysPoints, add=TRUE, col=scales::alpha('green', 0.4))
#' plot(mcpPolys, add=TRUE, col=scales::alpha('purple', 0.4))
#' points(redBelly, pch=16)
#' legend('bottomright', legend=c('Presences', '"Occupied" Faritras', 'MCP w/ polygons', 'MCP w/ polygons & points'), fill=c(NA, 'gray', scales::alpha('purple', 0.4), scales::alpha('green', 0.4)), pch=c(16, NA, NA, NA), border=c(NA, 'black', 'black', 'black'))
#' @export
mcpFromPolygons <- function(polys, pts=NULL) {

	### useful info

		# type of polygons
		polyIsDf <- ('SpatialPolygonsDataFrame' %in% class(polys))

		# number of polygons
		numPolys <- if (polyIsDf) {
			nrow(polys)
		} else {
			length(polys)
		}

	### focal centroid

		# polygon centroids
		polyCents <- rgeos::gCentroid(polys, byid=TRUE)

		# focal centroid
		center <- if (!is.null(pts)) {
			rgeos::gCentroid(pts)
		} else {
			rgeos::gCentroid(polyCents)
		}

	### find closest points to center

		# stores coordinates of intersections (and later maybe also reference coordinates)
		coords <- matrix(nrow=0, ncol=2)
		colnames(coords) <- c('longitude', 'latitude')

		# by polygon
		for (countPoly in 1:numPolys) {

			# get this polygon
			thisPoly <- if (polyIsDf) {
				polys[countPoly, ]
			} else {
				polys[countPoly]
			}

			# find closest point on this polygon to the focal centroid

			closestPoints <- rgeos::gNearestPoints(center, thisPoly)
			polyIntersectPoint <- closestPoints[2]
			polyIntersectPointCoords <- sp::coordinates(polyIntersectPoint)
			coords <- rbind(coords, polyIntersectPointCoords)

		} # next polygon

		# add record coordinates
		if (!is.null(pts)) {
			recordCoords <- sp::coordinates(pts)
			coords <- rbind(coords, recordCoords)
		}

		rownames(coords) <- 1:nrow(coords)

		# spatialize
		proj4 <- raster::projection(polys)
		proj4 <- sp::CRS(proj4)
		coords <- sp::SpatialPoints(coords, proj4string=proj4)

	### MCP

		minConvexPoly <- adehabitatHR::mcp(coords, 100, unin='m', unout='km2')
		minConvexPoly

}
