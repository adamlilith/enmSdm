#' Spatial designation of DawrinCore data (ie, GBIF)
#'
#' This function assigns a designation about the geolocation certainty of each record in a data frame in DarwinCore format (i.e., a GBIF download using the "DarwinCore" option). The designation is based on what data is available in the records and its comparison to a set of spatial polygons representing the study region over which records were collected. The default values for arguments assume the spatial objects were obtained from GADM. Each record is assigned a category indicating the kind of spatial uncertainty to use, based on whether or not it has values for state/province, county, coordinates and their precision, and coordinate uncertainty:
#' \itemize{
#' 	\item "certain/precise": The record represents a precise point occurrence with small uncertainty.
#' 	\item "uncertain/precise": The record represents a precise point occurrence with large uncertainty.
#' 	\item "certain/imprecise": The record represents an imprecise point occurrence with small uncertainty.
#' 	\item "uncertain/imprecise": The record represents an imprecise point occurrence with large uncertainty.
#' 	\item "county/precise": The record represents a precise point with no associated uncertainty, so can only be assigned to a county.
#' 	\item "county/imprecise": The record represents an imprecise point with no associated uncertainty, so can only be assigned to a county.
#' 	\item "state/imprecise": The record represents an very imprecise point with no associated uncertainty, so can only be assigned to a state.
#'  \item "county-only" and "state-only": The record has no coordinates so can only be assigned to a county or a state.
#'  \item "unusable": The record cannot be assigned to a unique location within the study region.
#' }
#' @param darwin Data frame in DarwinCore format. The data frame must contain at least these fields, although they can contain \code{NA} values: \code{stateProvince}, \code{county}, \code{decimalLatitude}, \code{decimalLongitude}, and \code{coordinateUncertaintyInMeters}. The function does strict matching of character values in the \code{stateProvince} and \code{county} fields (but ignores capitalization), so misspellings, diacritics, trailing or leading spaces, etc. can lead to a mismatch between a field in this data frame and the spatial objects.
#' @param geogCounty SpatialPolygonsDataFrame representing the area over which the records in \code{darwin} were assumed to have been collected with the finest level representing "counties" (parishes or similar "secondary-level" administrative units). The data frame component must contain fields named in \code{stateField} and \code{countyGeogField}.
#' @param geogState SpatialPolygonsDataFrame representing the area over which the records in \code{darwin} were assumed to have been collected with the finest level representing "states" or "provinces" ("primary-level" administrative units). The data frame component must contain fields named in \code{stateField}.
#' @param geogCountry SpatialPolygons or SpatialPolygonsDataFrame representing the area over which the records in \code{darwin} were assumed to have been collected with the finest level representing countries.
#' @param minCoordUncerPlusPrecision_m Numeric, smallest value of coordinate uncertainty (in meters) required for a record to be designated a "certain" designation (assuming other checks are OK). Generally, it assumed that records that pass all checks and have a value of coordinate uncertainty + imprecision equal to or less than this value can be represented by a single point (i.e., not a polygon).
#' @param minCoordPrecisionForceCounty_m Numeric. The smallest value of coordinate imprecision (in meters) required to declare the record as a "county" record. If the coordinate precision is larger than this value and no coordinate uncertainty is given, then the record is assumed to be best represented by the county in which it occurs.  If coordinate uncertainty is given, then the record is assigned to either a county or a circle with a radius equal to coordinate uncertainty and coordinate imprecision, whichever has the greater area. If this value is equal to or greater than \code{minCoordPrecisionForceState_m}, then a warning will be given.
#' @param minCoordPrecisionForceState_m Numeric. The smallest value of coordinate imprecision (in meters) required to declare the record as a "state" record. If the coordinate precision is larger than this value and no coordinate uncertainty is given, then the record is assumed to be best represented by the state/province in which it occurs.  If coordinate uncertainty is given, then the record is assigned to either a state/province or a circle with a radius equal to coordinate uncertainty and coordinate imprecision, whichever has the greater area.
#' @param calcDistToCentroids Logical, if \code{TRUE} (default), calculate distance from each record with coordinates to the nearest county, state, and country centroid.
#' @param coordSystemStringsDegMinSec Either \code{NULL} or character vector. This names the possible values in the DarwinCore field "verbatimCoordinateSystem" that correspond to a "degrees-minutes-seconds" reporting of the original coordinates. Coordinate imprecision due to using a degrees-minutes-seconds format will be calculated for records containing any of the strings in this argument. If the imprecision (in meters) is greater than that calculated from the coordinate precision in decimal format then this imprecision will be used to determine if the record can be confidently assigned to a given county or state/province. If \code{NULL}, no checking for imprecision due to degrees-minutes-seconds coordinates is performed.
#' @param coordSystemStringsDegMin Same as \code{coordSystemStringsDegMinSec}, but for records with a \code{verbatimCoordinateSystem} value indicating coordinates are to nearest degrees and minutes.
#' @param coordSystemStringsDeg Same as \code{coordSystemStringsDegMinSec}, but for records with a \code{verbatimCoordinateSystem} value indicating coordinates are to nearest degree.
#' @param eaProj PROJ4 string for an equal-area representation of the spatial polygons in \code{geogCounty} and \code{geogState}. This will look something like \code{'+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'} (which is Albers Equal-Area for North America).
#' @param countyGeogField Name of field (column) in \code{geog} with county names.
#' @param stateGeogField Name of field (column) in \code{geogCounty} and \code{geogState} with state/province names.
#' @param countryGeogField Name of field (column) in \code{geogState} with country names.
#' @param verbose Logical, if \code{TRUE} then display progress.
#' @param ... Other arguments, ignored for now.
#' @return A data frame with these fields:
#' \itemize{
#'		\item	\code{recordType}: Type of record.
#' 		\item	\code{uncerPrecisionArea_km2}: The area of uncertainty associated with the record in km2 (or \code{NA}).
#' 		\item	\code{uncerBasedOn}: A description of the basis for calculating the area of uncertainty.
#' 		\item	\code{adminMatch}: Logical or \code{NA}. If \code{TRUE}, then there is a match between the data in the \code{darwin} data frame and the geographies in \code{geogCounty} and/or \code{geogState}.
#' 		\item	\code{stateFromGeog} and \code{countyFromGeog}: Name of the state/province and county in which the record occurs obtained from the geographies in \code{geogCounty} and \code{geogState}. Note that \code{stateFromGeog} might be given but \code{countyFromGeog} \code{NA} if the record can only be assigned with confidence to a state/province.
#' 		\item	\code{coordUncer_m}: Coordinate uncertainty in meters (copied from the field \code{coordinateUncertaintyInMeters} in the records).
#' 		\item	\code{coordPrecision_m}: Coordinate precision in meters calculated from a best estimate of the precision in the decimal coordinates and/or coordinate reporting system (e.g., degrees-minutes). Note that this can be non-\code{NA} even if \code{coordUncer_m} is \code{NA} if coordinate are supplied but no coordinate uncertainty is given. This is \emph{not} a good replacement for a missing value of \code{coordinateUncertaintyInMeters}.
#'		\item \code{distToNearestCountyCentroid_m}, \code{distToNearestStateCentroid_m}, and \code{distToNearestCountryCentroid_m} (included only if \code{calcDistToCentroids} is \code{TRUE}): Distance between the record and the nearest county, state/province, and country centroids. Records are sometimes assigned to coordinates at or near the centroid of an administrative unit if they can only be located with confidence to that unit. As a result, they might appear precise, even if they are not. These values are useful for determining if that might be the case.
#' }
#' @export

darwinCoreSpatialAssign <- function(
	darwin,
	geogCounty,
	geogState,
	geogCountry,
	eaProj,
	minCoordUncerPlusPrecision_m,
	minCoordPrecisionForceCounty_m,
	minCoordPrecisionForceState_m,
	calcDistToCentroids = TRUE,
	coordSystemStringsDegMinSec = c('degrees minutes seconds', 'DMS', 'deg. min. sec.', 'degrÃ©es minutes secondes', 'grados minutos segundos'),
	coordSystemStringsDegMin = c('degrees minutes', 'deg. sec.', 'degrÃ©es minutes', 'degrÃ©s minutes', 'grados minutos'),
	coordSystemStringsDeg = c('degrees', 'degree', 'deg.', 'degrÃ©es', 'degrÃ©e', 'grados', 'grado'),
	countyGeogField = 'NAME_2',
	stateGeogField = 'NAME_1',
	countryGeogField = 'NAME_0',
	verbose = TRUE,
	...
) {

	### CONTENTS & BAUHAUS
	######################
	
		### CONTENTS
		### catch errors
		### constants
		### state/county names from records and geographies
		### state/county area
		### pre-allocate output
		### simplify input
		### distance to centroids
		### CLASSIFYING RECORDS
			### STATE present | COUNTY present | COORDS present | COORD UNCER present
			### STATE present | COUNTY present | COORDS present | COORD UNCER missing
			### STATE present | COUNTY present | COORDS missing | COORD UNCER (any)
			### STATE present | COUNTY missing | COORDS missing | COORD UNCER (any)
			### STATE missing | COUNTY present | COORDS missing | COORD UNCER (any)
			### STATE present | COUNTY missing | COORDS present | COORD UNCER missing
			### STATE present | COUNTY missing | COORDS present | COORD UNCER present
			### STATE missing | COUNTY present | COORDS present | COORD UNCER missing
			### STATE missing | COUNTY present | COORDS present | COORD UNCER present
			### STATE missing | COUNTY missing | COORDS present | COORD UNCER missing
			### STATE missing | COUNTY missing | COORDS present | COORD UNCER present
			### STATE missing | COUNTY missing | COORDS missing | COORD UNCER any
		
		# Each "classifying records" section begins by subsetting records that match the given conditions (indexed by "recs" which matches to the records data frame the user supplies).
		# If any records meet the criteria, they are subsetted to a data frame names "these".
		# It then generally determines if there is a match between the records' information and the geographies.
		# Records that have no match are assigned "unusable" (in a line at the end of each section).
		# Records that do match are assigned on the basis of whether or not state, county, coordinates, and maybe coordinate uncertainty information is present, and the size of coordinate uncertainty and precision of the coordinates. Depending on the case, several assignments may be made. In these cases, a sub-subset of matching records is defined (indexed by "subRecs" which matches to the records data frame supplied by the user). To facilitate assignments, a second index named "theseSubIndex" is also created which matches to "these" (the data frame that is a subset of the records supplied by the user.

	### catch errors
	################
	
		if (minCoordUncerPlusPrecision_m <= 0) warning('Argument "minCoordUncerPlusPrecision_m" is <= 0, so no records can be designated as "certain".')
		if (minCoordPrecisionForceCounty_m <= minCoordUncerPlusPrecision_m) warning('Argument "minCoordPrecisionForceCounty_m" is <= "minCoordUncerPlusPrecision_m". This seems illogical.')
		if (minCoordPrecisionForceState_m <= minCoordUncerPlusPrecision_m) warning('Argument "minCoordPrecisionForceState_m" is <= "minCoordUncerPlusPrecision_m". This seems illogical.')
		if (minCoordPrecisionForceState_m <= minCoordPrecisionForceCounty_m) warning('Argument "minCoordPrecisionForceState_m" is <= "minCoordPrecisionForceCounty_m". This seems illogical.')
		
		if (!(countyGeogField %in% names(geogCounty))) stop('The column specified in "countyGeogField" must appear in the "geogCounty" object.')
		if (!(stateGeogField %in% names(geogCounty))) stop('The column specified in "stateGeogField" must appear in the "geogCounty" object.')
		if (!(stateGeogField %in% names(geogState))) stop('The column specified in "stateGeogField" must appear in the "geogState" object.')
		
	### constants
	#############
		
		ll <- c('decimalLongitude', 'decimalLatitude')
		numRecords <- nrow(darwin)
		
		eaProj <- sp::CRS(eaProj)
		unprojProj <- sp::CRS(raster::projection(geogCounty))
		
	### state/county names from records and geographies
	###################################################
	
		## names in records
		stateOfDarwin <- tolower(darwin$stateProvince)
		stateCountyOfDarwin <- tolower(paste(darwin$stateProvince, darwin$county))

		## names in geography
		stateOfGeog <- unique(geogState@data[ , stateGeogField])
		countyOfGeog <- geogCounty@data[ , countyGeogField]
		stateCountyOfGeog <- unique(paste(geogCounty@data[ , stateGeogField], geogCounty@data[ , countyGeogField]))
		stateOfGeogLower <- tolower(stateOfGeog)
		stateCountyOfGeogLower <- tolower(stateCountyOfGeog)
			
	### state/county area
	#####################

		if (verbose) omnibus::say('Calculating administrative unit areas...')

		geogStateEa <- sp::spTransform(geogState, eaProj)
		geogCountyEa <- sp::spTransform(geogCounty, eaProj)
		
		stateArea_km2 <- rgeos::gArea(geogStateEa, byid=TRUE) / 1000^2
		countyArea_km2 <- rgeos::gArea(geogCountyEa, byid=TRUE) / 1000^2

		names(stateArea_km2) <- tolower(geogStateEa@data[ , stateGeogField])
		names(countyArea_km2) <- tolower(paste(geogCountyEa@data[ , stateGeogField], geogCountyEa@data[ , countyGeogField]))
		

	### pre-allocate output
	#######################
		
		out <- data.frame(
			# STATE = darwin$stateProvince,
			# COUNTY = darwin$county,
			# LONG = darwin$decimalLongitude,
			# LAT = darwin$decimalLatitude,
			# CU = darwin$coordinateUncertaintyInMeters,
			recordType = rep(NA, numRecords),					# character
			uncerPrecisionArea_km2 = rep(NA, numRecords),		# 0 or positive numeric/NA
			uncerBasedOn = rep(NA, numRecords),					# character/NA
			adminMatch = rep(NA, numRecords),					# TRUE/FALSE/NA
			stateFromGeog = rep(NA, numRecords),				# character/NA
			countyFromGeog = rep(NA, numRecords),				# character/NA
			coordUncer_m = darwin$coordinateUncertaintyInMeters,	# 0 or positive numeric/NA
			coordPrecision_m = rep(NA, numRecords)				# TRUE/FALSE/NA
		)
		
		if (calcDistToCentroids) {
		
			out <- cbind(
				out,
				data.frame(
					distToNearestCountyCentroid_m = rep(NA, numRecords),	# 0 or positive numeric/NA
					distToNearestStateCentroid_m = rep(NA, numRecords),	# 0 or positive numeric/NA
					distToNearestCountryCentroid_m = rep(NA, numRecords)	# 0 or positive numeric/NA
				)
			)
		
		}
				
		
		row.names(out) <- row.names(darwin)

	### simplify input
	##################
		
		darwin <- darwin[ , c('stateProvince', 'county', ll, 'coordinatePrecision', 'coordinateUncertaintyInMeters', 'verbatimCoordinateSystem')]
		geogState@data <- geogState@data[ , which(names(geogState@data) %in% stateGeogField), drop=FALSE]
		geogCounty@data <- geogCounty@data[ , c(which(names(geogCounty@data) %in% stateGeogField), which(names(geogCounty@data) %in% countyGeogField)), drop=FALSE]

	### distance to centroids
	#########################

		if (calcDistToCentroids) {
		
			recs <- which(!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude))
			
			if (length(recs) > 0) {
				
				if (verbose) omnibus::say('Calculating distance to county, state, and country centroids...')

				these <- darwin[recs, ]
				theseSp <- sp::SpatialPoints(these[ , ll], unprojProj)
				
				# distance to county centroids
				geogEa <- sp::spTransform(geogCounty, eaProj)
				centsSpEa <- rgeos::gCentroid(geogEa, byid=TRUE)
				centsSp <- sp::spTransform(centsSpEa, unprojProj)
				
				dists_m <- geosphere::distm(centsSp, theseSp)
				minDist_m <- apply(dists_m, 2, min, na.rm=TRUE)
				
				out$distToNearestCountyCentroid_m[recs] <- minDist_m

				# distance to state centroids
				geogEa <- sp::spTransform(geogState, eaProj)
				centsSpEa <- rgeos::gCentroid(geogEa, byid=TRUE)
				centsSp <- sp::spTransform(centsSpEa, unprojProj)
				
				dists_m <- geosphere::distm(centsSp, theseSp)
				minDist_m <- apply(dists_m, 2, min, na.rm=TRUE)
				
				out$distToNearestStateCentroid_m[recs] <- minDist_m

				# distance to country centroids
				geogEa <- sp::spTransform(geogCountry, eaProj)
				centsSpEa <- rgeos::gCentroid(geogEa, byid=TRUE)
				centsSp <- sp::spTransform(centsSpEa, unprojProj)
				
				dists_m <- geosphere::distm(centsSp, theseSp)
				minDist_m <- apply(dists_m, 2, min, na.rm=TRUE)
				
				out$distToNearestCountryCentroid_m[recs] <- minDist_m
				
			}
			
		}
				
	### CLASSIFYING RECORDS
	#######################
	
			
	### STATE present | COUNTY present | COORDS present | COORD UNCER present
	#########################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			!is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE present COUNTY present COORDINATES present COORD UNCER present: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			coordUncerPlusPrecision_m <- theseCoordPrecision_m + these$coordinateUncertaintyInMeters
			
			## distinguish matching/not matching
			doRecsMatch <- stateCountyOfDarwin[recs] %in% stateCountyExtractFromGeog
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: certain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'certain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
						
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'uncertain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(countyArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(countyArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE present | COUNTY present | COORDS present | COORD UNCER missing
	#########################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE present COUNTY present COORDINATES present COORD UNCER missing: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			
			## distinguish matching/not matching
			doRecsMatch <- stateCountyOfDarwin[recs] %in% stateCountyExtractFromGeog
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: precise at county
			subRecs <- recs[theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/precise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at county
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at state
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)

				uncerAreaOfMatches_km2 <- stateArea_km2[match(stateOfDarwin[subRecs], names(stateArea_km2))]

				out$recordType[subRecs] <- 'state/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'state area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE present | COUNTY present | COORDS missing | COORD UNCER (any)
	#######################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			(is.na(darwin$decimalLongitude) | is.na(darwin$decimalLatitude))
		)
		
		if (verbose) omnibus::say('STATE present COUNTY present COORDINATES missing COORD UNCER any: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]
			
			## distinguish matching/not matching
			doRecsMatch <- stateCountyOfDarwin[recs] %in% stateCountyOfGeogLower
			recsMatch <- recs[doRecsMatch]
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch
			
			# matches
			if (any(doRecsMatch)) {
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[recsMatch], names(countyArea_km2))]

				out$recordType[recsMatch] <- 'county-only'
				out$uncerPrecisionArea_km2[recsMatch] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[recsMatch] <- 'county area'
				out$stateFromGeog[recsMatch] <- stateOfGeog[match(stateOfDarwin[recsMatch], stateOfGeogLower)]
				out$countyFromGeog[recsMatch] <- countyOfGeog[match(stateCountyOfDarwin[recsMatch], stateCountyOfGeogLower)]
				
			}
			
			# mismatches
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE present | COUNTY missing | COORDS missing | COORD UNCER (any)
	#######################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			(is.na(darwin$decimalLongitude) | is.na(darwin$decimalLatitude))
		)
		
		if (verbose) omnibus::say('STATE present COUNTY missing COORDINATES missing COORD UNCER any: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]
			
			## distinguish matching/not matching
			doRecsMatch <- stateOfDarwin[recs] %in% stateOfGeogLower
			recsMatch <- recs[doRecsMatch]
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch
			
			# matches
			if (any(doRecsMatch)) {
			
				stateExtractFromGeog <- unique(geogState@data[ , stateGeogField])
				uncerAreaOfMatches_km2 <- stateArea_km2[match(stateOfDarwin[recsMatch], names(stateArea_km2))]

				out$recordType[recsMatch] <- 'state-only'
				out$uncerPrecisionArea_km2[recsMatch] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[recsMatch] <- 'state area'
				out$stateFromGeog[recsMatch] <- stateOfGeog[match(stateOfDarwin[recsMatch], stateOfGeogLower)]
				
			}
			
			# mismatches
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE missing | COUNTY present | COORDS missing | COORD UNCER (any)
	#######################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			(is.na(darwin$decimalLongitude) | is.na(darwin$decimalLatitude))
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY present COORDINATES missing COORD UNCER any: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]
			
			## distinguish matching/not matching
			countyFromGeog <- stateMatchesInGeog <- rep(NA, length(recs))
			
			for (i in seq_along(recs)) {
			
				matchingGeogCountyIndex <- which(tolower(these$county[i]) == tolower(geogCounty@data[ , countyGeogField]))
				thisState <- geogCounty@data[matchingGeogCountyIndex, stateGeogField]
				if (length(thisState) == 1) {
					stateMatchesInGeog[i] <- thisState
					countyFromGeog[i] <- geogCounty@data[matchingGeogCountyIndex, countyGeogField]
				}
			
			}
			
			doRecsMatch <- !is.na(stateMatchesInGeog)
			recsMatch <- recs[doRecsMatch]
			recsNoMatch <- recs[!doRecsMatch]

			stateFromGeogCountyFromGeog <- tolower(paste(stateMatchesInGeog[recs %in% recsMatch], countyFromGeog[recs %in% recsMatch]))
			
			## assign
			out$adminMatch[recs] <- doRecsMatch
			
			# matches
			if (any(doRecsMatch)) {
			
				theseSubIndex <- which(recs %in% recsMatch)
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateFromGeogCountyFromGeog, names(countyArea_km2))]

				out$recordType[recsMatch] <- 'county-only'
				out$uncerPrecisionArea_km2[recsMatch] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[recsMatch] <- 'county area'
				out$stateFromGeog[recsMatch] <- stateMatchesInGeog[theseSubIndex]
				out$countyFromGeog[recsMatch] <- countyFromGeog[theseSubIndex]
				
			}
			
			# mismatches
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE present | COUNTY missing | COORDS present | COORD UNCER missing
	#########################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE present COUNTY missing COORDINATES present COORD UNCER missing: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))
			stateFromDarwinCountyFromGeog <- tolower(paste(these$stateProvince, countyExtractFromGeog))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			
			## distinguish matching/not matching
			doRecsMatch <- stateFromDarwinCountyFromGeog %in% stateCountyOfGeogLower
			recsMatch <- recs[doRecsMatch]
			recsNoMatch <- recs[!doRecsMatch]
			
			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: precise at county
			subRecs <- recs[theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/precise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at county
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at state
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				uncerAreaOfMatches_km2 <- stateArea_km2[match(stateOfDarwin[subRecs], names(stateArea_km2))]

				theseSubIndex <- which(recs %in% subRecs)
				out$recordType[subRecs] <- 'state/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'state area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE present | COUNTY missing | COORDS present | COORD UNCER present
	#########################################################################

		recs <- which(
			!is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			!is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE present COUNTY missing COORDINATES present COORD UNCER present: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))
			stateFromDarwinCountyFromGeog <- tolower(paste(these$stateProvince, countyExtractFromGeog))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			coordUncerPlusPrecision_m <- theseCoordPrecision_m + these$coordinateUncertaintyInMeters
			
			## distinguish matching/not matching
			doRecsMatch <- stateFromDarwinCountyFromGeog %in% stateCountyOfGeogLower
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: certain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'certain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
						
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'uncertain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(countyArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}

	### STATE missing | COUNTY present | COORDS present | COORD UNCER missing
	#########################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY present COORDINATES present COORD UNCER missing: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			
			## distinguish matching/not matching
			doRecsMatch <- stateCountyOfDarwin[recs] %in% stateCountyExtractFromGeog
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: precise at county
			subRecs <- recs[theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/precise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at county
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyOfDarwin[subRecs], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at state
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)

				uncerAreaOfMatches_km2 <- stateArea_km2[match(stateOfDarwin[subRecs], names(stateArea_km2))]

				out$recordType[subRecs] <- 'state/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'state area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}

	### STATE missing | COUNTY present | COORDS present | COORD UNCER present
	#########################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			!is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			!is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY present COORDINATES present COORD UNCER present: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateFromGeogCountyFromDarwin <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], these$county))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			coordUncerPlusPrecision_m <- theseCoordPrecision_m + these$coordinateUncertaintyInMeters
			
			## distinguish matching/not matching
			doRecsMatch <- stateCountyOfDarwin[recs] %in% stateFromGeogCountyFromDarwin
			recsNoMatch <- recs[!doRecsMatch]

			## assign
			out$adminMatch[recs] <- doRecsMatch

			# matches: certain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'certain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
						
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'uncertain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(countyArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}
		
	### STATE missing | COUNTY missing | COORDS present | COORD UNCER missing
	#########################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY missing COORDINATES present COORD UNCER missing: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateCountyExtractFromGeog <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], stateCountyExtractFromGeog[ , countyGeogField]))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			
			## distinguish matching/not matching
			doRecsMatch <- !(is.na(stateExtractFromGeog) | is.na(countyExtractFromGeog))
			recsNoMatch <- recs[!doRecsMatch]

			# ## assign
			# out$adminMatch[recs] <- doRecsMatch

			# matches: precise at county
			subRecs <- recs[theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyExtractFromGeog[theseSubIndex], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/precise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at county
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				uncerAreaOfMatches_km2 <- countyArea_km2[match(stateCountyExtractFromGeog[theseSubIndex], names(countyArea_km2))]

				out$recordType[subRecs] <- 'county/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$uncerBasedOn[subRecs] <- 'county area'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: imprecise at state
			subRecs <- recs[theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)

				uncerAreaOfMatches_km2 <- stateArea_km2[match(tolower(stateExtractFromGeog[theseSubIndex]), names(stateArea_km2))]

				out$recordType[subRecs] <- 'state/imprecise'
				out$uncerPrecisionArea_km2[subRecs] <- uncerAreaOfMatches_km2
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$uncerBasedOn[subRecs] <- 'state area'
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}

	### STATE missing | COUNTY missing | COORDS present | COORD UNCER present
	#########################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			!is.na(darwin$decimalLongitude) & !is.na(darwin$decimalLatitude) &
			!is.na(darwin$coordinateUncertaintyInMeters)
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY missing COORDINATES present COORD UNCER present: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			these <- darwin[recs, ]

			## extract state/county from geography
			stateCountyExtractFromGeog <- raster::extract(geogCounty, these[ , ll])
			stateExtractFromGeog <- stateCountyExtractFromGeog[ , stateGeogField]
			countyExtractFromGeog <- stateCountyExtractFromGeog[ , countyGeogField]
			stateFromGeogCountyFromDarwin <- tolower(paste(stateCountyExtractFromGeog[ , stateGeogField], these$county))

			## coordinate precision
			theseCoordPrecision <- .subCoordPrecision(these, unprojProj, ll, coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg)
			theseCoordPrecision_m <- theseCoordPrecision$coordPrecision_m
			coordUncerPlusPrecision_m <- theseCoordPrecision_m + these$coordinateUncertaintyInMeters
			
			## distinguish matching/not matching
			doRecsMatch <- !(is.na(stateExtractFromGeog) | is.na(countyExtractFromGeog))
			recsNoMatch <- recs[!doRecsMatch]

			# ## assign
			# out$adminMatch[recs] <- doRecsMatch

			# matches: certain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'certain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/precise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m < minCoordPrecisionForceCounty_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
						
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				out$recordType[subRecs] <- 'uncertain/precise'
				out$uncerPrecisionArea_km2[subRecs] <- subTheseUncerSpEaArea_km2
				out$uncerBasedOn[subRecs] <- 'circle with radius = coordinate uncertainty + precision'
				out$coordPrecision_m[subRecs] <- theseCoordPrecision_m[theseSubIndex]
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(countyArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: certain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters <= minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'certain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'certain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'certain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'certain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or county/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceCounty_m & theseCoordPrecision_m < minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]

			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- countyArea_km2[match(stateCountyOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'county/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'county area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				out$countyFromGeog[subRecs] <- countyExtractFromGeog[theseSubIndex]
				
			}
			
			# matches: uncertain/imprecise or state/imprecise
			subRecs <- recs[these$coordinateUncertaintyInMeters > minCoordUncerPlusPrecision_m & theseCoordPrecision_m >= minCoordPrecisionForceState_m]
			subRecs <- subRecs[subRecs %in% recs[doRecsMatch]]
			
			if (length(subRecs) > 0) {
			
				theseSubIndex <- which(recs %in% subRecs)
			
				subThese <- these[theseSubIndex, ll]
				subTheseSp <- sp::SpatialPoints(subThese, proj4string=unprojProj)
				subTheseSpEa <- sp::spTransform(subTheseSp, eaProj)
				subTheseUncerSpEaBuffer <- rgeos::gBuffer(subTheseSpEa, width=coordUncerPlusPrecision_m[theseSubIndex], byid=TRUE)
				subTheseUncerSpEaArea_km2 <- rgeos::gArea(subTheseUncerSpEaBuffer, byid=TRUE) / 1000^2
				
				subAdminAreas_km2 <- stateArea_km2[match(stateOfDarwin[recs[theseSubIndex]], names(stateArea_km2))]
				
				designation <- ifelse(subTheseUncerSpEaArea_km2 > subAdminAreas_km2, 'uncertain/imprecise', 'state/imprecise')
				
				out$recordType[subRecs] <- designation
				out$uncerPrecisionArea_km2[subRecs] <- ifelse(designation == 'uncertain/imprecise', subTheseUncerSpEaArea_km2, subAdminAreas_km2)
				out$uncerBasedOn[subRecs] <- ifelse(designation == 'uncertain/imprecise', 'circle with radius = coordinate uncertainty + precision', 'state area')
				out$coordPrecision_m[subRecs] <- ifelse(designation == 'uncertain/imprecise', theseCoordPrecision_m[theseSubIndex], NA)
				out$stateFromGeog[subRecs] <- stateExtractFromGeog[theseSubIndex]
				
			}
			
			if (any(!doRecsMatch)) out$recordType[recsNoMatch] <- 'unusable'

		}

	### STATE missing | COUNTY missing | COORDS missing | COORD UNCER any
	#####################################################################

		recs <- which(
			is.na(darwin$stateProvince) &
			is.na(darwin$county) &
			(is.na(darwin$decimalLongitude) | !is.na(darwin$decimalLatitude))
		)
		
		if (verbose) omnibus::say('STATE missing COUNTY missing COORDINATES missing COORD UNCER any: ', length(recs), ' record(s).')
		
		if (length(recs) > 0) {
		
			out$recordType[recs] <- 'unusable'
			
		}

	### REPORTING
	if (verbose) {
	
		size <- nchar(numRecords)
		
		say('Final assignments (number of records):', pre=2)
		
		out <<- out
		
		say('certain/precise................ ', prefix(sum(out$recordType %==na% 'certain/precise'), size, pad=' '))
		say('uncertain/precise.............. ', prefix(sum(out$recordType %==na% 'uncertain/precise'), size, pad=' '))
		say('certain/imprecise.............. ', prefix(sum(out$recordType %==na% 'certain/imprecise'), size, pad=' '))
		say('uncertain/imprecise............ ', prefix(sum(out$recordType %==na% 'uncertain/imprecise'), size, pad=' '))
		say('county/precise................. ', prefix(sum(out$recordType %==na% 'county/precise'), size, pad=' '))
		say('county/imprecise............... ', prefix(sum(out$recordType %==na% 'county/imprecise'), size, pad=' '))
		say('state/imprecise................ ', prefix(sum(out$recordType %==na% 'state/imprecise'), size, pad=' '))
		say('county-only.................... ', prefix(sum(out$recordType %==na% 'county-only'), size, pad=' '))
		say('state-only..................... ', prefix(sum(out$recordType %==na% 'state-only'), size, pad=' '))
		say('unusable....................... ', prefix(sum(out$recordType %==na% 'unusable'), size, pad=' '))
		say('not assignable (should be 0)... ', prefix(sum(is.na(out$recordType)), size, pad=' '))
		say('TOTAL.......................... ', numRecords)
	
	}
		
	out
		
}

		
		
		
### calculate coordinate precision from coordinates and from verbatimCoordinateSystem
.subCoordPrecision <- function(
	records,
	unprojProj,
	ll = c('decimalLongitude', 'decimalLatitude'),
	coordSystemStringsDegMinSec = c('degrees minutes seconds', 'DMS', 'deg. min. sec.', 'degrÃ©es minutes secondes', 'grados minutos segundos'),
	coordSystemStringsDegMin = c('degrees minutes', 'deg. sec.', 'degrÃ©es minutes', 'degrÃ©s minutes', 'grados minutos'),
	coordSystemStringsDeg = c('degrees', 'degree', 'deg.', 'degrÃ©es', 'degrÃ©e', 'grados', 'grado')
) {

	# records			data frame with specimen records in DarwinCore format
	# unprojProj

	# ll				two-element character vector with names of longitude/latitude fields in "records"
	# coordSystemStringsDegMinSec, coordSystemStringsDegMin, coordSystemStringsDeg Either \code{NULL} or character vector. This names the possible values in the DarwinCore field "verbatimCoordinateSystem" that correspond to a "degrees-minutes-seconds", "degrees-minutes" or just "degrees" reporting of the original coordinates. Coordinate imprecision due to using a degrees-minutes-seconds format will be calculated for records containing any of the strings in this argument. If the imprecision (in meters) is greater than that calculated from the coordinate precision in decimal format then this imprecision will be used to determine if the record can be confidently assigned to a given county or state/province. If \code{NULL}, no checking for imprecision due to original coordinate system is performed.
	
	### coordinate precision from coordinates

		# number of decimal digits of coordinates (if NA then returns NA)
		numCoordDigitsLong <- omnibus::roundedSigDigits(records[ , ll[1]])
		numCoordDigitsLat <- omnibus::roundedSigDigits(records[ , ll[2]])
		
		numCoordDigitsLong <- -1 * numCoordDigitsLong
		numCoordDigitsLat <- -1 * numCoordDigitsLat
		
		numCoordDigits <- pmax(numCoordDigitsLong, numCoordDigitsLat)
		if (any(!is.na(numCoordDigits)) && any(numCoordDigits < 0)) numCoordDigits[!is.na(numCoordDigits) & numCoordDigits < 0] <- 0
		
		# coordinate precision ascertained directly from coordinates
		theseCoordDigits <- numCoordDigits
		records$decimalLongitude <- trunc(records$decimalLongitude * 10^theseCoordDigits) / 10^theseCoordDigits
		records$decimalLatitude <- trunc(records$decimalLatitude * 10^theseCoordDigits) / 10^theseCoordDigits

		theseSp <- sp::SpatialPoints(records[ , ll], proj4string=unprojProj)
		
		coordPrecision_m <- enmSdm::coordPrecision(theseSp)
		
		subOut <- data.frame(
			coordPrecision_digits = numCoordDigits,
			coordPrecisionFrom = 'coordinates',
			coordPrecision_m = coordPrecision_m
		)

	### coordinate precision ascertained from verbatim coordinate system

		# degrees-minutes-seconds
		if (!is.null(coordSystemStringsDegMinSec)) {
		
			recs <- which(records$verbatimCoordinateSystem %in% coordSystemStringsDegMinSec)
			
			if (length(recs) > 0) {

				coords <- records[recs, ll]
				coordsPlusHalfSec <- coords + (cbind(0.5, 0.5) / 3600)
				coordsMinusHalfSec <- coords - (cbind(0.5, 0.5) / 3600)
				theseSpPlusHalfSec <- sp::SpatialPoints(coordsPlusHalfSec, proj4string=unprojProj)
				theseSpMinusHalfSec <- sp::SpatialPoints(coordsMinusHalfSec, proj4string=unprojProj)
				
				coordPrecisionFromDms_m <- rep(NA, length(recs))
				for (i in seq_along(recs)) coordPrecisionFromDms_m[i] <- geosphere::distGeo(theseSpPlusHalfSec[i], theseSpMinusHalfSec[i])
				
				coordsPrecisionFromDms <- omnibus::whichPMax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE) == 2
				subOut$coordPrecisionFrom[recs] <- ifelse(coordsPrecisionFromDms, 'original coordinates in degrees-minutes-seconds', subOut$coordPrecisionFrom[recs])
				subOut$coordPrecision_m[recs] <- pmax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE)
			
			}
			
		}

		# degrees-minutes
		if (!is.null(coordSystemStringsDegMin)) {
		
			recs <- which(records$verbatimCoordinateSystem %in% coordSystemStringsDegMin & complete.cases(records[ , ll]))
			
			if (length(recs) > 0) {
			
				coords <- records[recs, ll]
				coordsPlusHalfSec <- coords + (cbind(0.5, 0.5) / 60)
				coordsMinusHalfSec <- coords - (cbind(0.5, 0.5) / 60)
				theseSpPlusHalfSec <- sp::SpatialPoints(coordsPlusHalfSec, proj4string=unprojProj)
				theseSpMinusHalfSec <- sp::SpatialPoints(coordsMinusHalfSec, proj4string=unprojProj)
				
				coordPrecisionFromDms_m <- rep(NA, length(recs))
				for (i in seq_along(recs)) coordPrecisionFromDms_m[i] <- geosphere::distGeo(theseSpPlusHalfSec[i], theseSpMinusHalfSec[i])
				
				coordsPrecisionFromDms <- omnibus::whichPMax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE) == 2
				subOut$coordPrecisionFrom[recs] <- ifelse(coordsPrecisionFromDms, 'original coordinates in degrees-minutes', subOut$coordPrecisionFrom[recs])
				subOut$coordPrecision_m[recs] <- pmax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE)
			
			}
			
		}
		
		# degrees
		if (!is.null(coordSystemStringsDeg)) {
		
			recs <- which(records$verbatimCoordinateSystem %in% coordSystemStringsDeg & complete.cases(records[ , ll]))
			
			if (length(recs) > 0) {
			
				coords <- records[recs, ll]
				coordsPlusHalfSec <- coords + cbind(0.5, 0.5)
				coordsMinusHalfSec <- coords - cbind(0.5, 0.5)
				theseSpPlusHalfSec <- sp::SpatialPoints(coordsPlusHalfSec, proj4string=unprojProj)
				theseSpMinusHalfSec <- sp::SpatialPoints(coordsMinusHalfSec, proj4string=unprojProj)
				
				coordPrecisionFromDms_m <- rep(NA, length(recs))
				for (i in seq_along(recs)) coordPrecisionFromDms_m[i] <- geosphere::distGeo(theseSpPlusHalfSec[i], theseSpMinusHalfSec[i])
				
				coordsPrecisionFromDms <- whichPMax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE) == 2
				subOut$coordPrecisionFrom[recs] <- ifelse(coordsPrecisionFromDms, 'original coordinates in degrees', subOut$coordPrecisionFrom[recs])
				subOut$coordPrecision_m[recs] <- pmax(subOut$coordPrecision_m[recs], coordPrecisionFromDms_m[recs], na.rm=TRUE)
			
			}
			
		}

	subOut
	
}
