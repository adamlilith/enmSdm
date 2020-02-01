# darwinCoreSpatialAssign DEVELOPMENT
# source('D:/ecology/Drive/R/enmSdm/working/darwinCoreSpatialAssign_development.r')


library(omnibus)
library(enmSdm)
library(dismo)
library(data.table)

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	gc()
	options(stringsAsFactors=FALSE)
	
	drive <- 'C:'
	# drive <- 'D:'
	
	setwd(paste0(drive, '/Ecology/Drive/Research/Vaguely Georeferenced Specimen Records'))

	### choice variables
		
		# minimum coordinate uncertainty in meters to be considered an "precise" record
		minCoordUncerForAccurate_m <- 5000
		
		# using records that were collected starting in...
		startYear <- 1970
		
		# number of digits to which to round coordinates to determine if they are spatial duplicates
		proximateDigits <- 3 # 3 digits corresponds to ~ 100 m
		
		# minimum number of precise records for a species to be included in the analysis
		minNumAccRecords <- 5
	
	### names

		llGbif <- c('decimalLongitude', 'decimalLatitude')
	
	### data objects
	
		if (file.exists('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 0.rda')) load('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 0.rda')
		if (file.exists('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 1.rda')) load('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 1.rda')
		if (file.exists('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 2.rda')) load('C:/ecology/!Scratch/GADM Ver 3pt6 North America Level 2.rda')

		records <- fread('./Data/GBIF Asclepias 2019-11-18/occurrence.txt', header=TRUE)
		records <- as.data.frame(records)

	say('### randomly re-order records in case any operations depend on order ###')
	###############################################################################
	
		set.seed(123)
		records <- records[sample(1:nrow(records), nrow(records)), ]
		
	say('### force empty strings to NA and rename countries ###')
	#############################################################
		
		records$species[records$species == ''] <- NA
		records$country[records$country == ''] <- NA
		records$stateProvince[records$stateProvince == ''] <- NA
		records$county[records$county == ''] <- NA

		records$country[records$country == 'CA'] <- 'Canada'
		records$country[records$country == 'US'] <- 'United States'
		records$country[records$country == 'MX'] <- 'Mexico'

		records$species <- trim(records$species)
		records$country <- trim(records$country)
		records$stateProvince <- trim(records$stateProvince)
		records$county <- trim(records$county)

	say('### clean state/province names for ALL records ###')
	#########################################################

		# # ### spelling corrections... these aren't caught below; I don't know why
		# this <- which(records$stateProvince == 'Chinuahua')
		# if (length(this) > 0) records$stateProvince[this] <- 'Chihuahua'
	
		### assign countries/states/counties to records with coordinates with NAs for country/state/county
		haveCoords <- which(complete.cases(records[ , llGbif]))
		admins <- raster::extract(nam2Sp, records[haveCoords, llGbif])

		for (countHaveCoords in seq_along(haveCoords)) {
		
			if (is.na(records$country[haveCoords[countHaveCoords]])) records$country[haveCoords[countHaveCoords]] <- admins$NAME_0[countHaveCoords]
			if (is.na(records$stateProvince[haveCoords[countHaveCoords]])) records$stateProvince[haveCoords[countHaveCoords]] <- admins$NAME_1[countHaveCoords]
			if (is.na(records$county[haveCoords[countHaveCoords]])) records$county[haveCoords[countHaveCoords]] <- admins$NAME_2[countHaveCoords]
		
		}
		
		# what state/province names don't match with GADM?
		gadmStateProv <- nam2Sp@data$NAME_0
		
		gadmCombined <- paste(gadmStateProv, nam2Sp@data$NAME_1)
		gadmCombined <- gadmCombined[!duplicated(gadmCombined)]
		gadmCombined <- sort(gadmCombined)
		gadmCombined <- tolower(gadmCombined)

		recordsCombined <- data.frame(
			country = records$country,
			stateProvince = records$stateProvince
		)
		
		recordsCombined$combined <- paste(recordsCombined$country, recordsCombined$stateProvince)
		recordsCombined$combined <- tolower(recordsCombined$combined)
		
		recordsCombined <- recordsCombined[!duplicated(recordsCombined$combined), ]
		recordsCombined <- recordsCombined[order(recordsCombined$combined), ]
		
		bads <- data.frame()
		for (i in 1:nrow(recordsCombined)) {
			if (!(recordsCombined$combined[i] %in% gadmCombined)) bads <- rbind(bads, recordsCombined[i, ])
		}
		
		# write.csv(bads, './Data/Species Record Cleaning/Mismatched State or Province vis-a-vis GADM.csv', row.names=FALSE)
		
		# MAKE CORRECTIONS using crosswalk created manually from previous file
		corrections <- read.csv('./Data/Species Record Cleaning/Mismatched State or Province vis-a-vis GADM - Corrections MANUALLY CREATED.csv')

		for (i in 1:nrow(corrections)) {
		
			recordsIndex <- which((records$country %in% corrections$country[i]) & (records$stateProvince %in% corrections$stateProvince[i]))
			records$stateProvince[recordsIndex] <- corrections$gadmStateProvince[i]
			
		}
		
	say('### clean county names for ALL records ###')
	#################################################
		
		# removes extraneous text
		# x is a character vector
		removeExtraText <- function(x) {
		
			bads <- c('(', ')', ' county', ' County', ' COUNTY', 'County of', 'Par.', 'Parish', 'Cty.', 'Cty', 'Co.', 'Municipio', 'Municipality', '[', ']', 'Co ')
			
			for (bad in bads) x <- gsub(x, pattern=bad, replacement='', fixed=TRUE)
			
			saints <- c('Ste. ', 'Ste ', 'St. ', 'St ')
			for (saint in saints) x <- gsub(x, pattern=saint, replacement='Saint ', fixed=TRUE)
			
			x <- raster::trim(x)
			x[x == ''] <- NA
			x
		}
		
		records$county <- removeExtraText(records$county)

		# what state/province names don't match with GADM?
		gadmCombined <- paste(nam2Sp@data$NAME_1, nam2Sp@data$NAME_2)
		gadmCombined <- gadmCombined[!duplicated(gadmCombined)]
		gadmCombined <- sort(gadmCombined)
		gadmCombined <- tolower(gadmCombined)

		recordsCombined <- data.frame(
			stateProvince = records$stateProvince,
			county = records$county
		)
		
		recordsCombined$combined <- paste(recordsCombined$stateProvince, recordsCombined$county)
		recordsCombined$combined <- tolower(recordsCombined$combined)
		
		recordsCombined <- recordsCombined[!duplicated(recordsCombined$combined), ]
		recordsCombined <- recordsCombined[order(recordsCombined$combined), ]
		
		bads <- data.frame()
		for (i in 1:nrow(recordsCombined)) {
			if (!(recordsCombined$combined[i] %in% gadmCombined)) bads <- rbind(bads, recordsCombined[i, ])
		}
		
		# write.csv(bads, './Data/Species Record Cleaning/Mismatched State-County with GADM.csv', row.names=FALSE)

		# MAKE CORRECTIONS using corrections created manually from previous file
		corrections <- read.csv('./Data/Species Record Cleaning/Mismatched State-County with GADM - Corrections MANUALLY CREATED.csv')
		
		records$notes <- NA

		for (crosswalkIndex in 1:nrow(corrections)) {

			# if (!is.na(corrections$notes[crosswalkIndex]) && corrections$notes[crosswalkIndex] != 'Possibly municipality-level record') {
			
				recordsIndex <- which((records$stateProvince %in% corrections$stateProvince[crosswalkIndex]) & (records$county %in% corrections$county[crosswalkIndex]))
				records$stateProvince[recordsIndex] <- corrections$gadmStateProvince[crosswalkIndex]
				records$county[recordsIndex] <- corrections$gadmCounty[crosswalkIndex]
				
				# notes
				records$notes[recordsIndex] <- paste0(records$notes[recordsIndex], '; ', corrections$notes[crosswalkIndex])
				
			# }
				
		}

darwin <- records
# darwin <- darwin[darwin$stateProvince %in% c('California', 'Nevada', 'Oregon', 'Utah'), ]

usa1 <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
usa2 <- getData('GADM', country='USA', level=2, path='C:/ecology/!Scratch')

west1 <- subset(usa1, NAME_1 %in% c('California', 'Nevada', 'Oregon'))
west2 <- subset(usa2, NAME_1 %in% c('California', 'Nevada', 'Oregon'))

# arguments
geogCounty <- usa2
geogState <- usa1
eaProj <- getCRS('albersNA')
minCoordUncerForPrecise_m <- 5000
precisionUncerForceCounty_m = 1000
precisionUncerForceState_m = 5000
countyGeogField = 'NAME_2'
stateGeogField = 'NAME_1'
verbose = TRUE

# NOW STEP THRU FUNCTION


