# darwinCoreSpatialAssign DEVELOPMENT
# source('C:/ecology/Drive/R/enmSdm/working/darwinCoreSpatialAssign_development.r')
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

darwin <- records
# darwin <- darwin[darwin$stateProvince %in% c('California', 'Nevada', 'Oregon', 'Utah'), ]

usa0 <- getData('GADM', country='USA', level=0, path='C:/ecology/!Scratch')
usa1 <- getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
usa2 <- getData('GADM', country='USA', level=2, path='C:/ecology/!Scratch')

usa1 <- usa1[usa1$NAME_1 != 'Alaska', ]
usa1 <- usa1[usa1$NAME_1 != 'Hawaii', ]

usa2 <- usa2[usa2$NAME_1 != 'Alaska', ]
usa2 <- usa2[usa2$NAME_1 != 'Hawaii', ]

mex0 <- getData('GADM', country='MEX', level=0, path='C:/ecology/!Scratch')
mex1 <- getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
mex2 <- getData('GADM', country='MEX', level=2, path='C:/ecology/!Scratch')

nam0 <- rbind(usa0, mex0)
nam1 <- rbind(usa1, mex1)
nam2 <- rbind(usa2, mex2)

# west1 <- subset(usa1, NAME_1 %in% c('California', 'Nevada', 'Oregon'))
# west2 <- subset(usa2, NAME_1 %in% c('California', 'Nevada', 'Oregon'))

# arguments
geogCounty <- nam2
geogState <- nam1
geogCountry <- nam0
eaProj <- getCRS('albersNA')
minCoordUncerPlusPrecision_m <- 1000
minCoordPrecisionForceCounty_m <- 2000
minCoordPrecisionForceState_m <- 5000
countyGeogField = 'NAME_2'
stateGeogField = 'NAME_1'
countryGeogField = 'NAME_0'
verbose = TRUE

coordSystemStringsDegMinSec <- c('degrees minutes seconds', 'DMS', 'deg. min. sec.', 'degrÃ©es minutes secondes', 'degrÃ©s minutes secondes', 'grados minutos segundos')
coordSystemStringsDegMin <- c('degrees minutes', 'deg. sec.', 'degrÃ©es minutes', 'degrÃ©s minutes', 'grados minutos')
coordSystemStringsDeg <- c('degrees', 'deg.', 'degrÃ©es', 'degrÃ©s', 'grados')

# NOW STEP THRU FUNCTION

source('C:/Ecology/Drive/R/enmSdm/working/darwinCoreSpatialAssign.r')

out <- darwinCoreSpatialAssign(
	darwin=records,
	geogCounty=nam2,
	geogState=nam1,
	geogCountry=nam0,
	eaProj=getCRS('albersNA'),
	minCoordUncerPlusPrecision_m = 1000,
	minCoordPrecisionForceCounty_m = 5000,
	minCoordPrecisionForceState_m = 10000,
	calcDistToCentroids = FALSE,
	coordSystemStringsDegMinSec = c('degrees minutes seconds', 'DMS', 'deg. min. sec.', 'degrÃ©es minutes secondes', 'grados minutos segundos'),
	coordSystemStringsDegMin = c('degrees minutes', 'deg. sec.', 'degrÃ©es minutes', 'degrÃ©s minutes', 'grados minutos'),
	coordSystemStringsDeg = c('degrees', 'degree', 'deg.', 'degrÃ©es', 'degrÃ©e', 'grados', 'grado'),
	countyGeogField = 'NAME_2',
	stateGeogField = 'NAME_1',
	countryGeogField = 'NAME_0',
	verbose = TRUE
)


