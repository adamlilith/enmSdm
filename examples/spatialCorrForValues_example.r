#' \donttest{
#' # get rasters for mean annual temperature and total precipitation
#' worldClim <- raster::getData('worldclim', var='bio', res=10)
#' worldClim <- raster::subset(worldClim, c(1, 12))
#' 
#' # crop to Madagascar
#' data(mad0)
#' madClim <- raster::crop(worldClim, mad0)
#' 
#' # remove non-Malagasy islands
#' madRast <- raster::rasterize(mad0, madClim)
#' madClim <- madClim * madRast
#' names(madClim) <- c('bio1', 'bio12')
#' madClim[['bio1']] <- madClim[['bio1']] / 10 # to deg C
#' 
#' ### spatial autocorrelation for raster (can take a long time!)
#' sacRast <- localSpatialCorrForValues(x=madClim, focal=NULL)
#' sacRast <- sacRast / 1000 # convert to km
#' 
#' # plot
#' par(mfrow=c(2, 2))
#' raster::plot(madClim[['bio1']], main='BIO 01\n(10 * deg C)')
#' raster::plot(madClim[['bio12']], main='BIO 12\n(mm)')
#' raster::plot(sacRast[['bio1']], main='BIO 01\nDistance (km)')
#' raster::plot(sacRast[['bio12']], main='BIO 12\nDistance (km)')
#' 
#' ### spatial autocorrelation for spatial points
#' # faster than rasters, but still takes a while
#' set.seed(123)
#' x <- dismo::randomPoints(madClim, 200)
#' env <- raster::extract(madClim, x)
#' x <- cbind(x, env)
#' breaks <- c(0, 100000, 10)
#' sacPoints <- localSpatialCorrForValues(x=x, breaks=breaks)
#' 
#' # plot: code point color by characteristic distance of spatial autocorrelation
#' maxSacBio1 <- max(sacPoints$sacDistMid_bio1, na.rm=TRUE)
#' maxSacBio12 <- max(sacPoints$sacDistMid_bio12, na.rm=TRUE)
#' sacBio1 <- 4 * sacPoints$sacDistMid_bio1 / maxSacBio1
#' sacBio12 <- 4 * sacPoints$sacDistMid_bio12 / maxSacBio12
#' 
#' # plot: notice points along transitional zones
#' # have smaller distances, as expected
#' par(mfrow=c(1, 2))
#' leg <- '\n(small: short dist, large: long dist)'
#' sp::plot(mad0, main=paste0('BIO 01', leg))
#' raster::plot(madClim[['bio1']], add=TRUE)
#' points(sacPoints, pch=1, cex=sacBio1, col='red')
#' sp::plot(mad0, main=paste0('BIO 12', leg))
#' raster::plot(madClim[['bio12']], add=TRUE)
#' points(sacPoints, pch=1, cex=sacBio12, col='blue')
#' #'
#' par(mfrow=c(1, 1))
#' maxDist <- max(c(sacPoints$sacDistMid_bio1, sacPoints$sacDistMid_bio12),
	#' na.rm=TRUE)
#' yMax <- max(
	#' hist(sacPoints$sacDistMid_bio1, breaks=10, plot=FALSE)$counts,
#' 	hist(sacPoints$sacDistMid_bio12, breaks=10, plot=FALSE)$counts
#' )
#' 
#' hist(sacPoints$sacDistMid_bio1, col='red', xlim=c(0, maxDist),
	#' ylim=c(0, yMax), breaks=10, xlab='Distance (m)')
#' hist(sacPoints$sacDistMid_bio12, border='blue', breaks=10, density=10,
	#' add=TRUE)
#' legend('topleft', inset=0.01, legend=c('BIO1', 'BIO12'), fill=c('red', NA),
	#' density=c(NA, 15), border=c('black', 'blue'))
#' 	
#' 
#' ### spatial autocorrelation for spatial points using raster as reference
#' data(lemurs)
#' lemur <- lemurs[lemurs$species == 'Eulemur rufifrons',
	#' c('longitude', 'latitude')]
#' 	
#' # remove record in water
#' lemurBio1 <- raster::extract(madClim[['bio1']], lemur)
#' nas <- which(is.na(lemurBio1))
#' lemur <- lemur[-nas, ]
#' 
#' sacLemur <- localSpatialCorrForValues(x=madClim, focal=lemur, breaks=10)
#' 
#' # plot: code point color by characteristic distance of spatial autocorrelation
#' bio1SacScaled <- 
	#' round(100 * omnibus::stretchMinMax(sacLemur$sacDistMid_bio1))
#' bio12SacScaled <-
	#' round(100 * omnibus::stretchMinMax(sacLemur$sacDistMid_bio12))
#' 
#' grayBio1 <- paste0('gray', bio1SacScaled)
#' grayBio12 <- paste0('gray', bio12SacScaled)
#' 
#' par(mfrow=c(1, 2))
#' leg <- '\n(dark: short dist, light: long dist)'
#' raster::plot(madClim[['bio1']], cex=1.2, main=paste0('BIO 01', leg))
#' points(lemur, pch=21, bg=grayBio1)
#' raster::plot(madClim[['bio12']], cex=1.2, main=paste0('BIO 12', leg))
#' points(lemur, pch=21, bg=grayBio12)
#' 
#' # plot buffers showing characteristic distance of spatial autocorrelation
#' # around each point
#' sacLemurEa <- sp::spTransform(sacLemur, getCRS('mollweide', TRUE))
#' buffsBio1Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
		#' width=sacLemurEa$sacDistMid_bio1)
#' buffsBio12Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
		#' width=sacLemurEa$sacDistMid_bio12)
#' buffsBio1 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' buffsBio12 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' par(mfrow=c(1, 2))
#' raster::plot(madClim[['bio1']], main=paste0('BIO 01', leg))
#' sp::plot(buffsBio1, add=TRUE)
#' points(lemur, pch=21, bg=grayBio1)
#' raster::plot(madClim[['bio12']], main=paste0('BIO 12', leg))
#' sp::plot(buffsBio12, add=TRUE)
#' points(lemur, pch=21, bg=grayBio12)
#' 
#' }
