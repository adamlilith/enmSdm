#' ### setup for all examples
#' ##########################
#' 
#' # simulate species starting at center of map
#' # will be using a spatially Gaussian distribution
#' # that uses only longitude and latitude as predictors
#' gauss <- function(x1, x2, mu1=0, mu2=0, sigma1=1, sigma2=0, rho=0) {
#' 
#'		first <- ((x1 - mu1) / sigma1)^2
#'		prod <- ((2 * rho * (x1 - mu1) * (x2 - mu2)) / (sigma1 * sigma2))
#'		second <- ((x2 - mu2) / sigma2)^2
#'		denom <- 2 * (1 - rho^2)
#' 
#'		inside <- first - prod + second
#'		inside <- (-1 * inside) / denom
#' 
#'		expo <- exp(inside)
#'		expo
#' 
#' }
#' 
#' r <- raster()
#' r[] <- 1
#' mollweide <- enmSdm::getCRS('mollweide', TRUE)
#' r <- raster::projectRaster(r, crs=mollweide)
#' ll <- enmSdm::longLatRasters(r, m=TRUE)
#' 
#' # simulated population trajectory:
#' # time 0 to 1: no change
#' # time 1 to 2: size doubles, no shift
#' # time 2 to 3: size halves, no shift
#' # time 3 to 4: size same, moves north
#' # time 4 to 5: size same, moves south
#' # time 5 to 6: size same, moves east
#' # time 6 to 7: size same, moves west
#' # times 7 to 8: size same, moves northwest
#' # times 8 to 10: size same, moves northwest, double time
#' # times 10 to 11: size doubles, moves northwest
#' # times 11 to 12: same size, moves northwest
#' x1 <- ll[['longitude']]
#' x2 <- ll[['latitude']]
#' 
#' s1 <- 2000000
#' s2 <- 1000000
#' shift <- 1500000
#' 
#' species_t0 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t1 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t2 <- gauss(x1=x1, x2=x2, sigma1=2 * s1, sigma2=2 * s2)
#' species_t3 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t4 <- gauss(x1=x1, x2=x2, mu2=shift, sigma1=s1, sigma2=s2)
#' species_t5 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t6 <- gauss(x1=x1, x2=x2, mu1=shift, sigma1=s1, sigma2=s2)
#' species_t7 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' species_t8 <- gauss(x1=x1, x2=x2, mu1=-1 * shift, mu2=1 * shift, sigma1=s1, sigma2=s2)
#' species_t10 <- gauss(x1=x1, x2=x2, mu1=-2 * shift, mu2=2 * shift, sigma1=s1, sigma2=s2)
#' species_t11 <- gauss(x1=x1, x2=x2, mu1=-3 * shift, mu2=3 * shift, sigma1=2 * s1, sigma2=s2)
#' species_t12 <- gauss(x1=x1, x2=x2, mu1=-4 * shift, mu2=4 * shift, sigma1=2 * s1, sigma2=s2,
#' 	rho=-0.5)
#' 	
#' rasts <- raster::stack(species_t0, species_t1, species_t2, species_t3,
#' species_t4, species_t5, species_t6, species_t7, species_t8,
#' 	species_t10, species_t11, species_t12)
#' times <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12)
#' names(rasts) <- paste0('t', times)
#' plot(rasts)
#' 
#' ### example with stationary land mass
#' #####################################
#' 
#' # across just start and end time periods
#' (vels <- bioticVelocity(rasts, times=times, atTimes=c(0, 12),
#' metrics='centroid'))
#' 
#' # across each time interval
#' (vels <- bioticVelocity(rasts, times=times, metrics='centroid'))
#' 
#' # plot movement of centroid
#' weightedLong <- species_t0 * ll[['longitude']]
#' weightedLat <- species_t0 * ll[['latitude']]
#' startLong <- cellStats(weightedLong, 'sum') / cellStats(species_t0, 'sum')
#' startLat <- cellStats(weightedLat, 'sum') / cellStats(species_t0, 'sum')
#' 
#' plot(species_t0)
#' 
#' for (i in 1:(nrow(vels) - 1)) {
#' 
#' x1pos <- vels$centroidLong[i]
#' y1pos <- vels$centroidLat[i]
#' x2pos <- vels$centroidLong[i + 1]
#' y2pos <- vels$centroidLat[i + 1]
#' 
#' 	move <- sqrt((x1pos - x2pos)^2 + (y1pos - y2pos)^2)
#' if (move > 10) arrows(x1pos, y1pos, x2pos, y2pos, angle=15, length=0.05)
#' 
#' }
#' 
#' # all metrics
#' (vels <- bioticVelocity(x=rasts, times=times))
#' 
#' ### example with shifting land mass
#' ###################################
#' 
#' species_t0 <- gauss(x1=x1, x2=x2, sigma1=s1, sigma2=s2)
#' rasts <- stack(species_t0, species_t0, species_t0, species_t0, species_t0)
#' # simulated population trajectory:
#' # time 0 to 1: species no change, lose eastward land
#' # time 1 to 2: species no change, lose southward land
#' # time 2 to 3: species no change, gain southward land
#' # time 3 to 4: species no change, reset land to t0
#' 
#' # create masks
#' ext <- extent(species_t0)
#' 
#' long <- ext@xmin + 0.55 * (ext@xmax - ext@xmin)
#' longMask <- extent(long, ext@xmax, ext@ymin, ext@ymax)
#' longMask <- as(longMask, 'SpatialPolygons')
#' projection(longMask) <- projection(species_t0)
#' 
#' lat <- ext@ymin + 0.45 * (ext@ymax - ext@ymin)
#' latMask <- extent(ext@xmin, ext@xmax, lat, ext@ymax)
#' latMask <- as(latMask, 'SpatialPolygons')
#' projection(latMask) <- projection(species_t0)
#' 
#' species_t1s <- mask(species_t0, longMask, inverse=TRUE) # lose eastward
#' species_t2s <- mask(species_t1s, latMask1) # lose southward
#' species_t3s <- species_t1s # gain southward
#' species_t4s <- species_t0s # gain eastward
#' 
#' rasts <- stack(species_t0s, species_t1s, species_t2s, species_t3s,
#' species_t4s)
#' 
#' plot(rasts)
#' 
#' # is this what you want? could be considered an artifact, could
#' # be considered "true" biotic velocity... range shift is extrinsic
#' (vels <- bioticVelocity(x=rasts, times=0:4, metrics=c('centroid', 'mean')))
#' 
#' # calculate velocity only for cells shared by start and end rasters
#' # in each time step
#' (vels <- bioticVelocity(x=rasts, times=0:4, metrics=c('centroid', 'mean'),
#' onlyInSharedCells=TRUE))
