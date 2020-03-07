#' ### model red-bellied lemurs
#' data(mad0)
#' data(lemurs)
#' 
#' # climate data
#' bios <- c(1, 5, 12, 15)
#' clim <- raster::getData('worldclim', var='bio', res=10)
#' clim <- raster::subset(clim, bios)
#' clim <- raster::crop(clim, mad0)
#' 
#' # occurrence data
#' occs <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' occsEnv <- raster::extract(clim, occs[ , c('longitude', 'latitude')])
#' 
#' # background sites
#' bg <- 2000 # too few cells to locate 10000 background points
#' bgSites <- dismo::randomPoints(clim, 2000)
#' bgEnv <- extract(clim, bgSites)
#' 
#' # collate
#' presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
#' env <- rbind(occsEnv, bgEnv)
#' env <- cbind(presBg, env)
#' env <- as.data.frame(env)
#' 
#' preds <- paste0('bio', bios)
#' 
#' regMult <- 1:3 # default values are probably better, but these will be faster
#' 
#' # calibrate MaxEnt model
#' ent <- trainMaxEnt(
#' 	data=env,
#' 	resp='presBg',
#' 	preds=preds,
#' 	regMult=regMult,
#' 	classes='lpq',
#' 	verbose=TRUE
#' )
#' 
#' # calibrate MaxNet model
#' net <- trainMaxNet(
#' 	data=env,
#' 	resp='presBg',
#' 	preds=preds,
#' 	regMult=regMult,
#' 	classes='lpq',
#' 	verbose=TRUE
#' )
#' 
#' # note the differences between the two models...
#' # this is because maxnet() (used by trainMaxNet())
#' # uses an approximation:
#' # (note maxnet() calculates hinges and thresholds differently
#' # so we will turn them off)
#' 
#' data(bradypus, package='dismo')
#' p <- bradypus$presence
#' data <- bradypus[ , 2:3] # easier to inspect betas
#' mn <- maxnet::maxnet(p, data, maxnet.formula(p, data, classes='lpq'))
#' mx <- dismo::maxent(data, p,
#' args=c('linear=true', 'product=true', 'quadratic=true', 'hinge=false',
#' 'threshold=false'))
#' 
#' predMx <- dismo::predict(mx, data)
#' predMn <- predict(mn, data, type='logistic')
#' 
#' plot(predMx, predMn)
#' abline(0, 1)

