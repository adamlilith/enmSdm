source('E:/Ecology/Drive/R/enmSdm/R/trainMaxEnt.r')

# ### model red-bellied lemurs
# data(mad0, package='enmSdm')
# data(lemurs, package='enmSdm')

# # climate data
# bios <- c(1, 5, 12, 15)
# clim <- raster::getData('worldclim', var='bio', res=10)
# clim <- raster::subset(clim, bios)
# clim <- raster::crop(clim, mad0)

# # occurrence data
# occs <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
# occsEnv <- raster::extract(clim, occs[ , c('longitude', 'latitude')])
# occsEnv <- as.data.frame(occsEnv) # need to do this for prediction later

# # background sites
# bg <- 2000 # too few cells to locate 10000 background points
# bgSites <- dismo::randomPoints(clim, 2000)
# bgEnv <- raster::extract(clim, bgSites)

# # collate
# presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
# env <- rbind(occsEnv, bgEnv)
# env <- cbind(presBg, env)
# env <- as.data.frame(env)

# preds <- paste0('bio', bios)

# regMult <- 1:3 # defa

ent <- trainMaxEnt(
	data=env,
	resp='presBg',
	preds=preds,
	regMult=regMult,
	classes='lpq',
	verbose=TRUE, 
	cores = 2
)
