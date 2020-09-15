
MAXENT AND MAXNET
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
#' occsEnv <- as.data.frame(occsEnv) # need to do this for prediction later
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
#' # prediction rasters
#' mapEnt <- predict(ent, clim, type='logistic')
#' mapNet <- predict(clim, net, type='logistic')
#'
#' par(mfrow=c(1, 2))
#' plot(mapEnt, main='MaxEnt')
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapNet, main='MaxNet')
#' points(occs[ , c('longitude', 'latitude')])
#'
#' # predictions to occurrences
#' (dismo::predict(ent, occsEnv, arguments=c('outputformat=logistic')))
#' (enmSdm::predictMaxEnt(ent, occsEnv, type='logistic'))
#' (c(predict(net, occsEnv, type='logistic')))
#' 
#' # note the differences between the tuning of the two models...
#' # this is because maxnet() (used by trainMaxNet())
#' # uses an approximation:
#' # (note maxnet() calculates hinges and thresholds differently
#' # so we will turn them off)
#' 
#' data(bradypus, package='maxnet')
#' p <- bradypus$presence
#' data <- bradypus[ , 2:3] # easier to inspect betas
#' mn <- maxnet::maxnet(p, data,
#' maxnet::maxnet.formula(p, data, classes='lpq'))
#' mx <- dismo::maxent(data, p,
#' args=c('linear=true', 'product=true', 'quadratic=true', 'hinge=false',
#' 'threshold=false'))
#' 
#' predMx <- dismo::predict(mx, data)
#' predMn <- predict(mn, data, type='logistic')
#' 
#' par(mfrow=c(1, 1))
#' plot(predMx, predMn)
#' abline(0, 1)

BRTs
#' \donttest{
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
#' # settings... defaults probably better, but these are faster
#' lr <- c(0.001, 0.1)
#' tc <- c(1, 3)
#' maxTrees <- 2000
#' set.seed(123)
#' model <- trainBrt(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' 	learningRate = lr,
#' 	treeComplexity = tc,
#' 	maxTrees = maxTrees,
#' 	verbose = TRUE
#' )
#' 
#' plot(model)
#' 
#' # prediction raster
#' nTrees <- model$gbm.call$n.trees
#' map <- predict(clim, model, type='response', n.trees=nTrees)
#' plot(map)
#' points(occs[ , c('longitude', 'latitude')])
#' 
#' }

RFs and CRFs
#' \dontest{
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
#' set.seed(123)
#'
#' # random forest
#' rf <- trainRf(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' )
#' 
#' # conditional random forest
#' crf <- trainCrf(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' )
#' 
#' plot(rf)
#' 
#' # prediction rasters
#' mapRf1 <- predict(clim, rf, type='prob') # opposite class!
#' mapRf2 <- 1 - predict(clim, rf, type='prob') # correct
#' pointsFx <- function() points(occs[ , c('longitude', 'latitude')])
#' plot(stack(mapRf1, mapRf2), addfun=pointsFx)
#'
#' # CRFs are tricky...
#' }

LARS
#' \dontest{
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
#' al <- c(0.01, 0.5, 1)
#' fit1 <- trainLars(data=data, penalty='cMCP', family='binomial',
#'    nfolds=3, alphas=al, quadratic=FALSE, cubic=FALSE, interaction=FALSE,
#'    interQuad=FALSE, verbose=TRUE)
#' fit2 <- trainLars(data=data, penalty='cMCP', family='binomial',
#'    nfolds=3, alphas=al, quadratic=TRUE, cubic=FALSE, interaction=FALSE,
#'    interQuad=FALSE, verbose=TRUE)
#' fit3 <- trainLars(data=data, penalty='cMCP', family='binomial',
#'    nfolds=3, alphas=al, quadratic=TRUE, cubic=TRUE, interaction=TRUE,
#'    interQuad=TRUE, verbose=TRUE)
#' 
#' summary(fit1)
#' summary(fit2)
#' summary(fit3)
#'
#' # predictions using all variables
#' pred1 <- predictLars(fit1, data, type='response')
#' pred2 <- predictLars(fit2, data, type='response')
#' pred3 <- predictLars(fit3, data, type='response')
#'
#' # partial predictions examining effect of just x1 (plus any interactions)
#' pred1bio1 <- predictLars(fit1, data, type='response', preds='bio1')
#' pred2bio1 <- predictLars(fit2, data, type='response', preds='bio1')
#' pred3bio1 <- predictLars(fit3, data, type='response', preds='bio1')
#'
#' par(mfrow=c(3, 3))
#' xlim <- c(0, 1)
#' breaks <- seq(0, 1, by=0.1)
#' plot(data$bio1, pred1bio1, ylim=c(0, 1))
#' points(data$bio1, pred2bio1, col='blue')
#' points(data$bio1, pred3bio1, col='red')
#' legend('topright', pch=1, col=c('black', 'blue', 'red'),
#' legend=c('linear-only', 'linear + quadratic', 'all terms'))
#'
#' # predictions using just bio1 and bio12
#' pred3bio1_12 <- predictLars(fit3, data, type='response', preds=c('bio1', 'bio12'))
#' plot(pred3, pred3bio1_12)
#' abline(0, 1)
#' }


#' # prediction rasters
#' mapGlm <- predict(clim, gl, type='response')
#' mapGam <- predict(clim, ga, type='response')
#' mapNs <- predict(clim, ga, type='response')
#'
#' par(mfrow=c(1, 3))
#' plot(mapGlm, main='GLM')
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapGam, main='GAM')
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapNs, main='NS')
#' points(occs[ , c('longitude', 'latitude')])
#' }

GLMs, GAMs, NSs
#' library(brglm2)
#'
#' \dontest{
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
#' # GLM
#' gl <- trainGlm(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # GAM
#' ga <- trainGam(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # NS
#' ns <- trainNs(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # prediction rasters
#' mapGlm <- predict(clim, gl, type='response')
#' mapGam <- predict(clim, ga, type='response')
#' mapNs <- predict(clim, ga, type='response')
#'
#' par(mfrow=c(1, 3))
#' plot(mapGlm, main='GLM')
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapGam, main='GAM')
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapNs, main='NS')
#' points(occs[ , c('longitude', 'latitude')])
#' }

