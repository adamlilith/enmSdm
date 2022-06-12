#' Generic predict function for LMs, GLMs, GAMs, RFs, BRTs, CRFs, Maxent, and more
#'
#' This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, Maxent, and more.
#' @param model  Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{MaxNet}, \code{prcomp}, \code{kde}, \code{gbm}, and possibly others (worth a try!).
#' @param newdata Data frame or matrix with data to which to predict
#' @param maxentFun This argument is only used if the \code{model} object is a MaxEnt model; otherwise, it is ignored. I takes a value of either \code{'dismo'}, in which case a MaxEnt model is predicted using the default \code{predict} function from the \pkg{dismo} package, or \code{'enmSdm'} in which case the function \code{\link[enmSdm]{predictMaxEnt}} function from the \pkg{enmSdm} package (this package) is used.
#' @param cores Number of cores to use. The default is 1. If >1 and \code{newdata} is a raster, only 1 core is used (i.e., basically, \code{cores} is ignored if you're writing to a raster... sorry!)
#' @param ... Arguments to pass to the algorithm-specific \code{predict} function.
#' @return Numeric.
#' @seealso \code{\link[stats]{predict}} from the \pkg{stats} package, \code{\link[dismo]{predict}} from the \pkg{dismo} package, \code{\link[raster]{predict}} from the \pkg{raster} package, \code{\link[terra]{predict}} from the \pkg{terra} package
#' @export

predictEnmSdm <- function(
	model,
	newdata,
	maxentFun='dismo',
	cores = 1,
	...
) {

	dataType <- if (inherits(newdata, c('matrix', 'data.frame', 'data.table'))) {
		'table'
	} else if (inherits(newdata, c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
		'rasterRaster'
	} else if (inherits(newdata, c('SpatRaster'))) {
		'terraRaster'
	}
	
	if (dataType %in% c('rasterRaster', 'terraRaster')) {
		cores <- 1L
	} else {
		cores <- min(cores, parallel::detectCores(logical=FALSE))
	}
	
	# multi-core
	if (cores > 1L) {

		`%makeWork%` <- foreach::`%dopar%`
		cl <- parallel::makePSOCKcluster(cores)
		doParallel::registerDoParallel(cl)

		paths <- .libPaths() # need to pass this to avoid "object '.doSnowGlobals' not found" error!!!
		mcOptions <- list(preschedule=TRUE, set.seed=FALSE, silent=FALSE)

		nPerJob <- floor(nrow(newdata) / cores)
		jobs <- rep(1L:cores, each=nPerJob)
		if (length(jobs) < nrow(newdata)) jobs <- c(jobs, rep(cores, nrow(newdata) - length(jobs)))

		combine <- if (inherits(model, 'prcomp')) {
			'rbind'
		} else {
			'c'
		}
		
		out <- foreach::foreach(i=1L:cores, .options.multicore=mcOptions, .combine=combine, .inorder=TRUE, .export=c('predictEnmSdm'),
		.packages = c('enmSdm')) %makeWork%
			predictEnmSdm(
				i = i,
				model = model,
				newdata = newdata,
				maxentFun = maxentFun,
				cores = 1L,
				...
			)

		if (cores > 1L) parallel::stopCluster(cl)

	# single-core
	} else {
		
		if (inherits(newdata, 'data.table')) newdata <- as.data.frame(newdata)

		# GAM
		if (inherits(model, 'gam')) {

			out <- if (dataType == 'table') {
				mgcv::predict.gam(model, newdata, type='response', ...)
			} else if (dataType == 'rasterRaster') {
				raster::predict(newdata, model, type='response', ...)
			} else if (dataType == 'terraRaster') {
				terra::predict(newdata, model, type='response', ...)
			}

		# GLM
		} else if (inherits(model, 'glm')) {

			out <- if (dataType == 'table') {
				stats::predict.glm(model, newdata, type='response', ...)
			} else if (dataType == 'rasterRaster') {
				raster::predict(newdata, model, type='response', ...)
			} else if (dataType == 'terraRaster') {
				terra::predict(newdata, model, type='response', ...)
			}

		# LM
		} else if (inherits(model, 'lm')) {

			out <- stats::predict.lm(model, newdata, ...)

		# BRT
		} else if (inherits(model, 'gbm')) {

			out <- gbm::predict.gbm(model, newdata, n.trees=model$gbm.call$n.trees, type='response', ...)

		# Maxent
		} else if (inherits(model, 'MaxEnt')) {

			out <- if (maxentFun == 'dismo') {
				dismo::predict(model, newdata, ...)
			} else if (maxentFun == 'enmSdm') {
				predictMaxEnt(model, newdata, ...)
			}

		# MaxNet
		} else if (inherits(model, 'maxnet')) {

			out <- predictMaxNet(object=model, newdata=newdata, ...)

		# random forest in party package
		} else if (inherits(model, 'RandomForest')) {

			out <- randomForest::predict.randomForest(model, newdata, type='prob', ...)
			out <- unlist(out)

		# anything else!
		} else {

			out <- do.call('predict', args=list(object=model, newdata=newdata, type='response', ...))

		}
		
	} # single-core

	out

}
