#' Generic predict function for LMs, GLMs, GAMs, RFs, BRTs, CRFs, Maxent, and more
#'
#' This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, Maxent, and more.
#' @param model  Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{MaxNet}, \code{prcomp}, \code{kde}, \code{gbm}, and possibly others (worth a try!).
#' @param newdata Data frame or matrix with data to which to predict
#' @param maxentFun This argument is only used if the \code{model} object is a MaxEnt model; otherwise, it is ignored. I takes a value of either \code{'dismo'}, in which case a MaxEnt model is predicted using the default \code{predict} function from the \pkg{dismo} package, or \code{'enmSdm'} in which case the function \code{\link[enmSdm]{predictMaxEnt}} function from the \pkg{enmSdm} package (this package) is used.
#' @param cores Number of cores to use. The default is 1. If >1 and \code{newdata} is a raster, only 1 core is used (i.e., basically, \code{cores} is ignored if you're writing to a raster... sorry!)
#' @param nrows Number of rows of \code{newdata} to predict at a time (assuming \code{newdata} is a \code{data.frame} or \code{matrix}). The default value is to predict all rows at once, but for very large data frames/matrices this can lead to memory issues in some cases. By setting the number of rows, \code{newdata} can be divided into chunks, and predictions made to each chunk, which may ease memory limitations. This can be combined with multi-coring (which will increase memory requirements). In this case, all cores combined will get \code{nrows} of data. How many rows are too many? You will have to decide depending on your data and the output. For example, predicting the outcome of a GLM on data with 10E6 rows ma be fine, but predicting a PCA (with multiple axes) to the data data may require too much memory. You can use \code{\link[omnibus]{memUse}} to see to help figure this out.
#' @param ... Arguments to pass to the algorithm-specific \code{predict} function.
#' @return Numeric.
#' @seealso \code{\link[stats]{predict}} from the \pkg{stats} package, \code{\link[dismo]{predict}} from the \pkg{dismo} package, \code{\link[raster]{predict}} from the \pkg{raster} package, \code{\link[terra]{predict}} from the \pkg{terra} package
#' @export

predictEnmSdm <- function(
	model,
	newdata,
	maxentFun = 'dismo',
	cores = 1,
	nrows = nrow(newdata),
	...
) {

	dataType <- if (inherits(newdata, c('matrix', 'data.frame', 'data.table'))) {
		'table'
	} else if (inherits(newdata, c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
		'rasterRaster'
	} else if (inherits(newdata, c('SpatRaster'))) {
		'terraRaster'
	}
	
	cores <- if (dataType %in% c('rasterRaster', 'terraRaster')) {
		1L
	} else {
		# min(cores, parallel::detectCores(logical=FALSE))
		min(cores, parallelly::availableCores(logical=FALSE))
	}
	
	### predict to each subset of rows
	##################################
	n <- nrow(newdata)
	if (dataType == 'table' & nrows < n) {

		# predict to each set of rows
		rowJobs <- ceiling(n / nrows)
		for (countRowJob in 1:rowJobs) {

			thisNewData <- newdata[(1 + (countRowJob - 1) * nrows):min(n, countRowJob * nrows), , drop=FALSE]
			thisOut <- predictEnmSdm(model = model, newdata = thisNewData, maxentFun = maxentFun, cores = cores, nrows = nrow(thisNewData), ...)
			
			if (exists('out', inherits = FALSE)) {
				if (inherits(out, c('data.frame', 'matrix'))) {
					out <- rbind(out, thisOut)
				} else {
					out <- c(out, thisOut)
				}
			} else {
				out <- thisOut
			}
		
		}
	
	### multi-core
	##############
	} else if (cores > 1L & nrows > n) {

		`%makeWork%` <- foreach::`%dopar%`
		options(future.globals.maxSize = +Inf)
		doFuture::registerDoFuture()
		future::plan(future::multisession, workers = cores)
		# cl <- parallel::makePSOCKcluster(cores, setup_strategy = 'sequential')
		# doParallel::registerDoParallel(cl)
		# parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths()) # can find non-standard paths

		mcOptions <- list(preschedule=TRUE, set.seed=FALSE, silent=FALSE)

		nPerJob <- floor(nrow(newdata) / cores)
		jobs <- rep(1:cores, each=nPerJob)
		if (length(jobs) < nrow(newdata)) jobs <- c(jobs, rep(cores, nrow(newdata) - length(jobs)))

		combine <- if (inherits(model, 'prcomp')) {
			'rbind'
		} else {
			'c'
		}
		
		# out <- foreach::foreach(i=1L:cores, .options.multicore=mcOptions, .combine=combine, .inorder=TRUE, .export=c('predictEnmSdm', 'predictMaxEnt', 'predictMaxNet'), .packages = c('enmSdm')) %makeWork%
		out <- foreach::foreach(i=1L:cores, .options.multicore=mcOptions, .combine=combine, .inorder=TRUE, .export=c('predictEnmSdm', 'predictMaxEnt', 'predictMaxNet'), .packages='enmSdm') %makeWork%
			predictEnmSdm(
				i = i,
				model = model,
				newdata = newdata[jobs == i, , drop = FALSE],
				maxentFun = maxentFun,
				cores = 1L,
				...
			)

		# if (cores > 1L) parallel::stopCluster(cl)

	### single-core
	###############
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
			
		# KDE from ks package
		} else if (inherits(model, 'kde')) {

			requireNamespace('ks', quietly = TRUE)
			out <- predict(model, x=as.matrix(newdata), ...)

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
