#' Generic predict function for LMs, GLMs, GAMs, RFs, BRTs, CRFs, Maxent, and more
#'
#' This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, Maxent, and more.
#' @param model  Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{MaxNet}, and possibly others (worth a try!).
#' @param newdata Data frame or matrix with data to which to predict
#' @param maxentFun Either \code{'dismo'} in which case a Maxent model is predicted using the default \code{predict} function from the \pkg{dismo} package, or \code{'enmSdm'} in which case the function \code{\link[enmSdm]{predictMaxEnt}} function from the \pkg{enmSdm} package is used.
#' @param ... Arguments to pass to the algorithm-specific \code{predict} function.
#' @return Numeric.
#' @seealso \code{\link[stats]{predict}} from the stats package, \code{\link[dismo]{predict}} from the dismo package, \code{\link[raster]{predict}} from the raster package
#' @export

predictEnmSdm <- function(
	model,
	newdata,
	maxentFun='dismo',
	...
) {

	dots <- list(...)

	dataType <- if (inherits(newdata, c('matrix', 'data.frame'))) {
		'table'
	} else if (inherits(newdata, c('RasterLayer', 'RasterStack', 'RasterBrick'))) {
		'rasterRaster'
	} else if (inherits(newdata, c('SpatRaster'))) {
		'terraRaster'
	}

	# GAM
	if (inherits(model, 'gam')) {

		out <- mgcv::predict.gam(model, newdata, type='response', ...)

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

	out

}
