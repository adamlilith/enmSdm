#' Generic predict function for LMs, GLMs, GAMs, RFs, BRTs, CRFs, Maxent, and more
#'
#' This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, Maxent, and more.
#' @param model  Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{MaxNet}, and possibly others (worth a try!).
#' @param newdata Data frame or matrix with data to which to predict
#' @param maxentFun Either \code{'dismo'} in which case a Maxent model is predicted using the default \code{predict} function from the \pkg{dismo} package, or \code{'enmSdm'} in which case the function \code{\link[enmSdm]{predictMaxEnt}} function from the \pkg{enmSdm} package is used.
#' @param ... Arguments to pass to the algorithm-specific \code{predict} function.
#' @return Numeric.
#' @seealso \code{\link[stats]{predict}} from the stats package, \code{\link[dismo]{predict}} from the dismo package, \code{\link[raster]{predict}} fromm the raster package
#' @examples
#' 
#' @export

predictEnmSdm <- function(
	model,
	newdata,
	maxentFun='dismo',
	...
) {

	cl <- class(model)

	# GAM/GLM
	if ('gam' %in% cl | 'glm' %in% cl) {
	
		out <- predict(model, newdata, type='response', ...)
		
	# LM
	} else if ('lm' %in% cl) {
	
		out <- predict(model, newdata, ...)
		
	# BRT
	} else if ('gbm' %in% cl) {
	
		out <- predict(model, newdata, n.trees=model$gbm.call$n.trees, ...)
	
	# Maxent
	} else if ('MaxEnt' %in% cl) {
	
		out <- if (maxentFun == 'dismo') {
			dismo::predict(model, newdata, ...)
		} else if (maxentFun == 'enmSdm') {
			predictMaxEnt(model, newdata, ...)
		}
	
	# MaxNet
	} else if ('maxnet' %in% cl) {
	
		out <- predict(model, newdata, ...)
		
	# random forest in party package
	} else if ('RandomForest' %in% cl) {
	
		out <- predict(model, newdata, type='prob', ...)
		out <- unlist(out)
		
	# anything else!
	} else {
	
		out <- predict(model, newdata, type='response', ...)
		
	}
	
	out

}
