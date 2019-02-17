#' Number of response data in a model object
#'
#' This function returns the number of response data used in a model (i.e., the sample size). If the data are binary it can return the number of 1s and 0s.
#' @param x A model object. This can be of many classes, including "gbm", "glm", "gam", "MaxEnt", and so on.
#' @param binary Logical, if \code{TRUE} (default) then the number of 1s and 0s in the response data is returned. If \code{FALSE} then the returned values is the total number of response data.
#' @param graceful Logical, if \code{TRUE} (default) then if the function cannot determine the sample size from the model object \code{NA} if returned. If \code{FALSE}, then the function exits with an error.
#' @return One or two integers named integers.
#' @examples
#' set.seed(123)
#' y <- runif(1:101)^2
#' yBinary <- as.integer(y > 0.6)
#' x <- data.frame(x1=1:101, x2=rnorm(101))
#' model <- lm(y ~ x1 + x2, data=x)
#' modelBinary <- glm(yBinary ~ x1 + x2, data=x, family='binomial')
#' modelSize(model, FALSE)
#' modelSize(model, TRUE) # not binary input... notice warning
#' modelSize(modelBinary)
#' modelSize(modelBinary, FALSE)
#' @export

modelSize <- function(
	x,
	binary = TRUE,
	graceful = TRUE
) {

	modelClass <- class(x)

	# LM/GAM/GLM
	samples <- if ('gam' %in% modelClass | 'glm' %in% modelClass | 'lm' %in% modelClass) {
	
		x$model[[1]]
		
	# BIOCLIM
	} else if ('Bioclim' %in% modelClass) {
	
		c(rep(1, nrow(x@presence)), rep(0, nrow(x@absence)))

	# BRT/GBM
	} else if ('gbm' %in% modelClass) {
	
		x$data$y
	
	# LARS
	} else if ('cv.grpregOverlap' %in% modelClass | 'cv.grpreg' %in% modelClass | 'larsModel' %in% modelClass) {
	
		x$lars$y
		
	# Maxent
	} else if ('MaxEnt' %in% modelClass) {
	
		c(rep(1, nrow(x@presence)), rep(0, nrow(x@absence)))
	
	# MaxNet
	} else if ('maxnet' %in% modelClass | 'lognet' %in% modelClass | 'glmnet' %in% modelClass) {
	
		rep(0.5, x$nobs)
		warning('Cannot determine number of presences and background sites for a MaxNet model.')
		
	# random forest in party package
	} else if ('randomForest' %in% modelClass) {
	
		if (binary) { (as.integer(x$y) == 2) } else { x$y }
		
	# random conditional forest in party package
	} else if ('RandomForest' %in% modelClass) {
	
		if (binary) { (as.integer(x@responses@variables$y) == 2) } else { x@responses@variables$y }
		
	# anything else!
	} else {
	
		NA
		
	}
	
	# sample size
	if (all(is.na(samples))) {
		if (graceful) {
			out <- NA
		} else {
			stop('Cannot extract sample size from model object.')
		}
	} else if (binary) {
		out <- c(sum(samples == 1), sum(samples == 0))
		names(out) <- c('num1s', 'num0s')
		if (out[1] == 0 & out[2] == 0) warning('Model does not seem to be using a binary response.', .immediate=TRUE)
	} else {
		out <- length(samples)
		names(out) <- 'sampleSize'
	}
	
	out

}
