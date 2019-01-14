#' Calibrate a random forest model
#'
#' This function trains a random forest model. It first finds the optimal value for \code{mtry} (number of variables sampled as candidates at each split). It then calls the function  \code{\link[randomForest]{randomForest}} from the \pkg{randomForest} package.
#' @param data Data frame.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Character. If "\code{binomial}" then the response is converted to a binary factor with levels 0 and 1. Otherwise, this argument has no effect.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) \emph{or} a numeric list of weights, one per class in \code{resp}. The default is to assign equal total weight to presences and contrast sites (\code{TRUE}).
#' @param verbose Logical. If \code{TRUE} then display progress for finding optimal value of \code{mtry}.
#' @param ... Arguments to pass to \code{\link[randomForest]{randomForest}}.
#' @return Object of class \code{\link[randomForest]{randomForest}}.
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{trainCrf}}
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(n = 6*100), ncol = 6)
#' # true variables will be #1, #2, #5, and #6, plus
#' # the squares of #1 and #6, plus
#' # interaction between #1 and #6
#' # the cube of #5
#' imp <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x1_pow2', 'x6_pow2', 'x1_by_x6', 'x5_pow3')
#' betas <- c(5, 2, 0, 0, 1, -1, 8, 1, 2, -4)
#' names(betas) <- imp
#' y <- 0.5 + x %*% betas[1:6] + betas[7] * x[ , 1] +
#' betas[8] * x[ , 6] + betas[9] * x[ , 1] * x[ , 6] + betas[10] * x[ , 5]^3
#' y <- as.integer(y > 10)
#' x <- cbind(y, x)
#' x <- as.data.frame(x)
#' names(x) <- c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' x$y <- as.factor(x$y)
#' model <- trainRf(x)
#' @export

trainRf <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	w = TRUE,
	verbose = FALSE,
	...
) {

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# weights
	classwt <- if (family == 'binomial') {
		if ('logical' %in% class(w)) {
			if (w) {
				c(1, 1)
			} else if (!w) {
				NULL
			}
		} else {
			w
		}
	}
		
	
	# binomial response
	if (family == 'binomial') data[ , resp] <- factor(data[ , resp], levels=0:1)

	model <- randomForest::tuneRF(x=data[ , preds, drop=FALSE], y=data[ , resp], trace=FALSE, plot=FALSE, doBest=TRUE, strata=data[ , resp], classwt=classwt, ...)
	model

}

