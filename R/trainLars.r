#' Calibrate a least angle regression (LARS) model
#'
#' This function calculates the least angle regression (LARS) using possibly overlapping grouped covariates. The model is fit using cross validation (the \code{cv.grpregOverlap} function). The cross-validation is calculated across values of the \code{alpha}, which controls the degree of ridge penalty (\code{alpha ~0} (bit not = 0) imposes the full ridge penalty and \code{alpha) = 1} imposes no ridge penalty). Higher-order terms are constructed (e.g., quadratic, 2-way interaction, etc.) and fitted in a manner that respects mrginality (i.e., all lower order terms will have non-zero coefficients if a high-order term is used).
#' @param data Data frame.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param alphas Numeric or numeric vector in the range (0, 1]. Degree of ridge penalty to impose (values close to 0 ==> full ridge penalty, while a value of 1 imposes no rifhe penalty).
#' @param scale Logical. If TRUE then scale values in \code{data[ , preds]} are rescaled to have mean of 0 and standard deviation of 1.
#' @param quadratic Logical. If TRUE then include quadratic terms in model construction stage for non-factor predictors. Quadratic columns will be named \code{<predictor name>_pow2}.
#' @param cubic Logical. If TRUE then include cubic terms in model construction stage for non-factor predictors. Cubic columns will be named \code{<predictor name>_pow3}.
#' @param interaction Logical. If TRUE then include 2-way interaction terms (including interactions between factor predictors). Interaction columns will be named \code{<predictor 1 name>_by_<predictor 2 name>}.
#' @param interQuad Logical. If TRUE then include all possible interactions of the form \code{x * y^2} unless \code{y} is a factor (linear-by-quadratic features). Linear-by-quadratic columns will be named \code{<predictor 1 name>_by_<predictor 2 name>_pow2}.
#' @param na.rm Logical. If TRUE then remove all rows of \code{data} in which there is at least one \code{NA} among \code{resp} or \code{preds}. The default is FALSE, which will cause an error if any row has an \code{NA}.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param ... Arguments to pass to \code{grpreg} \code{grpregOverlap}, and \code{cv.grpregOverlap}, especially \code{family} and \code{penalty}. \emph{Do not} include the \code{'group'} argument or \code{alpha} arguments.
#' @return Object of class \code{grpreg} and \code{grpregOverlap}.
#' @details If \code{scale} is \code{TRUE} then predictors with zero variance will be removed from the data before the model is trained.
#' @seealso \code{\link{makeLarsData}}, \code{\link{predictLars}}, \code{\link[grpreg]{grpreg}}, \code{\link[grpregOverlap]{grpregOverlap}}, \code{\link[grpregOverlap]{cv.grpregOverlap}}
#' @examples
#' \donttest{
#' set.seed(123)
#' X <- matrix(rnorm(n = 6*100), ncol = 6)
#' # true variables will be #1, #2, #5, and #6, plus
#' # the squares of #1 and #6, plus
#' # interaction between #1 and #6
#' # the cube of #5
#' imp <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x1_pow2', 'x6_pow2', 'x1_by_x6', 'x5_pow3')
#' betas <- c(5, 2, 0, 0, 1, -1, 8, 1, 2, -4)
#' names(betas) <- imp
#' y <- 0.5 + X %*% betas[1:6] + betas[7] * X[ , 1]^2 +
#'     betas[8] * X[ , 6]^2 + betas[9] * X[ , 1] * X[ , 6] + betas[10] * X[ , 5]^3
#' y <- as.integer(y >= 0)
#' X <- cbind(y, X)
#' X <- as.data.frame(X)
#' preds <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' names(X) <- c('y', preds)
#'
#' # NOTE: default is to use nfolds=10 and a finer sequence of alphas
#' fit1 <- trainLars(data=X, penalty='cMCP', family='binomial', nfolds=3, alphas=c(0.01, 0.5, 1),
#'    quadratic=FALSE, cubic=FALSE, interaction=FALSE, interQuad=FALSE)
#' fit2 <- trainLars(data=X, penalty='cMCP', family='binomial', nfolds=3, alphas=c(0.01, 0.5, 1),
#'     quadratic=TRUE, cubic=FALSE, interaction=FALSE, interQuad=FALSE)
#' fit3 <- trainLars(data=X, penalty='cMCP', family='binomial', nfolds=3, alphas=c(0.01, 0.5, 1),
#'     quadratic=TRUE, cubic=TRUE, interaction=TRUE, interQuad=TRUE)
#'
#' summary(fit1)
#' summary(fit2)
#' summary(fit3)
#'
#' XX <- matrix(rnorm(n = 6*100), ncol = 6)
#' XX <- as.data.frame(XX)
#' names(XX) <- preds
#'
#' # predictions using all variables
#' pred1 <- predictLars(fit1, XX, type='response')
#' pred2 <- predictLars(fit2, XX, type='response')
#' pred3 <- predictLars(fit3, XX, type='response')
#'
#' # partial predictions examining effect of just x1 (plus any interactions)
#' pred1x1 <- predictLars(fit1, XX, type='response', preds='x1')
#' pred2x1 <- predictLars(fit2, XX, type='response', preds='x1')
#' pred3x1 <- predictLars(fit3, XX, type='response', preds='x1')
#'
#' # partial predictions examining effect of just x4 (plus any interactions--
#' # note in reality x4 has no effect on y)
#' pred1x4 <- predictLars(fit1, XX, type='response', preds='x4')
#' pred2x4 <- predictLars(fit2, XX, type='response', preds='x4')
#' pred3x4 <- predictLars(fit3, XX, type='response', preds='x4')
#'
#' par(mfrow=c(3, 3))
#' xlim <- c(0, 1)
#' breaks <- seq(0, 1, by=0.1)
#' hist(pred1, main='linear only', xlim=xlim, breaks=breaks)
#' hist(pred2, main='linear + quadratic', xlim=xlim, breaks=breaks)
#' hist(pred3, main='all terms', xlim=xlim, breaks=breaks)
#'
#' hist(pred1x1, main='x1: linear only', xlim=xlim, breaks=breaks)
#' hist(pred2x1, main='x1: linear + quadratic', xlim=xlim, breaks=breaks)
#' hist(pred3x1, main='x1: all terms \n(linear, quad, cubic)', xlim=xlim, breaks=breaks)
#'
#' hist(pred1x4, main='x4: linear only', xlim=xlim, breaks=breaks)
#' hist(pred2x4, main='x4: linear + quadratic', xlim=xlim, breaks=breaks)
#' hist(pred3x4, main='x4: all terms \n(linear, quad, cubic)', xlim=xlim, breaks=breaks)
#'
#' # predictions using just x1 and x2
#' par(mfrow=c(1, 1))
#' pred3x1x2 <- predictLars(fit3, XX, type='response', preds=c('x1', 'x2'))
#' plot(pred3, pred3x1x2, xlim=c(0, 1), ylim=c(0, 1))
#' }
#' @export

trainLars <- function(
	data,
	resp = 1,
	preds = 2:ncol(data),
	alphas = c(0.01, seq(0.1, 1, by = 0.1)),
	scale = TRUE,
	quadratic = TRUE,
	cubic = TRUE,
	interaction = TRUE,
	interQuad = TRUE,
	na.rm = FALSE,
	verbose = FALSE,
	...
) {

	#############
	### setup ###
	#############

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# make LARS data object
	x <- makeLarsData(
		data=data,
		resp=resp,
		preds=preds,
		scale=scale,
		quadratic=quadratic,
		cubic=cubic,
		interaction=interaction,
		interQuad=interQuad,
		na.rm=na.rm
	)

	# parse x to just predictors
	y <- x$y
	X <- x$data

	# train model across ridge penalty gradient
	alphaTable <- data.frame(alpha=alphas, lambda=NA, cve=NA, cvse=NA)
	for (i in 1:nrow(alphaTable)) {

		if (verbose) omnibus::say('Training LARS model for alpha = ', alphaTable$alpha[i], '.')

		model <- grpregOverlap::cv.grpregOverlap(X=X, y=y, group=x$groups, alpha=alphaTable$alpha[i], ...)
		# model <- grpregOverlap::cv.grpregOverlap(X=X, y=y, group=x$groups, alpha=alphaTable$alpha[i], family='binomial', penalty='cMCP', nfolds=3)
		alphaTable$lambda[i] <- model$lambda.min
		alphaTable$cve[i] <- model$cve[model$min]
		alphaTable$cvse[i] <- model$cvse[model$min]

	}

	bestAlpha <- alphaTable$alpha[which.min(alphaTable$cve)]
	if (verbose) omnibus::say('Best alpha = ', bestAlpha, '.')

	# train final model
	model <- grpregOverlap::cv.grpregOverlap(X=X, y=y, group=x$groups, alpha=bestAlpha, ...)
	model$lars <- x
	model$alphas <- alphaTable

	if (verbose) print(summary(model))

	class(model) <- c(class(model), 'larsModel')
	model

}

