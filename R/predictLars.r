#' Predict using Least Angle Regression (LARS) model
#'
#' This function makes calculates predictions from a \code{grpregOverlap} object constructed using the \code{\link{trainLars}} function..
#' @param object Object of classes \code{larsModel} and \code{cv.grpregOverlapMulti}.
#' @param newdata Data frame.
#' @param type Character. Type of prediction to make: \code{link} produces predictions on the scale of the predictors, \code{response} produces predictions on the scale of the response (i.e., in the range [0, 1] if \code{family = 'binomial'} in \code{trainLars}).
#' @param preds \code{NULL} or a character list. If \code{NULL} then all predictors (and their higher-order terms) will be used in the prediction. If a character vector that is a subset of the predictors used in training the LARS model, then the function yields partial predictions using \emph{only} terms that include the specified variable(s) (plus any quadratic, cubic, interaction, and linear-by-quadratic terms, if the model was trained using them).
#' failOnMissing Logical. If \code{TRUE} then if no variable(s) stated in \code{pred} occurs in the model then cause an error. Variables can be missing either because the user did not include them in the original model or because the variable(s) had no variability and so were removed automatically in the model training process. Note that if \code{failOnMissing} is \code{TRUE} and there is at least one variable in \code{preds} that occurs in the model then this function will still return values (it will not cause an error). If \code{failOnMissing} is \code{FALSE} (default) and no variables in \code{preds} occurs in the model then a warning will be produced and the returned value will be \code{NULL}.
#' @param ... Arguments to pass to \code{\link[grpregOverlap]{predict.grpregOverlap}}.
#' @return Numeric.
#' @seealso \code{\link[grpreg]{grpreg}}, \code{\link[grpreg]{predict.grpreg}}, \code{\link[grpregOverlap]{grpregOverlap}}, \code{\link[grpregOverlap]{predict.grpregOverlap}}, \code{\link[enmSdm]{makeLarsData}}, \code{\link[enmSdm]{trainLars}}
#' @examples
#' # linear regression, a simulation demo
#' (from grpregOverlap() function in the grpregOverlap package)
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

predictLars <- function(
	object,
	newdata,
	type = c('link', 'response'),
	preds = NULL,
	failOnMissing = FALSE,
	...
) {

	#############
	### setup ###
	#############

	# get type of prediction
	type <- match.arg(type)

	# create larsData object
	x <- makeLarsData(
		data=newdata,
		resp=NULL,
		preds=object$lars$preds,
		scale=object$lars$scale,
		quadratic=object$lars$features$quadratic,
		cubic=object$lars$features$cubic,
		interaction=object$lars$features$interaction,
		interQuad=object$lars$features$interQuad
	)

	newdata <- x$data
	# if (class(newdata) == 'numeric') newdata <- t(as.matrix(newdata))

	# add intercept
	newdata <- cbind(rep(1, nrow(newdata)), newdata)
	newdata <- as.matrix(newdata)
	colnames(newdata)[1] <- '(Intercept)'

	newdata <- newdata[ , rownames(object$fit$beta), drop=FALSE]

	# # these functions return NA when object$lambda.min is outside object$beta$lambda
	# out <- grpregOverlap:::predict.cv.grpregOverlap(object=object, X=newdata, type=type, lambda=object$lambda.min, ...)
	# out <- grpregOverlap::predict(object=object, X=newdata, type=type, lambda=object$lambda.min, ...)

	### get betas for model with best lambda (or close to it)
	lambdas <- as.numeric(colnames(object$fit$beta))
	closestToBestLambda <- which.min(abs(lambdas - object$lambda.min))
	betas <- object$fit$beta[ , closestToBestLambda]

	### get just betas pertinent to the specified predictors
	if (!is.null(preds)) {

		wanted <- if (length(preds) == 1) {
			which(colnames(object$lars$inCol) %in% preds)
		} else {
			which(apply(object$lars$inCol[ , preds], 1, any))
		}

		if (length(wanted) == 0) {

			msg <- 'The selected variable(s) are not represented in the model\'s betas. These variable(s) either were not included in model training or were dropped in model training because they had no variability.'

			if (failOnMissing) {
				stop(msg)
			} else {
				warning(msg)
				out <- NULL
			}

		} else {

			# include intercept
			wanted <- c(which(names(betas) %in% '(Intercept)'), wanted + 1)

			betas <- betas[wanted, drop=FALSE]
			newdata <- newdata[ , wanted, drop=FALSE]
			out <- as.numeric(newdata %*% betas)

		}

	} else {
		out <- as.numeric(newdata %*% betas)
	}

	if (!is.null(out) && type == 'response' && object$fit$family == 'binomial') out <- statisfactory::probitAdj(out, 0)

	out

}

