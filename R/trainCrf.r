#' Calibrate a conditional random forest model.
#'
#' This function trains a  conditional random forest model. It is nearly identical to the \code{\link[party]{cforest }} function in the \pkg{party} package but is included for constistancy with \code{\link{trainGlm}}, \code{\link{trainGam}}, and similar functions.
#' @param data Data frame.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}). Default is to use the 'binomial' family.
#' @param w Either logical in which case TRUE causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param ... Arguments to pass to [party::cforest()] function.
#' @return Object of class \code{RandomForest}.
#' @seealso \code{cforest} (in the \pkg{party} package), \code{\link{trainRf}}
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
#' model <- trainCrf(x)
#' predict(model, newdata=x)
#' @export

trainCrf <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	w = ifelse(family == 'binomial', TRUE, FALSE),
	...
) {

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# get just desired columns
	data <- data[ , c(resp, preds)]

	# model weights
	if (length(w) == 1 && class(w) == 'logical') {
		w <- if (w & family == 'binomial') {
			c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
		} else {
			rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w]
	}

	# binomial response
	if (family == 'binomial') data[ , resp] <- factor(data[ , resp], levels=0:1)

	# formula
	form <- as.formula(paste(resp, '~ .'))

	# train model
	model <- party::cforest(
		formula=form,
		data=data,
		weights=w,
		...
	)

	model

}
