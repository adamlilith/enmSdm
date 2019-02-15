#' Calibrate a boosted regresion tree (generalied boosting machine) model
#'
#' This function is a wrapper for \code{gbm.step()}. It returns the model with best combination of learning rate, tree depth, and bag fraction based on cross-validated deviance. Models with >=1000 trees are preferred over those with lower deviance.
#' @param data data frame with first column being response
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Character. Name of error family.  See [dismo::gbm.step()].
#' @param learningRate Numeric. Learning rate at which model learns from successive trees (recommended range 0.0001 to 0.1).
#' @param treeComplexity Positive integer. Tree complexity: depth of branches in a single tree (2 to 16).
#' @param bagFraction Numeric in the range (0, 1). Bag fraction: proportion of data used for training in cross-validation (recommended range recommended 0.5 to 0.7).
#' @param maxTrees Positive integer. Maximum number of trees in model set (same as parameter \code{max.trees} in [dismo::gbm.step()]).
#' @param tries Integer > 0. Number of times to try to train a model with a particular set of tuning parameters. The function will stop training the first time a model converges (usually on the first attempt). Non-convergence seems to be related to the number of trees tried in each step.  So if non-convergence occurs then the function automatically increases the number of trees in the step size until \code{tries} is reached.
#' @param tryBy Character list. A list that contains one or more of \code{'learningRate'}, \code{'treeComplexity'}, \code{numTrees}, and/or \code{'stepSize'}. If a given combination of \code{learningRate}, \code{treeComplexity}, \code{numTrees}, \code{stepSize}, and \code{bagFraction} do not allow model convergence then then the function tries again but with alterations to any of the arguments named in \code{tryBy}:
#' * \code{learningRate}: Decrease the learning rate by a factor of 10.
#' * \code{treeComplexity}: Randomly increase/decrease tree complexity by 1 (minimum of 1).
#' * \code{maxTrees}: Increase number of trees by 20%.
#' * \code{stepSize}: Increase step size (argument \code{n.trees} in \code{gbm.step()}) by 50%.
#' If \code{tryBy} is NULL then the function attempts to train the model with the same parameters up to \code{tries} times.
#' @param w Either logical in which case TRUE causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character. Indicates type of value returned. If \code{model} (default) then returns an object of class \code{gbm}. If \code{tuning} then just return a data frame with tuning parameters and deviance of each model sorted by deviance. If both then return a 2-item list with the best model and the tuning table.
#' @param verbose Logical. If TRUE display progress.
#' @param ... Arguments to pass to `gbm.step`.
#' @return If \code{out = 'model'} this function returns an object of class \code{gbm}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and cross-validation deviance for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gbm} object and the data frame.
#' @seealso [dismo::gbm.step()]
#' @examples
#' \donttest{
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
#' out <- trainBrt(
#'     x,
#'     treeComplexity=c(1, 3),
#'     learningRate=c(0.01, 0.001),
#'     bagFraction=0.6,
#'     maxTrees=1000,
#'     out=c('tuning', 'model'),
#'     verbose=TRUE
#' )
#' plot(out$model)
#' out$tuning
#' }
#' @export

trainBrt <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'bernoulli',
	learningRate = c(0.0001, 0.001, 0.01),
	treeComplexity = c(9, 3, 1),
	bagFraction = 0.6,
	maxTrees = 4000,
	tries = 5,
	tryBy = c('learningRate', 'treeComplexity', 'maxTrees', 'stepSize'),
	w = TRUE,
	out = 'model',
	verbose = FALSE,
	...
) {

	#############
	### setup ###
	#############

	# # add dummy variable if using univariate model to avoid errors
	# if (ncol(data)==2) {
		# data$DUMMY <- 1
		# preds <- c(preds, 'DUMMY')
	# }

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# model weights
	if (class(w)[1] == 'logical') {
		w <- if (w & (family == 'binomial' | family == 'quasibinomial')) {
			c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
		} else {
			rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w]
	}

	w <- w / max(w)

	# initialize tuning table
	tuning <- data.frame()

	# initialize lowest cross-validation deviance
	lowestDeviance <- Inf

	bestModel <- NA
	
	############
	### MAIN ###
	############

	# by LEARNING RATE
	for (thisLr in learningRate) {

		# by TREE COMPLEXITY
		for (thisTc in treeComplexity) {

			tempTc <- thisTc

			# by BAG FRACTION
			for (thisBag in bagFraction) {

				# for each parameter combination, train model and get deviance, remember if model with lowest deviance so far
				for (thisMaxTrees in maxTrees) {

					# flag to indicate if model converged or not
					converged <- FALSE
					tempLr <- thisLr
					tempTc <- thisTc
					tempMaxTrees <- thisMaxTrees
					tempStepSize <- 50 # default for n.trees in gbm.step

					# by TRY
					numTries <- 0
					while (numTries <= tries & !converged) {

						numTries <- numTries + 1

						# try with different parameter combinations
						if (numTries > 1 && !is.null(tryBy)) {

							if ('learningRate' %in% tryBy) tempLr <- tempLr / 10
							if ('treeComplexity' %in% tryBy) tempTc <- max(1, tempTc + ifelse(runif(1) > 0.5, 1, -1))
							if ('maxTrees' %in% tryBy) tempMaxTrees <- round(1.2 * tempMaxTrees)
							if ('stepSize' %in% tryBy) tempStepSize <- round(0.8 * tempStepSize)

						}

						# display parameters
						if (verbose) omnibus::say('try: ', numTries, ' | max trees: ', tempMaxTrees, ' | step: ', tempStepSize, ' trees | learning: ', tempLr, ' | complexity: ', tempTc, ' | bag: ', thisBag, post=0)

						# train model... using tryCatch because model may not converge
						model <- tryCatch(
							model <- dismo::gbm.step(
								data=data,
								gbm.x=preds,
								gbm.y=resp,
								family=family,
								tree.complexity=tempTc,
								learning.rate=tempLr,
								bag.fraction=thisBag,
								max.trees=tempMaxTrees,
								n.trees=tempStepSize,
								plot.main=FALSE,
								plot.folds=FALSE,
								silent=TRUE,
								verbose=TRUE,
								site.weights=w,
								...
							),
							error=function(err) return(FALSE)
						)

						# if model training succeeded (model will be gbm object if training succeeded)
						if (class(model) == 'gbm') {

							converged <- TRUE
							dev <- model$cv.statistics$deviance.mean
							if (dev < lowestDeviance && model$gbm.call$best.trees >= 1000) {
								bestModel <- model
								lowestDeviance <- dev
								bestTc <- tempTc
								bestLr <- tempLr
								bestBF <- thisBag
								bestMaxTrees <- tempMaxTrees
								bestStepSize <- tempStepSize
							}

						} else {
							dev <- NA
						}

						if (verbose) omnibus::say('| CV deviance: ', sprintf('%.4f', dev))

						# save tuning table
						tuning <- rbind(
							tuning,
							data.frame(
								learningRate = tempLr,
								treeComplexity = tempTc,
								bagFraction = thisBag,
								maxTrees = tempMaxTrees,
								stepSize = tempStepSize,
								nTrees = ifelse(converged, model$gbm.call$best.trees, NA),
								converged = converged,
								deviance = dev
							)
						)

					} # while trying to train model

					tempLr <- thisLr
					tempTc <- thisTc
					tempMaxTrees <- thisMaxTrees
					tempStepSize <- 50

				} # next max number of trees

			} # next bag fraction

		} # next tree complexity

	} # next learning rate

	# accept models with >=1000 trees
	tuning <- tuning[order(tuning$dev), ]
	tuning$usable <- (tuning$nTrees >= 1000 & !is.na(tuning$deviance))
	if (any(tuning$usable) & any(!tuning$usable)) tuning <- rbind(tuning[tuning$usable, ], tuning[!tuning$usable, ])

	modelRank <- 1:nrow(tuning)
	modelRank[!tuning$converged] <- NA
	tuning <- cbind(modelRank, tuning)

	if (verbose) {
		omnibus::say('')
		print(tuning)
		omnibus::say('')
	}

	# ### get best model with >1000 trees
	# if (tuning$usable[1] & nrow(tuning) > 1) {

		# # if did not remember best model
		# if (tuning$treeComplexity[1] != bestTc | tuning$learningRate[1] != bestLr | tuning$bagFraction[1] != bestBF | tuning$maxTrees[1] != bestMaxTrees | tuning$stepSize[1] != bestStepSize) {

			# if (verbose) omnibus::say('Training best model...')

			# # train model... using tryCatch because model may not converge
			# model <- dismo::gbm.step(
				# data=data,
				# gbm.x=preds,
				# gbm.y=resp,
				# family=family,
				# tree.complexity=tuning$treeComplexity[1],
				# learning.rate=tuning$learningRate[1],
				# bag.fraction=tuning$bagFraction[1],
				# max.trees=tuning$maxTrees[1],
				# n.trees=tuning$stepSize[1],
				# plot.main=FALSE,
				# plot.folds=FALSE,
				# silent=TRUE,
				# verbose=TRUE,
				# site.weights=w,
				# ...
			# )

		# }

	# }

	# return
	if ('model' %in% out & 'tuning' %in% out) {
		out <- list()
		out$tuning <- tuning
		out$model <- bestModel
		out
	} else if ('tuning' %in% out) {
		tuning
	} else {
		bestModel
	}

}

