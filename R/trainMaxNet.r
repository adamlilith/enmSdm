#' Calibrate a Maxent (ver 3.4.0+ or "maxnet") model using AICc
#'
#' This function calculates the "best" Maxent model using AICc across all possible combinations of a set of master regularization parameters and feature classes.  See Warren, D.L. and S.N. Siefert.  2011.  Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria.  \emph{Ecological Applications} 21:335-342.  The function returns the best model and/or a data frame with AICc for each value of the multipler and combination of classes. Here, "best" is defined as the model with the lowest AICc, with ties broken by preferences for: lowest number of estimated parameters, and highest regularization multiplier.
#' @param data  Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh} or any combination of \code{lpqht}), where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If TRUE then test all possible combinations of classes (note that all tested models will at least have linear features). If FALSE then use the classes provided (these will not vary between models).
#' @param out Character. Indicates type of value returned. If \code{model} then returns an object of class \code{maxnet}. If \code{tuning} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param anyway Logical. If no model has fewer coefficients than predictors, return the model with the lowest AICc anyway.
#' @param verbose Logical. If TRUE report progress and AICc table.
#' @param ... Arguments to pass to \code{maxnet()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{maxnet}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log-likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{maxnet} object and the data frame.
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxEnt}}
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
#' model <- trainMaxNet(x, regMult=1:2, out=c('tuning', 'model'), verbose=TRUE)
#' model$tuning
#' summary(model$model)
#' @export

trainMaxNet <- function(
	data,
	resp=names(data)[1],
	preds=names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10),
	classes='default',
	testClasses=TRUE,
	out='model',
	anyway=TRUE,
	verbose=FALSE,
	...
) {

	###########
	## setup ##
	###########

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# get response and predictors
	presentBg <- data[ , resp]
	data <- data[ , preds, drop=FALSE]

	### get combinations of features to test for each regularization multiplier

	if (classes == 'default') {
		classesToTest <- if (ncol(data) > 1) {
			c('l', 'p', 'q', 'h')
		} else {
			c('l', 'q', 'h')
		}
	} else {
		classesToTest <- rep(NA, nchar(classes))
		for (i in 1:nchar(classes)) classesToTest[i] <- substr(classes, i, i)
	}

	# create df of 1/0 to indicate each combination of classes to test
	if (testClasses) {
		classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
		classGrid <- classGrid[-which(rowSums(classGrid) == 0), ]
		classGrid <- as.data.frame(classGrid)
	} else {
		classGrid <- data.frame(matrix(rep(1, length(classesToTest)), nrow=1))
	}

	names(classGrid) <- classesToTest
	if (any(classGrid$l == 0)) classGrid <- classGrid[-which(classGrid$l == 0), ]

	### collate presences and BG sites
	presences <- data[which(presentBg == 1), ]
	if (class(presences) != 'data.frame') presences <- as.data.frame(presences)
	names(presences) <- names(data) # names of data to which to predict

	bg <- data[which(presentBg == 0), ]
	if (class(bg) != 'data.frame') bg <- as.data.frame(bg)
	names(bg) <- names(data)

	##########
	## MAIN ##
	##########

	tuning <- data.frame()
	
	# for each regularization multiplier
	for (thisRegMult in regMult) {

		if (verbose) omnibus::say('Calculating AICc for multiplier ', thisRegMult, ' with features:', post=0)

		# for each combination of class features
		for (countCombo in 1:nrow(classGrid)) {

			if (verbose) omnibus::say(classesToTest[c(classGrid[countCombo, ]) == 1], post=0)

			theseClasses <- paste(classesToTest[as.logical(unlist(classGrid[countCombo, ]))], collapse='')

			# add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
			if (ncol(data) == 1 & theseClasses == 'l') {

				thisData <- data
				thisPresences <- presences
				thisBg <- bg

				thisData$DUMMY <- rep(1, nrow(thisData))
				thisPresences$DUMMY <- rep(1, nrow(presences))
				thisBg$DUMMY <- rep(1, nrow(bg))

			} else {

				thisData <- data
				thisPresences <- presences
				thisBg <- bg

			}

			# # train model
			# model <- maxnet::maxnet(
				# p=as.vector(presentBg),
				# data=thisData,
				# f=maxnet::maxnet.formula(p=as.vector(presentBg), data=thisData, classes=theseClasses),
				# regfun=maxnet.default.regularizationMOD,
				# regmult=thisRegMult,
				# ...
			# )

			# train model
			model <- maxnet::maxnet(
				p=as.vector(presentBg),
				data=thisData,
				f=maxnet::maxnet.formula(p=as.vector(presentBg), data=thisData, classes=theseClasses),
				regfun=maxnet::maxnet.default.regularization,
				regmult=thisRegMult,
				...
			)

			# predict presences
			predPres <- stats::predict(
				object=model,
				newdata=presences,
				type='exponential',
				...
			)

			# predict to background
			predBg <- stats::predict(
				object=model,
				newdata=bg,
				type='exponential',
				...
			)

			rawSum <- sum(c(predPres, predBg), na.rm=TRUE)

			## calculate log likelihood
			ll <- sum(log(predPres / rawSum), na.rm=TRUE)

			## number of parameters
			K <- length(model$betas)

			# AICc
			AICc <- -2 * ll + 2 * K + (2 * K * (K + 1)) / (sum(presentBg) - K - 1)

			# remember
			thisAicFrame <- data.frame(
				regMult=thisRegMult,
				n=sum(presentBg),
				classes=theseClasses,
				logLik=ll,
				K=K,
				AICc=AICc
			)

			tuning <- rbind(tuning, thisAicFrame)

		} # next combination of class features

		if (verbose) omnibus::say('')

	} # next reg mult

	# remove models with more parameters than data points that have more than 0 parameters
	tuning <- tuning[which(tuning$n >= tuning$K & tuning$K > 0), ]

	# re-order frame so one with lowest AICc, number of parameters, and reg mult are used (in that order, used to break ties)
	if (nrow(tuning) > 0) {

		tuning <- tuning[order(tuning$regMult, decreasing=TRUE), ]
		tuning <- tuning[order(tuning$AICc, tuning$K, tuning$regMult), ]

		tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
		tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
		tuning$aicWeight <- tuning$relLike / sum(tuning$relLike)

	}

	if (verbose) {

		omnibus::say('')
		print(tuning)
		omnibus::say('')

	}

	# if user wants best model returned
	if ('model' %in% out) {

		# train model
		if (nrow(tuning) > 0) {

			if (anyway) {

				# add dummy column if doing univariate model to avoid error in maxnet.default.regularizationMOD
				if (ncol(data) == 1 & theseClasses == 'l') {

					thisData <- data
					thisPresences <- presences
					thisBg <- bg

					thisData$DUMMY <- rep(1, nrow(thisData))
					thisPresences$DUMMY <- rep(1, nrow(presences))
					thisBg$DUMMY <- rep(1, nrow(bg))

				} else {

					thisData <- data
					thisPresences <- presences
					thisBg <- bg

				}

				model <- maxnet::maxnet(
					p=as.vector(presentBg),
					data=thisData,
					f=maxnet::maxnet.formula(p=as.vector(presentBg), data=thisData, classes=if (nrow(tuning) == 0) { 1 } else { tuning$classes[1] }),
					regfun=maxnet::maxnet.default.regularization,
					regmult=if (nrow(tuning) == 0) { 1 } else { tuning$regMult[1] },
					...
				)

				if (nrow(tuning) == 0) warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)

			} else {

				warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
				model <- 'No MAXENT model had number of parameters < number of training presences.'

			}

		}

	}

	# return stuff
	if ('model' %in% out & !('tuning' %in% out)) {
		model
	} else if (!('model' %in% out) & 'tuning' %in% out) {
		tuning
	} else {
		list(tuning=tuning, model=model)
	}

}

