#' Calibrate a Maxent (ver 3.3.3- or "maxent") model using AICc
#'
#' This function calculates the "best" Maxent model using AICc across all possible combinations of a set of master regularization parameters and feature classes.  See Warren, D.L. and S.N. Siefert.  2011.  Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria.  **Ecological Applications** 21:335-342.  The function returns the best model and/or a data frame with AICc for each value of the multipler and combination of classes.
#' @param data  Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh} or any combination of \code{lpqht}), where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quardatic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If TRUE then test all possible combinations of classes (note that all tested models will at least have linear features). If FALSE then use the classes provided (these will not varu between models).
#' @param scratchDir Character. Directory to which to write temporary files. Leave as NULL to create a temporary folder in the current working directory.
#' @param forceLinear Logical. If TRUE then require any tested models to include at least linear features.
#' @param jackknife Logical. If TRUE the the returned model will be also include jackknife testing of variable importance.
#' @param out Character. Indicates type of value returned. If \code{model} then returns an object of class \code{maxnet}. If \code{tuning} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param args Character list. Options to pass to \code{maxent()}'s \code{args} argument. (Do not include \code{l}, \code{p}, \code{q}, \code{h}, \code{t}, \code{betamultiplier}, or \code{jackknife}!)
#' @param anyway Logical. If no model has fewer coefficients than predictors, return the model with the lowest AICc anyway.
#' @param verbose Logical. If TRUE report progress and AICc table.
#' @param ... Arguments to pass to \code{maxent()} or \code{predict.maxent()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{MaxEnt}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{MaxEnt} object and the data frame.
#' @details This function is a wrapper for \code{maxent()}. That function relies on a maxent \code{jar} file being placed into the folder \code{./library/dismo/java}. See \code{\link[dismo]{maxent }}for more details. The \code{maxent()} function creates a series of files on disc for each model. This function assumes you do not want those files, so deletes most of them. However, there is one that cannot be deleted and the normal ways of changing its permissions in \code{R} do not work. So the function simply writes over that file (which is allowed) to make it smaller. Regardless, if you run many models your temporary directory (argument \code{scratchDir}) can fill up and require manual deletion.
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxNet}}
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
#' model <- trainMaxEnt(x, regMult=1:2, out=c('model', 'tuning'), verbose=TRUE)
#' model$tuning
#' model$model@lambdas
#' @export

trainMaxEnt <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 6:10, 12.5, 15, 17.5, 20),
	classes = 'default',
	testClasses = TRUE,
	forceLinear = TRUE,
	jackknife = TRUE,
	args = '',
	out = 'model',
	anyway = TRUE,
	scratchDir = NULL,
	verbose = FALSE,
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
		classesToTest <- c('l', 'p', 'q', 'h')
	} else {
		classesToTest <- rep(NA, nchar(classes))
		for (i in 1:nchar(classes)) classesToTest[i] <- substr(classes, i, i)
	}

	if (any('p' %in% classesToTest) & ncol(data) == 1) {
		product <- FALSE
		warning('Data has only one variable so forcing product features to FALSE.')
	}

	# create scratch directory
	scratchDir <- if (is.null(scratchDir)) {
		base::tempfile(pattern='/_maxentTempFiles/')
	} else {
		base::tempfile(pattern='/_maxentTempFiles/', tmpdir=scratchDir)
	}

	dirCreate(scratchDir)
	dirCreate(scratchDir, '/plots')

	## collate all presences
	allPres <- data[presentBg == 1, , drop=FALSE]

	## collate all background sites
	allBg <- data[presentBg == 0, , drop=FALSE]

	# create df of 1/0 to indicate each combination of classes to test
	if (testClasses) {
		classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
		classGrid <- classGrid[-which(rowSums(classGrid) == 0), ]
	} else {
		classGrid <- data.frame(matrix(rep(1, length(classesToTest)), nrow=1))
	}

	names(classGrid) <- classesToTest

	if (forceLinear & any(classGrid$l == 0)) classGrid <- classGrid[-which(classGrid$l == 0), ]

	##########
	## MAIN ##
	##########

	if (verbose) say('Testing models with regularization multiplier:', post=0)

	# by BETA
	for (thisRegMult in regMult) {

		if (verbose) omnibus::say(thisRegMult, post=0)

		# by FEATURE COMBINATION
		for (countParam in 1:nrow(classGrid)) {

			# get parameters
			params <- c(
				paste0('betamultiplier=', thisRegMult),
				paste0('linear=', ifelse('l' %in% names(classGrid) && classGrid$l[countParam] == 1, 'true', 'false')),
				paste0('product=', ifelse('p' %in% names(classGrid) && classGrid$p[countParam] == 1, 'true', 'false')),
				paste0('quadratic=', ifelse('q' %in% names(classGrid) && classGrid$q[countParam] == 1, 'true', 'false')),
				paste0('hinge=', ifelse('h' %in% names(classGrid) && classGrid$h[countParam] == 1, 'true', 'false')),
				paste0('threshold=', ifelse('t' %in% names(classGrid) && classGrid$t[countParam] == 1, 'true', 'false')),
				'jackknife=false'
			)

			if (args != '') params <- c(params, args)

			# train model
			trialModel <- dismo::maxent(
				x=data,
				p=as.vector(presentBg),
				path=scratchDir,
				args=params
			)

			## predict to training (and maybe test presences)
			predPres <- dismo::predict(
				object=trialModel,
				x=allPres,
				na.rm=TRUE,
				args='outputformat=raw',
				...
			)

			## predict to background
			predBg <- dismo::predict(
				object=trialModel,
				x=allBg,
				na.rm=TRUE,
				args='outputformat=raw',
				...
			)

			bgSum <- sum(predBg)

			## calculate log likelihood
			ll <- sum(log(predPres / bgSum), na.rm=TRUE)

			## calculate number of parameters
			K <- 0

			for (thisLambda in trialModel@lambdas) { # for each line in lambda object

				# if not a meta-data line
				if (!grepl(thisLambda, pattern='linearPredictorNormalizer') & !grepl(thisLambda, pattern='densityNormalizer') & !grepl(thisLambda, pattern='numBackgroundPoints') & !grepl(thisLambda, pattern='entropy')) {

					split <- strsplit(thisLambda, ', ')
					paramValue <- as.numeric(split[[1]][2])
					if (paramValue !=0) K <- K + 1 # increment number of parameters

				}

			}

			# AICc
			AICc <- -2 * ll + 2 * K + (2 * K * (K + 1)) / ( sum(presentBg) - K - 1)

			# remember
			thisAicFrame <- data.frame(
				regMult=thisRegMult,
				linear=as.logical(ifelse('l' %in% names(classGrid), classGrid$l[countParam], FALSE)),
				quadratic=as.logical(ifelse('q' %in% names(classGrid), classGrid$q[countParam], FALSE)),
				product=as.logical(ifelse('p' %in% names(classGrid), classGrid$p[countParam], FALSE)),
				hinge=as.logical(ifelse('h' %in% names(classGrid), classGrid$h[countParam], FALSE)),
				threshold=as.logical(ifelse('t' %in% names(classGrid), classGrid$t[countParam], FALSE)),
				numFeats=sum(classGrid[countParam, ]),
				n=sum(presentBg),
				logLik=ll,
				K=K,
				AICc=AICc
			)

			tuning <- if (exists('tuning', inherits=FALSE)) {
				rbind(tuning, thisAicFrame)
			} else {
				thisAicFrame
			}

		} # next set of parameters

	} # next regMult

	# remove models with more parameters than data points
	if (!anyway) tuning <- tuning[which(tuning$n >= tuning$K), ]

	# order by AIc then if ties reg multiplier then if ties by number of parameters
	if (nrow(tuning) > 0) {

		tuning <- tuning[order(tuning$K, decreasing=FALSE), ]
		tuning <- tuning[order(tuning$linear, tuning$quadratic, tuning$threshold, tuning$hinge, tuning$product, decreasing=TRUE), ]
		tuning <- tuning[order(tuning$numFeats, decreasing=FALSE), ]
		tuning <- tuning[order(tuning$regMult, decreasing=TRUE), ]
		tuning <- tuning[order(tuning$AICc, decreasing=FALSE), ]

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
		if (nrow(tuning) == 0 & !anyway) {

			warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
			model <- 'No MAXENT model had number of parameters < number of training presences.'

		} else {

			if (nrow(tuning) > 0) {

				params <- c(
					paste0('betamultiplier=', tuning$regMult[1]),
					paste0('linear=', tolower(tuning$linear[1])),
					paste0('quadratic=', tolower(tuning$quadratic[1])),
					paste0('product=', tolower(tuning$product[1])),
					paste0('hinge=', tolower(tuning$hinge[1])),
					paste0('threshold=', tolower(tuning$threshold[1])),
					paste0('jackknife=', tolower(jackknife))
				)

			} else if (anyway) {

				warning('Returning model with multipler = 1 even though no model had fewer coefficients than predictors.', immediate.=TRUE)

				params <- c(
					paste0('betamultiplier=1'),
					paste0('linear=', tolower(grepl(pattern='l', classes) | classses == 'default')),
					paste0('quadratic=', tolower(grepl(pattern='q', classes) | classses == 'default')),
					paste0('product=', tolower(grepl(pattern='p', classes) | classses == 'default')),
					paste0('hinge=', tolower(grepl(pattern='h', classes) | classses == 'default')),
					paste0('threshold=', tolower(grepl(pattern='t', classes))),
					paste0('jackknife=', tolower(jackknife))
				)

			}

			if (args != '') params <- c(params, args)

			model <- dismo::maxent(
				x=data,
				p=as.vector(presentBg),
				removeDuplicates=FALSE,
				path=scratchDir,
				args=params
			)

		}
	}

	# remove temporary files... note that "species.lambda" file cannot be removed unless R is closed, so we'll just make it smaller to reduce disk space usage
	write.csv(NULL, paste0(scratchDir, '/species.lambdas'))
	if (file.exists(paste0(scratchDir, '/presences'))) write.csv(NULL, paste0(scratchDir, '/presences'))
	if (file.exists(paste0(scratchDir, '/absences'))) write.csv(NULL, paste0(scratchDir, '/absences'))
	unlink(paste0(scratchDir, '/plots'), recursive=TRUE, force=TRUE)
	unlink(scratchDir, recursive=TRUE, force=TRUE)

	# return stuff
	if ('model' %in% out & !('tuning' %in% out)) {
		model
	} else if (!('model' %in% out) & 'tuning' %in% out) {
		tuning
	} else {
		list(tuning=tuning, model=model)
	}

}
