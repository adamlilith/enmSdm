#' Calibrate a Maxent (ver 3.3.3- or "maxent") model using AICc
#'
#' This function calculates the "best" Maxent model using AICc across all possible combinations of a set of master regularization parameters and feature classes.  See Warren, D.L. and S.N. Siefert.  2011.  Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria.  **Ecological Applications** 21:335-342.  The function returns the best model and/or a data frame with AICc for each value of the multiplier and combination of classes.
#' @param data  Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh} or any combination of \code{lpqht}), where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical.  If \code{TRUE} (default) then test all possible combinations of classes (note that all tested models will at least have linear features). If \code{FALSE} then use the classes provided (these will not vary between models).
#' @param scratchDir Character. Directory to which to write temporary files. Leave as NULL to create a temporary folder in the current working directory.
#' @param forceLinear Logical. If \code{TRUE} (default) then require any tested models to include at least linear features.
#' @param jackknife Logical. If \code{TRUE} (default) the the returned model will be also include jackknife testing of variable importance.
#' @param out Character. Indicates type of value returned. Values can be \code{'model'} (default; return model with lowest AICc), \code{'models'} (return a list of all models), and/or \code{'tuning'} (return a data frame with AICc for each model). If more than one value is specified, then the output will be a list with elements named "model", "models", and/or "tuning". If \code{'models'} is specified, they will only be produced if \code{select = TRUE}. The models will appear in the list in same order as they appear in the tuning table (i.e., model with the lowest AICc first, second-lowest next, etc.). If just one value is specified, the output will be either an object of class \code{glm}, a list with objects of class \code{glm}, or a data frame.
#' @param dropOverparam Logical, if \code{TRUE} (default), drop models if they have more coefficients than training occurrences. It is possible for no models to fulfill this criterion, in which case no models will be returned.
#' @param args Character list. Options to pass to \code{maxent()}'s \code{args} argument. (Do not include \code{l}, \code{p}, \code{q}, \code{h}, \code{t}, \code{betamultiplier}, or \code{jackknife}!)
#' @param anyway Logical. Same as \code{dropOverparam} (included for backwards compatibility. If \code{NULL} (default), then the value of \code{dropOverparam} will take precedence. If \code{TRUE} or \code{FALSE} then \code{anyway} will override the value of \code{dropOverparam}.
#' @param verbose Logical. If TRUE report progress and AICc table.
#' @param ... Arguments to pass to \code{maxent()} or \code{predict.maxent()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{MaxEnt}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{MaxEnt} object and the data frame.
#' @details This function is a wrapper for \code{maxent()}. That function relies on a maxent \code{jar} file being placed into the folder \code{./library/dismo/java}. See \code{\link[dismo]{maxent}} for more details. The \code{maxent()} function creates a series of files on disc for each model. This function assumes you do not want those files, so deletes most of them. However, there is one that cannot be deleted and the normal ways of changing its permissions in \code{R} do not work. So the function simply writes over that file (which is allowed) to make it smaller. Regardless, if you run many models your temporary directory (argument \code{scratchDir}) can fill up and require manual deletion.
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxNet}}
#' @examples
#' set.seed(123)
#' 
#' # contrived example
#' n <- 10000
#' x1 <- seq(-1, 1, length.out=n) + rnorm(n)
#' x2 <- seq(10, 0, length.out=n) + rnorm(n)
#' x3 <- rnorm(n)
#' y <- 2 * x1 + x1^2 - 10 * x2 - x1 * x2
#' 
#' y <- statisfactory::probitAdj(y, 0)
#' y <- y^3
#' hist(y)
#' 
#' presAbs <- runif(n) < y
#' 
#' trainData <- data.frame(presAbs=presAbs, x1=x1, x2=x2, x3=x3)
#' 
#' out <- trainMaxEnt(trainData, regMult=1:2,
#' 	out=c('models', 'model', 'tuning'))
#' str(out)
#' out$model@lambdas
#' out$tuning
#' 
#' predsLogistic <- raster::predict(out$model, trainData)
#' predsLogistic <- predictMaxEnt(out$model, trainData, type='logistic') # slow
#' predsCloglog <- predictMaxEnt(out$model, trainData)
#' plot(predsLogistic, predsCloglog, xlim=c(0, 1), ylim=c(0, 1))
#' abline(0, 1, col='gray')
#' @export

trainMaxEnt <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10),
	classes = 'default',
	testClasses = TRUE,
	forceLinear = TRUE,
	jackknife = TRUE,
	args = '',
	dropOverparam = TRUE,
	out = 'model',
	scratchDir = NULL,
	verbose = FALSE,
	anyway = TRUE,
	...
) {

	###########
	## setup ##
	###########

	if (!is.null(anyway)) dropOverparam <- anyway
	
	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# get response and predictors
	presentBg <- data[ , resp]
	data <- data[ , preds, drop=FALSE]

	### get combinations of features to test for each regularization multiplier
	classesToTest <- if (classes == 'default') {
		c('l', 'p', 'q', 'h')
	} else {
		unlist(strsplit(classes, ''))
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

	omnibus::dirCreate(scratchDir)
	omnibus::dirCreate(scratchDir, '/plots')

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

	if (verbose) omnibus::say('Testing models with regularization multiplier:', post=0)

	# remember all models and evaluation data
	models <- list()
	tuning <- data.frame()

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
				paste0('jackknife=', ifelse(jackknife, 'true', 'false'))
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

			rawSum <- sum(c(predPres, predBg), na.rm=TRUE)

			## log likelihood
			ll <- sum(log(predPres / rawSum), na.rm=TRUE)

			## number of parameters
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
			AICc <- -2 * ll + 2 * K + (2 * K * (K + 1)) / (sum(presentBg) - K - 1)

			# remember
			thisTuning <- data.frame(
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

			models[[length(models) + 1]] <- trialModel
			tuning <- rbind(tuning, thisTuning)

		} # next set of parameters

	} # next regMult
	
	# sort models by AICc
	modelOrder <- order(tuning$AICc, tuning$regMult, tuning$numFeats, tuning$linear, tuning$quadratic, tuning$threshold, tuning$hinge, tuning$product, tuning$K)
	tuning <- tuning[modelOrder, ]
	models <- models[modelOrder]

	# remove models with more parameters than data points
	bestModel <- models[[1]]
	bestTuning <- tuning[1, , drop=FALSE]
	overparamModels <- which(tuning$n < tuning$K)
	if (length(overparamModels) > 0 & !dropOverparam) {
		tuning <- tuning[-overparamModels, ]
		models <- rlist::list.remove(models, overparamModels)
	}

	if (dropOverparam & length(models) == 0) {
		tuning <- bestTuning
		models[[1]] <- bestModel
		model <- bestModel
	} else {
		model <- bestModel
	}
	
	# AICc weights
	if (nrow(tuning) > 0) {

		tuning$deltaAICc <- tuning$AICc - min(tuning$AICc)
		tuning$relLike <- exp(-0.5 * tuning$deltaAICc)
		tuning$aicWeight <- tuning$relLike / sum(tuning$relLike)

		rownames(tuning) <- 1:nrow(tuning)
	
	}

	if (verbose) {

		omnibus::say('')
		print(tuning)
		omnibus::say('')

	}

	# if user wants best model returned
	if ('model' %in% out & nrow(tuning) == 0 & (!anyway | !dropOverparam)) {

			warning('No models had fewer coefficients than predictors. No model returned.', immediate.=TRUE)
			model <- NA

	}

	# remove temporary files... note that "species.lambda" file cannot be removed unless R is closed, so we'll just make it smaller to reduce disk space usage
	write.csv(NULL, paste0(scratchDir, '/species.lambdas'))
	if (file.exists(paste0(scratchDir, '/presences'))) write.csv(NULL, paste0(scratchDir, '/presences'))
	if (file.exists(paste0(scratchDir, '/absences'))) write.csv(NULL, paste0(scratchDir, '/absences'))
	unlink(paste0(scratchDir, '/plots'), recursive=TRUE, force=TRUE)
	unlink(scratchDir, recursive=TRUE, force=TRUE)

	# return
	if (length(out) > 1) {
		output <- list()
		if ('models' %in% out) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuning
		output
	} else if ('models' %in% out) {
		models
	} else if ('model' %in% out) {
		model
	} else if ('tuning' %in% out) {
		tuning
	}

}
