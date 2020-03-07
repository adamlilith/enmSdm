#' Calibrate a Maxent (ver 3.4.0+ or "maxnet") model using AICc
#'
#' This function calculates the "best" Maxent model using AICc across all possible combinations of a set of master regularization parameters and feature classes. The "best" model has the lowest AICc, with ties broken by number of features (fewer is better), regularization multiplier (higher better), then finally the number of coefficients (fewer better). The function can return the best model (default), a list of models created using all possible combinations of feature classes and regularization multipliers, and/or a data frame with tuning statistics for each model. Models in the list and in the data frame are sorted from best to worst. See Warren, D.L. and S.N. Siefert. 2011. Ecological niche modeling in Maxent: The importance of model complexity and the performance of model selection criteria. \emph{Ecological Applications} 21:335-342. 
#' @param data  Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param regMult Numeric vector. Values of the master regularization parameters (called \code{beta} in some publications) to test.
#' @param classes Character list. Names of feature classes to use (either \code{default} to use \code{lpqh}) or any combination of \code{lpqht}, where \code{l} ==> linear features, \code{p} ==> product features, \code{q} ==> quadratic features, \code{h} ==> hinge features, and \code{t} ==> threshold features.
#' @param testClasses Logical. If \code{TRUE} (default) then test all possible combinations of classes (note that all tested models will at least have linear features). If \code{FALSE} then use the classes provided (these will not vary between models).
#' @param dropOverparam Logical, if \code{TRUE} (default), drop models if they have more coefficients than training occurrences. It is possible for no models to fulfill this criterion, in which case no models will be returned.
#' @param anyway Logical. If no model has fewer coefficients than predictors, return the model with the lowest AICc anyway. Default is \code{TRUE}.
#' @param out Character. Indicates type of value returned. If \code{model} then returns an object of class \code{maxnet}. If \code{tuning} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param verbose Logical. If \code{TRUE} report progress and AICc table.
#' @param ... Arguments to pass to \code{maxnet()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{maxnet}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters, log-likelihood, and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{maxnet} object and the data frame.
#' @details The function ranks models by AICc (first), then breaks any ties by sorting by number of features (fewer is better), regularization multiplier (higher is better), then finally the number of coefficients (fewer is better).
#' @seealso \code{\link[maxnet]{maxnet}}, \code{\link[dismo]{maxent}}, \code{\link{trainMaxEnt}}
#' @examples
#' ### model red-bellied lemurs
#' data(mad0)
#' data(lemurs)
#' 
#' # climate data
#' bios <- c(1, 5, 12, 15)
#' clim <- raster::getData('worldclim', var='bio', res=10)
#' clim <- raster::subset(clim, bios)
#' clim <- raster::crop(clim, mad0)
#' 
#' # occurrence data
#' occs <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' occsEnv <- raster::extract(clim, occs[ , c('longitude', 'latitude')])
#' 
#' # background sites
#' bg <- 2000 # too few cells to locate 10000 background points
#' bgSites <- dismo::randomPoints(clim, 2000)
#' bgEnv <- extract(clim, bgSites)
#' 
#' # collate
#' presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
#' env <- rbind(occsEnv, bgEnv)
#' env <- cbind(presBg, env)
#' env <- as.data.frame(env)
#' 
#' preds <- paste0('bio', bios)
#' 
#' regMult <- 1:3 # default values are probably better, but these will be faster
#' 
#' # calibrate MaxEnt model
#' ent <- trainMaxEnt(
#' 	data=env,
#' 	resp='presBg',
#' 	preds=preds,
#' 	regMult=regMult,
#' 	classes='lpq',
#' 	verbose=TRUE
#' )
#' 
#' # calibrate MaxNet model
#' net <- trainMaxNet(
#' 	data=env,
#' 	resp='presBg',
#' 	preds=preds,
#' 	regMult=regMult,
#' 	classes='lpq',
#' 	verbose=TRUE
#' )
#' 
#' # note the differences between the two models...
#' # this is because maxnet() (used by trainMaxNet())
#' # uses an approximation:
#' # (note maxnet() calculates hinges and thresholds differently
#' # so we will turn them off)
#' 
#' data(bradypus, package='dismo')
#' p <- bradypus$presence
#' data <- bradypus[ , 2:3] # easier to inspect betas
#' mn <- maxnet::maxnet(p, data, maxnet.formula(p, data, classes='lpq'))
#' mx <- dismo::maxent(data, p,
#' args=c('linear=true', 'product=true', 'quadratic=true', 'hinge=false',
#' 'threshold=false'))
#' 
#' predMx <- dismo::predict(mx, data)
#' predMn <- predict(mn, data, type='logistic')
#' 
#' plot(predMx, predMn)
#' abline(0, 1)
#' @export

trainMaxNet <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10),
	classes = 'default',
	testClasses = TRUE,
	dropOverparam = TRUE,
	anyway = TRUE,
	out = 'model',
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

	# remember all models and evaluation data
	models <- list()
	tuning <- data.frame()
	
	if (verbose) omnibus::say('Evaluating models with regularization multiplier:', post=0)

	# for each regularization multiplier
	for (thisRegMult in regMult) {

		if (verbose) omnibus::say(thisRegMult, post=0)

		# for each combination of class features
		for (countCombo in 1:nrow(classGrid)) {

			# if (verbose) omnibus::say(classesToTest[c(classGrid[countCombo, ]) == 1], post=0)

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
			numCoeff <- length(model$betas)

			# AICc
			AICc <- -2 * ll + 2 * numCoeff + (2 * numCoeff * (numCoeff + 1)) / (sum(presentBg) - numCoeff - 1)

			# remember
			thisAicFrame <- data.frame(
				regMult=thisRegMult,
				numPres=sum(presentBg),
				classes=theseClasses,
				numClasses=nchar(theseClasses),
				numCoeff=numCoeff,
				logLik=ll,
				AICc=AICc
			)

			models[[length(models) + 1]] <- model
			tuning <- rbind(tuning, thisAicFrame)

		} # next combination of class features

	} # next reg mult

	if (verbose) omnibus::say('')

	# sort from best to worst model
	modelOrder <- order(tuning$AICc, tuning$numClasses, tuning$regMult, tuning$numCoeff, decreasing=c(FALSE, FALSE, TRUE, FALSE))
	tuning <- tuning[modelOrder, ]
	models <- models[modelOrder]

	# remove models with more parameters than data points that have more than 0 parameters
	if (dropOverparam) {
	
		topModel <- models[[1]]
		topTuning <- tuning[1, , drop=FALSE]

		overparamModels <- which(tuning$n < tuning$numCoeff)
		if (length(overparamModels) > 0) {
			tuning <- tuning[-overparamModels, ]
			models <- models[-overparamModels]
		}

		if (length(models) == 0 & anyway) {
			warning('No models had fewer coefficients than predictors. Returning best model anyway.', immediate.=TRUE)
			model <- topModel
			tuning <- topTuning
		} else if (length(models) == 0 & !anyway) {
			warning('No models had fewer coefficients than predictors. No model(s) returned.', immediate.=TRUE)
			model <- NA
		} else {
			model <- models[[1]]
		}
			
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

	# return
	if (length(out) > 1) {
		output <- list()
		if ('models' %in% out) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuning
		output
	} else if (out == 'models') {
		models
	} else if (out == 'model') {
		model
	} else if (out == 'tuning') {
		tuning
	}

}

