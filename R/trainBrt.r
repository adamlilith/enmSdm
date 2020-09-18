#' Calibrate a boosted regression tree (generalized boosting machine) model
#'
#' This function is a wrapper for \code{\link[dismo]{gbm.step}}. It returns the model with best combination of learning rate, tree depth, and bag fraction based on cross-validated deviance. It can also return a table with deviance of different combinations of tuning parameters that were tested, and all of the models tested. See Elith, J., J.R. Leathwick, and T. Hastie. 2008. A working guide to boosted regression trees. \emph{Journal of Animal Ecology} 77:802-813.
#' @param data data frame with first column being response
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Character. Name of error family.  See \code{\link[dismo]{gbm.step}}.
#' @param learningRate Numeric. Learning rate at which model learns from successive trees (Elith et al. 2008 recommend 0.0001 to 0.1).
#' @param treeComplexity Positive integer. Tree complexity: depth of branches in a single tree (1 to 16).
#' @param bagFraction Numeric in the range [0, 1]. Bag fraction: proportion of data used for training in cross-validation (Elith et al. 2008 recommend 0.5 to 0.7).
#' @param minTrees Positive integer. Minimum number of trees to be scored as a "usable" model (Elith et al. 2008 recommend at least 1000). Default is 1000.
#' @param maxTrees Positive integer. Maximum number of trees in model set (same as parameter \code{max.trees} in \code{\link[dismo]{gbm.step}}).
#' @param tries Integer > 0. Number of times to try to train a model with a particular set of tuning parameters. The function will stop training the first time a model converges (usually on the first attempt). Non-convergence seems to be related to the number of trees tried in each step.  So if non-convergence occurs then the function automatically increases the number of trees in the step size until \code{tries} is reached.
#' @param tryBy Character list. A list that contains one or more of \code{'learningRate'}, \code{'treeComplexity'}, \code{numTrees}, and/or \code{'stepSize'}. If a given combination of \code{learningRate}, \code{treeComplexity}, \code{numTrees}, \code{stepSize}, and \code{bagFraction} do not allow model convergence then then the function tries again but with alterations to any of the arguments named in \code{tryBy}:
#' * \code{learningRate}: Decrease the learning rate by a factor of 10.
#' * \code{treeComplexity}: Randomly increase/decrease tree complexity by 1 (minimum of 1).
#' * \code{maxTrees}: Increase number of trees by 20%.
#' * \code{stepSize}: Increase step size (argument \code{n.trees} in \code{gbm.step()}) by 50%.
#' If \code{tryBy} is NULL then the function attempts to train the model with the same parameters up to \code{tries} times.
#' @param w Either logical in which case \code{TRUE} (default) causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) \emph{or} a numeric list of weights, one per row in \code{data} \emph{or} the name of the column in \code{data} that contains site weights. If \code{FALSE}, then each datum gets a weight of 1.
#' @param out Character. Indicates type of value returned. If \code{model} (default) then returns an object of class \code{gbm}. If \code{models} then all models that were trained are returned in a list in the order they appear in the tuning table (this may take a lot of memory!). If \code{tuning} then just return a data frame with tuning parameters and deviance of each model sorted by deviance. If both then return a 2-item list with the best model and the tuning table.
#' @param cores Integer >= 1. Number of cores to use when calculating multiple models. Default is 1.
#' @param verbose Logical. If \code{TRUE} display progress.
#' @param ... Arguments to pass to \code{\link[dismo]{gbm.step}}.
#' @return If \code{out = 'model'} this function returns an object of class \code{gbm}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and cross-validation deviance for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gbm} object and the data frame. Note that if a model does not converge or does not meet sufficiency criteria (i.e., the number of optimal trees is < \code{minTrees}, then the model is not returned (a \code{NULL} value is returned for \code{'model'} and models are simply missing from the \code{tuning} and \code{models} output.
#' @seealso \code{\link[dismo]{gbm.step}}
#' @examples
#' \donttest{
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
#' # settings... defaults probably better, but these are faster
#' lr <- c(0.001, 0.1)
#' tc <- c(1, 3)
#' maxTrees <- 2000
#' set.seed(123)
#' model <- trainBrt(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' 	learningRate = lr,
#' 	treeComplexity = tc,
#' 	maxTrees = maxTrees,
#' 	verbose = TRUE
#' )
#' 
#' plot(model)
#' 
#' # prediction raster
#' nTrees <- model$gbm.call$n.trees
#' map <- predict(clim, model, type='response', n.trees=nTrees)
#' plot(map)
#' points(occs[ , c('longitude', 'latitude')])
#' 
#' }
#' @export

trainBrt <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'bernoulli',
	learningRate = c(0.0001, 0.001, 0.01),
	treeComplexity = c(5, 3, 1),
	bagFraction = 0.6,
	minTrees = 1000,
	maxTrees = 8000,
	tries = 5,
	tryBy = c('learningRate', 'treeComplexity', 'maxTrees', 'stepSize'),
	w = TRUE,
	out = 'model',
	cores = 1,
	verbose = FALSE,
	...
) {
	
	### setup 
	#########

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
			w <- if (w) {
				c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
			} else {
				rep(1, nrow(data))
			}
		} else if (class(w) == 'character') {
			w <- data[ , w]
		}

		w <- w / max(w)

	### generate table of parameterizations
	#######################################
		
		params <- expand.grid(learningRate=learningRate, treeComplexity=treeComplexity, bagFraction=bagFraction, maxTrees=maxTrees)
		
	### MAIN
	########

		if (cores > 1) {
			`%makeWork%` <- foreach::`%dopar%`
			cl <- parallel::makeCluster(cores)
			doParallel::registerDoParallel(cl)
		} else {
			`%makeWork%` <- foreach::`%do%`
		}
		
		mcOptions <- list(preschedule=TRUE, set.seed=TRUE, silent=FALSE)
		
		# work <- foreach::foreach(i=1:nrow(tuning), .options.multicore=mcOptions, .combine='c', .inorder=FALSE, .export=c('.trainBrtWorker'), .packages = c('gbm')) %makeWork%
		work <- foreach::foreach(i=1:nrow(params), .options.multicore=mcOptions, .combine='c', .inorder=FALSE, .export=c('.trainBrtWorker')) %makeWork%
			.trainBrtWorker(
				i = i,
				params = params,
				data = data,
				preds = preds,
				resp = resp,
				family = family,
				learningRate = learningRate,
				treeComplexity = treeComplexity,
				bagFraction = bagFraction,
				minTrees = minTrees,
				maxTrees = maxTrees,
				tries = tries,
				tryBy = tryBy,
				w = w,
				...
			)
				
		if (cores > 1) parallel::stopCluster(cl)

	### collate models
	##################
		
		models <- list()
		tuning <- data.frame()

		for (i in seq_along(work)) {
		
			models[[i]] <- work[[i]]$model
			tuning <- rbind(tuning, work[[i]]$workerTuning)
		
		}
		
	### process models
	##################

		# remove non-converged models
		keeps <- which(tuning$converged)
		tuning <- tuning[keeps, , drop=FALSE]
		models <- models[keeps, drop=FALSE]

		if (length(models) > 0) {
		
			# remove models with fewer trees than required
			keeps <- which(omnibus::naCompare('>=', tuning$nTrees, minTrees))
			tuning <- tuning[keeps, , drop=FALSE]
			models <- models[keeps, drop=FALSE]
			
			if (length(models) > 0) {
			
				# sort from best to worst model
				modelOrder <- order(tuning$dev, decreasing=FALSE)
				tuning <- tuning[modelOrder, , drop=FALSE]
				models <- models[modelOrder]
				
				rownames(tuning) <- 1:nrow(tuning)
				
			}
				
		}
		
	### return
	##########
		
		if (verbose) {
			omnibus::say('')
			print(tuning, digits=4)
			omnibus::say('')
		}

		if (length(out) > 1) {
			output <- list()
			if ('models' %in% out) output$models <- models
			if ('model' %in% out) output$model <- models[[1]]
			if ('tuning' %in% out) output$tuning <- tuning
			output
		} else if (out == 'models') {
			models
		} else if (out == 'model') {
			print(str(models, 2))
			if (length(models) > 0) {
				models[[1]]
			} else {
				NULL
			}
		} else if (out == 'tuning') {
			tuning
		}

}


#######################
### worker function ###
#######################

.trainBrtWorker <- function(
	i,								# iterator
	params,							# parameterizations
	data,							# data frame
	preds,							# character
	resp,							# character
	family,							# character
	learningRate,					# learning rate
	treeComplexity,					# tree depth
	bagFraction,					# bag fraction
	minTrees,						# minimum number of trees in a model
	maxTrees,						# maximum number of trees in a model
	tries,							# number of times to try if non-convergence
	tryBy,							# one or more of c('learningRate', 'treeComplexity', 'maxTrees', 'stepSize')
	w,								# weights (numeric vector),
	...								# other (to pass to step.gbm)
) {

	# flag to indicate if model converged or not
	converged <- FALSE

	# starter values
	tempLr <- params$learningRate[i]
	tempTc <- params$treeComplexity[i]
	tempBf <- params$bagFraction[i]
	tempMaxTrees <- params$maxTrees[i]
	
	tempStepSize <- 50 # default for n.trees in gbm.step

	# tuning table
	workerTuning <- data.frame()
	
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

		# train model... using tryCatch because model may not converge
		model <- tryCatch(
			model <- dismo::gbm.step(
				data=data,
				gbm.x=preds,
				gbm.y=resp,
				family=family,
				tree.complexity=tempTc,
				learning.rate=tempLr,
				bag.fraction=tempBf,
				max.trees=tempMaxTrees,
				n.trees=tempStepSize,
				plot.main=FALSE,
				plot.folds=FALSE,
				silent=TRUE,
				verbose=TRUE,
				site.weights=w,
				...
			),
			error=function(err) return(NULL)
		)

		# if model training succeeded (model will be gbm object if training succeeded)
		if (!is.null(model)) {

			converged <- TRUE
			dev <- model$cv.statistics$deviance.mean

		} else {
			dev <- NA
		}

		# save tuning table
		workerTuning <- rbind(
			workerTuning,
			data.frame(
				learningRate = tempLr,
				treeComplexity = tempTc,
				bagFraction = tempBf,
				maxTrees = tempMaxTrees,
				stepSize = tempStepSize,
				nTrees = ifelse(converged, model$gbm.call$best.trees, NA),
				converged = converged,
				deviance = dev
			)
		)

	} # while trying to train model

	workerOut <- list(
		list(
			model=model,
			workerTuning=workerTuning
		)
	)
	
	workerOut
	
}
