#' Summarize distribution/niche model cross-validation object
#'
#' This function summarizes models calibrated using the \code{\link[enmSdm]{trainByCrossValid}} function. It returns aspects of the best models across k-folds (the particular aspects depends on the kind of models used).
#' @param x An object of class \code{crossValid} (which is also a list). Note that the object \emph{must} include a sublist named \code{tuning}.
#' @param trainFxName Character, name of function used to train the SDM (examples: \code{'trainGlm'}, \code{'trainMaxEnt'}, \code{'trainBrt'})
#' @param metric Metric by which to select the best model in each k-fold. This can be any of the columns that appear in the data frames in \code{x$tuning} (or any columns added manually), but typically is one of the following \emph{plus} either \code{Train}, \code{Test}, or \code{Delta} (e.g., \code{'logLossTrain'}, \code{'logLossTest'}, or \code{'logLossDelta'}):
#' \itemize{
#' 	\item \code{'logLoss'}: Log loss.
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdm]{contBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdm]{aucWeighted}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdm]{tssWeighted}}.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainSe95'} and/or \code{'trainSe90'}: Sensitivity at the threshold that ensures either 95% or 90% of all training presences are classified as "present" (training sensitivity = 0.95 or 0.9).
#' }
#' @param decreasing Logical, if \code{TRUE} (default), for each k-fold sort models by the value listed in \code{metric} in decreasing order (highest connotes "best", lowest "worst"). If \code{FALSE} use the lowest value of \code{metric}.
#' @return Data frame with statistics on the best set of models across k-folds. Depending on the model algorithm, this could be:
#' \itemize{
#' 	\item BRTs (boosted regression trees): Learning rate, tree complexity, and bag fraction.
#' 	\item GLMs (generalized linear models): Frequency of use of each term in the best models.
#' 	\item Maxent: Frequency of times each specific combination of feature classes was used in the best models plus mean master regularization multiplier for each feature set.
#' 	\item NSs (natural splines): Data frame, one row per fold and one column per predictor, with values representing the maximum degrees of freedom used for each variable in the best model of each fold.
#' }
#' @seealso \code{\link[enmSdm]{trainByCrossValid}}
#' @examples
#' \dontrun{
#' set.seed(123)
#' ### contrived example
#' # generate training/testing data
#' n <- 10000
#' x1 <- seq(-1, 1, length.out=n) + rnorm(n)
#' x2 <- seq(10, 0, length.out=n) + rnorm(n)
#' x3 <- rnorm(n)
#' y <- 2 * x1 + x1^2 - 10 * x2 - x1 * x2
#' y <- statisfactory::probitAdj(y, 0)
#' y <- y^3
#' presAbs <- runif(n) < y
#' data <- data.frame(presAbs=presAbs, x1=x1, x2=x2, x3=x3)
#'
#' model <- trainGlm(data)
#' summary(model)
#'
#' folds <- dismo::kfold(data, 3)
#' out <- trainByCrossValid(data, folds=folds, verbose=1)
#'
#' summaryByCrossValid(out)
#'
#' str(out, 1)
#' head(out$tuning[[1]])
#' head(out$tuning[[2]])
#' head(out$tuning[[3]])
#'
#' # can do following for each fold (3 of them)
#' lapply(out$models[[1]], coefficients)
#' sapply(out$models[[1]], logLik)
#' sapply(out$models[[1]], AIC)
#'
#' # select model for k = 1 with greatest CBI
#' top <- which.max(out$tuning[[1]]$cbiTest)
#' summary(out$models[[1]][[top]])
#'
#' # in fold k = 1, which models perform well but aren not overfit?
#' plot(out$tuning[[1]]$cbiTrain, out$tuning[[1]]$cbiTest, pch='.',
#' 		main='Model Numbers for k = 1')
#' abline(0, 1, col='red')
#' numModels <- nrow(out$tuning[[1]])
#' text(out$tuning[[1]]$cbiTrain, out$tuning[[1]]$cbiTest, labels=1:numModels)
#' usr <- par('usr')
#' x <- usr[1] + 0.9 * (usr[4] - usr[3])
#' y <- usr[3] + 0.1 * (usr[4] - usr[3])
#' text(x, y, labels='overfit', col='red', xpd=NA)
#' x <- usr[1] + 0.1 * (usr[4] - usr[3])
#' y <- usr[3] + 0.9 * (usr[4] - usr[3])
#' text(x, y, labels='suspicious', col='red', xpd=NA)
#' }
#' @export
summaryByCrossValid <- function(
	x,
	trainFxName = 'trainGlm',
	metric = 'cbiTest',
	decreasing = TRUE
) {

	trainFxName <- tolower(trainFxName)

	tuning <- x$tuning

	### order each tuning attempt by performance
	for (k in seq_along(tuning)) {

		tuningOrder <- order(tuning[[k]][ , metric], decreasing=decreasing)
		tuning[[k]] <- tuning[[k]][tuningOrder, ]

	}

	### BRT
	#######
	if (trainFxName == tolower('trainBrt')) {

		# get list of terms in best models
		params <- data.frame()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]
			if (thisTuning$converged[1]) {

				params <- rbind(
					params,
					data.frame(
						learningRate = thisTuning$learningRate[1],
						treeComplexity = thisTuning$treeComplexity[1],
						bagFraction = thisTuning$bagFraction[1],
						nTrees = thisTuning$nTrees[1]
					)
				)

			}

		}

		# summarize best models
		if (nrow(params) > 0) {

			out <- rbind(
				apply(params, 2, min, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.25, na.rm=TRUE),
				apply(params, 2, mean, na.rm=TRUE),
				apply(params, 2, stats::median, na.rm=TRUE),
				apply(params, 2, stats::quantile, 0.75, na.rm=TRUE),
				apply(params, 2, max, na.rm=TRUE)
			)

			out <- as.data.frame(out)
			out$treeComplexity <- max(1, round(out$treeComplexity))
			out$nTrees <- round(out$nTrees)

			value <- data.frame(value=c('min', '25th percent quantile', 'mean', 'median', '75th percent quantile', 'max'))
			out <- omnibus::insertCol(value, out, 1)

		} else {
			out <- data.frame()
		}


	} else if (trainFxName == tolower('trainGlm')) {

		# get list of terms in best models
		term <- character()
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1]
			thisTerm <- strsplit(thisTerm, ' ')[[1]]
			thisTerm <- thisTerm[3:length(thisTerm)]
			if (any(thisTerm == 'presBg')) thisTerm <- thisTerm[-which(thisTerm == 'presBg')]
			if (any(thisTerm == '~')) thisTerm <- thisTerm[-which(thisTerm == '~')]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			if (length(thisTerm) == 0) thisTerm <- 1

			term <- c(term, thisTerm)

		}

		term <- unique(term)
		out <- data.frame(term)
		out$frequency <- 0

		# tally usage of best term(s)
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			thisTerm <- thisTuning$model[1]
			thisTerm <- strsplit(thisTerm, ' ')[[1]]
			thisTerm <- thisTerm[3:length(thisTerm)]
			if (any(thisTerm == '1')) thisTerm <- thisTerm[-which(thisTerm == '1')]
			if (any(thisTerm == '+')) thisTerm <- thisTerm[-which(thisTerm == '+')]
			if (any(thisTerm == '-')) thisTerm <- thisTerm[-which(thisTerm == '-')]

			for (countTerm in seq_along(thisTerm)) {
				whichTerm <- which(out$term == thisTerm[countTerm])
				out$frequency[whichTerm] <- out$frequency[whichTerm] + 1
			}

		}

		out <- out[order(out$frequency, decreasing = TRUE), ]
		out$proportionOfModels <- out$frequency / length(x$tuning)

	### MAXENT
	##########

	} else if (trainFxName == tolower('trainMaxEnt')) {

		## tally frequency of feature class combinations across best models and mean regularization associated with each set of features

		# make vector of all possible feature combinations
		simpleFeats <- c('l', 'p', 'q', 'h')
		featCombos <- vector('list', length(simpleFeats))
		for (i in seq_along(simpleFeats)) {
			featCombos[[i]] <- utils::combn(simpleFeats, i, paste, collapse = '')
		}

		featCombos <- unique(unlist(featCombos))

		featFreqs <- featRegMultSums <- rep(0, length(featCombos))
		names(featFreqs) <- names(featRegMultSums) <- featCombos

		# tally number of times each feature combination used in top models
		targetFeatLengths <- nchar(featCombos)
		for (k in seq_along(tuning)) {

			thisModelFeats <- tuning[[k]]$classes[1]
			thisModelFeats <- strsplit(thisModelFeats, '')[[1]]
			modelFeatLength <- length(thisModelFeats)

			modelFeatInTarget <- rep(TRUE, length(featCombos)) # flags if potential features are *all* used in model
			names(modelFeatInTarget) <- featCombos
			for (countModelFeat in seq_along(thisModelFeats)) {

				thisModelFeatInTarget <- grepl(pattern=thisModelFeats[countModelFeat], featCombos)
				modelFeatInTarget <- (modelFeatInTarget & thisModelFeatInTarget)

			}

			modelFeatInTarget <- modelFeatInTarget & (modelFeatLength == targetFeatLengths)

			featFreqs[modelFeatInTarget] <- featFreqs[modelFeatInTarget] + 1
			featRegMultSums[modelFeatInTarget] <- featRegMultSums[modelFeatInTarget] + tuning[[k]]$regMult[1]

		}

		featRegMultMeans <- featRegMultSums / featFreqs
		# featFreqs <- featFreqs / sum(featFreqs)

		out <- data.frame(
			featureSet = featCombos,
			frequencyInBestModels = featFreqs,
			meanRegMult = featRegMultMeans
		)

	### NS
	######

	} else if (trainFxName == tolower('trainNs')) {

		# get predictor names
		strings <- colnames(tuning[[1]])
		strings <- strings[grepl(strings, pattern='splines::ns')]
		strings <- strsplit(strings, 'splines::ns')
		preds <- character()
		for (i in seq_along(strings)) {
			string <- strings[[i]][2]
			pred <- substr(string, 2, regexpr(string, pattern=' df') - 2)
			preds <- c(preds, pred)
		}

		preds <- sort(unique(preds))

		# create table to hold information on degrees of freedom for each top model
		subOut <- data.frame(DUMMY=NA)
		colnames(subOut) <- preds[1]
		if (length(preds) > 1) {
			for (i in 2:length(preds)) {
				subOut$DUMMY <- NA
				colnames(subOut)[ncol(subOut)] <- preds[i]
			}
		}
		subOut <- subOut[rep(1, length(tuning)), ]

		# find df for each predictor
		for (k in seq_along(tuning)) {

			thisTuning <- tuning[[k]]

			for (pred in preds) {

				colPattern <- paste0('splines::ns\\(', pred, ', df')
				cols <- colnames(thisTuning)[grepl(colnames(thisTuning), pattern=colPattern)]

				# get degrees of freedom used to evaluate this variable
				dfs <- numeric()
				for (countCol in seq_along(cols)) {
					before <- paste0('splines::ns\\(', pred, ', df')
					thisDf <- strsplit(cols[countCol], before)[[1]][2]
					thisDf <- strsplit(thisDf, ' ')[[1]][3]
					thisDf <- substr(thisDf, 1, nchar(thisDf) - 1)
					thisDf <- as.numeric(thisDf)
					dfs <- c(dfs, thisDf)
				}

				vals <- thisTuning[1, cols, drop=FALSE]
				bestDf <- c(!apply(vals, 1, FUN=is.na))
				bestDf <- if (all(!bestDf)) {
					NA
				} else {
					max(dfs[bestDf])
				}

				subOut[k, pred] <- bestDf

			} # next predictor

		} # next model

		for (pred in preds) {
			thisPredOut <- data.frame(
				term=pred,
				frequencyInBestModels=sum(!is.na(subOut[ , pred])),
				minDf=min(subOut[ , pred], na.rm=TRUE),
				quantile25thDf=stats::quantile(subOut[ , pred], 0.25, na.rm=TRUE),
				meanDf=mean(subOut[ , pred], na.rm=TRUE),
				medianDf=stats::median(subOut[ , pred], na.rm=TRUE),
				quantile7thDf=stats::quantile(subOut[ , pred], 0.75, na.rm=TRUE),
				maxDf=max(subOut[ , pred], na.rm=TRUE)
			)

			out <- if (exists('out', inherits=FALSE)) {
				out <- rbind(out, thisPredOut)
			} else {
				thisPredOut
			}
		}


	} # NSs

	rownames(out) <- NULL
	out

}
