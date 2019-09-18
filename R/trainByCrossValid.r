#' Calibrate a distribution/niche model using cross-validation
#'
#' This function is an extension of any of the "trainXYZ" functions for calibrating species distribution and ecological niche models. This function uses the "trainXYZ" function to calibrate and evaluate a suite of models using cross-validation. The models are evaluated against withheld data to determine the optimal settings for a "final" model using all available data.
#' @param data Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param folds Either a numeric vector, or matrix or data frame:
#' \itemize{
#' 	\item If a vector, there must be one value per row in \code{data}. If there are \emph{K} unique values in the vector, then \emph{K} unique models will be trained. Each model will use all of the data except for rows that match a particular value in the \code{folds} vector. For example, if \code{folds = c(1, 1, 1, 2, 2, 2, 3, 3, 3)}, then three models will be trained, one with all rows that match the 2s and 3s, one with all rows matching 1s and 2s, and one will all rows matching 1s and 3s. The models will be evaluated against the withheld data and against the training data. Use \code{NA} to exclude rows from all testing/training. The default is to construct 5 folds of roughly equal size.
#' \item If a matrix or data frame, there must be one row per row in \code{data}. Each column corresponds to a different model to be trained. For a given column there should be only two unique values (plus possibly \code{NA}s). When sorted, the lesser value will be used to identify the calibration data and the upper value the evaluation data. For example, a particular column could contain 1s, 2, and \code{NA}s. Data rows corresponding to 1s will be used as training data and rows corresponding to 2s as test data (after rows with \code{NA} are dropped). This option is useful for creating spatially-structured cross-validation folds where training and test sites are separated (spatially) by censored (ignored) data.
#' }
#' @param trainFx Function, name of the "trainXYZ" function to use (e.g., \code{\link[enmSdm]{trainGlm}} or \code{\link[enmSdm]{trainGam}}).
#' @param ... Arguments to pass to the "trainXYZ" function, as well as the \code{\link{predictEnmSdm}}, \code{\link{contBoyce}}, \code{\link{tssWeighted}}, and \code{\link{thresholdWeighted}} functions. See help for the appropriate function.
#' @param metrics Character vector, names of evaluation metrics to calculate. The default is to calculate all of:
#' \itemize{
#' 	\item \code{'logLoss'}: Log loss. Higher (less negative) values connote better fit.
#' 	\item \code{'logLossEqualWeight'}: Log loss with the total weight of presences and contrast sites equalized before calculation (i.e., weight of training presences equalized to weight of training contrast sites and weight of test presences equalized to weight of test contrast sites). Higher (less negative) values connote better fit. This version is better for cases where the number of presences and contrast sites is unequal because it gives both equal total weight.
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdm]{contBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdm]{aucWeighted}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdm]{tssWeighted}}.
#' 	\item \code{'fpb'}: Mean Fpb across thresholds from 0 to 1 in steps of 0.01. Calculated with \code{\link[enmSdm]{fpb}}. NOT USED AT PRESENT.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'orss'}: Maximum odds ratio skill score (Wunderlich et al. 2019).
#' 	\item \code{'sedi'}: Maximum of symmetrical extremal dependence index (Wunderlich et al. 2019).
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainSe95'} and/or \code{'trainSe90'}: Sensitivity at the threshold that ensures either 95% or 90% of all training presences are classified as "present" (training sensitivity = 0.95 or 0.9).
#' }
#' @param weightEvalTrain Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "train" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied weights they will be used for model calibration.
#' @param weightEvalTest Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "test" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied then weights will be used for model calibration.
#' @param na.rm Logical, if \code{TRUE} then remove \code{NA} predictions before calculating evaluation statistics. If \code{FALSE} (default), propagate \code{NA}s (meaning if predictions contain \code{NA}s, then the evaluation statistic will most likely also be \code{NA}.)
#' @param out Character. Indicates type of value returned. If \code{'models'} then returns a list of a list of candidate models (one sublist per fold. If \code{'tuning'} then just return the evaluation table for candidate models of each fold. If both then return a 2-item list with all candidate models and tuning tables. \emph{WARNING}: Depending on the type of model, using \code{'models'} may produce objects that are very large in memory.
#' @param verbose Numeric. If 0 show no progress updates. If > 0 then show minimal progress updates for this function only. If > 1 show detailed progress for this function. If > 2 show detailed progress plus detailed progress for the "trainXYZ" function.
#' @return A list object with several named elements:
#' \itemize{
#' 		\item \code{meta}: Meta-data on the model call.
#' 		\item \code{folds}: the \code{folds} object.
#' 		\item \code{models} (if \code{'models'} is in argument \code{out}): A list of model objects, one per k-fold
#'		\item \code{tuning} (if \code{'tuning'} is in argument \code{out}): One data frame per k-fold, each containing evaluation statistics for all candidate models in the fold.
#' }
#' @references Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation} 24:38-49.
#' @references Wunderlich, R.F., Lin, P-Y., Anthony, J., and Petway, J.R. 2019. Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models. Nature Conservation 35:97-116.
#' @seealso \code{\link[enmSdm]{trainBrt}}, \code{\link[enmSdm]{trainCrf}}, \code{\link[enmSdm]{trainGam}}, \code{\link[enmSdm]{trainGlm}}, \code{\link[enmSdm]{trainMaxEnt}}, \code{\link[enmSdm]{trainLars}}, \code{\link[enmSdm]{trainMaxNet}}, \code{\link[enmSdm]{trainRf}}
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
#' # in fold k = 1, which models perform well but are not overfit?
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
trainByCrossValid <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	folds = dismo::kfold(data),
	trainFx = enmSdm::trainGlm,
	...,
	metrics = c('logLoss', 'logLossEqualWeight', 'cbi', 'auc', 'fpb', 'tss', 'msss', 'mdss', 'orss', 'sedi', 'minTrainPres', 'trainSe95', 'trainSe90'),
	sensitivity = 0.9,
	weightEvalTrain = TRUE,
	weightEvalTest = TRUE,
	na.rm = FALSE,
	out = c('models', 'tuning'),
	verbose = 1
) {

	hasWeights <- ('w' %in% omnibus::ellipseNames(...))

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# total number of models
	foldsClass <- class(folds)
	if (any(c('matrix', 'data.frame') %in% foldsClass)) {
		foldsType <- 'custom'
		foldCodes <- na.omit(sort(unique(c(unlist(folds)))))
		numFolds <- ncol(folds)

		# codes for train/test sets
		trainCode <- foldCodes[1]
		testCode <- foldCodes[2]

	} else {
		foldsType <- 'simple'
		foldCodes <- na.omit(sort(unique(folds)))
		numFolds <- length(foldCodes)
	}

	# create list of threshold types to be analyzed
	# have to do this because thresholdStats uses slightly different names
	# since it does not know if a set of presences is test/train
	threshTypes <- character()
	if ('msss' %in% metrics) threshTypes <- c(threshTypes, 'msss')
	if ('mdss' %in% metrics) threshTypes <- c(threshTypes, 'mdss')
	if ('minTrainPres' %in% metrics) threshTypes <- c(threshTypes, 'minTrainPres')
	if ('trainSe95' %in% metrics) threshTypes <- c(threshTypes, 'trainSe95')
	if ('trainSe90' %in% metrics) threshTypes <- c(threshTypes, 'trainSe90')

	### by fold
	###########
	
	# models and their evaluation statistics
	tuning <- models <- list()

	for (k in 1:numFolds) {
	
		if (verbose > 0) omnibus::say('Modeling k = ', k, ' on ', date(), '...', post=ifelse(verbose > 1, 2, 1), pre=ifelse(verbose > 1, 2, 0))
	
		### get training/testing data according to type of folds being used
		###################################################################

			if (foldsType == 'simple') {

				# copy data
				thisData <- data
			
				# drop NAs
				nas <- which(is.na(folds))
				if (length(nas) > 0) {
					thisData <- thisData[-nas, , drop=FALSE]
					if (hasWeights) w <- w[-nas]
				}
			
				# train/test data
				testCode <- foldCodes[k]
				trainData <- data[folds != testCode, c(resp, preds), drop=FALSE]
				testData <- data[folds == testCode, c(resp, preds), drop=FALSE]
				if (hasWeights) {
					trainWeights <- w[folds != testCode]
					testWeights <- w[folds == testCode]
				} else {
					trainWeights <- rep(1, nrow(trainData))
					testWeights <- rep(1, nrow(testData))
				}
			
			} else {
			
				# copy data and get folds codes
				thisData <- data
				thisFolds <- folds[ , k, drop=TRUE]
				foldCodes <- sort(unique(na.omit(thisFolds)))
				if (length(foldCodes) != 2) stop('"folds" matrix/data frame must contain only two unique values (aside from NAs).')
			
				# drop NAs
				nas <- which(is.na(thisFolds))
				if (length(nas) > 0) {
					thisData <- thisData[-nas, , drop=FALSE]
					if (hasWeights) w <- w[-nas]
				}
			
				# train/test data
				trainData <- data[which(thisFolds == trainCode), c(resp, preds), drop=FALSE]
				testData <- data[which(thisFolds == testCode), c(resp, preds), drop=FALSE]
				if (hasWeights) {
					trainWeights <- w[which(thisFolds == trainCode)]
					testWeights <- w[which(thisFolds == testCode)]
				} else {
					trainWeights <- rep(1, nrow(trainData))
					testWeights <- rep(1, nrow(testData))
				}
			
			}

		### train model
		###############

			# thisOut <- trainFx(data=trainData, resp=resp, preds=preds, out=c('models', 'tuning'), w=trainWeights, verbose=verbose > 2, ...)
			thisOut <- trainFx(data=trainData, resp=resp, preds=preds, out=c('models', 'tuning'), w=trainWeights, verbose=verbose > 2)
			kModels <- thisOut$models
			kTuning <- thisOut$tuning

			rm(thisOut)

			kTuning <- omnibus::insertCol(data.frame(k = rep(k, nrow(kTuning))), kTuning, at=1)

		### indices and weights for this fold
		#####################################
		
			# training
			whichTrainPres <- which(trainData[ , resp] == 1)
			whichTrainContrast <- which(trainData[ , resp] == 0)

			# "training" evaluation weights
			if (hasWeights && weightEvalTrain) {
				trainPresWeights <- trainWeights[whichTrainPres]
				trainContrastWeights <- trainWeights[whichTrainContrast]
			} else {
				trainPresWeights <- rep(1, length(whichTrainPres))
				trainContrastWeights <- rep(1, length(whichTrainContrast))
			}	

			# testing
			whichTestPres <- which(testData[ , resp] == 1)
			whichTestContrast <- which(testData[ , resp] == 0)

			# "testing" evaluation weights
			if (hasWeights && weightEvalTest) {
				testPresWeights <- testWeights[whichTestPres]
				testContrastWeights <- testWeights[whichTestContrast]
			} else {
				testPresWeights <- rep(1, length(whichTestPres))
				testContrastWeights <- rep(1, length(whichTestContrast))
			}	

			insert <- data.frame(
				numTrainPres = if (na.rm) { length(na.omit(whichTrainPres)) } else { length(whichTrainPres) },
				numTrainContrast = if (na.rm) { length(na.omit(whichTrainContrast)) } else { length(whichTrainContrast) },
				numTestPres = if (na.rm) { length(na.omit(whichTestPres)) } else { length(whichTestPres) },
				numTestContrast = if (na.rm) { length(na.omit(whichTestContrast)) } else { length(whichTestContrast) },
				trainPresWeights = sum(trainPresWeights, na.rm=na.rm),
				trainContrastWeights = sum(trainContrastWeights, na.rm=na.rm),
				testPresWeights = sum(testPresWeights, na.rm=na.rm),
				testContrastWeights = sum(testContrastWeights, na.rm=na.rm)
			)
			
			insert <- insert[rep(1, nrow(kTuning)), ]
				
			kTuning <- omnibus::insertCol(insert, kTuning, at='k', before=FALSE)
				
		### evaluate model vs training data
		###################################

			for (countModel in seq_along(kModels)) {

				kModel <- kModels[[countModel]]
			
				# predict to training data
				predToTrain <- predictEnmSdm(model=kModel, newdata=trainData, ...)
				# predToTrain <- predictEnmSdm(model=kModel, newdata=trainData)
				predToTrainPres <- predToTrain[whichTrainPres]
				predToTrainContrast <- predToTrain[whichTrainContrast]
				
				# predict to testing data
				predToTest <- predictEnmSdm(model=kModel, newdata=testData, ...)
				# predToTest <- predictEnmSdm(model=kModel, newdata=testData)
				predToTestPres <- predToTest[whichTestPres]
				predToTestContrast <- predToTest[whichTestContrast]
				
				# log loss
				if ('logLoss' %in% metrics) {

					metricTrain <- mean(c(trainPresWeights * log(predToTrainPres), trainContrastWeights * log(1 - predToTrainContrast)), na.rm=na.rm)
					metricTest <- mean(c(testPresWeights * log(predToTestPres), testContrastWeights * log(1 - predToTestContrast)), na.rm=na.rm)
					
					if (countModel == 1) kTuning$logLossDelta <- kTuning$logLossTest <- kTuning$logLossTrain <- NA
					kTuning$logLossTrain[countModel] <- metricTrain
					kTuning$logLossTest[countModel] <- metricTest
					kTuning$logLossDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# log loss, equal total weight
				if ('logLossEqualWeight' %in% metrics) {

					# rescale weights
					totalTrainPresWeight <- sum(trainPresWeights, na.rm=na.rm)
					totalTestPresWeight <- sum(testPresWeights, na.rm=na.rm)

					totalTrainContrastWeight <- sum(trainContrastWeights, na.rm=na.rm)
					totalTestContrastWeight <- sum(testContrastWeights, na.rm=na.rm)

					thisTrainPresWeights <- trainPresWeights
					thisTrainContrastWeights <- trainContrastWeights

					thisTestPresWeights <- testPresWeights
					thisTestContrastWeights <- testContrastWeights

					if (totalTrainPresWeight > totalTrainContrastWeight) {
						thisTrainContrastWeights <- thisTrainContrastWeights * (totalTrainPresWeight / totalTrainContrastWeight)
					} else {
						thisTrainPresWeights <- thisTrainPresWeights * (totalTrainContrastWeight / totalTrainPresWeight)
					}

					if (totalTestPresWeight > totalTestContrastWeight) {
						thisTestContrastWeights <- thisTestContrastWeights * (totalTestPresWeight / totalTestContrastWeight)
					} else {
						thisTestPresWeights <- thisTestPresWeights * (totalTestContrastWeight / totalTestPresWeight)
					}

					maxWeight <- max(c(thisTrainPresWeights, thisTrainContrastWeights))
					thisTrainPresWeights <- thisTrainPresWeights / maxWeight
					thisTrainContrastWeights <- thisTrainContrastWeights / maxWeight

					# log loss
					metricTrain <- mean(c(thisTrainPresWeights * log(predToTrainPres), trainContrastWeights * log(1 - predToTrainContrast)), na.rm=na.rm)
					metricTest <- mean(c(thisTestPresWeights * log(predToTestPres), thisTestContrastWeights * log(1 - predToTestContrast)), na.rm=na.rm)
					
					if (countModel == 1) kTuning$logLossDeltaEqualWeight <- kTuning$logLossTestEqualWeight <- kTuning$logLossTrainEqualWeight <- NA
					kTuning$logLossTrainEqualWeight[countModel] <- metricTrain
					kTuning$logLossTestEqualWeight[countModel] <- metricTest
					kTuning$logLossDeltaEqualWeight[countModel] <- metricTrain - metricTest
					
				}
				
				# CBI
				if ('cbi' %in% metrics) {
				
					# metricTrain <- contBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- contBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)
					
					metricTrain <- contBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- contBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$cbiDelta <- kTuning$cbiTest <- kTuning$cbiTrain <- NA
					kTuning$cbiTrain[countModel] <- metricTrain
					kTuning$cbiTest[countModel] <- metricTest
					kTuning$cbiDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# AUC
				if ('auc' %in% metrics) {
				
					metricTrain <- aucWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- aucWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$aucDelta <- kTuning$aucTest <- kTuning$aucTrain <- NA
					kTuning$aucTrain[countModel] <- metricTrain
					kTuning$aucTest[countModel] <- metricTest
					kTuning$aucDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# # Fpb
				# if ('fpb' %in% metrics) {
				
					# metricTrain <- fpb(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- fpb(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)
					
					# metricTrain <- mean(metricTrain, na.rm=na.rm)
					# metricTest <- mean(metricTest, na.rm=na.rm)
					
					# if (countModel == 1) kTuning$fpbDelta <- kTuning$fpbTest <- kTuning$fpbTrain <- NA
					# kTuning$fpbTrain[countModel] <- metricTrain
					# kTuning$fpbTest[countModel] <- metricTest
					# kTuning$fpbDelta[countModel] <- metricTrain - metricTest
					
				# }
				
				# TSS
				if ('tss' %in% metrics) {
				
					# metricTrain <- tssWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- tssWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)

					metricTrain <- tssWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- tssWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$tssDelta <- kTuning$tssTest <- kTuning$tssTrain <- NA
					kTuning$tssTrain[countModel] <- metricTrain
					kTuning$tssTest[countModel] <- metricTest
					kTuning$tssDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# ORSS
				if ('orss' %in% metrics) {
				
					# metricTrain <- orssWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- orssWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)

					metricTrain <- orssWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- orssWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
					metricTrain <- max(metricTrain, na.rm=TRUE)
					metricTest <- max(metricTest, na.rm=TRUE)
					
					if (countModel == 1) kTuning$orssDelta <- kTuning$orssTest <- kTuning$orssTrain <- NA
					kTuning$orssTrain[countModel] <- metricTrain
					kTuning$orssTest[countModel] <- metricTest
					kTuning$orssDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# SEDI
				if ('sedi' %in% metrics) {
				
					# metricTrain <- sediWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- sediWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm, ...)

					metricTrain <- sediWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- sediWeighted(pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
					metricTrain <- max(metricTrain, na.rm=TRUE)
					metricTest <- max(metricTest, na.rm=TRUE)
					
					if (countModel == 1) kTuning$sediDelta <- kTuning$sediTest <- kTuning$sediTrain <- NA
					kTuning$sediTrain[countModel] <- metricTrain
					kTuning$sediTest[countModel] <- metricTest
					kTuning$sediDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# threshold-dependent
				if (length(threshTypes) > 0) {
					
					for (thisThreshType in threshTypes) {

						# threshold
						if (thisThreshType == 'minTrainPres') {
							threshCode <- 'minPres'
						} else if (thisThreshType == 'minTrainPres') {
							threshCode <- 'minPres'
							sensitivity <- 0
						} else if (thisThreshType == 'trainSe95') {
							threshCode <- 'sensitivity'
							sensitivity <- 0.95
						} else if (thisThreshType == 'trainSe90') {
							threshCode <- 'sensitivity'
							sensitivity <- 0.9
						} else {
							threshCode = thisThreshType
							sensitivity <- 0
						}
						
						# get thresholds
						# thresholds <- thresholdWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm, ...)
						thresholds <- thresholdWeighted(pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, sensitivity=sensitivity, na.rm=na.rm)

						thold <- thresholds[[threshCode]]

						# evaluate
						trainEval <- thresholdStats(thold, pres=predToTrainPres, contrast=predToTrainContrast, presWeight=trainPresWeights, contrastWeight=trainContrastWeights, na.rm=na.rm)
						testEval <- thresholdStats(thold, pres=predToTestPres, contrast=predToTestContrast, presWeight=testPresWeights, contrastWeight=testContrastWeights, na.rm=na.rm)
					
						# remember
						if (countModel == 1) {
						
							kTuning$TEMP1 <- kTuning$TEMP2 <- kTuning$TEMP3 <- kTuning$TEMP4 <- kTuning$TEMP5 <- kTuning$TEMP6 <- kTuning$TEMP7 <- NA
							names(kTuning)[(ncol(kTuning) - 6):ncol(kTuning)] <- paste0(thisThreshType, c('Thold', 'SeTrain', 'SeTest', 'SeDelta', 'SpTrain', 'SpTest', 'SpDelta'))
							
						}
						
						kTuning[countModel, paste0(thisThreshType, 'Thold')] <- thold
						kTuning[countModel, paste0(thisThreshType, 'SeTrain')] <- trainEval[['sensitivity']]
						kTuning[countModel, paste0(thisThreshType, 'SeTest')] <- trainEval[['sensitivity']] - testEval[['sensitivity']]
						kTuning[countModel, paste0(thisThreshType, 'SeDelta')] <- testEval[['sensitivity']]
						kTuning[countModel, paste0(thisThreshType, 'SpTrain')] <- trainEval[['specificity']]
						kTuning[countModel, paste0(thisThreshType, 'SpTest')] <- testEval[['specificity']]
						kTuning[countModel, paste0(thisThreshType, 'SpDelta')] <- trainEval[['specificity']] - testEval[['specificity']]
						
					} # next threshold-dependent evaluation metric
					
				} # threshold-specific metrics
				
			} # next model in this k-fold
		
		if ('tuning' %in% out) tuning[[k]] <- kTuning
		if ('models' %in% out) models[[k]] <- kModels
		gc()
		
	} # next fold
	
	### return
	##########
	output <- list()
	
	# meta data
	meta <- list(
		resp = resp,
		preds = preds,
		metrics = metrics,
		sensitivity = sensitivity,
		weightEvalTrain = weightEvalTrain,
		weightEvalTest = weightEvalTest,
		na.rm = na.rm,
		out = c('models', 'tuning')
	)
	
	ellipses <- list(...)
	if (length(ellipses) > 0) {
		dotNames <- omnibus::ellipseNames(...)
		for (i in seq_along(ellipses)) {
			meta <- c(meta, ellipses[i])
		}
	}
	
	output$meta <- meta
	output$folds <- folds
	
	# models and tuning
	if ('models' %in% out) output$models <- models
	if ('tuning' %in% out) output$tuning <- tuning

	class(output) <- c('crossValid', class(output))
	output
	
}
