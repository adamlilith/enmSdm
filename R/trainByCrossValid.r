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
#' 	\item \code{'logLoss'}: Log loss.
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdm]{contBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdm]{aucWeighted}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdm]{tssWeighted}}.
#' 	\item \code{'fpb'}: Mean Fpb across thresholds from 0 to 1 in steps of 0.01. Calculated with \code{\link[enmSdm]{fpb}}. NOT USED AT PRESENT.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainSe95'} and/or \code{'trainSe90'}: Sensitivity at the threshold that ensures either 95% or 90% of all training presences are classified as "present" (training sensitivity = 0.95 or 0.9).
#' }
#' @param weightEvalTrain Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "train" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied weights they will be used for model calibration.
#' @param weightEvalTest Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "test" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied then weights will be used for model calibration.
#' @param na.rm Logical, if \code{TRUE} then remove \code{NA} predictions before calculating evaluation statistics. If \code{FALSE} (default), propagate \code{NA}s (meaning if predictions contain \code{NA}s, then the evaluation statistic will most likely also be \code{NA}.)
#' @param out Character. Indicates type of value returned. If \code{'models'} then returns a list of a list of candidate models (one sublist per fold. If \code{'tuning'} then just return the evaluation table for candidate models of each fold. If both then return a 2-item list with all candidate models and tuning tables. \emph{WARNING}: Depending on the type of model, using \code{'models'} may produce objects that are very large in memory.
#' @param verbose Numeric. If 0 show no progress updates. If > 0 then show minimal progress updates for this function only. If > 1 show detailed progress for this function. If > 2 show detailed progress plus detailed progress for the "trainXYZ" function.
#' @return If \code{out = 'models'} this function returns a list with models. If If \code{out = 'tuning'} this function returns a data frame with performance statistics for the training and testing data set. If \code{out = c('model', 'tuning'} then it returns a list object with the models and tuning statistics.
#' @return A list object with one element named "models" which will have one set of models per k-fold, and/or another element named "tuning" which will have one data frame per k-fold.
#' @seealso \code{\link[enmSdm]{trainBrt}}, \code{\link[enmSdm]{trainCrf}}, \code{\link[enmSdm]{trainGam}}, \code{\link[enmSdm]{trainGlm}}, \code{\link[enmSdm]{trainMaxEnt}}, \code{\link[enmSdm]{trainLars}}, \code{\link[enmSdm]{trainMaxNet}}, \code{\link[enmSdm]{trainRf}}
#' @examples
#' set.seed(123)
#' 
#' ### contrived example
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
#' model <- trainGlm(trainData)
#' summary(model)
#' 
#' folds <- dismo::kfold(trainData, 3)
#' out <- trainByCrossValid(data=trainData, folds=folds, verbose=1)
#' str(out, 1)
#' head(out$tuning)
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
#' @export
trainByCrossValid <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	folds = dismo::kfold(data),
	trainFx = enmSdm::trainGlm,
	...,
	metrics = c('logLoss', 'cbi', 'auc', 'fpb', 'tss', 'msss', 'mdss', 'minTrainPres', 'trainSe95', 'trainSe90'),
	weightEvalTrain = TRUE,
	weightEvalTest = TRUE,
	na.rm = FALSE,
	out = c('models', 'tuning'),
	verbose = 1
) {

	hasWeights <- ('w' %in% omnibus::ellipseNames(...))

	# total number of models
	foldsClass <- class(folds)
	if (any(c('matrix', 'data.frame') %in% foldsClass)) {
		foldsType <- 'custom'
		foldCodes <- sort(unique(na.omit(folds)))
		numFolds <- ncol(folds)
	} else {
		foldsType <- 'simple'
		foldCodes <- sort(unique(na.omit(folds)))
		numFolds <- length(foldCodes)
	}

	# create list of threshold types to be analyzed
	# have to do this because sensSpecWeighted uses slightly different names
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
				thisFoldCode <- foldCodes[k]
				trainData <- data[foldCodes != thisFoldCode, , drop=FALSE]
				testData <- data[foldCodes == thisFoldCode, , drop=FALSE]
				if (hasWeights) {
					trainWeights <- w[foldCodes != thisFoldCode]
					testWeights <- w[foldCodes == thisFoldCode]
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
				thisFoldCode <- foldCodes[k]
				trainData <- data[foldCodes != thisFoldCode, , drop=FALSE]
				testData <- data[foldCodes == thisFoldCode, , drop=FALSE]
				if (hasWeights) {
					trainWeights <- w[foldCodes != thisFoldCode]
					testWeights <- w[foldCodes == thisFoldCode]
				} else {
					trainWeights <- rep(1, nrow(trainData))
					testWeights <- rep(1, nrow(testData))
				}
			
			}

		### train model
		###############

			thisOut <- trainFx(data=trainData, resp=resp, preds=preds, out=c('models', 'tuning'), verbose=verbose > 2, ...)
			# thisOut <- trainFx(data=trainData, resp=resp, preds=preds, out=c('models', 'tuning'), verbose=verbose > 2)
			kModels <- thisOut$models
			kTuning <- thisOut$tuning

			kTuning <- omnibus::insertCol(data.frame(k = rep(k, nrow(kTuning))), kTuning, at=1)
		
		### evaluate model vs training data
		###################################

			for (countModel in seq_along(kModels)) {

				kModel <- kModels[[countModel]]
			
				# predict to training data
				predToTrain <- predictEnmSdm(model=kModel, newdata=trainData, ...)
				# predToTrain <- predictEnmSdm(model=kModel, newdata=trainData)
				whichPres <- which(trainData[ , resp] == 1)
				whichContrast <- which(trainData[ , resp] == 0)
				predToTrainPres <- predToTrain[whichPres]
				predToTrainContrast <- predToTrain[whichContrast]
				
				# "training" evaluation weights
				if (hasWeights && weightEvalTrain) {
					trainPresWeights <- trainWeights[whichPres]
					trainContrastWeights <- trainWeights[whichContrast]
				} else {
					trainPresWeights <- rep(1, length(whichPres))
					trainContrastWeights <- rep(1, length(whichContrast))
				}	
				
				# predict to testing data
				predToTest <- predictEnmSdm(model=kModel, newdata=testData, ...)
				# predToTest <- predictEnmSdm(model=kModel, newdata=testData)
				whichPres <- which(testData[ , resp] == 1)
				whichContrast <- which(testData[ , resp] == 0)
				predToTestPres <- predToTest[whichPres]
				predToTestContrast <- predToTest[whichContrast]
				
				# "testing" evaluation weights
				if (hasWeights && weightEvalTest) {
					testPresWeights <- testWeights[whichPres]
					testContrastWeights <- testWeights[whichContrast]
				} else {
					testPresWeights <- rep(1, length(whichPres))
					testContrastWeights <- rep(1, length(whichContrast))
				}	

				# log loss

				if ('logLoss' %in% metrics) {

					metricTrain <- sum(log(predToTrainPres) + log(1 - predToTrainContrast))
					metricTest <- sum(log(predToTestPres) + log(1 - predToTestContrast))
					
					if (countModel == 1) kTuning$logLossDelta <- kTuning$logLossTest <- kTuning$logLossTrain <- NA
					kTuning$logLossTrain[countModel] <- metricTrain
					kTuning$logLossTest[countModel] <- metricTest
					kTuning$logLossDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# CBI
				if ('cbi' %in% metrics) {
				
					metricTrain <- contBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					metricTest <- contBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					# metricTrain <- contBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
					# metricTest <- contBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$cbiDelta <- kTuning$cbiTest <- kTuning$cbiTrain <- NA
					kTuning$cbiTrain[countModel] <- metricTrain
					kTuning$cbiTest[countModel] <- metricTest
					kTuning$cbiDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# AUC
				if ('auc' %in% metrics) {
				
					metricTrain <- aucWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
					metricTest <- aucWeighted(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$aucDelta <- kTuning$aucTest <- kTuning$aucTrain <- NA
					kTuning$aucTrain[countModel] <- metricTrain
					kTuning$aucTest[countModel] <- metricTest
					kTuning$aucDelta[countModel] <- metricTrain - metricTest
					
				}
				
				# # Fpb
				# if ('fpb' %in% metrics) {
				
					# metricTrain <- fpb(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					# metricTest <- fpb(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					# metricTrain <- mean(metricTrain, na.rm=na.rm)
					# metricTest <- mean(metricTest, na.rm=na.rm)
					
					# if (countModel == 1) kTuning$fpbDelta <- kTuning$fpbTest <- kTuning$fpbTrain <- NA
					# kTuning$fpbTrain[countModel] <- metricTrain
					# kTuning$fpbTest[countModel] <- metricTest
					# kTuning$fpbDelta[countModel] <- metricTrain - metricTest
					
				# }
				
				# TSS
				if ('tss' %in% metrics) {
				
					metricTrain <- tssWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					metricTest <- tssWeighted(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)

					# metricTrain <- tssWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
					# metricTest <- tssWeighted(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm)
					
					if (countModel == 1) kTuning$tssDelta <- kTuning$tssTest <- kTuning$tssTrain <- NA
					kTuning$tssTrain[countModel] <- metricTrain
					kTuning$tssTest[countModel] <- metricTest
					kTuning$tssDelta[countModel] <- metricTrain - metricTest
					
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
						thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
						# thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, sensitivity=sensitivity, na.rm=na.rm)

						thold <- thresholds[[threshCode]]

						# evaluate
						trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
						testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm)
					
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
		
	} # next fold
	
	# return
	if ('models' %in% out & 'tuning' %in% out) {
		output <- list()
		output$models <- models
		output$tuning <- tuning
		output
	} else if ('tuning' %in% out) {
		tuning
	} else if ('models' %in% out) {
		models
	}
	
}
