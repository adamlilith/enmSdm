#' Calibrate a distribution/niche model using cross-validation
#'
#' This function is an extension of any of the "trainXYZ" functions for calibrating species distribution and ecological niche models. The "trainXYZ" functions use various criteria to select a "best" model across a suite of models using different sets of variables and/or functional forms and parameterizations. This function uses the "trainXYZ" function to train a suite of models typically calibrated on somewhat non-overlapping data.  The models are evaluated against withheld data to determine the optimal settings for a "final" model using all available data.
#' @param data Data frame or matrix. Environmental predictors (and no other fields) for presences and background sites.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param folds Either a numeric vector, or matrix or data frame:
#' \itemize{
#' 	\item If a vector, there must be one value per row in \code{data}. If there are \emph{K} unique values in the vector, then \emph{K} unique models will be trained. Each model will use all of the data except for rows that match a particular value in the \code{folds} vector. For example, if \code{folds = c(1, 1, 1, 2, 2, 2, 3, 3, 3), then three models will be trained, one with all rows that match the 2s and 3s, one with all rows matching 1s and 2s, and one will all rows matching 1s and 3s. The models will be evaluated against the withheld data and against the training data. Use \code{NA} to exclude rows from all testing/training. The default is to construct 5 folds of roughly equal size.
#' \item If a matrix or data frame, there must be one row per row in \code{data}. Each column corresponds to a different model to be trained. For a given column there should be only two unique values (plus possibly \code{NA}s). When sorted, the lesser value will be used to identify the calibration data and the upper value the evaluation data. For example, a particular column could contain 1s, 2, and \code{NA}s. Data rows corresponding to 1s will be used as training data and rows corresponding to 2s as test data (after rows with \code{NA} are dropped). This option is useful for creating spatially-structured cross-validation folds where training and test sites are separated (spatially) by censored (ignored) data.
#' }
#' @param trainFx Function, name of the "trainXYZ" function to use (e.g., \code{\link[enmSdm]{trainGlm}} or \code{\link[enmSdm]{trainGam}}).
#' @param ... Arguments to pass to the "trainXYZ" function and \code{\link{predictEnmSdm}} functions. See help for the appropriate function.
#' @param metrics Character vector, names of evaluation metrics to calculate. The default is to calculate all of:
#' \itemize
#' 	\item \code{'cbi'}: Continuous Boyce Index (CBI). Calculated with \code{\link[enmSdm]{contBoyce}}.
#' 	\item \code{'auc'}: Area under the receiver-operator characteristic curve (AUC). Calculated with \code{\link[enmSdm]{aucWeighted}}.
#' 	\item \code{'tss'}: Maximum value of the True Skill Statistic. Calculated with \code{\link[enmSdm]{tssWeighted}}.
#' 	\item \code{'fpb'}: Mean Fpb across thresholds from 0 to 1 in steps of 0.01. Calculated with \code{\link[enmSdm]{fpb}}. NOT USED AT PRESENT.
#' 	\item \code{'msss'}: Sensitivity and specificity calculated at the threshold that maximizes sensitivity (true presence prediction rate) plus specificity (true absence prediction rate).
#' 	\item \code{'mdss'}: Sensitivity (se) and specificity (sp) calculated at the threshold that minimizes the difference between sensitivity and specificity.
#' 	\item \code{'minTrainPres'}: Sensitivity and specificity at the greatest threshold at which all training presences are classified as "present".
#' 	\item \code{'trainPres05'} and/or \code{'trainPres10'}: Sensitivity or specificity at the threshold that ensures the training sensitivity is either 0.95 or 0.9 (i.e., classifies 5% or 10% of all training presences as "absences/background").
#' }
#' @param weightEvalTrain Logical, if \code{TRUE} (default) and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "train" version of evaluation statistics}. If \code{FALSE}, there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied weights they will be used for model calibration.
#' @param weightEvalTest Logical, if \code{TRUE} and an argument named \code{w} is specified in \code{...}, then evaluation statistics that support weighting will use the weights specified by \code{w} \emph{for the "test" version of evaluation statistics}. If \code{FALSE} (default), there will be no weighting of test sites. Note that this applies \emph{only} to the calculation of evaluation statistics.  If \code{w} is supplied then weights will be used for model calibration.
#' @param na.rm Logical, if \code{TRUE} then remove \code{NA} predictions before calculating evaluation statistics. If \code{FALSE} (default), propagate \code{NA}s (meaning if predictions contain \code{NA}s, then the evaluation statistic will most likely also be \code{NA}.)
#' @param out Character. Indicates type of value returned. If \code{model} then returns a model object (e.g., of class \code{MaxEnt}, \code{glm}, etc.). If \code{tuning} then just return the evaluation table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the tuning table.
#' @param verbose Numeric. If 0 show no progress updates. If > 0 then show minimal progress updates for this function only. If > 1 show detailed progress for this function. If > 2 show detailed progress plus detailed progress for the "trainXYZ" function.
#' @return If \code{out = 'models'} this function returns a list with models. If If \code{out = 'tuning'} this function returns a data frame with performance statistics against the training and testing data set. If \code{out = c('model', 'tuning'} then it returns a list object with the models and tuuning statistics.
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
#' 
#' out <- trainByCrossValid(data=trainData, verbose=2)
#' out$tuning
#' lapply(out$models, coefficients)
#' lapply(out$models, logLik)
#' @export
trainByCrossValid <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	folds = dismo::kfold(data),
	trainFx = enmSdm::trainGlm,
	...,
	metrics = c('cbi', 'auc', 'fpb', 'tss', 'msss', 'mdss', 'minTrainPres', 'trainPres05', 'trainPres10'),
	w = NULL,
	weightEvalTrain = TRUE,
	weightEvalTest = TRUE,
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

	### by fold
	###########
	
	# evaluation data frame
	tuning <- data.frame()
	
	# models
	models <- list()
	
	for (k in 1:numFolds) {
	
		if (verbose > 0) omnibus::say('Modeling k = ', k, ' on ', date(), '...', post=ifelse(verbose > 1, 2, 1), pre=ifelse(verbose > 1, 2, 1))
	
		# evaluation statistics
		thisTuning <- data.frame(k=k)
			
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

			# kModel <- trainFx(data=trainData, resp=resp, preds=preds, ..., out='model', verbose=verbose >= 2)
			models[[k]] <- trainFx(data=trainData, resp=resp, preds=preds, out='model', verbose=verbose > 2)
		
		### evaluate model vs training data
		###################################
			
			# predict to training data
			# predToTrain <- predictEnmSdm(model=kModel, newdata=trainData, ...)
			predToTrain <- predictEnmSdm(model=kModel, newdata=trainData)
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
			# predToTest <- predictEnmSdm(model=kModel, newdata=testData, ...)
			predToTest <- predictEnmSdm(model=kModel, newdata=testData)
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
			
			# CBI
			if ('cbi' %in% metrics) {
				cbiTrain <- contBoyce(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
				cbiTest <- contBoyce(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
				
				this <- data.frame(cbiTrain=cbiTrain, cbiTest=cbiTest)
				thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
				
			}
			
			# AUC
			if ('auc' %in% metrics) {
				aucTrain <- aucWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
				aucTest <- aucWeighted(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm)
				
				this <- data.frame(aucTrain=aucTrain, aucTest=aucTest)
				thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
				
			}
			
			# # Fpb
			# if ('fpb' %in% metrics) {
				# fpbTrain <- fpb(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
				# fpbTest <- fpb(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
				
				# fpbTrain <- mean(fpbTrain, na.rm=na.rm)
				# fpbTest <- mean(fpbTest, na.rm=na.rm)
				
				# this <- data.frame(fpbTrain=fpbTrain, fpbTest=fpbTest)
				# thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
			# }
			
			# TSS
			if ('tss' %in% metrics) {
				tssTrain <- tssWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
				tssTest <- tssWeighted(pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
				
				this <- data.frame(tssTrain=tssTrain, tssTest=tssTest)
				thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
				
			}
			
			# threshold-dependent
			if (any('msss', 'mdss', 'minTrainPres', 'minTrain05', 'minTrain10') %in% metrics) {
				
				# thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
				thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm)
			
				# max Se + Sp
				if ('msss' %in% metrics) {

					thold <- thresholds[['msss']]

					trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					this <- data.frame(
						msssThold=thold,
						msssSeTrain=trainEval[['sensitivity']],
						msssSeTest=testEval[['sensitivity']],
						msssSpTrain=trainEval[['specificity']],
						msssSpTest=testEval[['specificity']]
					)
					
					thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
					
				}
			
				# equalize Se and Sp
				if ('mdss' %in% metrics) {

					thold <- thresholds[['mdss']]

					trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					this <- data.frame(
						mdssThold=thold,
						mdssSeTrain=trainEval[['sensitivity']],
						mdssSeTest=testEval[['sensitivity']],
						mdssSpTrain=trainEval[['specificity']],
						mdssSpTest=testEval[['specificity']]
					)
					
					thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
					
				}
			
				# minimum training presence
				if ('minTrainPres' %in% metrics) {

					thold <- thresholds[['minPres']]

					trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					this <- data.frame(
						minTrainPresThold=thold,
						minTrainPresSeTrain=trainEval[['sensitivity']],
						minTrainPresSeTest=testEval[['sensitivity']],
						minTrainPresSpTrain=trainEval[['specificity']],
						minTrainPresSpTest=testEval[['specificity']]
					)
					
					thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
					
				}
			
				# training sensitivity = 0.95
				if ('trainPres05' %in% metrics) {

					thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, at='sensitivity', sensitivity=0.95)
					thold <- thresholds[['sensitivity']]

					trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					this <- data.frame(
						trainPres05Thold=thold,
						trainPres05SeTrain=trainEval[['sensitivity']],
						trainPres05SeTest=testEval[['sensitivity']],
						trainPres05SpTrain=trainEval[['specificity']],
						trainPres05SpTest=testEval[['specificity']]
					)
					
					thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
						
				}
			
				# training sensitivity = 0.9
				if ('trainPres10' %in% metrics) {

					thresholds <- thresholdWeighted(pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, at='sensitivity', sensitivity=0.9)
					thold <- thresholds[['sensitivity']]

					trainEval <- sensSpecWeighted(thold, pres=predToTrainPres, bg=predToTrainContrast, presWeight=trainPresWeights, bgWeight=trainContrastWeights, na.rm=na.rm, ...)
					testEval <- sensSpecWeighted(thold, pres=predToTestPres, bg=predToTestContrast, presWeight=testPresWeights, bgWeight=testContrastWeights, na.rm=na.rm, ...)
					
					this <- data.frame(
						trainPres10Thold=thold,
						trainPres10SeTrain=trainEval[['sensitivity']],
						trainPres10SeTest=testEval[['sensitivity']],
						trainPres10SpTrain=trainEval[['specificity']],
						trainPres10SpTest=testEval[['specificity']]
					)
					
					thisTuning <- omnibus::insertCol(this, into=thisTuning, at=ncol(thisTuning), before=FALSE)
						
				}
				
			} # threshold-specific metrics
			
		if (verbose > 1) print(thisTuning)
		tuning <- rbind(tuning, thisTuning)
		
	} # next fold
	
	# return
	if ('models' %in% out & 'tuning' %in% out) {
		out <- list()
		out$models <- models
		out$tuning <- tuning
	} else if ('tuning' %in% out) {
		tuning
	} else {
		model
	}
	
}
