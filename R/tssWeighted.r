#' Weighted True Skill Statistic (TSS)
#'
#' This function calculates the True Skill Statistic (TSS).
#' @param pres Numeric vector. Predicted values at test presences
#' @param bg Numeric vector. Predicted values at background/absence sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param bgWeight Numeric vector same length as \code{bg}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param thresholds Numeric vector. Thresholds at which to calculate the sum of sensitivity and specificity. The default evaluates all values from 0 to 1 in steps of 0.01.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @return Numeric value.
#' @details This function calculates the maximum value of the True Skill Statistic (i.e., across all thresholds, the values that maximizes sensitivity plus specificity).
#' @references  See Allouche, O., Tsoar, A., and Kadmon, R. 2006. Assessing the accuracy of species distribution models: Prevalence, kappa and the true skill statistic (TSS). \emph{Journal of Applied Ecology} 43:1223-1232.
#' @seealso \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce}}, \code{link[enmSdm]{contBoyce2x}}
#' @examples
#' set.seed(123)
#' 
#' # set of bad and good predictions at presences
#' bad <- runif(100)^2
#' good <- runif(100)^0.1
#' hist(good, breaks=seq(0, 1, by=0.1), border='green', main='Presences')
#' hist(bad, breaks=seq(0, 1, by=0.1), border='red', add=TRUE)
#' pres <- c(bad, good)
#' bg <- runif(1000)
#' tssWeighted(pres, bg)
#' 
#' # upweight bad predictions
#' presWeight <- c(rep(1, 100), rep(0.1, 100))
#' tssWeighted(pres, bg, presWeight=presWeight)
#' 
#' # upweight good predictions
#' presWeight <- c(rep(0.1, 100), rep(1, 100))
#' tssWeighted(pres, bg, presWeight=presWeight)
#' 
#' e <- dismo::evaluate(pres, bg)
#' max(e@TPR + e@TNR) - 1
#' 
#' # why different from dismo's evaluate() function?
#' # because uses different thresholds based on values of
#' # presence and absence predictions
#' head(e@t)
#' @export

tssWeighted <- function(
	pres,
	bg,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	thresholds = seq(0, 1, by=0.01),
	na.rm = FALSE
) {

	# if all NAs
	if (all(is.na(pres)) | all(is.na(bg)) | all(is.na(presWeight)) | all(is.na(bgWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
	if (length(bgWeight) != length(bg)) stop('You must have the same number of absence/background predictions and absence/background weights ("bg" and "bgWeight").')
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleanedBg[[1]]
		bgWeight <- cleanedBg[[2]]

	}

	# stats
	sumPresWeights <- sum(presWeight)
	sumBgWeights <- sum(bgWeight)
	
	numPres <- length(pres)
	numBg <- length(bg)
	
	# TSS
	tss <- rep(NA, length(thresholds))
	
	# for each threshold
	for (i in seq_along(thresholds)) {
		
		thisThresh <- thresholds[i]
	
		# which presences/contrast sites are correctly predicted at this threshold
		whichCorrectPres <- which(pres >= thisThresh)
		whichCorrectBg <- which(bg < thisThresh)
		
		numCorrectPres <- length(whichCorrectPres)
		numCorrectBg <- length(whichCorrectBg)
		
		anyCorrectPres <- (numCorrectPres > 0)
		anyCorrectBg <- (numCorrectBg > 0)
		
		# weights of correctly predicted predictions
		correctPresWeights <- if (anyCorrectPres) {
			sum(presWeight[whichCorrectPres])
		} else {
			0
		}
		
		correctBgWeights <- if (anyCorrectBg) {
			sum(bgWeight[whichCorrectBg])
		} else {
			0
		}
		
		# true positive/negative rates
		tpr <- correctPresWeights / sumPresWeights
		tnr <- correctBgWeights / sumBgWeights
		tss[i] <- tpr + tnr - 1
	
	}
	
	thresholdMaxTss <- thresholds[which.max(tss)]
	tss <- max(tss)
	attr(tss, 'thresholdMaxTss') <- thresholdMaxTss
	tss

}
