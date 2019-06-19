#' Weighted thresholds for predictions
#'
#' This function is similar to the \code{\link[dismo]{threshold}} function in the \pkg{dismo} package, which calculates thresholds to create binary predictions from continuous values. However, unlike that function is allows the user to specify weights for presences and absence/background predictions. The output will thus be the threshold that best matches the specified criterion taking into account the relative weights of the inpyt values.
#' @param pres Numeric vector. Predicted values at test presences
#' @param bg Numeric vector. Predicted values at background/absence sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param bgWeight Numeric vector same length as \code{bg}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param at Character or character vector, name(s) of threshold(s) to calculate. The default is to calculate them all.
#' \itemize{
#' 		\item \code{'msss'}: Threshold that the maximums sum of sensitivity and specificity.
#' 		\item \code{'mdss'}: Threshold that minimizes the difference between sensitivity and specificity.
#' 		\item \code{'minPres'}: Minimum prediction across presences. This threshold is not weighted.
#' 		\item \code{'prevalence'}: Prevalence of presences (sum(presence weights) / sum(presence weights + background weights))'
#' 		\item \code{'sensitivity'}: Threshold that most closely returns the sensitivity specified by \code{sensitivity}.
#' }
#' @param sensitivity Value of specificity to match (used only if \code{at} contains \code{'sensitivity'}).
#' @param thresholds Numeric vector. Thresholds at which to calculate the sum of sensitivity and specificity. The default evaluates all values from 0 to 1 in steps of 0.01.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @return Named numeric vector.
#' @references Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation} 24:38-49.
#' @seealso \code{\link[dismo]{threshold}}, \code{\link[dismo]{evaluate}}
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
#' thresholdWeighted(pres, bg)
#' 
#' # upweight bad predictions
#' presWeight <- c(rep(1, 100), rep(0.1, 100))
#' thresholdWeighted(pres, bg, presWeight=presWeight)
#' 
#' # upweight good predictions
#' presWeight <- c(rep(0.1, 100), rep(1, 100))
#' thresholdWeighted(pres, bg, presWeight=presWeight)
#' @export

thresholdWeighted <- function(
	pres,
	bg,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	at = c('msss', 'mdss', 'minPres', 'prevalence', 'sensitivity'),
	sensitivity = 0.9,
	thresholds = seq(0, 1, by=0.01),
	na.rm = FALSE
) {

	# if all NAs
	if (all(is.na(pres)) | all(is.na(bg))) return(NA)

	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleanedBg[[1]]
		bgWeight <- cleanedBg[[2]]

	}

	### calculate thresholds
	########################

	# used for several thresholds
	sumPresWeights <- sum(presWeight)
	sumBgWeights <- sum(bgWeight)
		
	# output
	out <- numeric()
	
	### maximum Se + Sp and minimum difference between Se and Sp
	
	if (any(c('msss', 'mdss') %in% at)) {
	
		# first, calculate TPR and TNR across all possible thresholds
	
		# at
		numPres <- length(pres)
		numBg <- length(bg)
		
		# TPR and TNR
		tpr <- tnr <- rep(NA, length(thresholds))
		
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
			tpr[i] <- correctPresWeights / sumPresWeights
			tnr[i] <- correctBgWeights / sumBgWeights
		
		}
			
		# second, calculate threshold
		if ('msss' %in% at) {
			x <- thresholds[which.max(tpr + tnr)]
			out <- c(out, x)
			names(out)[length(out)] <- 'msss'
		}

		if ('mdss' %in% at) {
			x <- thresholds[which.min(abs(tpr - tnr))]
			out <- c(out, x)
			names(out)[length(out)] <- 'mdss'
		}
	
	}

	### minimum presence prediction (unweighted)
	
	if ('minPres' %in% at) {
	
		x <- min(pres)
		out <- c(out, x)
		names(out)[length(out)] <- 'minPres'
		
	}

	### prevalence (weighted)
	
	if ('prevalence' %in% at) {
	
		x <- sumPresWeights / (sumPresWeights + sumBgWeights)
		if (x > 1) warning('Prevalence threshold is >1. The sum of presence weights may be larger than the sum of absence/background weights.')
		out <- c(out, x)
		names(out)[length(out)] <- 'prevalence'
		
	}

	### fixed sensitivity (weighted)
	
	if ('sensitivity' %in% at) {
	
		presOrder <- order(pres)
		pres <- pres[presOrder]
		presWeight <- presWeight[presOrder]
		sens <- 1 - cumsum(presWeight) / sumPresWeights
		x <- pres[which.min(abs(sensitivity - sens))]
		out <- c(out, x)
		names(out)[length(out)] <- 'sensitivity'
		
	}

	out
	
}
