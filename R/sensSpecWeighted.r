#' Thresholded evaluation statistics
#'
#' This function calculates a series of evaluation statistics based on a threshold or thresholds used to convert continuous predictions to binary predictions.
#' @param thresholds Numeric or numeric vector. Threshold(s) at which to calculate sensitivity and specificity.
#' @param pres Numeric vector. Predicted values at test presences
#' @param bg Numeric vector. Predicted values at background/absence sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param bgWeight Numeric vector same length as \code{bg}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @return 8-column matrix with the following named columns. \emph{a} = weight of presences >= threshold, \emph{b} = weight of backgrounds >= threshold, \emph{c} = weight of presences < threshold, \emph{d} = weight of backgrounds < threshold, and \emph{N} = sum of presence and background weights.
#' \itemize{
#' 		\item \code{'threshold'}: Threshold
#' 		\item \code{'sensitivity'}: Sensitivity (\emph{a} / (\emph{a} + \emph{c}))
#' 		\item \code{'specificity'}: Specificity (\emph{d} / (\emph{d} + \emph{b}))
#' 		\item \code{'ccr'}: Correct classification rate ((\emph{a} + \emph{d}) / \emph{N})
#' 		\item \code{'ppp'}: Positive predictive power (\emph{a} / (\emph{a} + \emph{b}))
#' 		\item \code{'npp'}: Negative predictive power (\emph{d} / (\emph{c} + \emph{d}))
#' 		\item \code{'mr'}: Misclassification rate ((\emph{b} + \emph{c}) / \emph{N})
#' }
#' @references Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation} 24:38-49.
#' @seealso \code{\link[dismo]{threshold}}, \code{\link[enmSdm]{thresholdWeighted}}, \code{\link[dismo]{evaluate}}
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
#' thresholds <- c(0.1, 0.5, 0.9)
#' sensSpecWeighted(thresholds, pres, bg)
#' 
#' # upweight bad predictions
#' presWeight <- c(rep(1, 100), rep(0.1, 100))
#' sensSpecWeighted(thresholds, pres, bg, presWeight=presWeight)
#' 
#' # upweight good predictions
#' presWeight <- c(rep(0.1, 100), rep(1, 100))
#' sensSpecWeighted(thresholds, pres, bg, presWeight=presWeight)
#' @export

sensSpecWeighted <- function(
	thresholds,
	pres,
	bg,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
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

	# stats
	sumPresWeights <- sum(presWeight)
	sumBgWeights <- sum(bgWeight)
	sumWeights <- sumPresWeights + sumBgWeights
	
	# output
	out <- matrix(NA, ncol=7, nrow=length(thresholds))
	stats <- c('threshold', 'sensitivity', 'specificity', 'ccr', 'ppp', 'npp', 'mr')
	
	for (i in seq_along(thresholds)) {
	
		thisThold <- thresholds[i]
		
		sumWeightCorrectPres <- sum(presWeight[pres >= thisThold])
		sumWeightCorrectBg <- sum(bgWeight[bg < thisThold])
		
		sumWeightIncorrectPres <- sum(presWeight[pres < thisThold])
		sumWeightIncorrectBg <- sum(bgWeight[bg >= thisThold])
		
		sens <-  sumWeightCorrectPres / sumPresWeights
		spec <- sumWeightCorrectBg / sumBgWeights
		ccr <- (sumWeightCorrectPres + sumWeightCorrectBg) / sumWeights
		ppp <- sumWeightCorrectPres / (sumWeightCorrectPres + sumWeightIncorrectBg)
		npp <- sumWeightCorrectBg / (sumWeightIncorrectPres + sumWeightCorrectBg)
		mr <- (sumWeightIncorrectPres + sumWeightIncorrectBg) / sumWeights
		
		out[i, ] <- c(thisThold, sens, spec, ccr, ppp, npp, mr)
	
	}
	
	if (nrow(out) == 1) {
		out <- c(out)
		names(out) <- stats
	} else {
		colnames(out) <- stats
	}
	
	out
	
}
