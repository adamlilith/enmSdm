#' Thresholded evaluation statistics
#'
#' This function calculates a series of evaluation statistics based on a threshold or thresholds used to convert continuous predictions to binary predictions.
#' @param thresholds Numeric or numeric vector. Threshold(s) at which to calculate sensitivity and specificity.
#' @param pres Numeric vector. Predicted values at test presences
#' @param contrast Numeric vector. Predicted values at background/absence sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param delta Positive numeric >0 in the range [0, 1] and usually very small. This value is used only if calculating the SEDI threshold when any true positive rate or false negative rate is 0 or the false negative rate is 1. Since SEDI uses log(x) and log(1 - x), values of 0 and 1 will produce \code{NA}s. To obviate this, a small amount can be added to rates that equal 0 and subtracted from rates that equal 1.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (unused).
#' @return 8-column matrix with the following named columns. \emph{a} = weight of presences >= threshold, \emph{b} = weight of backgrounds >= threshold, \emph{c} = weight of presences < threshold, \emph{d} = weight of backgrounds < threshold, and \emph{N} = sum of presence and background weights.
#' \itemize{
#' 		\item \code{'threshold'}: Threshold
#' 		\item \code{'sensitivity'}: Sensitivity (\emph{a} / (\emph{a} + \emph{c}))
#' 		\item \code{'specificity'}: Specificity (\emph{d} / (\emph{d} + \emph{b}))
#' 		\item \code{'ccr'}: Correct classification rate ((\emph{a} + \emph{d}) / \emph{N})
#' 		\item \code{'ppp'}: Positive predictive power (\emph{a} / (\emph{a} + \emph{b}))
#' 		\item \code{'npp'}: Negative predictive power (\emph{d} / (\emph{c} + \emph{d}))
#' 		\item \code{'mr'}: Misclassification rate ((\emph{b} + \emph{c}) / \emph{N})
#' 		\item \code{'orss'}: Threshold that maximizes the odds ratio skill score (Wunderlich et al. 2019).
#' 		\item \code{'sedi'}: Threshold that maximizes the symmetrical extremal dependence index (Wunderlich et al. 2019).
#' }
#' @references Fielding, A.H. and J.F. Bell. 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. \emph{Environmental Conservation} 24:38-49.
#' Wunderlich, R.F., Lin, P-Y., Anthony, J., and Petway, J.R. 2019. Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribition models. Nature Conservation 35:97-116.
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
#' contrast <- runif(1000)
#' thresholds <- c(0.1, 0.5, 0.9)
#' thresholdStats(thresholds, pres, contrast)
#' 
#' # upweight bad predictions
#' presWeight <- c(rep(1, 100), rep(0.1, 100))
#' thresholdStats(thresholds, pres, contrast, presWeight=presWeight)
#' 
#' # upweight good predictions
#' presWeight <- c(rep(0.1, 100), rep(1, 100))
#' thresholdStats(thresholds, pres, contrast, presWeight=presWeight)
#' @export

thresholdStats <- function(
	thresholds,
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	delta = 0.001,
	na.rm = FALSE,
	bg = NULL,
	bgWeight = NULL,
	...
) {

	if (missing(contrast) & !is.null(bg)) contrast <- bg
	if (missing(contrastWeight) & !is.null(bgWeight)) contrastWeight <- bgWeight

	# names of statistics to calculate
	stats <- c('threshold', 'sensitivity', 'specificity', 'ccr', 'ppp', 'npp', 'mr', 'orss', 'sedi')

	# if all NAs
	if (all(is.na(pres)) | all(is.na(contrast))) return(NA)

	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
		contrast <- cleanedContrast[[1]]
		contrastWeight <- cleanedContrast[[2]]

	}

	# stats
	sumPresWeights <- sum(presWeight)
	sumContrastWeights <- sum(contrastWeight)
	sumWeights <- sumPresWeights + sumContrastWeights
	
	# output
	out <- matrix(NA, ncol=length(stats), nrow=length(thresholds))
	
	for (i in seq_along(thresholds)) {
	
		thisThold <- thresholds[i]
		
		sumWeightCorrectPres <- sum(presWeight[pres >= thisThold])
		sumWeightCorrectContrast <- sum(contrastWeight[contrast < thisThold])
		
		sumWeightIncorrectPres <- sum(presWeight[pres < thisThold])
		sumWeightIncorrectContrast <- sum(contrastWeight[contrast >= thisThold])
		
		sens <-  sumWeightCorrectPres / sumPresWeights
		spec <- sumWeightCorrectContrast / sumContrastWeights
		ccr <- (sumWeightCorrectPres + sumWeightCorrectContrast) / sumWeights
		ppp <- sumWeightCorrectPres / (sumWeightCorrectPres + sumWeightIncorrectContrast)
		npp <- sumWeightCorrectContrast / (sumWeightIncorrectPres + sumWeightCorrectContrast)
		mr <- (sumWeightIncorrectPres + sumWeightIncorrectContrast) / sumWeights
		orss <- orssWeighted(thresholds=thisThold, pres=pres, contrast=contrast, presWeight=presWeight, contrastWeight=contrastWeight, na.rm=na.rm)

		sedi <- sediWeighted(thresholds=thisThold, pres=pres, contrast=contrast, presWeight=presWeight, contrastWeight=contrastWeight, delta=delta, na.rm=na.rm)

		### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		### must be in same order at "stats" define at start of function!!!
		out[i, ] <- c(thisThold, sens, spec, ccr, ppp, npp, mr, orss, sedi)
	
	}
	
	if (nrow(out) == 1) {
		out <- c(out)
		names(out) <- stats
	} else {
		colnames(out) <- stats
	}
	
	out
	
}
