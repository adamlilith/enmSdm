#' Symmetric extremal dependence index (SEDI)
#'
#' This function calculates the symmetric extremal dependence index (SEDI) using a threshold applied to continuous data to demarcate "presence" from "contrast".
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param contrast Numeric vector. Predicted values at absence/background sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param thresholds Numeric vector, Values at which to threshold predictions for calculation of SEDI.
#' @param delta Positive numeric >0 in the range [0, 1] and usually very small. This value is used only if calculating the SEDI threshold when any true positive rate or false negative rate is 0 or the false negative rate is 1. Since SEDI uses log(x) and log(1 - x), values of 0 and 1 will produce \code{NA}s. To obviate this, a small amount can be added to rates that equal 0 and subtracted from rates that equal 1.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and absence/background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (unused).
#' @return Numeric value.
#' @references Wunderlich, R.F., Lin, Y-P., Anthony, J., and Petway, J.R.  2019.  Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models.  Nature Conservation 35:97-116.
#' @seealso \code{\link[stats]{cor}}, \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce}}, \code{link[enmSdm]{contBoyce2x}}, \code{link[enmSdm]{orssWeighted}}, \code{link[enmSdm]{thresholdWeighted}}, \code{link[enmSdm]{thresholdStats}}
#' @examples
#' set.seed(123)
#' pres <- sqrt(runif(100))
#' contrast <- runif(10000)
#' hist(contrast, col='gray', xlim=c(0, 1), breaks=20)
#' hist(pres, col='green', breaks=20, add=TRUE)
#' max(sediWeighted(pres, contrast), na.rm=TRUE)
#' # SEDI is fairly insensitive to weighting (as per Wunderlich et al. 2019)
#' presWeight <- c(rep(1, 50), rep(1000, 50))
#' max(sediWeighted(pres, contrast, presWeight=presWeight), na.rm=TRUE)
#' @export

sediWeighted <- function(
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	thresholds = seq(0, 1, by=0.01),
	delta = 0.001,
	na.rm = FALSE,
	bg = NULL,
	bgWeight = NULL,
	...
) {

	if (missing(contrast) & !is.null(bg)) contrast <- bg
	if (missing(contrastWeight) & !is.null(bgWeight)) contrastWeight <- bgWeight

	# if all NAs
	if (all(is.na(pres)) | all(is.na(contrast)) | all(is.na(presWeight)) | all(is.na(contrastWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
	if (length(contrastWeight) != length(contrast)) stop('You must have the same number of absence/background predictions and absence/background weights ("contrast" and "contrastWeight").')
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
		contrast <- cleanedContrast[[1]]
		contrastWeight <- cleanedContrast[[2]]

	}

	sumPresWeights <- sum(presWeight)
	sumContrastWeights <- sum(contrastWeight)
	
	# SEDI, true positive rate, true negative rate, false negative rate
	sedi <- tpr <- tnr <- fnr <- rep(NA, length(thresholds))
	
	for (i in seq_along(thresholds)) {
	
		thisThresh <- thresholds[i]
	
		# which presences/contrast sites are CORRECTLY predicted at this threshold
		whichCorrectPres <- which(pres >= thisThresh)
		whichCorrectContrast <- which(contrast < thisThresh)
		
		numCorrectPres <- length(whichCorrectPres)
		numCorrectContrast <- length(whichCorrectContrast)
		
		anyCorrectPres <- (numCorrectPres > 0)
		anyCorrectContrast <- (numCorrectContrast > 0)
		
		# which presences/contrast sites are INCORRECTLY predicted at this threshold
		whichIncorrectPres <- which(pres < thisThresh)
		whichIncorrectContrast <- which(contrast >= thisThresh)
		
		numIncorrectPres <- length(whichIncorrectPres)
		numIncorrectContrast <- length(whichIncorrectContrast)
		
		anyIncorrectPres <- (numIncorrectPres > 0)
		anyIncorrectContrast <- (numIncorrectContrast > 0)
		
		# weights of CORRECTLY predicted predictions
		correctPresWeights <- if (anyCorrectPres) {
			sum(presWeight[whichCorrectPres])
		} else {
			0
		}
		
		correctContrastWeights <- if (anyCorrectContrast) {
			sum(contrastWeight[whichCorrectContrast])
		} else {
			0
		}
		
		# weights of INCORRECTLY predicted predictions
		incorrectPresWeights <- if (anyIncorrectPres) {
			sum(presWeight[whichIncorrectPres])
		} else {
			0
		}
		
		incorrectContrastWeights <- if (anyIncorrectContrast) {
			sum(contrastWeight[whichIncorrectContrast])
		} else {
			0
		}
		
		# true positive/negative rates
		tpr[i] <- correctPresWeights / sumPresWeights
		tnr[i] <- correctContrastWeights / sumContrastWeights
	
		# false positive/negative rates
		fnr[i] <- incorrectContrastWeights / sumContrastWeights
		
	}
	
	# SEDI
	if (any(tpr == 0)) tpr[tpr == 0] <- delta
	if (any(fnr == 0)) fnr[fnr == 0] <- delta
	if (any(tpr == 1)) tpr[tpr == 1] <- 1 - delta
	if (any(fnr == 1)) fnr[fnr == 1] <- 1 - delta
	
	numer <- log10(fnr) - log10(tpr) - log10(1 - fnr) + log10(1 - tpr)
	denom <- log10(fnr) + log10(tpr) + log10(1 - fnr) + log10(1 - tpr)
	sedi <- numer / denom
		
	sedi
	
}
