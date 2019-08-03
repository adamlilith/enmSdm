#' Odds ratio skill score (ORSS)
#'
#' This function calculates the odds ratio skill score (ORSS) using a threshold applied to continuous data to demarcate "presence" from "contrast".
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param contrast Numeric vector. Predicted values at absence/background sites.
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param thresholds Numeric vector, Values at which to threshold predictions for calculation of ORSS.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and absence/background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param ... Other arguments (unused).
#' @return Numeric value.
#' @references Wunderlich, R.F., Lin, Y-P., Anthony, J., and Petway, J.R.  2019.  Two alternative evaluation metrics to replace the true skill statistic in the assessment of species distribution models.  Nature Conservation 35:97-116.
#' @seealso \code{\link[stats]{cor}}, \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce}}, \code{link[enmSdm]{contBoyce2x}}, \code{link[enmSdm]{sediWeighted}}, \code{link[enmSdm]{thresholdWeighted}}, \code{link[enmSdm]{thresholdStats}}
#' @examples
#' set.seed(123)
#' pres <- sqrt(runif(100))
#' contrast <- runif(10000)
#' hist(contrast, col='gray', xlim=c(0, 1), breaks=20)
#' hist(pres, col='green', breaks=20, add=TRUE)
#' max(orssWeighted(pres, contrast), na.rm=TRUE)
#' presWeight <- c(rep(1, 50), rep(0.5, 50))
#' max(orssWeighted(pres, contrast, presWeight=presWeight), na.rm=TRUE)
#' @export

orssWeighted <- function(
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	thresholds = seq(0, 1, by=0.01),
	na.rm = FALSE,
	bg = NULL,
	...
) {

	if (missing(contrast) & !is.null(bg)) contrast <- bg

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
	
	# ORSS, true positive rate, true negative rate, false negative rate
	orss <- tpr <- tnr <- fnr <- rep(NA, length(thresholds))
	
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
	
		# ORSS
		orss[i] <- (correctPresWeights * correctContrastWeights - incorrectPresWeights * incorrectContrastWeights) / (correctPresWeights * correctContrastWeights + incorrectPresWeights * incorrectContrastWeights)
		
	} # next threshold
	
	orss
	
}
