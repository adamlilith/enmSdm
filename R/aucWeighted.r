#' Weighted AUC
#'
#' This function calculates the area under the receiver-operator characteristic curve (AUC) following Mason & Graham (2002). Each case (presence/non-presence) can be assigned a weight, if desired.
#' @param pres Vector of predictions for "positive" cases (e.g., predictions at presence sites).
#' @param contrast Vector of predictions for "negative" cases (e.g., predictions at absence/background sites).
#' @param presWeight Weights of positive cases. The default is to assign each positive case a weight of 1.
#' @param contrastWeight Weights of contrast cases. The default is to assign each case a weight of 1.
#' @param na.rm Logical. If \code{TRUE} then remove any positive cases and associated weights and contrast predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (unused).
#' @return A Numeric value.
#'
#' @references Mason, S.J. and N.E. Graham.  2002.  Areas beneath the relative operating characteristics (ROC) and relative operating levels (ROL) curves: Statistical significance and interpretation.  \emph{Quarterly Journal of the Royal Meteorological Society} 128:2145-2166. \doi{10.1256/003590002320603584}
#' 
#' @seealso \code{\link[dismo]{evaluate}}, \code{\link{aucMultiWeighted}}
#' 
#' @examples
#' pres <- seq(0.5, 1, by=0.1)
#' contrast <- seq(0, 1, by=0.01)
#'
#' # unweighted
#' aucWeighted(pres, contrast)
#'
#' # weighted (weight presences with low predictions more)
#' presWeight <- c(1, 1, 1, 0.5, 0.5, 0.5)
#' aucWeighted(pres, contrast, presWeight=presWeight)
#'
#' # weighted (weight presences with high predictions more)
#' presWeight <- c(0.5, 0.5, 0.5, 1, 1, 1)
#' aucWeighted(pres, contrast, presWeight=presWeight)
#'
#' # weight presences and absences
#' contrastWeight <- sqrt(contrast)
#' aucWeighted(pres, contrast, presWeight=presWeight, contrastWeight=contrastWeight)
#' @export

aucWeighted <- function(
	pres,
	contrast,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	na.rm = FALSE,
	bg = NULL,
	bgWeight = NULL,
	...
) {

	if (missing(contrast) & !is.null(bg)) contrast <- bg
	if (missing(contrastWeight) & !is.null(bgWeight)) contrast <- bgWeight

	# if all NAs
	if (all(is.na(pres)) | all(is.na(contrast)) | all(is.na(presWeight)) | all(is.na(contrastWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of positive case predictions and positive case weights ("pres" and "presWeight").')
	if (length(contrastWeight) != length(contrast)) stop('You must have the same number of negative case predictions and negative case weights ("contrast" and "contrastWeight").')
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(contrast, contrastWeight)
		contrast <- cleanedBg[[1]]
		contrastWeight <- cleanedBg[[2]]

	}

	# calculate AUC
	numer <- 0
	for (i in seq_along(pres)) {
		numer <- numer + presWeight[i] * sum(contrastWeight * (pres[i] > contrast))
		numer <- numer + 0.5 * presWeight[i] * sum(contrastWeight * (pres[i] == contrast))
	}

	denom <- sum(presWeight) * sum(contrastWeight)

	auc <- numer / denom
	auc

}
