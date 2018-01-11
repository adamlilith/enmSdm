#' Calculate AUC with optional weighting
#'
#' This function calculates the area under the receiver-operator characteristic curve (AUC) following Mason, S.J. and N.E. Graham.  2002.  Areas beneath the relative operating characteristics (ROC) and relative operating levels (ROL) curves: Statistical significance and interpretation.  \emph{Quarterly Journal of the Royal Meteorological Society} 128:2145-2166. Positives and negatives values can be given weights.
#' @param pres Predictions at presence sites.
#' @param bg Predictions at absence sites.
#' @param presWeight Weights of presences.
#' @param bgWeight Weights of absences.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @return Numeric value.
#' @seealso \code{\link{Fpb}}, \code{\link{contBoyce}}, \code{\link[dismo]{evalulate}}
#' @examples
#' pres <- seq(0.5, 1, by=0.1)
#' bg <- seq(0, 1, by=0.01)
#'
#' # unweighted
#' aucWeighted(pres, bg)
#'
#' # weighted (weight presences with low predictions more)
#' presWeight <- c(1, 1, 1, 0.5, 0.5, 0.5)
#' aucWeighted(pres, bg, presWeight=presWeight)
#'
#' # weighted (weight presences with high predictions more)
#' presWeight <- c(0.5, 0.5, 0.5, 1, 1, 1)
#' aucWeighted(pres, bg, presWeight=presWeight)
#'
#' # weight presences and absences
#' bgWeight <- sqrt(bg)
#' aucWeighted(pres, bg, presWeight=presWeight, bgWeight=bgWeight)
#' @export

aucWeighted <- function(
	pres,
	bg,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	na.rm = FALSE
) {

	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleanedBg[[1]]
		bgWeight <- cleanedBg[[2]]

	}

	# calculate AUC
	numer <- 0
	for (i in seq_along(pres)) {
		numer <- numer + presWeight[i] * sum(bgWeight * (pres[i] > bg))
		numer <- numer + 0.5 * presWeight[i] * sum(bgWeight * (pres[i] == bg))
	}

	denom <- sum(presWeight) * sum(bgWeight)

	auc <- numer / denom
	auc

}
