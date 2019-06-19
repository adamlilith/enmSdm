#' Calculate the Fpb measure of model accuracy for a single threshold
#'
#' @param threshold Numeric value or numeric vector within the range [0, 1]. Threshold value(s) at which to calculate Fpb.
#' @param pres Numeric vector. Predictions at presences.
#' @param bg Numeric vector. Predictions at background sites.
#' @param presWeight Numeric same length as \code{pres}. Weights for presence predictions.
#' @param bgWeight Numeric same length as \code{bg}. Weights for background predictions.
#' @param ttr Same as \code{threshold} but deprecated. Included for backwards compatability.
#' @return Numeric.
#' @export

fpbOneThreshold <- function(
	threshold,
	pres,
	bg,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	tr = NULL
) {

	if (!is.null(tr)) threshold <- tr

	out <- (2 * (sum(presWeight * (pres >= tr)) / sum(presWeight))) / (1 + (sum(bgWeight * (bg >= tr)) / sum(bgWeight)))
	out

}

