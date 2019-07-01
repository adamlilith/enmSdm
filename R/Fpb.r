#' Fpb measure of model accuracy with weighting
#'
#' This function calculates Li and Guo's Fpb measure of model accuracy for test presences and randomly located background sites.
#' @param pres Numeric vector. Predictions at presence sites.
#' @param bg Numeric vector. Predictions at absence/background sites.
#' @param thresholds Numeric value or numeric vector within the range [0, 1]. Threshold value(s) at which to calculate Fpb.
#' @param presWeight Numeric same length as \code{pres}. Weights for presence predictions. The default is to assign each presence a weight of 1.
#' @param bgWeight Numeric same length as \code{bg}. Weights for background predictions. The default is to assign each presence a weight of 1.
#' @param na.rm Logical. If \code{TRUE} remove any \code{NA}s in predictions at presences and/or absences before calculation. If \code{NA}s occur but are not removed then the output will be \code{NA}.
#' @param tr Same as \code{thresholds}. Deprecated, but included for backwards compatibility.
#' @param ... Other arguments (unused).
#' @return Numeric.
#' @references Li, W. and Guo, Q.  2013.  How to assess the prediction accuracy of species presence-absence models without absence data? \emph{Ecography} 36:788-799.
#' @seealso \code{\link{aucWeighted}}, \code{\link{contBoyce}}, \code{\link[dismo]{evaluate}}
#' @examples
#' pres <- seq(0.5, 1, by=0.1)
#' bg <- seq(0, 1, by=0.01)
#' thresholds <- seq(0, 1, by=0.1)
#'
#' # unweighted
#' f1 <- fpb(pres, bg, thresholds)
#'
#' # weighted (weight presences with low predictions more)
#' presWeight <- c(1, 1, 1, 0.5, 0.5, 0.5)
#' f2 <- fpb(pres, bg, thresholds, presWeight=presWeight)
#'
#' # weighted (weight presences with high predictions more)
#' presWeight <- c(0.5, 0.5, 0.5, 1, 1, 1)
#' f3 <- fpb(pres, bg, thresholds, presWeight=presWeight)
#'
#' # weight presences and absences
#' bgWeight <- sqrt(bg)
#' f4 <- fpb(pres, bg, thresholds, presWeight=presWeight, bgWeight=bgWeight)
#'
#' plot(thresholds, f1, type='b', xlab='Threshold', ylab='fpb', ylim=c(0, 1.5))
#' points(thresholds, f2, type='b', pch=2)
#' points(thresholds, f3, type='b', pch=3)
#' points(thresholds, f4, type='b', pch=4)
#' legend('topright', inset=0.01,
#'    legend=c('no weights', 'high presences upweighted',
#'       'low presences upweighted', 'pres and bg weighted'),
#'    pch=1:4)
#' @export

fpb <- function(
	pres,
	bg,
	thresholds = seq(0, 1, by=0.01),
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	na.rm = FALSE,
	tr = NULL,
	...
) {

	# if all NAs
	if (all(is.na(pres)) | all(is.na(bg)) | all(is.na(presWeight)) | all(is.na(bgWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
	if (length(bgWeight) != length(bg)) stop('You must have the same number of absence/background predictions and absence/background weights ("bg" and "bgWeight").')
	
	# for backwards compatibility
	if (!is.null(tr)) thresholds <- tr
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleanedBg[[1]]
		bgWeight <- cleanedBg[[2]]

		thresholds <- na.omit(thresholds)
		
	}

	out <- sapply(X=thresholds, FUN=fpbOneThreshold, pres=pres, bg=bg, presWeight=presWeight, bgWeight=bgWeight)
	out

}
