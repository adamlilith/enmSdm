#' Calculate the Fpb measure of model accuracy
#'
#' This function calculates Li and Guo's Fpb measure of model accuracy for test presences and randomly located background sites.
#' @param tr Numeric value or numeric vector within the range [0, 1]. Threshold value(s) at which to calculate Fpb.
#' @param pres Numeric vector. Predictions at presences.
#' @param bg Numeric vector. Predictions at background sites.
#' @param presWeight Numeric same length as \code{pres}. Weights for presence predictions.
#' @param bgWeight Numeric same length as \code{bg}. Weights for background predictions.
#' @param na.rm Logical. If TRUE remove any \code{NA}s in predictions at preseneces and/or absences before calculation. If \code{NA}s occur but are not removed then the outout will be \code{NA}.
#' @return Numeric.
#' @references Li, W. and Guo, Q.  2013.  How to assess the prediction accuracy of species presence-absence models without absence data? \emph{Ecography} 36:788-799.
#' @seealso \code{\link{aucWeighted}}, \code{\link{contBoyce}}, \code{\link[dismo]{evaluate}}
#' @examples
#' pres <- seq(0.5, 1, by=0.1)
#' bg <- seq(0, 1, by=0.01)
#' tr <- seq(0, 1, by=0.1)
#'
#' # unweighted
#' f1 <- Fpb(pres, bg, tr)
#'
#' # weighted (weight presences with low predictions more)
#' presWeight <- c(1, 1, 1, 0.5, 0.5, 0.5)
#' f2 <- Fpb(pres, bg, tr, presWeight=presWeight)
#'
#' # weighted (weight presences with high predictions more)
#' presWeight <- c(0.5, 0.5, 0.5, 1, 1, 1)
#' f3 <- Fpb(pres, bg, tr, presWeight=presWeight)
#'
#' # weight presences and absences
#' bgWeight <- sqrt(bg)
#' f4 <- Fpb(pres, bg, tr, presWeight=presWeight, bgWeight=bgWeight)
#'
#' plot(tr, f1, type='b', xlab='Threshold', ylab='Fpb', ylim=c(0, 1.5))
#' points(tr, f2, type='b', pch=2)
#' points(tr, f3, type='b', pch=3)
#' points(tr, f4, type='b', pch=4)
#' legend('topright', inset=0.01,
#'    legend=c('no weights', 'high presences upweighted',
#'       'low presences upweighted', 'pres and bg weighted'),
#'    pch=1:4)
#' @export

Fpb <- function(
	pres,
	bg,
	tr = seq(0, 1, by = 0.1),
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	na.rm = FALSE
) {

	if (na.rm) {
		cleaned <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleaned[[1]]
		presWeight <- cleaned[[2]]

		cleaned <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleaned[[1]]
		bgWeight <- cleaned[[2]]

		tr <- na.omit(tr)
	}

	out <- sapply(X=tr, FUN=FpbOneThreshold, pres=pres, bg=bg, presWeight=presWeight, bgWeight=bgWeight)
	out

}
