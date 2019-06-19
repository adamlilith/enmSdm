#' Continuous Boyce Index (CBI) with weighting (2X-coverage version)
#'
#' This function calculates the Continuous Boyce Index (CBI), a measure of model accuracy for presence-only test data. Overlapping bins are placed such that any given point (prediction) along [0, 1] is covered by at most 2 bins. See the function \code{\link[enmSdm]{contBoyce}} for a version that allows coverage by 2 or more bins.
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param bg Numeric vector.  Predicted values at absence/background sites.
#' @param numClasses Positive integer. Number of classes into which to divide predictions at background sites. Hirzel et al. suggest using 10 (default).
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param bgWeight Numeric vector same length as \code{bg}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param upweightTails Logical. \code{TRUE} ==> weights of presences and background sites that occur in the first half of the lowest bin or in the second half of the last bin  have their weights multiplied by 2.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param autoWindow Logical. If TRUE then calculate bin boundaries starting at minimum predicted value and ending at maximum predicted value (default). If FALSE calculate bin boundaries starting at 0 and ending at 1.
#' @param method Character. Type of correlation to calculate. The default is \code{'spearman'}, the Spearman rank correlation coefficient used by Boyce et al. (2002) and Hirzel et al. (2006), which is the traditional CBI. In contrast, \code{'pearson'} or \code{'kendall'} can be used instead.  See [stats::cor()] for more details.
#' @param graph Logical. If TRUE then plot P vs E and P/E versus bin.
#' @return Numeric value.
#' @details The CBI is the Spearman rank correlation coefficient between the proportion of sites in each prediction class and the expected proportion of predictions in each prediction class based on the proportion of the landscape that is in that class. Values >0 indicate the model's output is positively correlated with the true probability of presence.  Values <0 indicate it is negatively correlated with the true probabilty of presence.
#' @references Boyce, M.S., Vernier, P.R., Nielsen, S.E., and Schmiegelow, F.K.A.  2002.  Evaluating resource selection functions.  \emph{Ecological Modeling} 157:281-300)
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randon, C., and Guisan, A.  2006.  Evaluating the ability of habitat suitability models to predict species presences.  \emph{Ecological Modeling} 199:142-152.
#' @seealso \code{\link[enmSdm]{contBoyce}}
#' @examples
#' set.seed(57)
#' pres <- seq(0.5, 1, length.out=100)
#' bg <- runif(1000)
#' contBoyce2x(pres, bg)
#' contBoyce(pres, bg)
#' presWeight <- c(rep(1, 50), rep(0.5, 50))
#' contBoyce2x(pres, bg, presWeight=presWeight)
#' contBoyce(pres, bg, presWeight=presWeight)
#' @export
contBoyce2x <- function(
	pres,
	bg,
	numClasses = 10,
	presWeight = rep(1, length(pres)),
	bgWeight = rep(1, length(bg)),
	upweightTails = TRUE,
	na.rm = FALSE,
	autoWindow = TRUE,
	method = 'spearman',
	graph = FALSE
) {

	# if all NAs
	if (all(is.na(pres)) | all(is.na(bg)) | all(is.na(presWeight)) | all(is.na(bgWeight))) return(NA)

	# catch errors
	if (length(presWeight) != length(pres)) stop('You must have the same number of presence predictions and presence weights ("pres" and "presWeight").')
	if (length(bgWeight) != length(bg)) stop('You must have the same number of absence/background predictions and absence/background weights ("bg" and "bgWeight").')
	
	# remove NAs
	if (na.rm) {

		cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		pres <- cleanedPres[[1]]
		presWeight <- cleanedPres[[2]]

		cleanedBg <- omnibus::naOmitMulti(bg, bgWeight)
		bg <- cleanedBg[[1]]
		bgWeight <- cleanedBg[[2]]

	}

	## range of predicted values across which to calculate bins
	rangeOfPredVals <- if (autoWindow) {
		(max(c(pres, bg), na.rm=T) - min(c(pres, bg), na.rm=T))
	} else {
		1
	}

	# right hand side of each class (assumes max value is >0)
	lowestBound <- if (autoWindow) { min(c(pres, bg), na.rm=T) } else { 0 }
	highestBound <- if (autoWindow) { max(c(pres, bg), na.rm=T) + .Machine$double.eps } else { 1 + omnibus::eps() }
	classBound <- seq(from=lowestBound, to=highestBound, length.out=numClasses + 2)

	# up-weight tails
	if (upweightTails) {

		# lowest values
		presWeight[pres < classBound[2]] <- 2 * presWeight[pres < classBound[2]]
		bgWeight[bg < classBound[2]] <- 2 * bgWeight[bg < classBound[2]]

		# highest values
		presWeight[pres >= classBound[numClasses - 1]] <- 2 * presWeight[pres >= classBound[numClasses - 1]]
		bgWeight[bg >= classBound[numClasses - 1]] <- 2 * bgWeight[bg >= classBound[numClasses - 1]]

	}

	##########
	## MAIN ##
	##########

	## initiate variables to store predicted/expected (P/E) values
	freqPres <- freqBg <- numeric(numClasses)

	### tally proportion of test presences/background sites in each class
	for (countClass in 1:numClasses) {

		# number of presence predictions in this class
		presInBin <- pres >= classBound[countClass] & pres < classBound[countClass + 2]
		presInBin <- presInBin * presWeight
		freqPres[countClass] <- sum(presInBin)

		# number of background predictions in this class
		bgInBin <- bg >= classBound[countClass] & bg < classBound[countClass + 2]
		bgInBin <- bgInBin * bgWeight
		freqBg[countClass] <- sum(bgInBin)

	} # next predicted value class

	# add small number to each background frequency class to avoid division by 0 ("small" number is "half" a presence)
	freqBg <- freqBg + min(0.5 * c(presWeight[presWeight > 0], bgWeight[bgWeight > 0]))

	P <- freqPres / sum(presWeight)
	E <- freqBg / sum(bgWeight)
	PE <- P / E

	# plot
	if (graph) {
		par(mfrow=c(1, 2))
		plot(E, P, col='white', xlab='Expected', ylab='Predicted', main='P/E\n(letters from lowest to highest class)')
		text(E, P, labels=letters[1:numClasses], col='red')
		plot(1:numClasses, PE, type='l', xlab='Bin', ylab='P/E Ratio', main='CBI\n(letters from lowest to highest class)')
		text(1:numClasses, PE, labels=letters[1:numClasses], col='red')
	}

	# calculate continuous Boyce index (cbi)
	cbi <- stats::cor(x=1:length(PE), y=P/E, method=method)

	#####################
	## post-processing ##
	#####################

	cbi

}
