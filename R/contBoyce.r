#' Continuous Boyce Index (CBI) with weighting
#'
#' This function calculates the Continuous Boyce Index (CBI), a measure of model accuracy for presence-only test data.
#' @param pres Numeric vector. Predicted values at presence sites.
#' @param contrast Numeric vector. Predicted values at absence/background sites.
#' @param numBins Positive integer. Number of (overlapping) bins into which to divide predictions.
#' @param binWidth Positive numeric. Size of a bin. Each bin will be \code{binWidth * (max - min)} where \code{max} and \code{min}. If \code{autoWindow} is \code{FALSE} (the default) then \code{min} is 0 and \code{max} is 1. If \code{autoWindow} is \code{TRUE} then \code{min} and \code{max} are the maximum and minimum value of all predictions in the background and presence sets (i.e., not necessarily 0 and 1).
#' @param presWeight Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.
#' @param contrastWeight Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.
#' @param autoWindow Logical. If \code{FALSE} calculate bin boundaries starting at 0 and ending at 1 + epsilon (where epsilon is a very small number to assure inclusion of cases that equal 1 exactly). If \code{TRUE} (default) then calculate bin boundaries starting at minimum predicted value and ending at maximum predicted value.
#' @param method Character. Type of correlation to calculate. The default is \code{'spearman'}, the Spearman rank correlation coefficient used by Boyce et al. (2002) and Hirzel et al. (2006), which is the traditional CBI. In contrast, \code{'pearson'} or \code{'kendall'} can be used instead.  See [stats::cor()] for more details.
#' @param dropZeros Logical. If TRUE then drop all bins in which the frequency of presences is 0.
#' @param graph Logical. If TRUE then plot P vs E and P/E versus bin.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @param bg Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.
#' @param bgWeight Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.
#' @param ... Other arguments (not used).
#' @return Numeric value.
#' @details CBI is the Spearman rank correlation coefficient between the proportion of sites in each prediction class and the expected proportion of predictions in each prediction class based on the proportion of the landscape that is in that class. Values >0 indicate the model's output is positively correlated with the true probability of presence.  Values <0 indicate it is negatively correlated with the true probability of presence.
#' @references Boyce, M.S., Vernier, P.R., Nielsen, S.E., and Schmiegelow, F.K.A.  2002.  Evaluating resource selection functions.  \emph{Ecological Modeling} 157:281-300)
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randon, C., and Guisan, A.  2006.  Evaluating the ability of habitat suitability models to predict species presences.  \emph{Ecological Modeling} 199:142-152.
#' @seealso \code{\link[stats]{cor}}, \code{\link{fpb}}, \code{\link{aucWeighted}}, \code{link[enmSdm]{contBoyce2x}}
#' @examples
#' set.seed(123)
#' pres <- sqrt(runif(100))
#' contrast <- runif(1000)
#' contBoyce(pres, contrast)
#' contBoyce2x(pres, contrast)
#' presWeight <- c(rep(1, 10), rep(0.5, 90))
#' contBoyce(pres, contrast, presWeight=presWeight)
#' contBoyce2x(pres, contrast, presWeight=presWeight)
#' \donttest{
#' # compare stability of CBI calculated with ecospat.boyce() in ecospat package
#' library(ecospat)
#' set.seed(123)
#' results <- data.frame()
#' for (perform in c(1, 1.5, 2)) {
#' 	for (i in 1:30) {
#'
#'    pres <- runif(100)^(1 / perform)
#'    contrast <- runif(1000)
#'
#'    cbi_enmSdm <- contBoyce(pres, contrast)
#'    cbi_ecospat <- ecospat.boyce(contrast, pres, PEplot=FALSE)$Spearman.cor
#'
#'    results <- rbind(
#'      results,
#'      data.frame(
#'        performance = rep(perform, 2),
#'        method = c('enmSdm', 'ecospat'),
#'        cbi = c(cbi_enmSdm, cbi_ecospat)
#'      )
#'    )
#'
#' 	}
#'
#' }
#'
#' results$performance[results$performance == 1] <- 'poor'
#' results$performance[results$performance == 1.5] <- 'OK'
#' results$performance[results$performance == 2] <- 'good'
#'
#' results$category <- paste0(results$method, '\n', results$performance)
#'
#' par(mfrow=c(1, 2))
#' boxplot(cbi ~ category,
#' 	data=results,
#' 	ylab='CBI',
#' 	main='CBI of poor, OK, and good models',
#' 	border=c(rep('darkred', 3),
#' 	rep('darkblue', 3))
#' )
#' plot(results$cbi,
#' 	pch=rep(c(21, 22, 23, 24), each=2),
#' 	contrast=ifelse(results$method == 'ecospat', 'darkred', 'cornflowerblue'),
#' 	main='Pairs of CBIs',
#' 	ylab='CBI'
#' )
#' legend('bottomright', fill=c('darkred', 'cornflowerblue'), legend=c('ecospat', 'enmSdm'))
#' }
#' @export

contBoyce <- function(
	pres,
	contrast,
	numBins = 101,
	binWidth = 0.1,
	presWeight = rep(1, length(pres)),
	contrastWeight = rep(1, length(contrast)),
	autoWindow = TRUE,
	method = 'spearman',
	dropZeros = TRUE,
	graph = FALSE,
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
	
	# # remove NAs
	# if (na.rm) {

		# cleanedPres <- omnibus::naOmitMulti(pres, presWeight)
		# pres <- cleanedPres[[1]]
		# presWeight <- cleanedPres[[2]]

		# cleanedContrast <- omnibus::naOmitMulti(contrast, contrastWeight)
		# contrast <- cleanedContrast[[1]]
		# contrastWeight <- cleanedContrast[[2]]

	# }
	
	# right hand side of each class (assumes max value is >0)
	lowest <- if (autoWindow) { min(c(pres, contrast), na.rm=na.rm) } else { 0 }
	highest <- if (autoWindow) { max(c(pres, contrast), na.rm=na.rm) + omnibus::eps() } else { 1 + omnibus::eps() }

	windowWidth <- binWidth * (highest - lowest)

	lows <- seq(lowest, highest - windowWidth, length.out=numBins)
	highs <- seq(lowest + windowWidth + omnibus::eps(), highest, length.out=numBins)

	##########
	## MAIN ##
	##########

	## initiate variables to store predicted/expected (P/E) values
	freqPres <- freqContrast <- rep(NA, length(numBins))

	### tally proportion of test presences/background sites in each class
	for (countClass in 1:numBins) {

		# number of presence predictions in this class
		presInBin <- pres >= lows[countClass] & pres < highs[countClass]
		presInBin <- presInBin * presWeight
		freqPres[countClass] <- sum(presInBin, na.rm=na.rm)

		# number of background predictions in this class
		bgInBin <- contrast >= lows[countClass] & contrast < highs[countClass]
		bgInBin <- bgInBin * contrastWeight
		freqContrast[countClass] <- sum(bgInBin, na.rm=na.rm)

	} # next predicted value class

	# mean bin prediction
	meanPred <- rowMeans(cbind(lows, highs))

	# add small number to each bin that has 0 background frequency but does have a presence frequency > 0
	if (any(freqPres > 0 & freqContrast == 0)) {
		smallValue <- min(0.5 * c(presWeight[presWeight > 0], contrastWeight[contrastWeight > 0]))
		freqContrast[freqPres > 0 & freqContrast == 0] <- smallValue
	}

	# remove classes with 0 presence frequency
	if (dropZeros && 0 %in% freqPres) {
		zeros <- which(freqPres == 0)
		meanPred[zeros] <- NA
		freqPres[zeros] <- NA
		freqContrast[zeros] <- NA
	}

	# remove classes with 0 background frequency
	if (any(0 %in% freqContrast)) {
		zeros <- which(freqContrast == 0)
		meanPred[zeros] <- NA
		freqPres[zeros] <- NA
		freqContrast[zeros] <- NA
	}

	P <- freqPres / sum(presWeight, na.rm=TRUE)
	E <- freqContrast / sum(contrastWeight, na.rm=TRUE)
	PE <- P / E

	# plot
	if (graph) {
		par(mfrow=c(1, 2))
		lims <- c(0, max(P, E, na.rm=TRUE))
		plot(E, P, col='white', xlab='Expected', ylab='Predicted', main='P/E\nNumbered from lowest to highest class', xlim=lims, ylim=lims)
		text(E, P, labels=1:numBins, col=1:20)
		plot(meanPred, PE, type='l', xlab='Mean Prediction in Bin', ylab='P/E Ratio', main='CBI\nNumbered from lowest to highest class')
		text(meanPred, PE, labels=1:numBins, col='blue')
	}

	# remove NAs
	out <- omnibus::naOmitMulti(meanPred, PE)
	meanPred <- out[[1]]
	PE <- out[[2]]

	# calculate continuous Boyce index (cbi)
	cbi <- stats::cor(x=meanPred, y=PE, method=method)
	cbi

}
