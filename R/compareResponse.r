#' Compare two response curves along one or more predictors
#'
#' This function calculates a suite of metrics reflecting of niche overlap for two response curves. Response curves are predicted responses of a uni- or multivariate model along a single variable. Depending on the user-specified settings the function calculates these values either at each pair of values of \code{pred1} and \code{pred2} \emph{or} along a smoothed version of \code{pred1} and \code{pred2}.
#' @param pred1 Numeric list. Predictions from first model along \code{data} (one value per row in \code{data}).
#' @param pred2 Numeric list. Predictions from second model along \code{data} (one value per row in \code{data}).
#' @param data Data frame or matrix corresponding to \code{pred1} and \code{pred2}.
#' @param predictor Character list. Name(s) of predictor(s) for which to calculate comparisons. These must appear as column names in \code{data}.
#' @param adjust Logical. If \code{TRUE} then subtract the mean of \code{pred1} from \code{pred1} and the mean of \code{pred2} from \code{pred2} before analysis. Useful for comparing the shapes of curves while controlling for different elevations (intercepts).
#' @param gap Numeric >0. Proportion of range of predictor variable across which to assume a gap exists. Calculation of \code{areaAbsDiff} will  ignore gaps wide than this. To ensure the entire range of the data is included set this equal to \code{Inf} (default).
#' @param smooth Logical. If TRUE then the responses are first smoothed using loess() then compared at \code{smoothN} values along each predictor. If FALSE then comparisons are conducted at the raw values \code{pred1} and \code{pred2}.
#' @param smoothN \code{NULL} or positive integer. Number of values along "pred" at which to calculate comparisons. Only used if \code{smooth} is \code{TRUE}. If \code{NULL}, then comparisons are calculated at each value in data. If a number, then comparisons are calculated at \code{smoothN} values of \code{data[ , pred]} that cover the range of \code{data[ , pred]}.
#' @param smoothRange 2-element numeric list or \code{NULL}. If \code{smooth} is TRUE, then force loess'ed predictions < \code{smoothRange[1]} to equal \code{smoothRange[1]} and predictions > \code{smoothRange[2]} to equal \code{smoothRange[2]}. Ignored if \code{NULL}.
#' @param graph Logical. If \code{TRUE} then plot predictions.
#' @param ... Arguments to pass to functions like \code{sum()} (for example, \code{na.rm=TRUE}) and to \code{overlap()} (for example, \code{w} for weights). Note that if \code{smooth} = TRUE then passing an argument called \code{w} will likely cause a warning and make results circumspect \emph{unless} weights are pre-calculated for each of the \code{smoothN} points along a particular predictor.
#' @return Either a data frame (if \code{smooth = FALSE} or a list object with the smooth model plus a data frame (if \code{smooth = TRUE}) . The data frame represents metrics comparing response curves of \code{pred1} and \code{pred2}:
#' \itemize{
#' \item \code{predictor} Predictor for which comparison was made
#' \item \code{n} Number of values of predictor at which comparison was calculated
#' \item \code{adjust} \code{adjust} argument.
#' \item \code{smooth} \code{smooth} argument.
#' \item \code{meanDiff} Mean difference between predictions of \code{pred1} and \code{pred2} (higher ==> more different).
#' \item \code{meanAbsDiff} Mean absolute value of difference between predictions of \code{pred1} and \code{pred2} (higher ==> more different).
#' \item \code{areaAbsDiff} Sum of the area between curves predicted by \code{pred1} and \code{pred2}, standardized by total potential area between the two curves (i.e., the area available between the minimum and maximum prediction along the minimum and maximum values of the predictor) (higher ==> more different).
#' \item \code{d} Schoener's \emph{D}
#' \item \code{i} Hellinger's \emph{I} (adjusted to have a range [0, 1])
#' \item \code{esp} Godsoe's ESP
#' \item \code{cor} Pearson correlation between predictions of \code{pred1} and \code{pred2}.
#' \item \code{rankCor} Spearman rank correlation between predictions of \code{pred1} and \code{pred2}.
#' }
#' @references Warren, D.L., Glor, R.E., and Turelli, M.  2008.  Environmental niche equivalency versus conservatism: Quantitative approaches to niche evolution.  Evolution 62:2868-2883.
#' @references Warren, D.L., Glor, R.E., and Turelli, M.  2008.  Erratum.  Evolution 62:2868-2883.
#' @references Godsoe, W.  2014.  Inferring the similarity of species distributions using Species’ Distribution Models.  Ecography 37:130-136.
#' @seealso \code{\link[enmSdm]{compareNiches}}
#' @examples
#' set.seed(123)
#' data <- data.frame(
#' 	x1=seq(-1, 1, length.out=100),
#' 	x2=seq(-1, 1, length.out=100) + rnorm(100, 0, 0.3)
#' )
#' 
#' pred1 <- 1 / (1 + exp(-(0.3 + 2 * (data$x1 - 0.2) -0.3 * data$x2)))
#' pred2 <- 1 / (1 + exp(-(-0 + 0.1 * data$x1 - 4 * data$x1^2 + 0.4 * data$x2)))
#' 
#' compareResponse(pred1, pred2, data, graph=TRUE)
#' compareResponse(pred1, pred2, data, smooth=TRUE, graph=TRUE)
#' compareResponse(pred1, pred2, data, adjust=TRUE, graph=TRUE)
#' @export

compareResponse <- function(
	pred1, pred2,
	data,
	predictor = names(data),
	adjust = FALSE,
	gap = Inf,
	smooth = FALSE,
	smoothN = 1000,
	smoothRange = c(0, 1),
	graph = FALSE,
	...
) {

	args <- list(...)

	# for each predictor
	for (thisPred in predictor) {

		# sort predictor and predictions by predictor
		x <- data[ , thisPred]
		orderOfX <- order(x)
		x <- x[orderOfX]
		pred1 <- pred1[orderOfX]
		pred2 <- pred2[orderOfX]

		# remember original predictions (for plotting)
		origPred1 <- pred1
		origPred2 <- pred2

		### if adjusting elevation
		if (adjust) {

			adjust1 <- mean(pred1, na.rm=TRUE)
			adjust2 <- mean(pred2, na.rm=TRUE)

			pred1 <- pred1 - adjust1
			pred2 <- pred2 - adjust2

			minPred <- min(min(pred1, na.rm=TRUE), min(pred2, na.rm=TRUE)) - omnibus::eps()

			pred1 <- pred1 - minPred
			pred2 <- pred2 - minPred

		}

		### if using smooth to smooth predictions
		if (smooth) {

			# smooth
			smooth1 <- stats::loess(pred1 ~ x, data=data.frame(pred1=pred1, x=x), na.rm=TRUE)
			smooth2 <- stats::loess(pred2 ~ x, data=data.frame(pred2=pred2, x=x), na.rm=TRUE)

			# remember original X and re-define X
			if (!is.null(smoothN)) {
				origX <- x
				x <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=smoothN)
			}

			pred1 <- stats::predict(smooth1, newdata=x)
			pred2 <- stats::predict(smooth2, newdata=x)

			# force to user-specified range
			if (!is.null(smoothRange)) {

				pred1[pred1 < smoothRange[1]] <- smoothRange[1]
				pred1[pred1 > smoothRange[2]] <- smoothRange[2]

				pred2[pred2 < smoothRange[1]] <- smoothRange[1]
				pred2[pred2 > smoothRange[2]] <- smoothRange[2]

			}

		}

		### calculate comparisons

		# basic comparisons
		sim <- compareNiches(x1=pred1, x2=pred2, logit=!adjust, na.rm=TRUE)

		## differences in area under each curve (curve1 - curve2)
		gapRange <- gap * (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))

		areaAbsDiff <- possibleArea <- 0

		# for each set of 4 points defining a quadrilateral (two values of predictor and two predictions)
		for (count in 2:length(pred1)) {

			# if distance between the two value of the predictor aren't too big (this is not a gap)
			if (x[count] - x[count - 1] <= gapRange) {

				# absolute areal difference between pred1 and pred2 curves
				thisArea <- 0.5 * (x[count] - x[count - 1]) * abs(pred1[count - 1] - pred2[count - 1]) * abs(pred1[count] - pred2[count])

				areaWeight <- if (exists('w', where=args)) {
					mean(args$w[c(count, count - 1)], ...)
				} else {
					1
				}

				areaAbsDiff <- areaAbsDiff + thisArea * areaWeight
				possibleArea <- possibleArea + (x[count] - x[count - 1]) * areaWeight

			}

		}

		areaAbsDiff <- areaAbsDiff / possibleArea

		### remember
		thisOut <- data.frame(
			predictor=thisPred,
			n=length(x),
			adjust=adjust,
			smooth=smooth,
			smoothN=if (smooth & !is.null(smoothN)) { smoothN } else { NA },
			smoothRange=if (smooth & !is.null(smoothRange)) { paste(smoothRange, collapse=' ') } else { NA },
			gap=gap,
			meanDiff=sim['meanDiff'],
			meanAbsDiff=sim['meanAbsDiff'],
			areaAbsDiff=areaAbsDiff,
			d=sim['d'],
			i=sim['i'],
			esp=sim['esp'],
			cor=sim['cor'],
			rankCor=sim['rankCor']
		)

		out <- if (exists('out', inherits=FALSE)) {
			rbind(out, thisOut)
		} else {
			thisOut
		}

		### plot
		if (graph) {

			graphics::par(mfrow=c(1, 2), ask=TRUE)

			lims <- c(min(pred1, pred2, origPred1, origPred2, na.rm=TRUE), max(pred1, pred2, origPred1, origPred2, na.rm=TRUE))

			### predictions #1 vs predictions #2
			plot(pred1, pred2, col='white', main=paste0('Model 1 vs Model 2 for ', thisPred), xlab=paste('Model 1 Prediction', ifelse(adjust, '(Adjusted)', '')), ylab=paste('Model 2 Prediction', ifelse(adjust, '(Adjusted)', '')), xlim=lims, ylim=lims)
			graphics::abline(0, 1, col='gray')

			# original predictions
			if (smooth) graphics::points(origPred1, origPred2, pch=16, cex=0.8, col=scales::alpha('black', 0.2))

			# pred1 vs pred2
			graphics::points(pred1, pred2, col='darkgreen')

			# fake legend (stats)
			graphics::legend('bottomright', inset=0.01, bty='n', lwd=NA, pch=NA,
				legend=c(
					paste0('cor = ', format(round(sim['cor'], 2), nsmall=2)),
					paste0('rankCor = ', format(round(sim['rankCor'], 2), nsmall=2))
				)
			)

			### predictions #1 and #2 vs this predictor
			plot(x, pred1, col='white', main=paste('Responses versus', thisPred), xlab=thisPred, ylab=paste('Prediction', ifelse(adjust, '(Adjusted)', '')), ylim=lims)

			# plot original predictions if using smooth
			if (smooth) {
				# origPred1 <- origPred1[order(x)]
				# origPred2 <- origPred2[order(x)]
				graphics::points(origX, origPred1, pch=16, cex=0.8, col=scales::alpha('blue', 0.5))
				graphics::points(origX, origPred2, pch=16, cex=0.8, col=scales::alpha('red', 0.3))
			}

			# plot areal difference between curves
			graphics::polygon(x=c(x, rev(x)), y=c(pred1, rev(pred2)), border=NA, col=scales::alpha('black', 0.3))

			graphics::lines(x, pred1, type='l', col='blue')
			graphics::lines(x, pred2, type='l', col='red')
			graphics::rug(x, ticksize=0.01)
			graphics::rug(origPred1, ticksize=0.01, side=2, col='blue')
			graphics::rug(origPred2, ticksize=0.01, side=4, col='red')

			# real legend
			graphics::legend('topright', inset=0.01, bty='n', pch=NA, col=c('blue', 'red', NA), fill=c(NA, NA, 'gray70'), border=c(NA, NA, 'black'), legend=c('model 1', 'model 2', 'areaAbsDiff'), lwd=2)

			# fake legend (stats)
			graphics::legend('bottomright', inset=0.01, bty='n', lwd=NA, pch=NA,
				legend=c(
					paste0('meanDiff = ', format(round(sim['meanDiff'], 2), nsmall=2)),
					paste0('meanAbsDiff = ', format(round(sim['meanAbsDiff'], 2), nsmall=2)),
					paste0('areaAbsDiff = ', format(round(areaAbsDiff, 2), nsmall=2)),
					paste0('d = ', format(round(sim['d'], 2), nsmall=2)),
					paste0('i = ', format(round(sim['i'], 2), nsmall=2)),
					paste0('esp = ', format(round(sim['esp'], 2), nsmall=2))
				)
			)

		}

	} # next predictor

	#####################
	## POST-PROCESSING ##
	#####################

	if (!smooth) {
		out
	} else {
		list(out, smooth)
	}

}
