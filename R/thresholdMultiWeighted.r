#' Weighted thresholds for multi-class predictions
#'
#' This function calculates thresholds for cases where predictions are intended to be indicative of the state of three or more ordered categories. For example, a model might be intended to differentiate between sites where a species is present, sites where the species has been absent for a short time, and sites where the species has been absent for a long time. In this case, two thresholds would be needed to separate presences from short-term absences, and short-term absences from long-term absences. It is possible to simply calculate thresholds for each of these pairwise cases while ignoring the third case. However, this is not guaranteed to achieve any particular thresholding rule across all categories (e.g., maximization of the sum of sensitivity and specificity). Note that it is possible to have multiple combinations of thresholds that do equally well under any given rule for finding the "best" set of thresholds.
#'
#' @param ... A set of two or more objects, each of which can be any of a 1- or 2-column matrix or data frame or a numeric vector. The objects must be listed in \emph{reverse} order of \emph{expected} probability. For example, you might have a set of predictions for objects you expect to have a high predicted probability of occurrence (e.g., presences of an animal), a set that you expect to have middle levels of probability of occurrence (e.g., sites that were recently vacated), and a set for which you expect a low level of predicted probability (e.g., sites that have been long-vacated). In this example you should list the cases in order: present-cases, short-term absence cases, long-term absence cases. If a 2-column matrix or data frame is supplied, then the first column is assumed to represent predictions and the second assumed to represent weights.  If you do not specify weights then every case will be assumed to have the same weight.
#' @param thresholds Numeric vector. Thresholds at which to calculate the sum of sensitivity and specificity. The default evaluates all values from 0 to 1 in steps of 0.01.
#' @param na.rm Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.
#' @return A data frame is returned, with one row per valid combination of thresholds separating category A from B, B from C, and so on. Columns are as:
#' \itemize
#' 		\item The true positive rate ("tpr", also known as sensitivity) and true negative rate ("tnr", also known as specificity) are reported for each category.
#'		\item Column \code{meanMinDiffTprTnr} is the mean difference between each category's TPR and TNR from the mean across all TPRs and TNRs, and column \code{minDiffTprTnr} is \code{TRUE} for the row(s) that minimize this difference.
#'		\item Column \code{sumSqDiffTprTnr} is the sum of the squared mean difference between each category's TPR and TNR from the mean across all TPRs and TNRs, and column \code{minSqDiffTprTnr} is \code{TRUE} for the row(s) that minimize this value.
#'		\item Column \code{sumTprTnr} is the sum all TPRs and TNRs, and column \code{maxSumTprTnr} is \code{TRUE} for the row(s) that maximize this value.
#' }
#'
#' @seealso \code{\link{thresholdWeighted}}, \code{\link[dismo]{threshold}}, \code{\link[dismo]{evaluate}}
#' @examples
#' set.seed(123)
#' par(mfrow=c(1, 3))
#' 
#' thresholds <- seq(0, 1, by=0.05)
#' 
#' ### 3 cases of equal size, no weights assigned so all are equal
#' A <- runif(20, 0.4, 1)
#' B <- runif(20, 0.3, 0.8)
#' C <- runif(20, 0, 0.6)
#' 
#' r <- thresholdMultiWeighted(A, B, C, thresholds=thresholds)
#' 
#' # size of circle indicates sum of sensitivity and specificity for each
#' # combination of thresholds
#' plot(r$thold_A_over_B, r$thold_B_over_C, cex=r$sumTprTnr,
#' 	xlab='Threshold separating A and B',
#' 	ylab='Threshold separating B and C',
#' 	main='Equal cases, equal weighting',
#' 	xlim=c(0, 1), ylim=c(0, 1)
#' )
#' 
#' # pluses indicate which combination of thresholds maximize the
#' # sum of sensitivities and specificities
#' maxs <- which(r$maxSumTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=3)
#' 
#' # down-triangles indicate which combination of thresholds minimize
#' # differences mean differences between sensitivities and specificities
#' # and their mean
#' maxs <- which(r$minDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=6)
#' 
#' # up-triangles indicate which combination of thresholds minimize
#' # differences sum of squared differences between sensitivities and
#' # specificities and their mean
#' # (in this example, they are the the same as previous thresholds)
#' maxs <- which(r$minSqDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=2)
#' 
#' 
#' ### 3 cases of unequal size, no weights assigned so all are equal
#' A <- runif(20, 0.4, 1)
#' B <- runif(200, 0.3, 0.8)
#' C <- runif(20, 0, 0.6)
#' 
#' r <- thresholdMultiWeighted(A, B, C, thresholds=thresholds)
#' 
#' # size of circle indicates sum of sensitivity and specificity for each
#' # combination of thresholds
#' plot(r$thold_A_over_B, r$thold_B_over_C, cex=r$sumTprTnr,
#' 	xlab='Threshold separating A and B',
#' 	ylab='Threshold separating B and C',
#' 	main='Case B most numerous',
#' 	xlim=c(0, 1), ylim=c(0, 1)
#' )
#' 
#' # pluses indicate which combination of thresholds maximize the
#' # sum of sensitivities and specificities
#' maxs <- which(r$maxSumTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=3)
#' 
#' # down-triangles indicate which combination of thresholds minimize
#' # differences mean differences between sensitivities and specificities
#' # and their mean
#' maxs <- which(r$minDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=6)
#' 
#' # up-triangles indicate which combination of thresholds minimize
#' # differences sum of squared differences between sensitivities and
#' # specificities and their mean
#' # (in this example, they are the the same as previous thresholds)
#' maxs <- which(r$minSqDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=2)
#' 
#' ### 3 cases of equal size, weights assigned so getting case A right
#' ### is "worth" more
#' A <- runif(20, 0.4, 1)
#' B <- runif(20, 0.3, 0.8)
#' C <- runif(20, 0, 0.6)
#' 
#' A <- cbind(A, rep(1, 20))
#' B <- cbind(B, rep(0.5, 20))
#' C <- cbind(C, rep(0.1, 20))
#' 
#' r <- thresholdMultiWeighted(A, B, C, thresholds=thresholds)
#' 
#' # size of circle indicates sum of sensitivity and specificity for each
#' # combination of thresholds
#' plot(r$thold_A_over_B, r$thold_B_over_C, cex=r$sumTprTnr,
#' 	xlab='Threshold separating A and B',
#' 	ylab='Threshold separating B and C',
#' 	main='Weight of case A weight > case B > case C',
#' 	xlim=c(0, 1), ylim=c(0, 1)
#' )
#' 
#' # pluses indicate which combination of thresholds maximize the
#' # sum of sensitivities and specificities
#' maxs <- which(r$maxSumTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=3)
#' 
#' # down-triangles indicate which combination of thresholds minimize
#' # differences mean differences between sensitivities and specificities
#' # and their mean
#' maxs <- which(r$minDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=6)
#' 
#' # up-triangles indicate which combination of thresholds minimize
#' # differences sum of squared differences between sensitivities and
#' # specificities and their mean
#' # (in this example, they are the the same as previous thresholds)
#' maxs <- which(r$minSqDiffTprTnr)
#' points(r$thold_A_over_B[maxs], r$thold_B_over_C[maxs], pch=2)
#' 
#' @export
thresholdMultiWeighted <- function(
	...,
	thresholds = seq(0, 1, by=0.01),
	na.rm = FALSE
) {


	cases <- list(...)
	numCases <- length(cases)
	names(cases) <- if (numCases > 26) {
		c(LETTERS, letters[1:(numCases - 25)])
	} else {
		LETTERS[1:numCases]
	}
	caseNames <- names(cases)

	# if input is a data frame or matrix with just one column, add another
	# to represent weights (all weights = 1), if input is a vector then
	# convert to two-column matrix with second column being weights of 1
	for (i in seq_along(cases)) {
		this <- cases[[i]]
		if (!any(c('data.frame', 'matrix') %in% class(this))) {
			cases[[i]] <- cbind(this, w=rep(1, length(this)))
		} else {
			if (ncol(this) == 1) {
				cases[[i]] <- cbind(this, w=rep(1, length(this)))
			} else if (ncol(this) > 2) {
				stop('Category #', i, ' has >2 columns. Values must be a vector or a 1- or 2-column matrix or data frame.')
			}
		}
	}
	
	# remove NAs
	if (na.rm) {

		for (i in seq_along(cases)) {
			
			cases[[i]] <- cases[[i]][complete.cases(cases[[i]]), ]
			
		}

	}

	### calculate thresholds
	########################

	# construct matrix with all valid thresholds and TPR and TNR values
	tholds <- list(thresholds)
	for (i in 2:(numCases - 1)) tholds <- c(tholds, list(thresholds))
	names(tholds) <- paste0('thold_', caseNames[1:(numCases - 1)], '_over_', caseNames[2:numCases])
	
	tholds <- expand.grid(tholds)
	for (i in 2:ncol(tholds)) tholds <- tholds[which(tholds[ , i] <= tholds[ , i - 1]), ]
	
	tprTnr <- matrix(NA, ncol=2 * numCases, nrow=nrow(tholds))
	colnames(tprTnr) <- c(paste0('tpr', caseNames), paste0('tnr', caseNames))
	tprTnr <- tprTnr[ , -which(colnames(tprTnr) %in% c(paste0('tpr', tail(caseNames, 1)), 'tnrA'))]
	
	out <- cbind(tholds, tprTnr)
	out <- as.data.frame(out)
	
	### calculate TPR and TNR for all valid combinations
	for (i in 1:nrow(out)) {
	
		for (one in 1:(numCases - 1)) {
	
			two <- one + 1
			tholdName <- paste0('thold_', caseNames[one], '_over_', caseNames[two])
			thisThold <- out[i, tholdName]
	
			w <- cases[[one]][ , 2]
			W <- sum(w)
			out[i, paste0('tpr', caseNames[one])] <- (sum(w * (cases[[one]][ , 1] >= thisThold)) / W)
			
			w <- cases[[two]][ , 2]
			W <- sum(w)
			out[i, paste0('tnr', caseNames[two])] <- (sum(w * (cases[[two]][ , 1] < thisThold)) / W)
	
		}
	
	}
		
	# minimum differences between TPRs and TNRs
	tprsTnrs <- out[ , c(paste0('tpr', caseNames[1:(numCases - 1)]), paste0('tnr', caseNames[2:numCases]))]
	tprTnrMeans <- rowMeans(tprsTnrs)
	
	diffs <- abs(tprsTnrs - cbind(tprTnrMeans))
	diffs <- rowMeans(diffs)
	out$meanMinDiffTprTnr <- diffs
	out$minDiffTprTnr <- FALSE
	out$minDiffTprTnr[diffs == min(diffs)] <- TRUE
	
	# minimum SQUARED differences between TPRs and TNRs
	sqDiffs <- (tprsTnrs - cbind(tprTnrMeans))^2
	sqDiffs <- rowSums(sqDiffs)
	out$sumSqDiffTprTnr <- sqDiffs
	out$minSqDiffTprTnr <- FALSE
	out$minSqDiffTprTnr[sqDiffs == min(sqDiffs)] <- TRUE

	# maximum TPR + TNR
	sumTprTnr <- rowSums(tprsTnrs)
	out$sumTprTnr <- sumTprTnr
	out$maxSumTprTnr <- FALSE
	out$maxSumTprTnr[sumTprTnr == max(sumTprTnr)] <- TRUE
	
	
	
	out
	
}
