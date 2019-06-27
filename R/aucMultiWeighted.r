#' Calculate multivariate weighted AUC
#'
#' This function calculates a multivariate version of the area under the receiver-operator characteristic curve (AUC). The multivariate version is simply the mean AUC across all possible pairwise AUCs for all cases (based on Hand, DJ and Till, RJ. 2001. A simple generalisation of the area under the ROC curve for multiple class classification problems. \emph{Machine Learning} 45:171-186.).
#' @param ... A set of two or more objects, each of which can be any of a 1- or 2-column matrix or data frame or a numeric vector. The objects must be listed in order of \emph{expected} probability. For example, you might have a set of predictions for objects you expect to have a low predicted probability (e.g., long-term absences of an animal), a set that you expect to have middle levels of probability (e.g., sites that were recently vacated), and a set for which you expect a high level of predicted probability (e.g., sites that are currently occupied). In this case you should list the cases in order: low, middle, high. If a 2-column matrix or data frame is supplied, then the first column is assumed to represent predictions and the second assumed to represent weights.
#' @param weightBySize Logical, if \code{FALSE} then the multivariate measure of AUC will treat all comparisons as equal (e.g., low versus middle will weigh as much as middle versus high), and so will simply be the mean AUC across all possible comparisons. If \code{TRUE} then multivariate AUC is the weighted mean across all possible comparisons where weights are the number of comparisons between each of the two cases. For example, if a set of "low" predictions ("low") has 10 data points, "middle" has 10, and "high" has 20, then the multivariate AUC will be (10 * low + 10 * middle + 20 * high) / (10 + 10 + 20).
#' @param na.rm Logical. If \code{TRUE} then remove any cases in \code{...} that are \code{NA}.
#' @return Named numeric vector.
#' @seealso \code{\link{fpb}}, \code{\link{contBoyce}}, \code{\link[dismo]{evaluate}}, \code{link[enmSdm]{aucWeighted}}
#' @examples
#' set.seed(123)
#' 
#' # no weights
#' low <- runif(10)^2
#' middle <- runif(10)
#' high <- sqrt(runif(20))
#' 
#' aucMultiWeighted(low, middle, high)
#' 
#' # equal weights
#' low <- matrix(c(low, rep(1, length(low))), ncol=2)
#' middle <- matrix(c(middle, rep(1, length(middle))), ncol=2)
#' high <- matrix(c(high, rep(1, length(high))), ncol=2)
#' aucMultiWeighted(low, middle, high)
#' 
#' # equal weights with weighting by number of comparisons
#' aucMultiWeighted(low, middle, high, weightBySize=TRUE)
#' 
#' # unequal weights
#' middle[ , 2] <- ifelse(middle[ , 1] > 0.5, 0.1, 1)
#' aucMultiWeighted(low, middle, high)
#' 
#' # unequal weights with weighting by number of comparisons
#' aucMultiWeighted(low, middle, high, weightBySize=TRUE)
#' @export

aucMultiWeighted <- function(
	...,
	weightBySize = FALSE,
	na.rm = FALSE
) {

	cases <- list(...)
	namae <- ellipseNames(...)

	# if input is a data frame or matrix with just one column, add another
	# to represent weights (all weights = 1), if input is a vector then
	# convert to two-column matrix with second column being weights of 1
	for (i in seq_along(cases)) {
		this <- cases[[i]]
		if (class(this) == 'data.frame') this <- as.matrix(this)
		if (class(this) == 'matrix') {
			if (ncol(this) == 1) {
				cases[[i]] <- cbind(this, matrix(1, nrow=nrow(this)))
			}
		} else {
			cases[[i]] <- matrix(
				c(this, rep(1, length(this))),
				ncol=2
			)
		}
	}
	
	# remove NAs
	if (na.rm) {

		for (i in seq_along(cases)) {
			
			this <- cases[[i]]
			nas <- omnibus::naRows(this)
			if (length(nas) > 0) cases[[i]] <- this[-nas, , drop=FALSE]
			
		}

	}

	# AUC and number of cases
	aucs <- numCases <- numeric()
	
	# calculate AUC
	aucs <- numeric()
	for (one in 1:(length(cases) - 1)) {
		for (two in (one + 1):length(cases)) {

			# neatify
			contrast <- cases[[one]][ , 1]
			pres <- cases[[two]][ , 1]
			
			contrastWeight <- cases[[one]][ , 2]
			presWeight <- cases[[two]][ , 2]
			
			# pairwise AUC
			thisAuc <- aucWeighted(
				pres=pres,
				contrast=contrast,
				presWeight=presWeight,
				contrastWeight=contrastWeight
			)
			
			# remember and assign name
			aucs <- c(aucs, thisAuc)
			name1 <- namae[one]
			name2 <- namae[two]
			name <- paste0(name2, '_over_', name1)
			names(aucs)[length(aucs)] <- name
			
			# remember number of cases
			thisNumCases <- length(pres) * length(contrast)
			numCases <- c(numCases, thisNumCases)
			names(numCases)[length(numCases)] <- name
			
		}
	}
	
	# multivariate AUC
	aucMulti <- if (weightBySize) {
		sum(aucs * numCases) / sum(numCases)
	} else {
		mean(aucs)
	}
		
	aucs <- c(aucs, aucMulti)
	names(aucs)[length(aucs)] <- 'multivariate'
	
	aucs

}
