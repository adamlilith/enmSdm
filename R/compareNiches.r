#' Metrics of niche overlap
#'
#' This function calculates several metrics of niche overlap based on predictions for two species (or for the same species but different models) at the same sites.
#' @param x1 Numeric. Vector or matrix of predictions from a model.
#' @param x2 Numeric. Vector or matrix of predictions from another model.
#' @param method Character vector, indicates type of metric to calculate:
#' \itemize{
#' \item \code{meanDiff} mean difference between \code{x1} and \code{x2}
#' \item \code{meanAbsDiff} mean absolute difference between \code{x1} and \code{x2} (ie, \code{sum(abs(x1 - x2))})
#' \item \code{d} Schoener's \emph{D}
#' \item \code{i} Warren's \emph{I}
#' \item \code{esp} Godsoe's \emph{ESP}
#' \item \code{rho} Correlation between \code{x1} and \code{x2} (will apply \code{logitAdj()} first unless logit=FALSE).
#' \item \code{rankCor}  Pearson rank correlation.
#' }
#' @param w Numeric list. Weights of predictions in \code{x1} and \code{x2}.
#' @param na.rm Logical.  If T\code{TRUE} then remove elements in \code{x1} and \code{2} that are \code{NA} in \emph{either} \code{x1} or \code{x2}.
#' @param ... Other arguments (not used).
#' @return List object with one element per value specified by the argument in \code{method}.
#' @seealso \code{\link{compareResponse}}
#' @examples
#' x1 <- seq(0, 1, length.out=100)
#' x2 <- x1^2
#' compareNiches(x1, x2)
#' @export

compareNiches <- function(
	x1,
	x2,
	method = c('meanDiff', 'meanAbsDiff', 'd', 'i', 'esp', 'rho', 'rankCor'),
	w = rep(1, length(x1)),
	na.rm = FALSE,
	...
) {

	x1 <- c(x1)
	x2 <- c(x2)

	# remove NAs
	if (na.rm) {
		out <- omnibus::naOmitMulti(x1, x2, w)
		x1 <- out[[1]]
		x2 <- out[[2]]
		w <- out[[3]]
	}

	sim <- numeric()

	# mean difference
	if ('meanDiff' %in% method) sim <- c(sim, sum(w * (x1 - x2)) / sum(w))

	# mean absolute difference
	if ('meanAbsDiff' %in% method) sim <- c(sim, sum(w * abs(x1 - x2)) / sum(w))

	# calculate Schoener's D
	if ('d' %in% method) {

		x1Stand <- (x1 * w) / sum(x1 * w)
		x2Stand <- (x2 * w) / sum(x2 * w)

		sim <- c(sim, 1 - 0.5 * sum(abs(x1Stand - x2Stand)))

	}

	# calculate Warren's I
	if ('i' %in% method) {

		x1Stand <- (x1 * w) / sum(x1 * w)
		x2Stand <- (x2 * w) / sum(x2 * w)

		sim <- c(sim, 1 - 0.5 * sqrt(sum((sqrt(x1Stand) - sqrt(x2Stand))^2)))

	}

	# Godsoe's ESP
	if ('esp' %in% method) sim <- c(sim, 2 * sum(w * x1 * x2) / sum(w * c(x1 + x2)))

	# correlation
	if ('rho' %in% method) sim <- c(sim, boot::corr(cbind(x1, x2), w=w))

	# rank correlation
	if ('rankCor' %in% method) sim <- c(sim, stats::cor(w * x1, w * x2, method='spearman'))

	names(sim) <- method
	sim

}
