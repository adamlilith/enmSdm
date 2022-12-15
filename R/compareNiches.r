#' Metrics of niche overlap
#'
#' This function calculates several metrics of niche overlap based on predictions for two species (or for the same species but different models) at the same sites.
#' @param x1 Numeric. Vector of predictions from a model.
#' @param x2 Numeric. Vector of predictions from another model.
#' @param method Character vector, indicates type of metric to calculate:
#' \itemize{
#' 
#' \item \code{meanDiff}: Average difference
#' \item \code{meanAbsDiff}: Average of absolute values of difference
#' \item \code{rmsd}: Root-mean square deviation
#' \item \code{d}: Schoener's \emph{D}
#' \item \code{i}: Warren's \emph{I}
#' \item \code{esp}: Godsoe's \emph{ESP}
#' \item \code{cor}: Pearson correlation between \code{x1} and \code{x2} (will apply \code{logitAdj()} first unless logit=FALSE).
#' \item \code{rankCor}: Spearman rank correlation.
#' }
#' @param w Numeric list. Weights of predictions in \code{x1} and \code{x2}.
#' @param na.rm Logical.  If T\code{TRUE} then remove elements in \code{x1} and \code{2} that are \code{NA} in \emph{either} \code{x1} or \code{x2}.
#' @param ... Other arguments (not used).
#' 
#' @return List object with one element per value specified by the argument in \code{method}.
#' @references Warren, D.L., Glor, R.E., and Turelli, M.  2008.  Environmental niche equivalency versus conservatism: Quantitative approaches to niche evolution. \emph{Evolution} 62:2868-2883. \doi{10.1111/j.1558-5646.2008.00482.x}
#' @references Warren, D.L., Glor, R.E., and Turelli, M.  2008.  Erratum.  \emph{Evolution} 62:2868-2883. \doi{10.1111/j.1558-5646.2010.01204.x}
#' @references Godsoe, W.  2014.  Inferring the similarity of species distributions using Species’ Distribution Models.  \emph{Ecography} 37:130-136. \doi{10.1111/j.1600-0587.2013.00403.x}
#' 
#' @seealso \code{\link[enmSdm]{compareResponse}}
#' @examples
#'
#' x1 <- seq(0, 1, length.out=100)
#' x2 <- x1^2
#' compareNiches(x1, x2)
#'
#' @export

compareNiches <- function(
	x1,
	x2,
	method = c('meanDiff', 'meanAbsDiff', 'rmsd', 'd', 'i', 'esp', 'cor', 'rankCor'),
	w = rep(1, length(x1)),
	na.rm = FALSE,
	...
) {

	x1 <- c(x1)
	x2 <- c(x2)
	
	# remove NAs
	if (na.rm) {
		out <- omnibus::naOmitMulti(x1, x2, w)
		x1 <- out[[1L]]
		x2 <- out[[2L]]
		w <- out[[3L]]
	}

	# weight
	x1w <- w * x1
	x2w <- w * x2

	# standardize weighted values to sum to 1
	wSum <- sum(w)
	x1wStand <- x1w / sum(x1w)
	x2wStand <- x2w / sum(x2w)
	
	sim <- numeric()

	# mean difference
	if ('meanDiff' %in% method) sim <- c(sim, sum(x1w - x2w) / wSum)

	# mean absolute difference
	if ('meanAbsDiff' %in% method) sim <- c(sim, sum(abs(x1w - x2w)) / wSum)

	# RMSD
	if ('rmsd' %in% method) sim <- c(sim, sqrt(sum((x1w - x2w)^2) / wSum))

	# Schoener's D
	if ('d' %in% method) sim <- c(sim, 1 - 0.5 * sum(abs(x1wStand - x2wStand)))

	# calculate Warren's I
	if ('i' %in% method) sim <- c(sim, 1 - 0.5 * sum((sqrt(x1wStand) - sqrt(x2wStand))^2))

	# Godsoe's ESP
	if ('esp' %in% method) sim <- c(sim, 2 * sum(w * x1 * x2) / (sum(w * (x1 + x2))))

	# correlation
	if ('cor' %in% method) sim <- c(sim, boot::corr(cbind(x1, x2, w=w)))

	# rank correlation
	if ('rankCor' %in% method) sim <- c(sim, stats::cor(w * x1, w * x2, method='spearman'))

	names(sim) <- method
	sim

}
