#' Calculates biotic velocity
#'
#' XXX
#' @param x XXX
#' @param ... Other arguments (not used).
#' @return XXXX
#' @seealso
#' @examples
#' @export

bioticVelocity <- function(
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
