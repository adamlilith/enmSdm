#' Niche overlap for a set of iterated "randPointsRespecting~" functions
#'
#' This function is called using a list object typically generated using the \code{\link[enmSdm]{randPontsMaster}} function (plus sometimes followed by the \code{\link[enmSdm]{randPointsExtract}} and \code{\link[enmSdm]{randPointsSampled}} functions). It calculates niche overlap for each set of points. Essentially it is a wrapper for \code{\link[enmSdm]{nicheOverlap}}.
#' @param rands A list object typically generated using the \code{\link[enmSdm]{randPontsMaster}} function.
#' @param env Either a data frame, matrix, or any object that can be coerced to a data frame containing environmental data at available background sites, \emp{or} an object of class \code{princomp} representing a principal components analysis generated using the \code{\link[stats]{princomp}} function with argument \code{scores = TRUE}.
#' @param vars Either a character list naming columns in \code{x1}, \code{x2}, and \code{x3} to be used as environmental data, \emph{or} positive integers indexing the columns to be used as environmental data.
#' @param x Either \code{NULL} (default) or a data frame, matrix, SpatialPointsDataFrame, or other object that can be coerced to a data frame. If supplied then the objects (usually species) represented in \code{rands} are compared to this set of environmental values. This argument \emph{must} be supplied if \code{rands} was generated using \code{\link[enmSdm]{randPointsRespectingSelf}} or \code{\link[enmSdm]{randPointsRespectingSelfOther1}}.
#' @param bins Number of bins into which to divide the environmental space (default is 100 on each side).
#' @param cor Logical, if \code{TRUE} (default), then the PCA used to construct the environmental space will use the correlation matrix (this is highly recommended if the variables are on different scales). This is ignored if \code{env} is an object of class \code{princomp}.

randPointsBatchNicheOverlap <- function(
	rands,
	env,
	vars,
	x = NULL,
	bins = 100,
	cor = TRUE
) {
	
	if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
		if (is.null(x)) stop('Argument "x" must be specified if argument "rands" was generated using either "randPointsRespectingSelf" or "randPointsRespectingSelfOther1".')
	}
	
	out <- data.frame()
	
	for (i in seq_along(rands)) {
		
		if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
			x1 <- rands[[i]]
			x2 <- x
		} else {
			x1 <- rands[[i]]$x1rand
			x2 <- rands[[i]]$x2rand
		}
	
		thisOverlap <- nicheOverlap(x1=x1, x2=x2, env=env, vars=vars, bins=bins, cor=cor)
		thisOverlap <- as.data.frame(thisOverlap)
		out <- rbind(out, thisOverlap)
	
	}
	
	out <- cbind(
		data.frame(iter = 1:nrow(out)),
		out
	)
	
	out
	
}
