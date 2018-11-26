#' Calculate niche overlap as per Broennimann et al. (2012)
#' 
#' This function calculates niche overlap between two species.
#' @param x1 Data frame, matrix, or any object that can be coerced to a data frame containing environmental data at occurrence sites of a species.
#' @param x2 Data frame, matrix, or any object that can be coerced to a data frame containing environmental data at occurrence sites of another species.
#' @param env Either a data frame, matrix, or any object that can be coerced to a data frame containing environmental data at available background sites, \emp{or} an object of class \code{princomp} representing a principal components analysis generated using the \code{\link[stats]{princomp}} function with argument \code{scores = TRUE}.
#' @param vars Either a character list naming columns in \code{x1}, \code{x2}, and \code{x3} to be used as environmental data, \emph{or} positive integers indexing the columns to be used as environmental data.
#' @param bins Number of bins into which to divide the environmental space (default is 100 on each side).
#' @param cor Logical, if \code{TRUE} (default), then the PCA used to construct the environmental space will use the correlation matrix (this is highly recommended if the variables are on different scales). This is ignored if \code{env} is an object of class \code{princomp}.
#' @return List with these named elements:
#' \itemize{
#' \item \code{meanDiff} mean difference between binned, standardized densities of \code{x1} and \code{x2} in environmental space.
#' \item \code{meanAbsDiff} mean absolute difference between binned, standardized densities of \code{x1} and \code{x2} (ie, \code{sum(abs(x1 - x2))}) in environmental space.
#' \item \code{d} Schoener's \emph{D}.
#' \item \code{i} Warren's \emph{I}.
#' \item \code{esp} Godsoe's \emph{ESP}.
#' \item \code{rho} Correlation between binned, standardized densities of \code{x1} and \code{x2} in environmental space.
#' \item \code{rankCor}  Pearson rank correlation between binned, standardized densities of \code{x1} and \code{x2}.
#' }
#' @details This function replicates the procedure presented in Broennimann, O., Fitzpatrick, M.C., Pearman, P.B., Petitpierre, B., Pellissier, L., Yoccoz, N.G., Thuiller, W., Fortin, M-J., Randin, C., Zimmermann, N.E., Graham, C.H., and Guisan, A.  2012.  Measuring ecological niche overlap from occurrence and spatial environmental data.  Global Ecology and Biogeography 21:481-497.
#' @examples
#' data(lemur)
#' @export

nicheOverlap <- compiler::cmpfun(function(
	x1,
	x2,
	env,
	vars,
	bins = 100,
	cor = TRUE
) {
	
	### construct PCA
	if (class(env) != 'princomp') {
		env <- as.data.frame(env)
		env <- env[ , vars]
		pca <- princomp(env, cor=TRUE)
	} else {
		pca <- env
	}
	
	env <- pca$scores[ , 1:2]
	colnames(env) <- paste0('pc', 1:2)
	
	### parse input
	
	x1 <- as.data.frame(x1)
	x2 <- as.data.frame(x2)
	x1 <- x1[ , vars]
	x2 <- x2[ , vars]
	
	# convert environment to PC axes
	x1 <- predict(pca, x1)
	x2 <- predict(pca, x2)
	
	x1 <- x1[ , 1:2]
	x2 <- x2[ , 1:2]
	colnames(x1) <- paste0('pc', 1:2)
	colnames(x2) <- paste0('pc', 1:2)

	### construct kernel density estimator
	xlims <- range(x1[ , 1], x2[ , 1], env[ , 1])
	ylims <- range(x1[ , 2], x2[ , 2], env[ , 2])
	
	# add extra "empty" bin to accommodate difference in interpretation of bins in kde2d and hist2d functions
	xlims[1] <- xlims[1] - omnibus::eps()
	ylims[1] <- ylims[1] - omnibus::eps()
	
	xBinEdges <- seq.int(xlims[1L], xlims[2L], length.out=1L + bins[1L])
	yBinEdges <- seq.int(ylims[1L], ylims[2L], length.out=1L + bins[1L])
	
	xBinWidth <- xBinEdges[length(xBinEdges)] - xBinEdges[length(xBinEdges) - 1L]
	yBinWidth <- yBinEdges[length(yBinEdges)] - yBinEdges[length(yBinEdges) - 1L]
	xlims[2] <- xlims[2] + xBinWidth
	ylims[2] <- ylims[2] + yBinWidth
	
	kde <- MASS::kde2d(x=env[ , 1], y=env[ , 2], n=bins + 1, lims=c(xlims, ylims))
	
	# calculate frequency of occupancy
	envDens <- kde$z
	envDens <- envDens[1:(nrow(envDens) - 1), 1:(ncol(envDens) - 1)]
	envDens <- envDens / sum(envDens)
	
	breaks1 <- kde$x
	breaks2 <- kde$y
	
	x1hist <- statisfactory::hist2d(x1, breaks1=breaks1, breaks2=breaks2)
	x2hist <- statisfactory::hist2d(x2, breaks1=breaks1, breaks2=breaks2)
	
	freqOcc1 <- x1hist / envDens
	freqOcc2 <- x2hist / envDens

	freqOcc1 <- freqOcc1 / sum(freqOcc1)
	freqOcc2 <- freqOcc2 / sum(freqOcc2)
	
	freqOcc1 <- c(freqOcc1)
	freqOcc2 <- c(freqOcc2)
	
	# calculate niche overlap statistics
	out <- enmSdm::compareNiches(freqOcc1, freqOcc2)
	out
	
})
