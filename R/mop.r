#' Identify sites in two sets based on mobility-oriented parity (MOP).
#'
#' This function takes two data frames or matrices that represent the environment of sites measured in univariate or multivariate space and returns the \emph{x}-th percentile of sites in each set that are closest to each other in Euclidean space. Note that in many cases it is advisable to first transform the raw environmental values using, for example, PCA.  MOP was first presented formally in Owens, H.L., Campbell, L.P., Dornak, L.L., Saupe, E.E., Barve, N., Soberón, Ingenloff, K., Lira-Noriega, A., Hensz, C., Myers, C.E., and Peterson, A.T.  2013.  Constraints on interpretation of ecological niche models by limited environmental ranges on calibration area.  \emph{Ecological Modeling} 263:10-18.
#' @param set1 Data frame or matrix one or more columns wide.
#' @param set2 Data frame or matrix one or more columns wide.
#' @param p Numeric value(s) in the range [0, 1]. The \emph{p}-th percentile of sites in \code{set1} and \code{set2} that are closest to one another are returned.  Note that if \code{p} = 1 then all sites in \code{set1} and \code{set2} are returned.
#' @param index Logical, if \code{TRUE} then return the indices of the rows in \code{set1} and \code{set2} that correspond to each value of \code{p}. If \code{FALSE} then return data frames or matrices (depending on the class of \code{set1} and \code{set2}).
#' @param na.rm Logical, if \code{TRUE} then any rows in \code{set1} or \code{set2} with at least one \code{NA} are removed first.
#' @return List with three elements. The first two elements correspond to \code{set1} and \code{set2}. Each of these elements is a list the same length of \code{p}, with each data frame/matrix coresponding to a value of \code{p}. The third element is a matrix of statistics reporting the statistics pertaining to the environmental distances between each subset of \code{set1} and \code{set2}.
#' @seealso \code{\link[dismo]{mess}}
#' @examples
#' set1 <- data.frame(x1=1:20, x2=round(100 * rnorm(20)))
#' set2 <- data.frame(x1=sample(1:30, 30), x2=sort(round(100 * rnorm(30))))
#' # return data frames that are subsets of set1 and set2
#' out <- mop(set1, set2, p=c(0.1, 0.5))
#' out
#' # return indices of subsets of set1 and set2
#' out <- mop(set1, set2, p=c(0.1, 0.5), index=TRUE)
#' out
#' @export

mop <- function(
	set1,
	set2,
	p,
	index = FALSE,
	na.rm = FALSE
) {

	# calculate Euclidean distances
	dists <- matrix(NA, nrow=nrow(set1), ncol=nrow(set2))
	for (i in 1:nrow(set1)) {
		for (j in 1:nrow(set2)) {
			dists[i, j] <-  sqrt(rowSums((set1[i, , drop=FALSE] - set2[j, , drop=FALSE])^2))
		}
	}

	# find nearest set of points in each set
	near1 <- apply(dists, 1, min, na.rm=na.rm)
	near2 <- apply(dists, 2, min, na.rm=na.rm)
	
	nearRank1 <- rank(near1)
	nearRank2 <- rank(near2)
	
	out <- list()
	out$set1 <- list()
	out$set2 <- list()
	
	# get nearest set for each value of p
	for (i in seq_along(p)) {
	
		thisNearest1 <- which(nearRank1 <= round(p[i] * nrow(set1)))
		thisNearest2 <- which(nearRank2 <= round(p[i] * nrow(set2)))
		out$set1[[i]] <- if (index) { thisNearest1 } else { set1[thisNearest1, , drop=FALSE] }
		out$set2[[i]] <- if (index) { thisNearest2 } else { set2[thisNearest2, , drop=FALSE] }
		attr(out$set1[[i]], 'p') <- p[i]
		attr(out$set2[[i]], 'p') <- p[i]
		
		theseDists <- c(
			apply(dists[thisNearest1, , drop=FALSE], 1, mean),
			apply(dists[ , thisNearest2, drop=FALSE], 2, mean)
		)
		
		thisStats <- matrix(c(p[i], min(theseDists), quantile(theseDists, 0.025), mean(theseDists), median(theseDists), quantile(theseDists, 0.975), max(theseDists), sd(theseDists), length(thisNearest1), length(thisNearest2)), ncol=10)
		if (exists('stats', where=out, inherits=FALSE)) {
			out$stats <- rbind(out$stats, thisStats)
		} else {
			out$stats <- thisStats
		}
		
	}

	colnames(out$stats) <- c('p', 'minDist', 'quantile0p025', 'meanDist', 'medianDist', 'quantile0p975', 'maxDist', 'sdDist', 'numSitesSet1', 'numSitesSet2')
	out
	
}
