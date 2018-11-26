#' Automatically calculate tolerance for null model point randomization
#'
#' This function takes a triangular or rectangular distance matrix representing pairwise distances between points and calculates a tolerance value for spatial randomizations of points meant to reflect the spatial autocorrelation inherent in the observed points. It is meant as a helper function for the family of \code{randPoints---} functions.
#' @param dists Either a lower-triangular square matrix (with the diagonal set to \code{NA}) or a rectangular matrix. Values represent pairwise distances between points.
#' @param removeTopRow Logical, if \code{TRUE} then remove the top row of the \code{dists} matrix before calculating summary statistics. If the matrix is lower-triangular and the diagonal is \code{NA} then retaining the top row will cause a warning. This value should be \code{FALSE} if \code{dists} is a rectangular matrix or otherwise not triangular!
#' @details Tolerance is calculated as half of the value of minimum pairwise distances that are >0.
#' @return Numeric.
#' @keywords internal
.calcTol <- function(dists, removeTopRow) {

	if (removeTopRow) dists <- dists[2:nrow(dists), ]
	
	minDists <- apply(dists, 1, min, na.rm=TRUE)
	minPosDists <- minDists[minDists > 0]
	tol <- 0.5 * min(minPosDists)
	tol
}
