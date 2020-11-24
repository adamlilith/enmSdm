#' Calculate niche overlap as per Broennimann et al. (2012)
#'
#' This function calculates niche overlap between two species.
#' @param x1 Data frame, matrix, or any object that can be coerced to a data frame containing environmental data at occurrence sites of a species.
#' @param x2 Data frame, matrix, or any object that can be coerced to a data frame containing environmental data at occurrence sites of another species.
#' @param env Either a data frame, matrix, or any object that can be coerced to a data frame containing environmental data at available background sites, \emph{or} an object of class \code{princomp} representing a principal components analysis generated using the \code{\link[stats]{princomp}} function with argument \code{scores = TRUE}.
#' @param vars Either a character list naming columns in \code{x1}, \code{x2}, and \code{x3} to be used as environmental data, \emph{or} positive integers indexing the columns to be used as environmental data.
#' @param bins Number of bins into which to divide the environmental space (default is 100 on each side).
#' @param cor Logical, if \code{TRUE} (default), then the PCA used to construct the environmental space will use the correlation matrix (this is highly recommended if the variables are on different scales). This is ignored if \code{env} is an object of class \code{princomp}.
#' @param densities Logical. If \code{TRUE}, then return not only metrics of niche similarity but also the density matrices (environment, species #1, and species #2).
#' @return If \code{densities} is \code{FALSE} (default), return a vector these named elements:
#' \itemize{
#' \item \code{meanDiff}: Mean difference between binned, standardized densities of \code{x1} and \code{x2} in environmental space.
#' \item \code{meanAbsDiff}: Mean absolute difference between binned, standardized densities of \code{x1} and \code{x2} (ie, \code{sum(abs(x1 - x2))}) in environmental space.
#' \item \code{rmsd}: Root mean squared difference.
#' \item \code{d}: Schoener's \emph{D}.
#' \item \code{i}: Warren's \emph{I}.
#' \item \code{esp}: Godsoe's \emph{ESP}.
#' \item \code{rho}: Correlation between binned, standardized densities of \code{x1} and \code{x2} in environmental space.
#' \item \code{rankCor}: Pearson rank correlation between binned, standardized densities of \code{x1} and \code{x2}.
#' }
#' If \code{densities} is \code{TRUE}, then return a list with a vector of metrics of niche overlap as above, plus three matrices:
#' \itemize{
#' \item \code{environDens}: Density of available environment.
#' \item \code{x1density}: Density of species #1 in environmental space, normalized to sum to 1 but not normalized by available environment.
#' \item \code{x2density}: Density of species #2 in environmental space, normalized to sum to 1 but not normalized by available environment.
#' }
#' @references This function replicates the procedure presented in Broennimann, O., Fitzpatrick, M.C., Pearman, P.B., Petitpierre, B., Pellissier, L., Yoccoz, N.G., Thuiller, W., Fortin, M-J., Randin, C., Zimmermann, N.E., Graham, C.H., and Guisan, A.  2012.  Measuring ecological niche overlap from occurrence and spatial environmental data.  Global Ecology and Biogeography 21:481-497.
#' @seealso \code{\link{compareNiches}}
#' @examples
#' # comparing niches between the common brown leumr (Eulemur fulvus)
#' # and the red-bellied lemur (Eulemur rubriventer)
#' 
#' data(mad0)
#' data(lemurs)
#' 
#' # climate data
#' bios <- c(1, 5, 12, 15)
#' clim <- raster::getData('worldclim', var='bio', res=10)
#' clim <- raster::subset(clim, bios)
#' clim <- raster::crop(clim, mad0)
#' 
#' # occurrence data
#' occs1 <- lemurs[lemurs$species == 'Eulemur fulvus', ]
#' occs2 <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' 
#' ll <- c('longitude', 'latitude')
#' plot(mad0)
#' points(occs1[ , ll])
#' points(occs2[ , ll], col='red', pch=3)
#' 
#' occsEnv1 <- raster::extract(clim, occs1[ , ll])
#' occsEnv2 <- raster::extract(clim, occs2[ , ll])
#' 
#' # background sites
#' bg <- 2000 # too few cells to locate 10000 background points
#' bgSites <- dismo::randomPoints(clim, 2000)
#' bgEnv <- extract(clim, bgSites)
#' 
#' vars <- paste0('bio', bios)
#' nicheOverlap(occsEnv1, occsEnv2, env=bgEnv, vars=vars)
#' 
#' @export

nicheOverlap <- compiler::cmpfun(function(
	x1,
	x2,
	env,
	vars,
	bins = 100,
	cor = TRUE,
	densities = FALSE
) {

	### construct PCA
	if (!('princomp' %in% class(env))) {
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

	### get limits of environmental space
	pc1lims <- range(x1[ , 1], x2[ , 1], env[ , 1])
	pc2lims <- range(x1[ , 2], x2[ , 2], env[ , 2])

	# add an extra "half" bin to each side
	pc1range <- diff(pc1lims)
	pc2range <- diff(pc2lims)
	
	pc1binWidth <- 0.5 * 0.01 * pc1range
	pc2binWidth <- 0.5 * 0.01 * pc2range
	
	pc1inc <- pc1binWidth / 2
	pc2inc <- pc2binWidth / 2
	
	pc1lims[1] <- pc1lims[1] - pc1inc
	pc2lims[1] <- pc2lims[1] - pc2inc
	
	pc1lims[2] <- pc1lims[2] + pc1inc
	pc2lims[2] <- pc2lims[2] + pc2inc
	
	### construct kernel density estimator
	kdeEnv <- MASS::kde2d(x=env[ , 1], y=env[ , 2], n=bins, lims=c(pc1lims, pc2lims))
	x1kde <- MASS::kde2d(x=x1[ , 1], y=x1[ , 2], n=bins, lims=c(pc1lims, pc2lims))
	x2kde <- MASS::kde2d(x=x2[ , 1], y=x2[ , 2], n=bins, lims=c(pc1lims, pc2lims))

	### calculate frequency of occupancy
	envDens <- kdeEnv$z
	x1dens <- x1kde$z
	x2dens <- x2kde$z
	
	envDens <- envDens / sum(envDens)
	x1dens <- x1dens / sum(x1dens)
	x2dens <- x2dens / sum(x2dens)

	freqOcc1 <- x1dens / envDens
	freqOcc2 <- x2dens / envDens

	# calculate niche overlap statistics
	out <- compareNiches(freqOcc1, freqOcc2)
	
	if (densities) {
		out <- list(similarity=out, environDensity=envDens, x1density=x1dens, x2density=x2dens)
	}
	out

})
