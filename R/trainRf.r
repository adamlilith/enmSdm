#' Calibrate a random forest model
#'
#' This function trains a random forest model. It first finds the optimal value for \code{mtry} (number of variables sampled as candidates at each split). It then calls the function  \code{\link[randomForest]{randomForest}} from the \pkg{randomForest} package.
#' @param data Data frame.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Character. If "\code{binomial}" then the response is converted to a binary factor with levels 0 and 1. Otherwise, this argument has no effect.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) \emph{or} a numeric list of weights, one per class in \code{resp}. The default is to assign equal total weight to presences and contrast sites (\code{TRUE}).
#' @param verbose Logical. If \code{TRUE} then display progress for finding optimal value of \code{mtry}.
#' @param ... Arguments to pass to \code{\link[randomForest]{randomForest}}.
#' @return Object of class \code{\link[randomForest]{randomForest}}.
#' @seealso \code{\link[randomForest]{randomForest}}, \code{\link{trainCrf}}
#' @examples
#' \dontrun{
#' ### model red-bellied lemurs
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
#' occs <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
#' occsEnv <- raster::extract(clim, occs[ , c('longitude', 'latitude')])
#' 
#' # background sites
#' bg <- 2000 # too few cells to locate 10000 background points
#' bgSites <- dismo::randomPoints(clim, 2000)
#' bgEnv <- raster::extract(clim, bgSites)
#' 
#' # collate
#' presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
#' env <- rbind(occsEnv, bgEnv)
#' env <- cbind(presBg, env)
#' env <- as.data.frame(env)
#' 
#' preds <- paste0('bio', bios)
#' 
#' set.seed(123)
#'
#' # random forest
#' rf <- trainRf(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' )
#' 
#' # conditional random forest
#' crf <- trainCrf(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#' )
#' 
#' plot(rf)
#' 
#' # prediction rasters
#' mapRf1 <- predict(clim, rf, type='prob') # opposite class!
#' mapRf2 <- 1 - predict(clim, rf, type='prob') # correct
#' pointsFx <- function() points(occs[ , c('longitude', 'latitude')])
#' plot(stack(mapRf1, mapRf2), addfun=pointsFx)
#'
#' # CRFs are tricky...
#' }
#' @export

trainRf <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	w = TRUE,
	verbose = FALSE,
	...
) {

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# weights
	classwt <- if (family == 'binomial') {
		if ('logical' %in% class(w)) {
			if (w) {
				c(1, 1)
			} else if (!w) {
				NULL
			}
		} else {
			w
		}
	}
	
	# binomial response
	if (family == 'binomial') data[ , resp] <- factor(data[ , resp], levels=0:1)

	model <- randomForest::tuneRF(x=data[ , preds, drop=FALSE], y=data[ , resp], trace=FALSE, plot=FALSE, doBest=TRUE, strata=data[ , resp], classwt=classwt, ...)
	model

}

