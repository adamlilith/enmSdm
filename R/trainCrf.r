#' Calibrate a conditional random forest model
#'
#' This function trains a  conditional random forest model. It is nearly identical to the \code{\link[party]{cforest}} function in the \pkg{party} package but is included for consistency with \code{\link{trainGlm}}, \code{\link{trainGam}}, and similar functions.
#' @param data Data frame.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}). Default is to use the 'binomial' family.
#' @param w Either logical in which case TRUE causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param ... Arguments to pass to \code{\link[party]{cforest}} function.
#' @return Object of class \code{RandomForest}.
#' @seealso \code{\link[party]{cforest}}, \code{\link[enmSdm]{trainRf}}
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

trainCrf <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	w = ifelse(family == 'binomial', TRUE, FALSE),
	...
) {

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# model weights
	if (length(w) == 1 && class(w) == 'logical') {
		w <- if (w & family == 'binomial') {
			c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
		} else {
			rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w]
	}

	# get just desired columns
	data <- data[ , c(resp, preds)]

	# binomial response
	if (family == 'binomial') {
		data[ , resp] <- if (data[1, resp] == 0) {
			factor(data[ , resp], levels=0:1)
		} else if (data[ , resp] == 1) {
			factor(data[ , resp], levels=1:0)
		}
	}

	# formula
	form <- stats::as.formula(paste(resp, '~ .'))

	# train model
	model <- party::cforest(
		formula=form,
		data=data,
		weights=w
	)

	model

}
