#' Interpolate a stack of rasters
#'
#' This function returns a stack of rasters interpolated from a stack of rasters. For example, the input might represent rasters of a process measured at times t, t + 1, and t + 4. The rasters at t + 2 and t + 3 could be interpolated based on the values in the other rasters. Note that this function can take a lot of time and memory, even for relatively small rasters.
#' @param rasts A "stack" of \code{raster}s.
#' @param interpFrom Numeric vector, one value per raster in \code{rasts}. Values represent "distance" along the set of rasters rasters (e.g., time).
#' @param interpTo Numeric vector, values of "distances" at which to interpolate the rasters.
#' @param type Character. The type of model used to do the interpolation. Note that some of these (the first few) are guaranteed to go through every point being interpolated from. The second set, however, are effectively regressions so are not guaranteed to do through \emph{any} of the points. Note that some methods cannot handle cases where at least some series of cells have < a given number of non-\code{NA} values (e.g., smooth splines will not work if there are < 4 cells with non-\code{NA} values).
#' \itemize{
#' \item \code{linear}: A model based on linear segments "fastened" at each value of \code{interpFrom}. The segments will intersect each value being interpolated from.
#' \item \code{spline}: A natural splines-based model. Splines will intersect each value being interpolated from.
#' \item \code{gam}: A generalized additive model. Note that the GAM is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link[mgcv]{gam}} can be supplied via \code{...}. Especially note the \code{family} argument! You can use the \code{onFail} argument with this method since in some cases \code{\link[mgcv]{gam}} if there are too few data points.
#' \item \code{glm}: A generalized linear model. Note that the GLM is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link[mgcv]{gam}} can be supplied via \code{...}. Especially note the \code{family} argument (the main reason for why you would use a GLM versus just linear interpolation)! You can use the \code{onFail} argument with this method since in some cases \code{\link[stats]{glm}} if there are too few data points.
#' \item \code{ns}: A natural splines model. Note that the NS is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link{trainNs}} can be supplied via \code{...}. Especially note the \code{family} argument and the \code{df} argument! If \code{df} is not supplied, then the number of splines attempted will be equal to \code{1:(length(interpFrom) - 1)}. You can use the \code{onFail} argument with this method.
#' \item \code{poly}: A polynomial model. This method constructs an \emph{n}-degree polynomial where \emph{n} = \code{length(interpFrom) - 1}. The most parsimonious model is then selected from all possible subsets of models (including an intercept-only model) using AICc. This method is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link{glm}} can be supplied via \code{...}. Especially note the \code{family} argument! If \code{family} is not supplied, then the response is assumed to have a Gaussian distribution. You can use the \code{onFail} argument with this method.
#' \item \code{bs}: A basis-spline model. This method constructs a series of models with \emph{n}-degree basis-spline model where \emph{n} ranges from 3 to \code{length(interpFrom) - 1}. The most parsimonious model is then selected from all possible subsets of models (including an intercept-only model) using AICc. This method is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link{glm}} can be supplied via \code{...}. Especially note the \code{family} argument! If \code{family} is not supplied, then the response is assumed to have a Gaussian distribution. You can use the \code{onFail} argument with this method.
#' \item \code{smooth.spline}: A smooth-spline model (see \code{\link{smooth.spline}}). This method is \emph{not} guaranteed to intersect each value being interpolated from. Arguments to \code{\link{smooth.spline}} can be supplied via \code{...}. Unlike some other methods, a \code{family} cannot be specified (Gaussian is assumed)! You can use the \code{onFail} argument with this method.
#' }
#' @param onFail Either \code{NA} (default) or any one of \code{'linear'}, \code{'spline'}, or \code{'poly'}. If a method specified by \code{type} fails (i.e., because there are fewer than the required number of values to interpolate from), this method is used in its place. If this is \code{NA} and the method fails, then an error occurs.
#' @param useRasts Logical. If \code{FALSE} (default), then the calculations are done using arrays. This can be substantially faster than using rasters (when \code{useRasts = TRUE}), but also run into memory issues.
#' @param na.rm Logical, if \code{TRUE} (default) then ignore cases where all values in the same cells across rasters from which interpolations are made are \code{NA} (i.e., do not throw an error). If \code{FALSE}, then throw an error when this occurs.
#' @param verbose Logical. If \code{TRUE} (default), display progress.
#' @param ... Other arguments passed to \code{approx} or \code{spline} (\emph{do not} include any of these arguments: \code{x}, \code{y}, or \code{xout}), or to \code{\link{glm}}, \code{\link[mgcv]{gam}}, or \code{\link{smooth.spline}}.
#' @return A raster stack with one layer per element in \code{interpTo}.
#' @details This function can be very memory-intensive for large rasters.  It may speed things up (and make them possible) to do interpolations piece by piece (e.g., instead of interpolating between times t0, t1, t2, t3, ..., interpolate between t0 and t1, then t1 and t2, etc. This may give results that differ from using the entire set, however. Note that using linear and splines will often yield very similar results except that in a small number of cases splines may produce very extreme interpolated values.
#' @seealso \code{\link[raster]{approxNA}}, \code{\link[stats]{approxfun}}, \code{\link[stats]{splinefun}}, \code{\link{trainNs}}, \code{\link{glm}}, , \code{\link[splines]{bs}}, \code{\link{smooth.spline}}.
#' @examples
#' \dontrun{
#' interpFrom <- c(1, 3, 4, 8, 10, 11, 15)
#' interpTo <- 1:15
#' rx <- rast(nrows=10, ncols=10)
#' r1 <- setValues(rx, rnorm(100, 1))
#' r3 <- setValues(rx, rnorm(100, 3))
#' r4 <- setValues(rx, rnorm(100, 5))
#' r8 <- setValues(rx, rnorm(100, 11))
#' r10 <- setValues(rx, rnorm(100, 3))
#' r11 <- setValues(rx, rnorm(100, 5))
#' r15 <- setValues(rx, rnorm(100, 13))
#' rasts <- c(r1, r3, r4, r8, r10, r11, r15)
#' names(rasts) <- paste0('rasts', interpFrom)
#' 
#' linear <- interpolateRasters(rasts, interpFrom, interpTo)
#' spline <- interpolateRasters(rasts, interpFrom, interpTo, type='spline')
#' gam <- interpolateRasters(rasts, interpFrom, interpTo, type='gam', onFail='linear')
#' ns <- interpolateRasters(rasts, interpFrom, interpTo, type='ns', onFail='linear', verbose=FALSE)
#' poly <- interpolateRasters(rasts, interpFrom, interpTo, type='poly', onFail='linear')
#' bs <- interpolateRasters(rasts, interpFrom, interpTo, type='bs', onFail='linear')
#' ss <- interpolateRasters(rasts, interpFrom, interpTo, type='smooth.spline', onFail='linear',
#' verbose=FALSE)
#' 
#' # examine trends for a particular point on the landscape
#' pts <- rbind(c(-9, 13))
#' linearExt <- unlist(raster::extract(linear, pts))
#' splineExt <- unlist(raster::extract(spline, pts))
#' gamExt <- unlist(raster::extract(gam, pts))
#' nsExt <- unlist(raster::extract(ns, pts))
#' polyExt <- unlist(raster::extract(poly, pts))
#' bsExt <- unlist(raster::extract(bs, pts))
#' ssExt <- unlist(raster::extract(ss, pts))
#' 
#' mins <- min(linearExt, splineExt, gamExt, nsExt, polyExt, bsExt, ssExt)
#' maxs <- max(linearExt, splineExt, gamExt, nsExt, polyExt, bsExt, ssExt)
#' 
#' plot(interpTo, linearExt, type='l', lwd=2, ylim=c(mins, maxs), ylab='Value')
#' lines(interpTo, splineExt, col='blue')
#' lines(interpTo, gamExt, col='green')
#' lines(interpTo, nsExt, col='orange')
#' lines(interpTo, polyExt, col='gray')
#' lines(interpTo, bsExt, col='magenta')
#' lines(interpTo, ssExt, col='cyan')
#' 
#' ext <- c(extract(rasts, pts))
#' points(interpFrom, ext)
#' 
#' legend('topleft', inset=0.01, lty=c(rep(1, 7), NA),
#' legend=c('linear', 'spline', 'GAM', 'NS', 'polynomial', 'B-spline',
#' 'Smooth spline', 'Observed'), col=c('black', 'blue', 'green',
#' 'orange', 'gray', 'magenta', 'cyan'), pch=c(rep(NA, 7), 1))
#'
#' }
#' @export

interpolateRasters <- function(
	rasts,
	interpFrom,
	interpTo,
	type = 'linear',
	onFail = NA,
	useRasts = FALSE,
	na.rm = TRUE,
	verbose = TRUE,
	...
) {

	### check for errors
	####################
		
		if (!any(c('linear', 'spline', 'gam', 'ns', 'poly', 'bs', 'smooth.spline') %in% type)) stop('Argument "type" is not a valid value.')
		if (raster::nlayers(rasts) < 2) stop('Argument "rasts" must have >1 raster layer.')
		if (length(interpFrom) != raster::nlayers(rasts)) stop('Argument "interpFrom" must have same length as number of rasters in argument "rasts".')
		
	### reserve blank array for output
	##################################

		rows <- nrow(rasts)
		cols <- ncol(rasts)
		numInterps <- length(interpTo)
		
		xmin <- raster::xmin(rasts)
		xmax <- raster::xmax(rasts)
		ymin <- raster::ymin(rasts)
		ymax <- raster::ymax(rasts)
		proj4 <- raster::projection(rasts)

		thisOut <- array(NA, dim=c(rows, cols, numInterps))
		if (useRasts) {
			out <- raster::raster(thisOut, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj4)
			rm(thisOut); gc()
		}
		
	### interpolate each cell
	#########################
		
		if (verbose) progress <- utils::txtProgressBar(min=1, max=raster::ncell(rasts), style=3, width=min(64, getOption('width')))

		countCell <- 1
		for (countRow in 1:nrow(rasts)) {
		
			for (countCol in 1:ncol(rasts)) {
			
				y <- if (useRasts) {
					rasts[countCell]
				} else {
					rasts[countRow, countCol, ]
				}
				
				y <- c(y)
					
				if (sum(!is.na(y)) < 2) {
					cellInterpol <- rep(NA, numInterps)
				} else if (type == 'linear') {
					cellInterpol <- stats::approx(x=interpFrom, y=y, xout=interpTo, ...)$y
					# cellInterpol <- stats::approx(x=interpFrom, y=y, xout=interpTo)$y
				} else if (type == 'spline') {
					# stats::spline(x=interpFrom, y=y, xout=interpTo, ...)$y
					cellInterpol <- stats::spline(x=interpFrom, y=y, xout=interpTo, method='natural', ...)$y
					# stats::spline(x=interpFrom, y=y, xout=interpTo, method='natural')$y
				} else {
					
					dataFrom <- data.frame(y=y, x=interpFrom)
					dataFrom <- dataFrom[stats::complete.cases(dataFrom), ]
				
					dataTo <- data.frame(x=interpTo)
					
					dots <- list(...)
					args <- c(list(data=dataFrom), dots)
					# args <- list(data=dataFrom)
					
					if (type == 'gam') {

						k <- length(y)
						form <- y ~ s(x, bs='cs', k=k)
						thisArgs <- c(args, formula=form)
						interpModel <- tryCatch(do.call(mgcv::gam, thisArgs), error=function(err) FALSE)
						
					} else if (type == 'glm') {

						form <- y ~ x
						thisArgs <- c(args, formula=form)
						interpModel <- tryCatch(do.call(stats::glm, thisArgs), error=function(err) FALSE)
						
					} else if (type == 'ns') {
					
						thisArgs <- c(args, resp='y', preds='x')
						if (!('df' %in% names(dots))) {
							df <- 1:(length(interpFrom) - 1)
							df <- list(df)
							thisArgs <- c(thisArgs, df=df)
						}

						if (!('family' %in% names(dots))) {
							family <- 'gaussian'
							thisArgs <- c(thisArgs, family=family)
						}

						interpModel <- tryCatch(do.call(trainNs, thisArgs), error=function(err) FALSE)

						
					} else if (type == 'poly') {

						degrees <- 1:nrow(dataFrom)
						
						models <- list()
						k <- integer()
						for (degree in degrees) {
						
							form <- paste0('y ~ poly(x, degree=', degree, ')')
							form <- stats::as.formula(form)
							thisArgs <- c(args, formula=form)
							tryModel <- tryCatch(do.call('glm', thisArgs), error=function(err) FALSE)
							if (!is.logical(tryModel)) {
								models[[length(models) + 1]] <- tryModel
								k <- c(k, length(tryModel$coefficients))
							}
						
						}
						
						form <- y ~ 1
						thisArgs <- c(args, form=form)
						model <- do.call('glm', args=thisArgs)
						models[[length(models) + 1]] <- model
						k <- c(k, length(model$coefficients))
						
						aic <- sapply(models, MuMIn::AICc)
						aicc <- aic + (2 * k^2 + 2 * k) / (nrow(dataFrom) - k - 1)
						interpModel <- models[[which.min(aicc)]]

					} else if (type == 'bs') {

						degrees <- 3:max(3, nrow(dataFrom))

						# try all possible models
						models <- list()
						k <- integer()
						for (i in seq_along(degrees)) {
						
							degree <- degrees[i]
						
							form <- paste0('y ~ splines::bs(x, df=', degree, ')')
							form <- stats::as.formula(form)
							thisArgs <- c(args, formula=form)
							tryModel <- tryCatch(do.call('glm', thisArgs), error=function(err) FALSE)
							if (!is.logical(tryModel)) {
								models[[length(models) + 1]] <- tryModel
								k <- c(k, length(tryModel$coefficients))
							}
						
						}
						
						form <- y ~ 1
						thisArgs <- c(args, form=form)
						model <- do.call('glm', args=thisArgs)
						models[[length(models) + 1]] <- model
						k <- c(k, length(model$coefficients))
						
						aic <- sapply(models, MuMIn::AICc)
						aicc <- aic + (2 * k^2 + 2 * k) / (nrow(dataFrom) - k - 1)
						interpModel <- models[[which.min(aicc)]]

					} else if (type == 'smooth.spline') {

						if (nrow(dataFrom) < 4 & is.na(onFail)) stop('At least one set of cells has fewer than 4 unique values. Cannot use interpolation function "smooth.spline".')
					
						# interpModel <- tryCatch(stats::smooth.spline(x=dataFrom$x, y=dataFrom$y, keep.data=FALSE, cv=TRUE, control.spar=list(trace=FALSE), ...), error=function(err) FALSE)
						interpModel <- tryCatch(stats::smooth.spline(x=dataFrom$x, y=dataFrom$y, keep.data=FALSE, cv=TRUE, control.spar=list(trace=FALSE)), error=function(err) FALSE)

					}

					### if model succeeded, remember values
					#######################################
					
					if (!is.logical(interpModel)) {
					
						cellInterpol <- do.call('predict', args=list(object=interpModel, newdata=dataTo, type='response'))
						if (type == 'smooth.spline') cellInterpol <- cellInterpol$y$x
					
					### if model failed, try simpler models
					} else if (!is.na(onFail)) {
					
						if (onFail == 'linear') {
						
							cellInterpol <- stats::approx(x=interpFrom, y=y, xout=interpTo)$y
							
						} else if (onFail == 'spline') {

							cellInterpol <- stats::spline(x=interpFrom, y=y, xout=interpTo, method='natural')$y
							
						} else if (onFail == 'poly') {
						
							degrees <- 1:nrow(dataFrom)
							
							models <- list()
							k <- integer()
							for (degree in degrees) {
							
								form <- paste0('y ~ poly(x, degree=', degree, ')')
								form <- stats::as.formula(form)
								thisArgs <- c(args, formula=form)
								tryModel <- tryCatch(do.call('glm', thisArgs), error=function(err) FALSE)
								if (!is.logical(tryModel)) {
									models[[length(models) + 1]] <- tryModel
									k <- c(k, length(tryModel$coefficients))
								}
							
							}
							
							form <- y ~ 1
							thisArgs <- c(args, form=form)
							model <- do.call('glm', args=thisArgs)
							models[[length(models) + 1]] <- model
							k <- c(k, length(model$coefficients))
							
							aic <- sapply(models, MuMIn::AICc)
							aicc <- aic + (2 * k^2 + 2 * k) / (nrow(dataFrom) - k - 1)
							interpModel <- models[[which.min(aicc)]]
							
						}
					
					# interpolation method failed
					} else {
					
						stop('Interpolation method failed. You could try using the "onFail" argument.')
						
					}
					
				}
				
				if (useRasts) {
					out[countCell] <- cellInterpol
				} else {
					thisOut[countRow, countCol, ] <- cellInterpol
				}
				
				countCell <- countCell + 1
				
				if (verbose) utils::setTxtProgressBar(progress, value=countCell)
				
			} # next column

		} # next row

		if (verbose) close(progress)

	### reconfigure back to raster format
	#####################################

		if (!useRasts) {
			if (verbose) omnibus::say('Compiling rasters...')
			out <- raster::raster(thisOut[ , , 1], xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj4)
			dims <- dim(thisOut)[3]
			if (dims > 1) {
				for (i in 2:dims) {
					out <- stack(
						out,
						raster::raster(thisOut[ , , i], xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj4)
					)
				}
			}
		}

		interpToNames <- gsub(interpTo, pattern='-', replacement='Neg')
		interpToNames <- gsub(interpToNames, pattern='\\.', replacement='p')
		names(out) <- paste0('interpTo', interpToNames)
		out
	
}
