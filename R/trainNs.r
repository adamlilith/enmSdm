#' Calibrate a natural splines model
#'
#' This function constructs a natural-spline model piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param df Integer > 0 \emph{or} vector of integers > 0. Sets flexibility of model fit. See documentation for \code{\link[splines]{ns}}.  If \code{construct} is \code{TRUE}, then univariate models for each term will be evaluated using each value in \code{df}. Note that \code{NULL} is also valid, but it can create problems when used with other functions in this package (and usually defaults to \code{df=3} anyway).
#' @param construct Logical. If TRUE then construct model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{initialTerms}.
#' @param select Logical. If TRUE then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is TRUE.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{\link[MuMIn]{dredge}} is 31, so if this number is >31 and \code{select} is \code{TRUE} then it is forced to 31 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character or character vector. Indicates type of value returned. Values can be \code{'model'} (default; return model with lowest AICc), \code{'models'} (return a list of all models), and/or \code{'tuning'} (return a data frame with AICc for each model). If more than one value is specified, then the output will be a list with elements named "model", "models", and/or "tuning". If \code{'models'} is specified, they will only be produced if \code{select = TRUE}. The models will appear in the list in same order as they appear in the tuning table (i.e., model with the lowest AICc first, second-lowest next, etc.). If just one value is specified, the output will be either an object of class \code{glm}, a list with objects of class \code{glm}, or a data frame.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device. Default is \code{FALSE}.
#' @param ... Arguments to send to \code{gam()} or \code{dredge()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{gam}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gam} object and the data frame.
#' @seealso \code{\link[splines]{ns}}, \code{\link[mgcv]{gam}}, \code{\link{trainGam}}
#' @examples
#' \donttest{
#' library(brglm2)
#'
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
#' bgEnv <- extract(clim, bgSites)
#' 
#' # collate
#' presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
#' env <- rbind(occsEnv, bgEnv)
#' env <- cbind(presBg, env)
#' env <- as.data.frame(env)
#' 
#' preds <- paste0('bio', bios)
#' 
#' # GLM
#' gl <- trainGlm(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # GAM
#' ga <- trainGam(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # NS
#' ns <- trainNs(
#' 	data = env,
#' 	resp = 'presBg',
#' 	preds = preds,
#'  verbose = TRUE
#' )
#' 
#' # prediction rasters
#' mapGlm <- predict(clim, gl, type='response')
#' mapGam <- predict(clim, ga, type='response')
#' mapNs <- predict(clim, ga, type='response')
#'
#' par(mfrow=c(1, 3))
#' plot(mapGlm, main='GLM')
#' plot(mad0, add=TRUE)
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapGam, main='GAM')
#' plot(mad0, add=TRUE)
#' points(occs[ , c('longitude', 'latitude')])
#' plot(mapNs, main='NS')
#' plot(mad0, add=TRUE)
#' points(occs[ , c('longitude', 'latitude')])
#' }
#' @export

trainNs <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	df = 1:3,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 20,
	initialTerms = 8,
	w = TRUE,
	out = 'model',
	verbose = FALSE,
	...
) {

	ellipses <- list(...)

	# force number of starting terms to 31 or less
	if (select & initialTerms > 31) {
		initialTerms <- 31
		warning('initialTerms must be 31 or less. Forcing to 31.')
	}

	# degrees of freedom
	if (is.null(df)) df <- 'NULL'

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	#############
	## weights ##
	#############

	# model weights
	if (is.logical(w)) {
		if (w[1L] & (family %in% c('binomial', 'quasibinomial'))) {
			posCases <- sum(data[ , resp, drop=TRUE] == 1)
			negCases <- sum(data[ , resp, drop=TRUE] == 0)
			w <- c(rep(1, posCases), rep(posCases / (posCases + negCases), negCases))
		} else {
			w <- rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w, drop=TRUE]
	}
	w <<- w / max(w) # declare to global because dredge() and pdredge() have problems if it is not

	################################
	## initial model construction ##
	################################

	# create starting formula
	form <- paste0(resp, ' ~ 1')

	if (construct) {

		### SINGLE-variable terms
		for (thisPred in preds) { # for each predictor test single-variable terms

			for (thisDf in df) {
			
				term <- if (class(data[ , thisPred]) != 'factor') {
					paste0('splines::ns(', thisPred, ', df=', thisDf, ')')
				} else {
					thisPred
				}

				thisThisForm <- paste0(form, ' + ', term)

				thisModel <- stats::glm(stats::as.formula(thisThisForm), family=family, data=data, weights=w, method='brglmFit', ...)
				thisAic <- AIC(thisModel)

				# remember
				tuning <- if (exists('tuning', inherits=FALSE)) {
					rbind(tuning, data.frame(term=term, AIC=thisAic, df=thisDf))
				} else {
					data.frame(term=term, AIC=thisAic, df=thisDf)
				}
				
			}

		} # next single-variable term

		# sort by AIC
		tuning <- tuning[order(tuning$AIC), ]

		# print AICc frame
		if (verbose) {

			omnibus::say('Model construction results for each term tested:')
			print(tuning)
			omnibus::say('')

		}

		## construct final model
		form <- paste0(form, ' + ', tuning$term[1]) # add first term

		# for each set of presences > min num required, add a term
		if (floor(sum(data[ , resp]) / presPerTermInitial ) - 1 > 1 & initialTerms > 1) {

			termsToAdd <- 2:max(2, min(initialTerms, c(floor(sum(data[ , resp]) / presPerTermInitial ) - 1, nrow(tuning) - 1)))

			form <- paste0(form, ' + ', paste0(tuning$term[termsToAdd], collapse=' + '))

		} # if there are sufficient presences for additional terms beyond first

	# NO AUTOMATED MODEL CONSTRUCTION
	# use all single-variable terms and two-variable terms
	} else {

		if (length(df) > 1) warning('Multiple values of "df" assigned. Using the first one.')
	
		# single terms
		for (thisPred in preds) {

			if (class(data[ , thisPred]) != 'factor') {
				form <- paste0(form, ' + splines::ns(', thisPred, ', df=', df[1], ')')
			} else {
				form <- paste0(form, ' + ', thisPred)
			}

		}

	} # if not doing automated model construction

	###########################################################################
	## train model ############################################################
	## while model hasn't converged and while gamma is <= tolerance value... ##
	###########################################################################

	# get GAM model... using automated scale selection with weights so influence of absences equals influence of presences... using tryCatch because sometimes for variables with too little variation the default df of the basis is too high, in which case it is reduced and attempted again (for univariate and bivariate models only)
	model <- stats::glm(stats::as.formula(form), family=family, data=data, weights=w, na.action='na.fail', method='brglmFit', ...)

	if (verbose) {

		omnibus::say('Starting full model:')
		print(summary(model))
		omnibus::say('')

	}

	#######################################################################################
	## if doing model construction, evaluate all possible models using AIC then get best ##
	#######################################################################################

	if (select) {

		if (verbose) omnibus::say('Calculating AICc across all possible models...')

		# calculate all possible models and rank by AIC
		lims <- c(0, max(1, min(c(floor(sum(data[ , resp]) / presPerTermFinal), initialTerms, nrow(tuning)))))

		tuningModels <- MuMIn::dredge(
			global.model=model,
			rank='AICc',
			m.lim=lims,
			trace=FALSE,
			...
		)

		# get model with best AIC
		model <- MuMIn::get.models(tuningModels, subset = 1)[[1]]
		if ('models' %in% out) models <- MuMIn::get.models(tuningModels, subset=TRUE)

	} # if model selection

	# return
	if (length(out) > 1) {
		output <- list()
		if ('models' %in% out) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuningModels
		output
	} else if (out == 'models') {
		models
	} else if (out == 'model') {
		model
	} else if (out == 'tuning') {
		tuningModels
	}
	
}
