#' Calibrate a generalized linear model (GLM)
#'
#' This function constructs a GLM piece-by-piece by first calculating AICc for all models with univariate, quadratic, and 2-way-interaction terms. It then creates a "full" model with the highest-ranked uni/bivariate terms. Finally, it implements an all-subsets model selection routine using AICc. Its output is any or all of: a table with AICc for all possible models, all possible models (after model construction), and/or the model with the lowest AICc.
#' @param data Data frame. Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}). Default is to use the 'binomial' family.
#' @param tooBig Numeric. Used to catch errors when fitting a model fit with the \code{brglmFit} function in the \pkg{brglm2} package. In some cases fitted coefficients are unstable and tend toward very high values, even if training data is standardized. Models with such coefficients will be discarded if any one coefficient is \code{> tooBig}. Set equal to \code{Inf} to keep all models.
#' @param anyway Logical. If \code{FALSE} (default), then during model construction, if no univariate models have valid coefficients (< \code{tooBog}), then do not proceed and return \code{NULL}. If \code{TRUE}, then proceed with instable models (with a warning), but if teh final "best" model has unstable coefficients, then return \code{NULL} for the best model.
#' @param construct Logical. If \code{TRUE} (default) then construct model from individual terms entered in order from lowest to highest AICc up to limits set by \code{presPerTermInitial} or \code{initialTerms} is met. If \code{FALSE} then the "full" model consists of all terms allowed by \code{quadratic} and \code{interaction}.
#' @param select Logical. If \code{TRUE} (default) then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param quadratic Logical. Used only if \code{construct} is \code{TRUE}. If \code{TRUE} (default) then include quadratic terms in model construction stage for non-factor predictors.
#' @param interaction Logical. Used only if \code{construct} is \code{TRUE}. If \code{TRUE} (default) then include 2-way interaction terms (including interactions between factor predictors).
#' @param verboten Either \code{NULL} (default) in which case \code{forms} is returned without any manipulation. Alternatively, this is a character list of terms that are not allowed to appear in any model in \code{forms}. Models with these terms are removed from \code{forms}. Note that the order of variables in interaction terms does not matter (e.g., \code{x1:x2} will cause the removal of models with this term verbatim as well as \code{x2:x1}). All possible permutations of three-way interaction terms are treated similarly.
#' @param verbotenCombos Either \code{NULL} or a list of lists. This argument allows excluding particular combinations of variables using exact matches (i.e., a variable appears exactly as stated) or general matches (i.e., a variable appears in any term). Please see the \emph{Details} section of \code{\link[statisfactory]{makeFormulae}} for more information on how to use this argument. The default is \code{NULL} in which case any combination of variables is allowed.
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is TRUE.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model. Used only if \code{select} is \code{TRUE}.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is \code{TRUE}.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign equal total weights to presences and contrast sites (\code{TRUE}).
#' @param method Character, name of function used to solve. This can be \code{'glm.fit'} (default), \code{'brglmFit'} (from the \pkg{brglm2} package), or another function.
#' @param out Character. Indicates type of value returned. Values can be \code{'model'} (default; return model with lowest AICc), \code{'models'} (return a list of all models), and/or \code{'tuning'} (return a data frame with AICc for each model). If more than one value is specified, then the output will be a list with elements named "model", "models", and/or "tuning". The models will appear in the list in same order as they appear in the tuning table (i.e., model with the lowest AICc first, second-lowest next, etc.). If just one value is specified, the output will be either an object of class \code{MaxEnt}, a list with objects of class \code{MaxEnt}, or a data frame.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device.
#' @param ... Arguments to pass to \code{glm}.
#' @examples
#' \dontrun{
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
trainGlm <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	construct = TRUE,
	select = TRUE,
	anyway = FALSE,
	quadratic = TRUE,
	interaction = TRUE,
	verboten = NULL,
	verbotenCombos = NULL,
	presPerTermInitial = 10,
	presPerTermFinal = 10,
	initialTerms = 10,
	w = TRUE,
	method = 'glm.fit',
	out = 'model',
	tooBig = 10E6,
	verbose = FALSE,
	...
) {

	#####################
	### for debugging ###
	#####################

	if (FALSE) {
	
		family <- 'binomial'
		construct <- TRUE
		select <- TRUE
		anyway <- FALSE
		quadratic <- TRUE
		interaction <- TRUE
		verboten <- NULL
		verbotenCombos <- NULL
		presPerTermInitial <- 10
		presPerTermFinal <- 20
		initialTerms <- 10
		w <- TRUE
		method <- 'glm.fit'
		out <- 'model'
		tooBig <- 10E6
		verbose <- TRUE
		
	}

	#############
	### setup ###
	#############

	# force number of starting terms to 31 or less
	if (select & initialTerms > 30) {
		initialTerms <- 30
		warning('initialTerms must be 30 or fewer. Forcing to 30.')
	}

	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# number of data
	sampleSize <- if (family=='binomial') {
		sum(data[ , resp])
	} else {
		nrow(data)
	}

	# model weights
	if (class(w)[1] == 'logical') {
		w <- if (w & (family == 'binomial' | family == 'quasibinomial')) {
			c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
		} else {
			rep(1, nrow(data))
		}
	} else if (class(w) == 'character') {
		w <- data[ , w]
	}

	w <- w / max(w) # declare to global because dredge() has problems if it is not

	## MODEL CONSTRUCTION
	#####################

	# create starting formula
	form <- paste0(resp, ' ~ 1')

	if (construct) {

		tuning <- data.frame()

		## UNIVARIATE terms
		for (thisPred in preds) { # for each predictor test single-variable terms

			# train model
			# thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', thisPred)), family=family, data=data, weights=w, method=method)
			thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', thisPred)), family=family, data=data, weights=w, method=method, ...)

			# get AICc
			thisAic <- MuMIn::AICc(thisModel)
			k <- length(thisModel$coefficients)

			aicc <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

			# remember if coefficients were stable
			if (all(!is.na(stats::coef(thisModel)))) {
				if (all(abs(stats::coef(thisModel)) < tooBig)) {

					tuning <- rbind(tuning, data.frame(type='linear', term=thisPred, aicc=aicc, terms=1))
					
				}
			}

		} # next single-variable term
		
		if (nrow(tuning) == 0 & !anyway) {
			warning('No univariate models were stable (at least one covariate > "tooBig").')
			return(NULL)
		}

		## QUADRATIC terms
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (quadratic & ((sampleSize / 2 >= presPerTermInitial & initialTerms >= 2) | (sampleSize / 2 >= presPerTermInitial & initialTerms >= 2))) {

			for (thisPred in preds) { # for each predictor test single-variable terms

				if (class(data[ , thisPred]) != 'factor') {

					term <- paste0(thisPred, ' + I(', thisPred, '^2)')

					# train model
					# thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method)
					thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# get AICc
					thisAic <- MuMIn::AICc(thisModel)
					k <- length(thisModel$coefficients)

					aicc <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# remember if coefficients were stable
					if (all(!is.na(stats::coef(thisModel)))) {
						if (all(abs(stats::coef(thisModel)) < tooBig | anyway)) {
							tuning <- rbind(tuning, data.frame(type='quadratic', term=term, aicc=aicc, terms=2))
						}
					}

				}

			} # next quadratic term

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		# ## CUBIC TERMS
		# # if there are more than desired number of presences per term and initial model can have more than 1 term

		# if (cubic & ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			# for (thisPred in preds) { # for each predictor test cubic terms

				# if (class(data[ , thisPred]) != 'factor') {

					# term <- paste0(thisPred, ' + I(', thisPred, '^2) + I(', thisPred, '^3)')

					# # train model
					# thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# # get AICc
					# thisAic <- MuMIn::AICc(thisModel)
					# k <- length(thisModel$coefficients)

					# aicc <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# # remember if coefficients were stable
					# if (all(!is.na(stats::coef(thisModel)))) {
						# if (all(abs(stats::coef(thisModel)) < tooBig | anyway)) {
							# tuning <- rbind(tuning, data.frame(type='cubic', term=term, aicc=aicc, terms=3))
						# }
					# }

				# }

			# } # next cubic term

		# } # if there are more than desired number of presences per term and initial model can have more than 1 term

		## 2-WAY INTERACTION TERMS
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (interaction & ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			for (countPred1 in 1:(length(preds) - 1)) { # for each predictor test two-variable terms

				for (countPred2 in (countPred1 + 1):length(preds)) { # for each second predictor test two-variable terms

					thisPred <- preds[countPred1]
					thatPred <- preds[countPred2]

					term <- paste0(thisPred, ' + ', thatPred, ' + ', thisPred, ':', thatPred)

					# train model
					# thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method)
					thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# get AICc
					thisAic <- MuMIn::AICc(thisModel)
					k <- length(thisModel$coefficients)

					aicc <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# remember if coefficients were stable
					if (all(!is.na(stats::coef(thisModel)))) {
						if (all(abs(stats::coef(thisModel)) < tooBig | anyway)) {
							tuning <- rbind(tuning, data.frame(type='interaction', term=term, aicc=aicc, terms=3))
						}
					}

				} # for each second predictor test interaction terms

			} # for each predictor test interaction terms

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		# ## INTERACTION-QUADRATIC terms
		# # if there are more than desired number of presences per term and initial model can have more than 1 term
		# if (interQuad & ((sampleSize / 5 >= presPerTermInitial & initialTerms >= 5) | (sampleSize / 5 >= presPerTermInitial & initialTerms >= 5))) {

			# for (thisPred in preds) { # for each predictor

				# for (thatPred in preds[!(preds %in% thisPred)]) { # for each second predictor

					# term <- if (class(data[ , thatPred]) != 'factor') {

						# paste0(thisPred, ' + ', thatPred, ' + I(', thatPred, '^2) + ', thisPred, ':', thatPred, ' + ', thisPred, ':I(', thatPred, '^2)')

					# } else { NA }

					# if (!is.na(term)) {

						# # train model
						# thisModel <- stats::glm(formula=stats::as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

						# # get aicc
						# thisAic <- MuMIn::AICc(thisModel)
						# k <- length(thisModel$coefficients)

						# thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

						# # remember if coefficients were stable
						# if (all(!is.na(stats::coef(thisModel)))) {
							# if (all(abs(stats::coef(thisModel)) < tooBig | anyway)) {
								# tuning <- rbind(tuning, data.frame(type='interaction-quadratic', term=term, AICc=thisAic, terms=4))
							# }
						# }

					# }

				# } # for each second predictor test interaction terms

			# } # for each predictor test interaction terms

		# } # if there are more than desired number of presences per term and initial model can have more than 1 term

		# sort by AIC
		tuning <- tuning[order(tuning$aicc), ]

		if (verbose) {
			omnibus::say('GLM construction results for each term tested:');
			print(tuning)
			omnibus::say('')
		}

		### train all possible models and select best by AIC
		####################################################

		## construct final model
		form <- paste0(resp, ' ~ 1 + ', tuning$term[1]) # add first term

		numTerms <- length(colnames(attr(stats::terms(stats::as.formula(form)), 'factors')))

		# if there are more presence sites than required per term in model and if terms in model are fewer than specified limit
		if ((sampleSize / numTerms) > presPerTermInitial & numTerms < initialTerms) {

			# initialize number of candidate term
			glmFrameRow <- 2

			# add terms
			while ((sampleSize / numTerms) > presPerTermInitial & numTerms < initialTerms & glmFrameRow <= nrow(tuning)) {

				# make trial formula
				trialForm <- stats::as.formula(paste0(form, ' + ', tuning$term[glmFrameRow]))
				termsInTrial <- length(colnames(attr(stats::terms(stats::as.formula(trialForm)), 'factors')))

				# update formula if there are enough presences per term
				if (sampleSize / termsInTrial >= presPerTermInitial & termsInTrial <= initialTerms) {
					form <- paste0(form, ' + ', tuning$term[glmFrameRow])
				}

				# get number of unique terms that would be added
				numTerms <- length(colnames(attr(stats::terms(stats::as.formula(form)), 'factors')))

				# look at next row of tuning
				glmFrameRow <- glmFrameRow + 1

			} # next term

		} # if there are sufficient presences for additional terms beyond first

	} else { # use all terms (no stepwise model construction)

		# univariate terms
		for (thisPred in preds) form <- paste0( form, ' + ', thisPred)
		if (quadratic) for (thisPred in preds) form <- paste0( form, ' + I(', thisPred, '^2)')
		# if (cubic) for (thisPred in preds) form <- paste0( form, ' + I(', thisPred, '^3)')

		# interactions
		if (length(preds) > 1) {

			# 2-WAY INTERACTIONS
			if (interaction) {
				for (thisPred in preds[1:(length(preds) - 1)]) { # for each initial predictor
					for (thatPred in preds[2:length(preds)]) {
						form <- paste0(form, ' + ', thisPred, ':', thatPred)
					}
				}
			}

			# # INTERACTION with QUADRATIC
			# if (interQuad) {
				# for (thisPred in preds) {
					# for (thatPred in preds[!(preds %in% thisPred)]) {
						# if (class(data[ , thatPred]) != 'factor') {
							# if (!quadratic) form <- paste0(form, ' + I(', thatPred, '^2)')
							# if (!interaction) form <- paste0(form, ' + ', thisPred, ':', thatPred)
							# form <- paste0(form, ' + ', thisPred , ':I(', thatPred, '^2)')
						# }
					# }
				# }
			# }

		}

	} # if not doing automated model construction

	# convert to formula
	form <- stats::as.formula(form)

	## MODEL SELECTION
	##################

	if (!select) {

		# train (starting) GLM model
		model <- stats::glm(form, family=family, data=data, weights=w, method=method, ...)

		if (verbose) {
			omnibus::say('Full model:', pre=1);
			print(summary(model))
			utils::flush.console()
		}

	} else {

		# store *all* models
		models <- list()

		# maximum number of terms to accommodate in any model
		maxTerms <- if (family == 'binomial') {
			max(1, 1 + floor(sampleSize / presPerTermFinal))
		} else {
			Inf
		}

		# make all possible formula
		forms <- statisfactory::makeFormulae(
			form,
			maxTerms=maxTerms,
			intercept=TRUE,
			interceptOnly=TRUE,
			linearOnly=TRUE,
			quad=quadratic,
			ia=interaction,
			verboten=verboten,
			verbotenCombos=verbotenCombos,
			returnFx=as.character
		)

		tuning <- data.frame()

		# evaluate each model
		for (i in seq_along(forms)) {

			form <- stats::as.formula(forms[[i]])

			# models[[i]] <- stats::glm(form, family=family, data=data, weights=w, method=method)
			models[[i]] <- stats::glm(form, family=family, data=data, weights=w, method=method, ...)
			ll <- stats::logLik(models[[i]])
			aicc <- MuMIn::AICc(models[[i]])
			models[[i]]$aicc <- aicc

			tuning <- rbind(
				tuning,
				data.frame(
					model = forms[[i]],
					logLik = ll,
					aicc = aicc
				)
			)


		}

		models <- rlist::list.sort(models, aicc)
		tuning <- tuning[order(tuning$aicc), ]

		if (verbose) {
			omnibus::say('')
			print(tuning)
			omnibus::say('')
		}

		# get best model
		if ('model' %in% out) {

			model <- models[[1]]
			
			if (all(abs(stats::coef(model)) < tooBig)) {

				if (verbose) {
					omnibus::say('Final model:', pre=1);
					print(summary(model))
					omnibus::say('')
				}
				
			} else {
				warning('The best model had unstable coefficients (> "tooBig").', immediate.=TRUE)
				model <- NULL
			}

		}

	} # model selection

	# return
	if (length(out) > 1) {
		output <- list()
		if ('models' %in% out & select) output$models <- models
		if ('model' %in% out) output$model <- model
		if ('tuning' %in% out) output$tuning <- tuning
		output
	} else if ('models' %in% out) {
		models
	} else if ('model' %in% out) {
		model
	} else if ('tuning' %in% out) {
		tuning
	}

}
