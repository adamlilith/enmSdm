#' Calibrate a generalized linear model (GLM)
#'
#' This function constructs a GLM piece-by-piece by first calculating AICc for all models with univariate, quadratic, and 2-way-interaction terms. It then creates a "full" model with the highest-ranked uni/bivariate terms. Finally, it implements an all-subsets model selection routine using AICc. Its output is any or all of: a table with AICc for all possible models, all possible models (after model construction), and/or the model with the lowest AICc.
#' @param data Data frame. Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}). Default is to use the 'binomial' family.
#' @param tooBig Numeric. Used to catch errors when fitting a model fit with the \code{brglmFit} function in the \pkg{brglm2} package. In some cases fitted coefficients are unstable and tend toward very high values, even if training data is standardized. Models with such coefficients will be discarded if any one coefficient is \code{> tooBig}. Set equal to \code{Inf} to keep all models.
#' @param construct Logical. If \code{TRUE} (default) then construct model from individual terms entered in order from lowest to highest AICc up to limits set by \code{presPerTermInitial} or \code{initialTerms} is met. If \code{FALSE} then the "full" model consists of all terms allowed by \code{quadratic} and \code{interaction}.
#' @param select Logical. If \code{TRUE} (default) then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param quadratic Logical. Used only if \code{construct} is \code{TRUE}. If \code{TRUE} (default) then include quadratic terms in model construction stage for non-factor predictors.
#' @param interaction Logical. Used only if \code{construct} is \code{TRUE}. If \code{TRUE} (default) then include 2-way interaction terms (including interactions between factor predictors).
#' @param verboten Either \code{NULL} (default) in which case \code{forms} is returned without any manipulation. Alternatively, this is a character list of terms that are not allowed to appear in any model in \code{forms}. Models with these terms are removed from \code{forms}. Note that the order of variables in interaction terms does not matter (e.g., \code{x1:x2} will cause the removal of models with this term verbatim as well as \code{x2:x1}). All possible permutations of three-way interaction terms are treated similarly.
#' @param verbotenCombos Either \code{NULL} or a list of lists. This argument allows excluding particular combinations of variables using exact matches (i.e., a variable appears exactly as stated) or general matches (i.e., a variable appears in any term). Please see the \emph{Details} section of \code{\link[statisfactory]{makeFormulae}} for more information on how to use this argument. The default is \code{NULL} in which case any combination of variables is allowed.
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is TRUE.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model. Used only if \code{select} is TRUE.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is \code{TRUE}.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign equal total weights to presences and contrast sites (\code{TRUE}).
#' @param method Character, name of function used to solve. This can be \code{'glm.fit'} (default), \code{'brglmFit'} (from the \pkg{brglm2} package), or another function.
#' @param out Character. Indicates type of value returned. Values can be \code{'model'} (default; return model with lowest AICc), \code{'models'} (return a list of all models), and/or \code{'tuning'} (return a data frame with AICc for each model). If more than one value is specified, then the output will be a list with elements named "model", "models", and/or "tuning". The models will appear in the list in same order as they appear in the tuning table (i.e., model with the lowest AICc first, second-lowest next, etc.). If just one value is specified, the output will be either an object of class \code{MaxEnt}, a list with objects of class \code{MaxEnt}, or a data frame.
#' @param verbose Logical. If \code{TRUE} then display intermediate results on the display device.
#' @param ... Arguments to pass to \code{glm}.
#' @examples
#' set.seed(123)
#' x <- matrix(rnorm(n = 6*100), ncol = 6)
#' # true variables will be #1, #2, #5, and #6, plus
#' # the squares of #1 and #6, plus
#' # interaction between #1 and #6
#' # the cube of #5
#' imp <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x1_pow2', 'x6_pow2', 'x1_by_x6', 'x5_pow3')
#' betas <- c(5, 2, 0, 0, 1, -1, 8, 1, 2, -4)
#' names(betas) <- imp
#' y <- 0.5 + x %*% betas[1:6] + betas[7] * x[ , 1] +
#' betas[8] * x[ , 6] + betas[9] * x[ , 1] * x[ , 6] + betas[10] * x[ , 5]^3
#' y <- as.integer(y > 10)
#' x <- cbind(y, x)
#' x <- as.data.frame(x)
#' names(x) <- c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' model <- trainGlmDredge(x, verbose=TRUE)
#' @seealso \code{\link[enmSdm]{trainGlmDredge}}, \code{\link[stats]{glm}} in the \pkg{stats} package, \code{\link[brglm2]{brglmFit}} in the \pkg{brglm2} package
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- matrix(rnorm(n = 6*100), ncol = 6)
#' # true variables will be #1, #2, #5, and #6, plus
#' # the squares of #1 and #6, plus
#' # interaction between #1 and #6
#' # the cube of #5
#' imp <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x1_pow2', 'x6_pow2', 'x1_by_x6', 'x5_pow3')
#' betas <- c(5, 2, 0, 0, 1, -1, 8, 1, 2, -4)
#' names(betas) <- imp
#' y <- 0.5 + x %*% betas[1:6] + betas[7] * x[ , 1] +
#' betas[8] * x[ , 6] + betas[9] * x[ , 1] * x[ , 6] + betas[10] * x[ , 5]^3
#' y <- as.integer(y > 10)
#' x <- cbind(y, x)
#' x <- as.data.frame(x)
#' names(x) <- c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' model <- trainGlm(x, verbose=TRUE)
#' }
#' @export
trainGlm <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	tooBig = 10E6,
	construct = TRUE,
	select = TRUE,
	quadratic = TRUE,
	interaction = TRUE,
	verboten = NULL,
	verbotenCombos = NULL,
	presPerTermInitial = 10,
	presPerTermFinal = 20,
	initialTerms = 10,
	w = TRUE,
	method = 'glm.fit',
	out = 'model',
	verbose = FALSE,
	...
) {

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

	w <<- w / max(w) # declare to global because dredge() has problems if it is not

	## MODEL CONSTRUCTION
	#####################

	# create starting formula
	form <- paste0(resp, ' ~ 1')

	if (construct) {

		## UNIVARIATE terms
		for (thisPred in preds) { # for each predictor test single-variable terms

			# train model
			# thisModel <- glm(formula=as.formula(paste0(form, ' + ', thisPred)), family=family, data=data, weights=w, method=method)
			thisModel <- glm(formula=as.formula(paste0(form, ' + ', thisPred)), family=family, data=data, weights=w, method=method, ...)

			# get AICc
			thisAic <- AIC(thisModel)
			k <- length(thisModel$coefficients)

			thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

			# remember if coefficients were stable
			if (all(!is.na(stats::coef(thisModel)))) {
				if (all(abs(stats::coef(thisModel)) < tooBig)) {

					tuning <- if (exists('tuning', inherits=FALSE)) {
						rbind(tuning, data.frame(type='linear', term=thisPred, AICc=thisAic, terms=1))
					} else {
						data.frame(type='linear', term=thisPred, AICc=thisAic, terms=1)
					}

				}
			}

		} # next single-variable term

		## QUADRATIC terms
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (quadratic && ((sampleSize / 2 >= presPerTermInitial & initialTerms >= 2) | (sampleSize / 2 >= presPerTermInitial & initialTerms >= 2))) {

			for (thisPred in preds) { # for each predictor test single-variable terms

				if (class(data[ , thisPred]) != 'factor') {

					term <- paste0(thisPred, ' + I(', thisPred, '^2)')

					# train model
					# thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method)
					thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# get AICc
					thisAic <- AIC(thisModel)
					k <- length(thisModel$coefficients)

					thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# remember if coefficients were stable
					if (all(!is.na(coef(thisModel)))) {
						if (all(abs(coef(thisModel)) < tooBig)) {
							tuning <- rbind(tuning, data.frame(type='quadratic', term=term, AICc=thisAic, terms=2))
						}
					}

				}

			} # next quadratic term

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		# ## CUBIC TERMS
		# # if there are more than desired number of presences per term and initial model can have more than 1 term

		# if (cubic && ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			# for (thisPred in preds) { # for each predictor test cubic terms

				# if (class(data[ , thisPred]) != 'factor') {

					# term <- paste0(thisPred, ' + I(', thisPred, '^2) + I(', thisPred, '^3)')

					# # train model
					# thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# # get AICc
					# thisAic <- AIC(thisModel)
					# k <- length(thisModel$coefficients)

					# thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# # remember if coefficients were stable
					# if (all(!is.na(coef(thisModel)))) {
						# if (all(abs(coef(thisModel)) < tooBig)) {
							# tuning <- rbind(tuning, data.frame(type='cubic', term=term, AICc=thisAic, terms=3))
						# }
					# }

				# }

			# } # next cubic term

		# } # if there are more than desired number of presences per term and initial model can have more than 1 term

		## 2-WAY INTERACTION TERMS
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (interaction && ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			for (countPred1 in 1:(length(preds) - 1)) { # for each predictor test two-variable terms

				for (countPred2 in (countPred1 + 1):length(preds)) { # for each second predictor test two-variable terms

					thisPred <- preds[countPred1]
					thatPred <- preds[countPred2]

					term <- paste0(thisPred, ' + ', thatPred, ' + ', thisPred, ':', thatPred)

					# train model
					# thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method)
					thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# get AICc
					thisAic <- AIC(thisModel)
					k <- length(thisModel$coefficients)

					thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# remember if coefficients were stable
					if (all(!is.na(coef(thisModel)))) {
						if (all(abs(coef(thisModel)) < tooBig)) {
							tuning <- rbind(tuning, data.frame(type='interaction', term=term, AICc=thisAic, terms=3))
						}
					}

				} # for each second predictor test interaction terms

			} # for each predictor test interaction terms

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		# ## INTERACTION-QUADRATIC terms
		# # if there are more than desired number of presences per term and initial model can have more than 1 term
		# if (interQuad && ((sampleSize / 5 >= presPerTermInitial & initialTerms >= 5) | (sampleSize / 5 >= presPerTermInitial & initialTerms >= 5))) {

			# for (thisPred in preds) { # for each predictor

				# for (thatPred in preds[!(preds %in% thisPred)]) { # for each second predictor

					# term <- if (class(data[ , thatPred]) != 'factor') {

						# paste0(thisPred, ' + ', thatPred, ' + I(', thatPred, '^2) + ', thisPred, ':', thatPred, ' + ', thisPred, ':I(', thatPred, '^2)')

					# } else { NA }

					# if (!is.na(term)) {

						# # train model
						# thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

						# # get AICc
						# thisAic <- AIC(thisModel)
						# k <- length(thisModel$coefficients)

						# thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

						# # remember if coefficients were stable
						# if (all(!is.na(coef(thisModel)))) {
							# if (all(abs(coef(thisModel)) < tooBig)) {
								# tuning <- rbind(tuning, data.frame(type='interaction-quadratic', term=term, AICc=thisAic, terms=4))
							# }
						# }

					# }

				# } # for each second predictor test interaction terms

			# } # for each predictor test interaction terms

		# } # if there are more than desired number of presences per term and initial model can have more than 1 term

		# sort by AIC
		tuning <- tuning[order(tuning$AIC, tuning$terms), ]

		if (verbose) {
			omnibus::say('GLM construction results for each term tested:');
			print(tuning)
			omnibus::say('')
		}

		### train all possible models and select best by AIC
		####################################################

		## construct final model
		form <- paste0(resp, ' ~ 1 + ', tuning$term[1]) # add first term

		numTerms <- length(colnames(attr(terms(as.formula(form)), 'factors')))

		# if there are more presence sites than required per term in model and if terms in model are fewer than specified limit
		if ((sampleSize / numTerms) > presPerTermInitial & numTerms < initialTerms) {

			# initialize number of candidate term
			glmFrameRow <- 2

			# add terms
			while ((sampleSize / numTerms) > presPerTermInitial & numTerms < initialTerms & glmFrameRow <= nrow(tuning)) {

				# make trial formula
				trialForm <- as.formula(paste0(form, ' + ', tuning$term[glmFrameRow]))
				termsInTrial <- length(colnames(attr(terms(as.formula(trialForm)), 'factors')))

				# update formula if there are enough presences per term
				if (sampleSize / termsInTrial >= presPerTermInitial & termsInTrial <= initialTerms) {
					form <- paste0(form, ' + ', tuning$term[glmFrameRow])
				}

				# get number of unique terms that would be added
				numTerms <- length(colnames(attr(terms(as.formula(form)), 'factors')))

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
	form <- as.formula(form)

	## MODEL SELECTION
	##################

	if (!select) {

		# train (starting) GLM model
		model <- glm(form, family=family, data=data, weights=w, method=method, ...)

		if (verbose) {
			omnibus::say('Full model:', pre=1);
			print(summary(model))
			flush.console()
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

			form <- as.formula(forms[[i]])

			models[[i]] <- glm(form, family=family, data=data, weights=w, method=method, ...)
			ll <- logLik(models[[i]])
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

			if (verbose) {
				omnibus::say('Final model:', pre=1);
				print(summary(model))
				omnibus::say('')
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
