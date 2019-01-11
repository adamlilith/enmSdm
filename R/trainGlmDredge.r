#' Calibrate a generalized linear model (GLM)
#'
#' This is a pseudo-deprecated function to construct a GLM piece-by-piece by first calculating AICc for all models with univariate, quadratic, cubic, 2-way-interaction, and linear-by-quadratic terms. It then creates a "full" model with the highest-ranked uni/bivariate terms. Finally, it implements an all-subsets model selection routine using AICc. Its output is a table with AICc for all possible models (resulting from the "full" model) and/or the model of these with the lowest AICc. The procedure uses Firth's penalized likelihood to address issues related to seperability, small sample size, and bias. This function uses \code{\link[MuMIn]{dredge}} in the \pkg{MuMIn} package to cycle through all possible models.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}). Default is to use the 'binomial' family.
#' @param tooBig Numeric. Used to catch errors when fitting a model fit with the \code{brglmFit} function in the \pkg{brglm2} package. In some cases fitted coefficients are unstable and tend toward very high values, even if training data is standardized. Models with such coefficients will be discarded if any one coefficient is \code{> tooBig}. Set equal to \code{Inf} to keep all models.
#' @param construct Logical. If TRUE then construct model from individual terms entered in order from lowest to highest AICc up to limits set by \code{presPerTermInitial} or \code{initialTerms} is met. If \code{FALSE} then the "full" model consists of all terms allowed by \code{quadratic}, \code{cubic}, \code{interaction}, and \code{interQuad}.
#' @param select Logical. If TRUE then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param quadratic Logical. Used only if \code{construct} is TRUE. If TRUE then include quadratic terms in model construction stage for non-factor predictors.
#' @param cubic Logical. Used only if \code{construct} is TRUE. If TRUE then include cubic terms in model construction stage for non-factor predictors.
#' @param interaction Logical. Used only if \code{construct} is TRUE. If TRUE then include 2-way interaction terms (including interactions between factor predictors).
#' @param interQuad Logical. Used only if \code{construct} is TRUE. If TRUE then include all possible interactions of the form 'x * y^2' unless 'y' is a factor.
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is TRUE.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model. Used only if \code{select} is TRUE.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{dredge()} is 30, so if this number is >30 and \code{select} is \code{TRUE} then it is forced to 30 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign equal total weights to presences and contrast sites (\code{TRUE}).
#' @param method Character, name of function used to solve. This can be \code{'glm.fit'} (default), \code{'brglmFit'} (from the \pkg{brglm2} package), or another function.
#' @param out Character. Indicates type of value returned. If \code{model} (default) then returns an object of class \code{brglm2}/\code{glm}. If \code{table} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param verbose Logical. If TRUE then display intermediate results on the display device.
#' @param ... Arguments to pass to \code{brglm()} or \code{dredge()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{glm}. If \code{out = 'table'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'table'} then it returns a list object with the \code{glm} object and the data frame.
#' @seealso \code{\link[enmSdm]{trainGlm}}, \code{\link[stats]{glm}} in the \pkg{stats} package, \code{\link[brglm2]{brglmFit}} in the \pkg{brglm2} package
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
#' model <- trainGlmDredge(x, verbose=TRUE)
#' }
#' @export

trainGlmDredge <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	tooBig = 10E6,
	construct = TRUE,
	select = TRUE,
	quadratic = TRUE,
	cubic = TRUE,
	interaction = TRUE,
	interQuad = TRUE,
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

	################################
	## initial model construction ##
	################################

	# create starting formula
	form <- paste0(resp, ' ~ 1')

	if (construct) {

		## UNIVARIATE terms
		for (thisPred in preds) { # for each predictor test single-variable terms

			# train model
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

		## CUBIC TERMS
		# if there are more than desired number of presences per term and initial model can have more than 1 term

		if (cubic && ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			for (thisPred in preds) { # for each predictor test cubic terms

				if (class(data[ , thisPred]) != 'factor') {

					term <- paste0(thisPred, ' + I(', thisPred, '^2) + I(', thisPred, '^3)')

					# train model
					thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

					# get AICc
					thisAic <- AIC(thisModel)
					k <- length(thisModel$coefficients)

					thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

					# remember if coefficients were stable
					if (all(!is.na(coef(thisModel)))) {
						if (all(abs(coef(thisModel)) < tooBig)) {
							tuning <- rbind(tuning, data.frame(type='cubic', term=term, AICc=thisAic, terms=3))
						}
					}

				}

			} # next cubic term

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		## 2-WAY INTERACTION TERMS
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (interaction && ((sampleSize / 3 >= presPerTermInitial & initialTerms >= 3) | (sampleSize / 3 >= presPerTermInitial & initialTerms >= 3))) {

			for (countPred1 in 1:(length(preds) - 1)) { # for each predictor test two-variable terms

				for (countPred2 in (countPred1 + 1):length(preds)) { # for each second predictor test two-variable terms

					thisPred <- preds[countPred1]
					thatPred <- preds[countPred2]

					term <- paste0(thisPred, ' + ', thatPred, ' + ', thisPred, ':', thatPred)

					# train model
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

				}  # for each second predictor test interaction terms

			} # for each predictor test interaction terms

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

		## INTERACTION-QUADRATIC terms
		# if there are more than desired number of presences per term and initial model can have more than 1 term
		if (interQuad && ((sampleSize / 5 >= presPerTermInitial & initialTerms >= 5) | (sampleSize / 5 >= presPerTermInitial & initialTerms >= 5))) {

			for (thisPred in preds) { # for each predictor

				for (thatPred in preds[!(preds %in% thisPred)]) { # for each second predictor

					term <- if (class(data[ , thatPred]) != 'factor') {

						paste0(thisPred, ' + ', thatPred, ' + I(', thatPred, '^2) + ', thisPred, ':', thatPred, ' + ', thisPred, ':I(', thatPred, '^2)')

					} else { NA }

					if (!is.na(term)) {

						# train model
						thisModel <- glm(formula=as.formula(paste0(form, ' + ', term)), family=family, data=data, weights=w, method=method, ...)

						# get AICc
						thisAic <- AIC(thisModel)
						k <- length(thisModel$coefficients)

						thisAic <- thisAic + (2 * k * (k + 1)) / (sampleSize - k - 1)

						# remember if coefficients were stable
						if (all(!is.na(coef(thisModel)))) {
							if (all(abs(coef(thisModel)) < tooBig)) {
								tuning <- rbind(tuning, data.frame(type='interaction-quadratic', term=term, AICc=thisAic, terms=4))
							}
						}

					}

				}  # for each second predictor test interaction terms

			} # for each predictor test interaction terms

		} # if there are more than desired number of presences per term and initial model can have more than 1 term

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
		if (cubic) for (thisPred in preds) form <- paste0( form, ' + I(', thisPred, '^3)')

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

			# INTERACTION with QUADRATIC
			if (interQuad) {
				for (thisPred in preds) {
					for (thatPred in preds[!(preds %in% thisPred)]) {
						if (class(data[ , thatPred]) != 'factor') {
							if (!quadratic) form <- paste0(form, ' + I(', thatPred, '^2)')
							if (!interaction) form <- paste0(form, ' + ', thisPred, ':', thatPred)
							form <- paste0(form, ' + ', thisPred , ':I(', thatPred, '^2)')
						}
					}
				}
			}

		}

	} # if not doing automated model construction

	# convert to formula
	form <- as.formula(form)

	##################
	## train model ###
	##################

	# train (starting) GLM model
	model <- glm(formula=form, family=family, data=data, weights=w, na.action=stats::na.fail, method=method, ...)

	if (verbose) {
		omnibus::say('Full model:', pre=1);
		print(summary(model))
		flush.console()
	}

	########################################################################################
	## if doing model construction, evaluate all possible models using AICc then get best ##
	########################################################################################

	if (select) {

		sizeGlmFrame <- if (exists('tuning', inherits=FALSE)) { nrow(tuning) } else { Inf }
		lims <- c(0, max(1, min(floor(sampleSize / presPerTermFinal), initialTerms, sizeGlmFrame)))

		# calculate all possible models and rank by AIC
		tuning <- MuMIn::dredge(
			global.model=model,
			rank='AICc',
			m.lim=lims,
			trace=FALSE,
			...
		)

		### remove any models with wild coefficients
		############################################

		coeffOk <- abs(coefficients(tuning)) < tooBig
		coeffOk <- apply(coeffOk, 1, function(y) all(y, na.rm=TRUE))
		tuning <- subset(tuning, coeffOk, recalc.weights=TRUE, recalc.delta=TRUE)

		### remove models that ignore marginality for polynomial terms
		##############################################################

		allModelsDf <- as.data.frame(tuning)

		modTerms <- terms(tuning)
		modTerms <- sort(modTerms[!(modTerms %in% '(Intercept)')])

		for (countTerm in seq_along(modTerms)) {

			modTerm <- modTerms[countTerm]

			# ensure QUADRATIC marginality
			if (substr(modTerm, nchar(modTerm) - 2, nchar(modTerm)) == '^2)' & !grepl(modTerm, pattern=':')) {

				# find column with linear modTerm that matches this quadratic modTerm
				linearTerm <- substr(modTerm, 3, nchar(modTerm) - 3)

				badModel <- which(!is.na(allModelsDf[ , modTerm]) & is.na(allModelsDf[ , linearTerm]))
				goodModel <- which(!({1:nrow(allModelsDf)} %in% badModel))

				if (length(goodModel) > 0) {
					tuning <- subset(tuning, goodModel, recalc.weights=TRUE, recalc.delta=TRUE)
					allModelsDf <- as.data.frame(tuning)
				}

			}

			# ensure QUADRATIC with INTERACTION marginality
			quadFirst <- (substr(modTerm, 1, 2) == 'I(' & grepl(modTerm, pattern='\\^2)') & grepl(modTerm, pattern=':'))
			quadSecond <- (substr(modTerm, 1, 2) != 'I(' & substr(modTerm, nchar(modTerm) - 2, nchar(modTerm)) == '^2)' & grepl(modTerm, pattern=':'))
			if (quadFirst | quadSecond) {

				# find column with linear modTerm that matches this quadratic modTerm
				linearTerm <- substr(modTerm, regexpr('[(]', modTerm) + 1, regexpr('\\^2)', modTerm) - 1)
				secondTerm <- if (quadFirst) {
					substr(modTerm, regexpr(':', modTerm) + 1, nchar(modTerm))
				} else {
					substr(modTerm, 1, regexpr(':', modTerm) - 1)
				}
				quadTerm <- paste0('I(', linearTerm, '^2)')
				iaTerm1 <- paste0(linearTerm, ':', secondTerm)
				iaTerm2 <- paste0(secondTerm, ':', linearTerm)
				iaTerm <- if (iaTerm1 %in% modTerms) { iaTerm1 } else { iaTerm2 }

				badModel <- which(!is.na(allModelsDf[ , modTerm]) & (is.na(allModelsDf[ , linearTerm]) | is.na(allModelsDf[ , quadTerm]) | is.na(allModelsDf[ , iaTerm])))
				goodModel <- which(!({1:nrow(allModelsDf)} %in% badModel))

				if (length(goodModel) > 0) {
					tuning <- subset(tuning, goodModel, recalc.weights=TRUE, recalc.delta=TRUE)
					allModelsDf <- as.data.frame(tuning)
				}

			}

			# ensure CUBIC marginality
			if (substr(modTerm, nchar(modTerm) - 2, nchar(modTerm)) == '^3)' & !grepl(modTerm, pattern=':')) {

				# find column with linear modTerm that matches this quadratic modTerm
				linearTerm <- substr(modTerm, 3, nchar(modTerm) - 3)
				quadraticTerm <- paste('I(', linearTerm, '^2)', sep='')

				badModel <- which(!is.na(allModelsDf[ , modTerm]) & (is.na(allModelsDf[ , linearTerm]) | is.na(allModelsDf[ , quadraticTerm])))
				goodModel <- which(!({1:nrow(allModelsDf)} %in% badModel))

				if (length(goodModel) > 0) {
					tuning <- subset(tuning, goodModel, recalc.weights=TRUE, recalc.delta=TRUE)
					allModelsDf <- as.data.frame(tuning)
				}

			}

			# ensure CUBIC with INTERACTION marginality
			cubicFirst <- (substr(modTerm, 1, 2) == 'I(' & grepl(modTerm, pattern='\\^3)') & grepl(modTerm, pattern=':'))
			cubicSecond <- (substr(modTerm, 1, 2) != 'I(' & substr(modTerm, nchar(modTerm) - 2, nchar(modTerm)) == '^3)' & grepl(modTerm, pattern=':'))
			if (cubicFirst | cubicSecond) {

				# find linear version of this term
				linearTerm <- substr(modTerm, regexpr('[(]', modTerm) + 1, regexpr('\\^3)', modTerm) - 1)
				secondTerm <- if (cubicFirst) {
					substr(modTerm, regexpr(':', modTerm) + 1, nchar(modTerm))
				} else {
					substr(modTerm, 1, regexpr(':', modTerm) - 1)
				}
				quadTerm <- paste0('I(', linearTerm, '^2)')
				cubicTerm <- paste0('I(', linearTerm, '^3)')

				# 2-way interacton term
				iaTerm1 <- paste0(linearTerm, ':', secondTerm)
				iaTerm2 <- paste0(secondTerm, ':', linearTerm)
				iaTerm <- if (iaTerm1 %in% modTerms) { iaTerm1 } else { iaTerm2 }

				# interaction-quadratic term
				quadTermIa1 <- paste0(secondTerm, ':I(', linearTerm, '^2)')
				quadTermIa2 <- paste0('I(', linearTerm, '^2):', secondTerm)
				quadTermIa <- if (quadTermIa1 %in% modTerms) { quadTermIa1 } else { quadTermIa2 }

				badModel <- which(!is.na(allModelsDf[ , modTerm]) & (is.na(allModelsDf[ , linearTerm]) | is.na(allModelsDf[ , secondTerm]) | is.na(allModelsDf[ , quadTerm]) | is.na(allModelsDf[ , cubicTerm]) |  is.na(allModelsDf[ , iaTerm]) |  is.na(allModelsDf[ , quadTermIa])))
				goodModel <- which(!({1:nrow(allModelsDf)} %in% badModel))

				if (length(goodModel) > 0) {
					tuning <- subset(tuning, goodModel, recalc.weights=TRUE, recalc.delta=TRUE)
					allModelsDf <- as.data.frame(tuning)
				}

			}

		} # next term

		if (any(is.na(tuning[ , '(Intercept)']))) tuning <- tuning[-which(is.na(tuning[ , '(Intercept)'])), ]

		if (verbose) {
			omnibus::say('')
			print(tuning)
			omnibus::say('')
		}

		# retrain best model
		if ('model' %in% out) {

			model <- MuMIn::get.models(tuning, subset = 1)[[1]]

			if (verbose) {
				omnibus::say('Final model:', pre=1);
				print(summary(model))
				omnibus::say('')
			}

		}

	} # if want model construction

	gc()

	# return
	if ('model' %in% out & 'table' %in% out) {
		out <- list()
		out$table <- tuning
		out$model <- model
		out
	} else if ('table' %in% out) {
		tuning
	} else {
		model
	}

}
