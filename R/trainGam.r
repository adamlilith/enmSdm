#' Calibrate a generalized additive model (GAM)
#'
#' This function constructs a GAM piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{?family}).
#' @param gamma Initial penalty to degrees of freedom to use (larger ==> smoother fits).
#' @param construct Logical. If TRUE then construct model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{initialTerms}.
#' @param select Logical. If TRUE then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only if \code{construct} is \code{TRUE}.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is TRUE.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{dredge()} is 31, so if this number is >31 and \code{select} is \code{TRUE} then it is forced to 31 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param interaction Character or \code{NULL}. Type of interaction term to use (\code{te}, \code{ts}, \code{s}, etc.). See \code{?te} (for example) for help on any one of these. If \code{NULL} then interactions are not used.
#' @param w Either logical in which case TRUE causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character. Indicates type of value returned. If \code{model} (default) then returns an object of class \code{brglm} or \code{glm} (depending on the value of \code{use}). If \code{tuning} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param verbose Logical. If TRUE then display intermediate results on the display device.
#' @param ... Extra arguments (not used).
#' @return If \code{out = 'model'} this function returns an object of class \code{gam}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gam} object and the data frame.
#' @seealso \code{\link[mgcv]{gam}}
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
#' model <- trainGam(x, verbose=TRUE)
#' model$tuning
#' summary(model$model)
#' }
#' @export

trainGam <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	gamma = 1,
	construct = TRUE,
	select = TRUE,
	presPerTermInitial = 10,
	presPerTermFinal = 20,
	initialTerms = 8,
	interaction = 'te',
	w = TRUE,
	out = 'model',
	verbose = FALSE,
	...
) {

	###########
	## setup ##
	###########

		# ellipses <- list(...)

		# force number of starting terms to 31 or less
		if (select & initialTerms > 31) {
			initialTerms <- 31
			warning('initialTerms must be 31 or less. Forcing to 31.')
		}

		# response and predictors
		if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
		if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	#############
	## weights ##
	#############

		# model weights
		if (class(w)[1] == 'logical') {
			w <- if (w && (family %in% c('binomial', 'quasibinomial'))) {
				c(rep(1, sum(data[ , resp])), rep(sum(data[ , resp]) / sum(data[ , resp] == 0), sum(data[ , resp] == 0)))
			} else {
				rep(1, nrow(data))
			}
		} else if (class(w) == 'character') {
			w <- data[ , w]
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

				term <- if (class(data[ , thisPred]) != 'factor') {
					paste0('s(', thisPred, ', bs=\'cs\')')
				} else {
					thisPred
				}

				thisThisForm <- paste0(form, ' + ', term)

				thisAic <- stats::AIC(
					mgcv::gam(
						formula=stats::as.formula(thisThisForm),
						family=family,
						data=data,
						method='ML',
						optimizer=c('outer', 'newton'),
						scale=-1,
						select=TRUE,
						gamma=gamma,
						weights=w
					)
				)

				# remember
				gamFrame <- if (exists('gamFrame')) {
					rbind(gamFrame, data.frame(term=term, AIC=thisAic))
				} else {
					data.frame(term=term, AIC=thisAic)
				}

			} # next single-variable term

			### TWO-variable terms
			if (length(preds) > 1 & !is.null(interaction)) {
				
				for (thisPred in preds[1:(length(preds)-1)]) { # for each predictor test two-variable terms

					for (thatPred in preds[ (which(preds==thisPred) + 1):length(preds) ]) { # for each second predictor test two-variable terms

						# create term
						term <- if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) != 'factor') {

							term <- paste0(interaction, '(', thisPred, ', ', thatPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) != 'factor') {

							paste0(interaction, '(', thatPred, ', by=', thisPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) == 'factor') {

							paste0(interaction, '(', thisPred, ', by=', thatPred, ', bs=\'cs\')')

						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) == 'factor') {

							paste0(thisPred, ' * ', thatPred)

						}

						thisAic <- AIC(
							mgcv::gam(
								formula=as.formula(paste0(form, ' + ', term)),
								family=family,
								data=data,
								method='ML',
								optimizer=c('outer', 'newton'),
								scale=-1,
								select=TRUE,
								gamma=gamma,
								weights=w
							)
						)

						# remember
						gamFrame <- if (exists('gamFrame')) {
							rbind(gamFrame, data.frame(term=term, AIC=thisAic))
						} else {
							data.frame(term=term, AIC=thisAic)
						}

					}  # for each second predictor test two-variable terms

				} # for each predictor test two-variable terms
				
			} # if interactions

			# sort by AIC
			gamFrame <- gamFrame[order(gamFrame$AIC), ]

			# print AICc frame
			if (verbose) {

				omnibus::say('GAM construction results for each term tested:', level=1)
				print(gamFrame)
				omnibus::say('')

			}

			## construct final model
			form <- paste0(form, ' + ', gamFrame$term[1]) # add first term

			# for each set of presences > min num required, add a term
			if (floor(sum(data[ , resp]) / presPerTermInitial ) - 1 > 1 & initialTerms > 1) {

				termsToAdd <- 2:min(initialTerms, c(floor(sum(data[ , resp]) / presPerTermInitial ) - 1, nrow(gamFrame) - 1))

				form <- paste0(form, ' + ', paste0(gamFrame$term[termsToAdd], collapse=' + '))

			} # if there are sufficient presences for additional terms beyond first

		# NO AUTOMATED MODEL CONSTRUCTION
		# use all single-variable terms and two-variable terms
		} else {

			# single terms
			for (thisPred in preds) {

				if (class(data[ , thisPred]) != 'factor') {
					# form <- paste0(form, ' + s(', thisPred, ', bs=\'cs\', k=basisK)')
					form <- paste0(form, ' + s(', thisPred, ')')
				} else {
					form <- paste0(form, ' + ', thisPred)
				}

			}

			# interaction terms
			if (length(preds) > 1 & !is.null(interaction)) {

				numPreds <- length(preds)
			
				for (countPred1 in 1:(numPreds - 1)) { # for each initial predictor

					thisPred <- preds[countPred1]
				
					for (countPred2 in (countPred1 + 1):numPreds) {
					
						thatPred <- preds[countPred2]

						form <- if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) != 'factor') {
							paste0( form, ' + ', interaction, '(', thisPred, ', ', thatPred, ')')
						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) != 'factor') {
							paste0( form, ' + ', interaction, '(', thatPred, ', by=', thisPred, ')')
						} else if (class(data[ , thisPred]) != 'factor' & class(data[ , thatPred]) == 'factor') {
							paste0( form, ' + ', interaction, '(', thisPred, ', by=', thatPred, ')')
						} else if (class(data[ , thisPred]) == 'factor' & class(data[ , thatPred]) == 'factor') {
							paste0(form, ' + ', thisPred, ' * ', thatPred)
						}

					}

				}
				
			}

		} # if not doing automated model construction


	###########################################################################
	## train model ############################################################
	## while model hasn't converged and while gamma is <= tolerance value... ##
	###########################################################################

	# # initialize flag to indicate if model converged
	# modelConverged <- FALSE
	# thisGamma <- gamma

	# while (modelConverged==FALSE & thisGamma <= 10) {

		# # if (verbose) omnibus::say('Basis K = ', basisK, ' | gamma = ', thisGamma)

		# # get GAM model... using automated scale selection with weights so influence of absences equals influence of presences... using tryCatch because sometimes for variables with too little variation the default df of the basis is too high, in which case it is reduced and attempted again (for univariate and bivariate models only)
		# model <- tryCatch(
			# mgcv::gam(
				# formula=as.formula(form),
				# family=family,
				# data=data,
				# method='ML',
				# optimizer=c('outer', 'newton'),
				# scale=-1,
				# select=TRUE,
				# gamma=thisGamma,
				# weights=w,
				# na.action='na.fail'
			# ),
			# error=function(err) return(FALSE)
		# )

		# # # if basis starting degrees of freedom was too large, reduce it by 1 and try again
		# # if (class(model)[1]=='logical' & basisK > 2) {

			# # basisK <- basisK - 1

		# # if basisK = 1 and model hasn't converged yet, then try increasing gamma and restart with initial basisK
		# # } else if (class(model)[1]=='logical' & basisK==2) {
		# if (class(model)[1]=='logical') {

			# thisGamma <- thisGamma + 0.4
			# # basisK <- startBasisK # reset basis to starting value if too low

		# # if model worked and converged, trip flag to say so
		# } else if (class(model)[1]!='logical') {

			# if (model$converged) modelConverged <- TRUE
			# # if (model$outer.info$conv=='full convergence' & model$converged==TRUE) modelConverged <- TRUE
			# if (modelConverged & verbose) say('Model converged.')

		# }

	# } # while model hasn't converged and while gamma is <= tolerance value

	# train FULL model
	model <- mgcv::gam(
		formula=as.formula(form),
		family=family,
		data=data,
		method='ML',
		optimizer=c('outer', 'newton'),
		scale=-1,
		select=TRUE,
		gamma=gamma,
		weights=w,
		na.action='na.fail'
	)

	if (verbose) {

		omnibus::say('Full model:', level=1);
		print(summary(model))
		flush.console()

	}

	#######################################################################################
	## if doing model construction, evaluate all possible models using AIC then get best ##
	#######################################################################################

	if (select & 'model' %in% out) {

		if (verbose) { omnibus::say('Calculating AICc across all possible models...', pre=1) }

		# calculate all possible models and rank by AIC
		sizeGamFrame <- if (exists('gamFrame', inherits=FALSE)) { nrow(gamFrame) } else { Inf }
		lims <- c(0, max(1, min(c(floor(sum(data[ , resp]) / presPerTermFinal), initialTerms, sizeGamFrame))))

		tuning <- MuMIn::dredge(
			global.model=model,
			rank='AICc',
			m.lim=lims,
			trace=FALSE
		)

		if ('model' %in% out) {

			if (verbose) omnibus::say('Final model:', level=1)

			# get model with best AIC
			thisFormula <- paste0(resp, ' ~ 1 + ', paste(names(tuning)[which(tuning[1, ] == '+')], sep=' + ', collapse=' + '))
			if (thisFormula == paste0(resp, ' ~ 1 + ')) thisFormula <- paste0(resp, ' ~ 1')

			# # initialize penalty for effective degrees of freedom... incrementing this until GAM converges
			# thisGamma <- initialGamma

			# initialize flag to indicate if model converged
			modelConverged <- FALSE

			# ## compute final GAM
			# while (modelConverged==FALSE & thisGamma <= 10) {

				# # using tryCatch because sometimes for variables with too little variation the default df of the basis is too high, in which case it is reduced and attempted again
				# model <- tryCatch(
					# mgcv::gam(
						# formula=as.formula(thisFormula),
						# family=family,
						# data=data,
						# method='REML',
						# optimizer=c('outer', 'newton'),
						# scale=-1,
						# select=TRUE,
						# gamma=thisGamma,
						# weights=w,
						# ...
					# ),
					# error=function(err) return(FALSE)
				# )

				# # if model hasn't converged yet, then try increasing gamma
				# if (class(model)[1]=='logical' ){
					# thisGamma <- thisGamma + 0.4
				# # if model worked and converged, trip flag to say so
				# } else {
					# if (model$converged) {
						# modelConverged <- TRUE
					# } else {
						# thisGamma <- thisGamma + 0.4
					# }
				# }

			# } # while model hasn't converged and while gamma is <= tolerance value

			## compute final GAM
			model <- mgcv::gam(
				formula=as.formula(thisFormula),
				family=family,
				data=data,
				method='REML',
				optimizer=c('outer', 'newton'),
				scale=-1,
				select=TRUE,
				gamma=gamma,
				weights=w,
				...
			)

			# # get best model (simple way... throws an error)
			# model <- MuMIn::get.models(tuning, subset = 1)[[1]]
			if (verbose) print(summary(model))

		} # if wanting model as output

	} # if doing model construction, evaluate all possible models using AIC then get best

	rm(w)

	# return
	if ('model' %in% out & 'tuning' %in% out) {
		out <- list()
		out$tuning <- tuning
		out$model <- model
		out
	} else if ('tuning' %in% out) {
		tuning
	} else {
		model
	}

}
