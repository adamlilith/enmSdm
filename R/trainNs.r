#' Calibrate a natural splines model
#'
#' This function constructs a natural-spline model piece-by-piece by first calculating AICc for all models with univariate and bivariate (interaction) terms. It then creates a "full" model with the highest-ranked uni/bivariate terms then implements an all-subsets model selection routine.
#' @param data Data frame.  Must contain fields with same names as in \code{preds} object.
#' @param resp Character or integer. Name or column index of response variable. Default is to use the first column in \code{data}.
#' @param preds Character list or integer list. Names of columns or column indices of predictors. Default is to use the second and subsequent columns in \code{data}.
#' @param family Name of family for data error structure (see \code{\link[stats]{family}}).
#' @param df Integer > 0 OR \code{NULL}. Sets flexibility of model fit. See documentation for \code{ns()}.
#' @param construct Logical. If TRUE then construct model by computing AICc for all univariate and bivariate models. Then add terms up to maximum set by \code{presPerTermInitial} and \code{initialTerms}.
#' @param select Logical. If TRUE then calculate AICc for all possible subsets of models and return the model with the lowest AICc of these. This step if performed \emph{after} model construction (if any).
#' @param presPerTermInitial Positive integer. Minimum number of presences needed per model term for a term to be included in the model construction stage. Used only is \code{construct} is TRUE.
#' @param presPerTermFinal Positive integer. Minimum number of presence sites per term in initial starting model; used only if \code{select} is TRUE.
#' @param initialTerms Positive integer. Maximum number of terms to be used in an initial model. Used only if \code{construct} is TRUE. The maximum that can be handled by \code{\link[MuMIn]{dredge }}is 31, so if this number is >31 and \code{select} is \code{TRUE} then it is forced to 31 with a warning. Note that the number of coefficients for factors is not calculated correctly, so if the predictors contain factors then this number might have to be reduced even more.
#' @param w Either logical in which case \code{TRUE} causes the total weight of presences to equal the total weight of absences (if \code{family='binomial'}) OR a numeric list of weights, one per row in \code{data} OR the name of the column in \code{data} that contains site weights. The default is to assign a weight of 1 to each datum.
#' @param out Character. Indicates type of value returned. If \code{model} (default) then returns an object of class \code{gam}. If \code{tuning} then just return the AICc table for each kind of model term used in model construction. If both then return a 2-item list with the best model and the AICc table.
#' @param verbose Logical. If TRUE then display intermediate results on the display device.
#' @param ... Arguments to send to \code{gam()} or \code{dredge()}.
#' @return If \code{out = 'model'} this function returns an object of class \code{gam}. If \code{out = 'tuning'} this function returns a data frame with tuning parameters and AICc for each model tried. If \code{out = c('model', 'tuning'} then it returns a list object with the \code{gam} object and the data frame.
#' @seealso \code{\link[splines]{ns}}, \code{\link[mgcv]{gam}}, \code{\link{trainGam}}
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
#' y <- as.integer(y > 0)
#' x <- cbind(y, x)
#' x <- as.data.frame(x)
#' names(x) <- c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' model <- trainNs(x, out=c('model', 'tuning'), verbose=TRUE)
#' model$tuning
#' summary(model$model)
#' }
#' @export

trainNs <- function(
	data,
	resp = names(data)[1],
	preds = names(data)[2:ncol(data)],
	family = 'binomial',
	df = NULL,
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
	if (class(w)[1] == 'logical') {
		w <- if (w & (family %in% c('binomial', 'quasibinomial'))) {
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
				paste0('splines::ns(', thisPred, ', df=', df, ')')
			} else {
				thisPred
			}

			thisThisForm <- paste0(form, ' + ', term)

			thisModel <- stats::glm(stats::as.formula(thisThisForm), family=family, data=data, weights=w, method='brglmFit', ...)
			thisAic <- AIC(thisModel)

			# remember
			tuning <- if (exists('tuning')) {
				rbind(tuning, data.frame(term=term, AIC=thisAic))
			} else {
				data.frame(term=term, AIC=thisAic)
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

		# single terms
		for (thisPred in preds) {

			if (class(data[ , thisPred]) != 'factor') {
				form <- paste0(form, ' + splines::ns(', thisPred, ', df=', df, ')')
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
		say('')

	}

	#######################################################################################
	## if doing model construction, evaluate all possible models using AIC then get best ##
	#######################################################################################

	if (select) {

		if (verbose) { omnibus::say('Calculating AICc across all possible models...') }

		# calculate all possible models and rank by AIC
		lims <- c(0, max(1, min(c(floor(sum(data[ , resp]) / presPerTermFinal), initialTerms, nrow(tuning)))))

		tuning <- MuMIn::dredge(
			global.model=model,
			rank='AICc',
			m.lim=lims,
			trace=FALSE,
			...
		)

		# get model with best AIC
		model <- MuMIn::get.models(tuning, subset = 1)[[1]]

	} # if model selection

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
