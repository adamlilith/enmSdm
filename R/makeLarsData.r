#' Add columns to a data matrix to represent polynomial and interaction terms.
#'
#' This function adds columns to a data matrix representing quadratic, cubic, 2-way-interaction, and linear:quadratic interactions. It is especially useful for preparing a data matrix for \code{\link{trainLars}} or \code{\link{predictLars}}.
#' @param data Data frame.
#' @param resp Character or integer or \code{NULL}. Name or index of colum in \code{data} that represents the response variable. If \code{NULL} then it is assumed there is no response column in \code{data}.
#' @param preds Character or integer. Names or indices of columns to use as predictors.
#' @param scale Logical \emph{or} a list. If TRUE then scale values in \code{data[ , preds]} to have mean of 0 and unit variance. Note that scaling is done before adding terms. If a list, then this is the same as, for example, the two attributes from \code{attributes(scale(data[ , preds]))} named \code{`scaled:center`}and \code{`scaled:scale`}. Ergo, they each are list of the centers and scales (usually means and standard deviation) of each column of \code{data[ , preds]}, and each has names given by \code{preds}.
#' @param quadratic Logical. If TRUE then include quadratic terms in model construction stage for non-factor predictors. Quadratic columns will be named \code{<predictor name>_pow2}.
#' @param cubic Logical. If TRUE then include cubic terms in model construction stage for non-factor predictors. Cubic columns will be named \code{<predictor name>_pow3}.
#' @param interaction Logical. If TRUE then include 2-way interaction terms (including interactions between factor predictors). Interaction columns will be named \code{<predictor 1 name>_by_<predictor 2 name>}.
#' @param interQuad Logical. If TRUE then include all possible interactions of the form \code{x * y^2} unless \code{y} is a factor (linear-by-quadratic features). Linear-by-quadratic columns will be named \code{<predictor 1 name>_by_<predictor 2 name>_pow2}.
#' @param na.rm Logical. If TRUE then remove all rows of \code{data} in which there is at least one \code{NA} among \code{resp} or \code{preds}. The default is FALSE, which will cause an error if any row has an \code{NA}.
#' @return An object of class \code{larsData} (which is also a \code{list}) with six elements:
#' * A character named "\code{resp}" indicating the name of the column of that contains the response variable.
#' * A character list named "\code{preds}" indicating the name of the columns of that contain the original predictors.
#' * A data frame named "\code{data}" containing \code{data} but with extra columns representing the added terms;
#' * A list object named "\code{scales}" representing the scale (mean and standard derviation) used to center and rescale values in the data frame; and
#' * A list named "\code{groups}" with groups (names of predictors that should be considered together based on marginality).
#' * A list named "\code{features}" indicating what kind of features were added to \code{data} (e.g., quadratic, cubic, etc.).
#' @details If \code{scale} is \code{TRUE} then predictors with zero variance will be removed from the data before creating higher-order terms.
#' @seealso \code{\link{trainLars}}, \code{\link{predictLars}}
#' @examples
#' set.seed(123)
#' x <- data.frame(y=c(rep(1, 10), rep(0, 10)), x1=1:10, x2=runif(20) * 1:20, x3=rnorm(20) - 1:20)
#' out <- makeLarsData(x, resp='y', preds=c('x1', 'x2', 'x3'))
#' str(out)

makeLarsData <- function(
	data,
	resp,
	preds,
	scale,
	quadratic,
	cubic,
	interaction,
	interQuad,
	na.rm=FALSE
) {

	### illogical feature requests
	if (cubic & !quadratic) {
		quadratic <- TRUE
		warning('Forcing quadratic features because cubic features are requested.', .immediate=TRUE)
	}

	if (interQuad & !quadratic) {
		quadratic <- TRUE
		warning('Forcing quadratic features because linear-by-quadratic features are requested.', .immediate=TRUE)
	}

	### get just columns with response and predictors
	if (!is.null(resp)) {
		if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
		y <- data[ , resp]
	} else {
		y <- NULL
	}
	
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]
	data <- data[ , preds, drop=FALSE]

	### remove NAs
	if (na.rm) {
	
		yNa <- if (is.null(y)) {
			NULL
		} else {
			which(is.na(y))
		}
		
		predsNa <- naRows(data)
		
		if (length(yNa) > 0) {
			y <- y[-yNa]
			data <- data[-yNa, ]
		}
	
		if (length(predsNa) > 0) {
			y <- y[-predsNa]
			data <- data[-predsNa, ]
		}
		
	}
	
	### rescale
	scales <- list()
	
	# if user supplies scale values, scale data by these
	if (class(scale) == 'list') {

		data <- data[ , names(scale$`scaled:center`)]
		data <- as.data.frame(base::scale(data, center=scale$`scaled:center`, scale=scale$`scaled:scale`))
		scales$`scaled:center` <- scale$`scaled:center`
		scales$`scaled:scale` <- scale$`scaled:scale`
	
	} else if (class(scale) == 'logical') {
	
		# automatically scale
		if (scale) {

			data <- base::scale(data)
			scales$`scaled:center` <- attributes(data)$`scaled:center`
			scales$`scaled:scale` <- attributes(data)$`scaled:scale`
			data <- as.data.frame(data)
			
			# remove variables with no variation
			if (any(scales$`scaled:scale` == 0)) {
				
				zeroVar <- names(scales$`scaled:scale`[scales$`scaled:scale` == 0])
				
				scales$`scaled:center` <- scales$`scaled:center`[-which(names(scales$`scaled:center`) %in% zeroVar), drop=FALSE]
				scales$`scaled:scale` <- scales$`scaled:scale`[-which(names(scales$`scaled:scale`) %in% zeroVar), drop=FALSE]
				data <- data[ , -which(names(data) %in% zeroVar), drop=FALSE]
				preds <- preds[-which(preds %in% zeroVar)]
				
				if (ncol(data) == 1 & interaction) interaction <- FALSE
				if (ncol(data) == 1 & interQuad) interQuad <- FALSE
				
				warning(paste0('Removing variables from predictor set because they have 0 variance: ', paste(zeroVar, collapse=' ')))
				
			}
			
		# no scaling
		} else {
			scales$`scaled:center` <- rep(0, ncol(data))
			scales$`scaled:scale` <- rep(1, ncol(data))
			names(scales$`scaled:center`) <- names(scales$`scaled:scale`) <- names(data)
		}
	}

	### matrix to record if a predictor is in a given column of the design matrix
	inCol <- matrix(FALSE, ncol=length(preds), nrow=length(preds))
	rownames(inCol) <- 1:nrow(inCol)
	colnames(inCol) <- preds
	diag(inCol) <- TRUE
	
	### define linear groups (one variable per group)
	groups <- list()
	for (i in seq_along(preds)) groups[[i]] <- preds[i]

	#############################
	## construct design matrix ##
	#############################

	# QUADRATIC TERMS
	if (quadratic) {

		for (thisPred in preds) {

			# add new data column(s)
			add <- as.data.frame(data[ , thisPred]^2)
			newPred <- paste0(thisPred, '_pow2')
			names(add) <- newPred
			data <- cbind(data, add)
			
			# record which predictor(s) occur in the new data column
			thisInCol <- matrix(preds %in% thisPred, ncol=length(preds))
			rownames(thisInCol) <- ncol(data)
			inCol <- rbind(inCol, thisInCol)

			# remember original and synthetic variables in this group
			vars <- c(
				thisPred,
				newPred
			)

			groups[[length(groups) + 1]] <- vars

		}

	}

	# CUBIC TERMS
	if (cubic) {

		for (thisPred in preds) {

			# add new data column(s)
			add <- as.data.frame(data[ , thisPred]^3)
			newPred <- paste0(thisPred, '_pow3')
			names(add) <- newPred
			data <- cbind(data, add)

			# record which predictor(s) occur in the new data column
			thisInCol <- matrix(preds %in% thisPred, ncol=length(preds))
			rownames(thisInCol) <- ncol(data)
			inCol <- rbind(inCol, thisInCol)

			# remember original and synthetic variables in this group
			vars <- c(
				thisPred,
				paste0(thisPred, '_pow2'),
				newPred
			)

			groups[[length(groups) + 1]] <- vars

		}

	}

	# 2-WAY INTERACTION TERMS
	if (interaction) {

		for (countPred1 in 1:(length(preds) - 1)) {
			for (countPred2 in (countPred1 + 1):length(preds)) {

				thisPred <- preds[countPred1]
				thatPred <- preds[countPred2]

				newPred <- paste0(thisPred, '_by_', thatPred)

				# add new data column(s)
				add <- as.data.frame(data[ , thisPred] * data[ , thatPred])
				names(add) <- newPred
				data <- cbind(data, add)

				# record which predictor(s) occur in the new data column
				thisInCol <- matrix(preds %in% thisPred | preds %in% thatPred, ncol=length(preds))
				rownames(thisInCol) <- ncol(data)
				inCol <- rbind(inCol, thisInCol)

				# remember original and synthetic variables in this group
				vars <- c(
					thisPred,
					thatPred,
					newPred
				)

				groups[[length(groups) + 1]] <- vars

			}
		}
	}

	# LINEAR-by-QUADRATIC INTERACTION TERM
	if (interQuad) {

		for (thisPred in preds) {
			for (thatPred in preds[!(preds %in% thisPred)]) {

				# add new data column(s)
				add <- as.data.frame(data[ , thisPred] * (data[ , thatPred])^2)
				names(add) <- paste0(thisPred, '_by_', thatPred, '_pow2')
				data <- cbind(data, add)

				# record which predictor(s) occur in the new data column
				thisInCol <- matrix(preds %in% thisPred | preds %in% thatPred, ncol=length(preds))
				rownames(thisInCol) <- ncol(data)
				inCol <- rbind(inCol, thisInCol)

				# remember original and synthetic variables in this group
				vars <- c(
					thisPred,
					thatPred,
					paste0(thisPred, '_by_', thatPred),
					paste0(thatPred, '_by_', thisPred),
					paste0(thatPred, '_pow2'),
					paste0(thisPred, '_by_', thatPred, '_pow2')
				)

				vars <- names(data)[which(names(data) %in% vars)]
				groups[[length(groups) + 1]] <- vars

			}
		}
	}
	
	out <- list()
	out$resp <- resp
	out$preds <- preds
	out$y <- y
	out$data <- data
	out$scales <- scales
	out$groups <- groups
	out$inCol <- inCol
	out$features$quadratic <- quadratic
	out$features$cubic <- cubic
	out$features$interaction <- interaction
	out$features$interQuad <- interQuad
	class(out) <- c(class(out), 'larsData')
	out
	
}
