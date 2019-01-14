#' Predict a MaxEnt model object (with optional feature-level permutation)
#'
#' Takes a MaxEnt \code{lambda} object or a MaxEnt object and returns raw or logistic predictions.  Its output is the same as the function \code{raster::predict(maxentModelObject, dataFrame)} or \code{raster::predict(maxentModelObject, dataFrame, args='outputformat=raw')} (to within rounding error), and in fact those functions should be faster.  However, this function does allow custom manipulations that those functions do not allow (e.g., permuting product features while leaving other features with the same variables intact).  This function is based on Peter D. Wilson's document "Guidelines for computing MaxEnt model output values from a lambdas file".
#' @param x Either a Maxent lambda object or a Maxent model object
#' @param data Data frame. Data to which to make predictions
#' @param type Charcater.  One of: \code{'raw'} ==> produce Maxent raw values; \code{'logistic'} ==> Maxent logistic values; or \code{'cloglog'} ==> complementary log-log output (as per version 3.4.0+ of maxent--called "\code{maxnet()}" in the package of the same name)
#' @param perm Character list. Name(s) of variable to permute before calculating predictions. This permutes the variables for \emph{all} features in which they occur.  If a variable is named here, it overrides permutation settings for each feature featType.  Note that for product features the variable is permuted before the product is taken. This permutation is performed before any subsequent permutations (i.e., so if both variables in a product feature are included in \code{perms}, then this is equivalent to using the \code{'before'} rule for \code{permProdRule}). Ignored if \code{NULL}.
#' @param permLinear Character list. Names(s) of variables to permute in linear features before calculating predictions.  Ignored if \code{NULL}.
#' @param permQuad Names(s) of variables to permute in quadratic features before calculating predictions.  Ignored if \code{NULL}.
#' @param permHinge Character list. Names(s) of variables to permute in forward/reverse hinge features before calculating predictions.  Ignored if \code{NULL}.
#' @param permThresh Character list. Names(s) of variables to permute in threshold features before calculating predictions.  Ignored if \code{NULL}.
#' @param permProd Character list. A list object of \code{n} elements, each of which has two character elements naming the variables to permute if they occur in a product feature.  Depending on the value of \code{permProdRule}, the function will either permute the individual variables then calculate their product or calculate their product, then permute the product across observations.  Any other features containing the variables will produce values as normal.  Example: \code{permProd=list(c('precipWinter', 'tempWinter'), c('tempSummer', 'precipFall'))}.  The order of the variables in each element of \code{permProd} doesn't matter, so \code{permProd=list(c('temp', 'precip'))} is the same as \code{permProd=list(c('precip', 'temp'))}.  Ignored if \code{NULL}.
#' @param permProdRule Character. Rule for how permutation of product features is applied: \code{'before'} ==> Permute individual variable values then calculate product; \code{'after'} ==> calculate product then permute across these values. Ignored if \code{permProd} is \code{NULL}.
#' @param ... Extra arguments (not used).
#' @return Numeric.
#' @seealso \code{\link[dismo]{maxent}}
#' @examples
#' \donttest{
#' set.seed(123)
#' data <- matrix(rnorm(n = 6*1000), ncol = 6)
#' # true variables will be #1, #2, #5, and #6, plus
#' # the squares of #1 and #6, plus
#' # interaction between #1 and #6
#' # the cube of #5
#' imp <- c('x1', 'x2', 'x3', 'x4', 'x5', 'x6', 'x1_pow2', 'x6_pow2', 'x1_by_x6', 'x5_pow3')
#' betas <- c(5, 2, 0, 0, 1, -1, 8, 1, 2, -4)
#' names(betas) <- imp
#' y <- 0.5 + data %*% betas[1:6] + betas[7] * data[ , 1] +
#' betas[8] * data[ , 6] + betas[9] * data[ , 1] * data[ , 6] + betas[10] * data[ , 5]^3
#' y <- as.integer(y > 10)
#' data <- cbind(y, data)
#' data <- as.data.frame(data)
#' names(data) <- c('y', 'x1', 'x2', 'x3', 'x4', 'x5', 'x6')
#' 
#' x <- maxent(data[ , -1], data$y)
#' 
#' pred <- predict(x, data)
#' predNew <- predictMaxEnt(x, data, type='logistic')
#' plot(pred, predNew)
#' 
#' predPermAll <- predictMaxEnt(
#'     x=x,
#'     data=data,
#'     perm=c('x1', 'x2', 'x3', 'x4', 'x5', 'x6'),
#' 	type='logistic'
#' )
#' 
#' predPermLinear <- predictMaxEnt(x=x, data=data, permLinear='x1', type='logistic')
#' predPermQuad <- predictMaxEnt(x=x, data=data, permQuad='x1', type='logistic')
#' predPermHinge <- predictMaxEnt(x=x, data=data, permHinge='x1', type='logistic')
#' predPermThresh <- predictMaxEnt(x=x, data=data, permThresh='x1', type='logistic')
#' predPermProdBefore <- predictMaxEnt(x=x, data=data, permProd=list(c('x1', 'x6')),
#' 	permProdRule='before', type='logistic')
#' predPermProdAfter <- predictMaxEnt(x=x, data=data, permProd=list(c('x1', 'x6')),
#' 	permProdRule='after', type='logistic')
#' predPermX1 <- predictMaxEnt(x=x, data=data, perm='x1', type='logistic')
#' 
#' predPermAll <- predPermAll[order(pred)]
#' predPermLinear <- predPermLinear[order(pred)]
#' predPermQuad <- predPermQuad[order(pred)]
#' predPermHinge <- predPermHinge[order(pred)]
#' predPermThresh <- predPermThresh[order(pred)]
#' predPermProdBefore <- predPermProdBefore[order(pred)]
#' predPermProdAfter <- predPermProdAfter[order(pred)]
#' predPermX1 <- predPermX1[order(pred)]
#' predNew <- predNew[order(pred)]
#' pred <- pred[order(pred)]
#' 
#' plot(pred, ylab='Prediction', xlab='Index (ordered from low to high unpermuted prediction)')
#' points(predNew, pch=3, col=2, cex=2)
#' points(predPermLinear, pch=2, col=3)
#' points(predPermQuad, pch=4, col=4)
#' points(predPermHinge, pch=5, col=5)
#' points(predPermThresh, pch=6, col=6)
#' points(predPermProdBefore, pch=7, col=7)
#' points(predPermProdAfter, pch=8, col=8)
#' points(predPermAll, pch=9, col=9)
#' points(predPermX1, pch=10, col=10)
#' 
#' legend('topleft',
#' 	inset=0.01,
#' 	pch=c(1, 3, 2, 4, 5, 6, 7, 8, 9, 10),
#' 	col=1:10,
#' 	bg='white',
#' 	legend=c('raster::predict()',
#' 		'enmSdm::predictMaxEnt()',
#' 		'permuted x1 linear',
#' 		'permuted x1 quadratic',
#' 		'permuted x1 hinge',
#' 		'permuted x1 threshold',
#' 		'permuted x1 product (before)',
#' 		'permuted x1 product (after)',
#' 		'permuted (all variables/features)',
#' 		'permuted x1 all features'
#' 	)
#' )
#' }
#' @export
predictMaxEnt <- function(
	x,
	data,
	type='cloglog',
	perm=NULL,
	permLinear=NULL,
	permQuad=NULL,
	permHinge=NULL,
	permThresh=NULL,
	permProd=NULL,
	permProdRule=NULL,
	...
) {

	####################
	## PRE-PROCESSING ##
	####################

	# get lambdas is object is a model object
	if (class(x)=='MaxEnt') x <- x@lambdas

	##########
	## MAIN ##
	##########

	### create meta value data frame
	################################
	metaLines <- x[which(grepl(x, pattern='linearPredictorNormalizer')):length(x)]
	meta <- data.frame(
		item=substr(metaLines, 1, unlist(gregexpr(metaLines, pattern=',')) - 1),
		value=as.numeric(substr(metaLines, unlist(gregexpr(metaLines, pattern=',')) + 2, nchar(metaLines)))
	)

	### create data frame of feature functions
	##########################################
	feats <- data.frame(
		feature=x[1:(which(grepl(x, pattern='linearPredictorNormalizer')) - 1)],
		featType='linear',
		var1=NA,
		var2=NA,
		crit=NA,
		lambda=NA,
		beta1=NA,
		beta2=NA
	)

	feats$feature <- as.character(feats$feature)
	feats$featType <- as.character(feats$featType)
	feats$crit <- as.numeric(feats$crit)
	feats$lambda <- as.numeric(feats$lambda)
	feats$beta1 <- as.numeric(feats$beta1)
	feats$beta2 <- as.numeric(feats$beta2)

	## identify feature types
	feats$featType[which(grepl(feats$feature, pattern='[(]') & grepl(feats$feature, pattern='[)]'))] <- 'threshold'
	feats$featType[which(grepl(feats$feature, pattern='\\^'))] <- 'quadratic'
	feats$featType[which(grepl(feats$feature, pattern='[*]'))] <- 'product'
	feats$featType[which(grepl(feats$feature, pattern='\''))] <- 'forward hinge'
	feats$featType[which(grepl(feats$feature, pattern='`'))] <- 'reverse hinge'

	## identify (first) variable used
	# see which features have this variable string... note that there is a risk that the parts of a "feature" that have the *substring* including the variable name will falsely imply the variable is in the feature... sorting by length should obviate this
	names <- sort(names(data))
	names <- names[order(nchar(names))] # sort variable names by length
	for (i in 1:length(names)) {

		varIn <- grepl(feats$feature, pattern=names[i]) & is.na(feats$var1)
		feats$var1[which(varIn)] <- names[i]
		
	}

	## identify second variable used (if any)... same risks as above
	if ('product' %in% feats$featType) {
		indexProd <- which(feats$featType=='product')
		for (i in indexProd) {
			for (j in (which(feats$var1[i] == names) + 1):length(names)) {
				secondVarIn <- grepl(feats$feature[i], pattern=names[j])
				if (secondVarIn) feats$var2[i] <- names[j]
			}
		}
	}

	## extract "betas"
	commaPos <- unlist(gregexpr(pattern=',', feats$feature))
	feats$beta1 <- as.numeric(substr(feats$feature, start=commaPos[seq(2, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(3, 3 * nrow(feats), by=3)] - 1))
	feats$beta2 <- as.numeric(substr(feats$feature, start=commaPos[seq(3, 3 * nrow(feats), by=3)] + 2, stop=nchar(feats$feature)))

	## extract critical values (hinge/threshold)
	feats$crit[which(feats$featType=='threshold')] <- as.numeric(substr(feats$feature[which(feats$featType=='threshold')], start=2, stop=unlist(gregexpr(feats$feature[which(feats$featType=='threshold')], pattern='<')) - 1))
		
	feats$crit[which(feats$featType=='forward hinge')] <- feats$beta1[which(feats$featType=='forward hinge')]
	feats$crit[which(feats$featType=='reverse hinge')] <- feats$beta2[which(feats$featType=='reverse hinge')]

		
	## extract lambdas
	feats$lambda <- as.numeric(substr(feats$feature, start=commaPos[seq(1, 3 * nrow(feats), by=3)] + 2, stop=commaPos[seq(2, 3 * nrow(feats), by=3)] - 1))

	## get values
	value1 <- t(data[ , feats$var1])
	rownames(value1) <- paste0('feat', 1:nrow(feats))
	value2 <- value1 * NA
		
	if ('product' %in% feats$featType) {
		value2Compact <- t(data[ , feats$var2[which(feats$featType=='product')]])
		value2[which(feats$featType=='product'), ] <- value2Compact
	}

	# value1 <<- value1
	# value2 <<- value2
	# feats <<- feats
	# print(feats)

	### permutations
	################

	## permute values REGARDLESS OF FEATURE TYPE
	if (!is.null(perm)) {

		for (i in 1:nrow(feats)) {
			if (any(perm %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
			if (any(perm %in% feats$var2[i])) value2[i, ] <- sample(value2[i, ], ncol(value2))
		}
			
	}

	## permute values in LINEAR features
	if (!is.null(permLinear) & 'linear' %in% feats$featType) {
		for (i in which(feats$featType %in% 'linear')) {
			if (any(permLinear %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in QUADRATIC features
	if (!is.null(permQuad) & 'quadratic' %in% feats$featType) {
		for (i in which(feats$featType %in% 'quadratic')) {
			if (any(permQuad %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in HINGE features
	if (!is.null(permHinge) & ('forward hinge' %in% feats$featType | 'reverse hinge' %in% feats$featType)) {
		for (i in sort(unique(c(which(feats$featType %in% 'forward hinge'), which(feats$featType %in% 'reverse hinge'))))) {
			if (any(permHinge %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in THRESHOLD features
	if (!is.null(permThresh) & 'threshold' %in% feats$featType) {
		for (i in which(feats$featType %in% 'threshold')) {
			if (any(permThresh %in% feats$var1[i])) value1[i, ] <- sample(value1[i, ], ncol(value1))
		}
	}

	## permute values in PRODUCT features
	if (!is.null(permProd) & 'product' %in% feats$featType) {

		# permute variables then calculate product
		if (permProdRule=='before') {
		
			value1Perm <- value1
			value2Perm <- value2
		
			for (i in which(feats$featType %in% 'product')) {

				# see if permutation is wanted for this particular var1 and var2 of this product feature
				wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
				
				if (sum(wantPairPerm)==2) {
					value1Perm[i, ] <- sample(value1Perm[i, ], ncol(value1Perm))
					value2Perm[i, ] <- sample(value2Perm[i, ], ncol(value2Perm))
				}

			}

			# calculate product
			prodMat <- value1Perm * value2Perm
		
		# calculate product then permute product
		} else if (permProdRule=='after') {
		
			# calculate product
			prodMat <- value1 * value2

			for (i in which(feats$featType %in% 'product')) {

				# see if permutation is wanted for this particular var1 and var2 of this product feature
				wantPairPerm <- c(sapply(X=permProd, FUN=function(x, y) {x %in% y}, c(feats$var1[i], feats$var2[i])))
				if (sum(wantPairPerm)==2) prodMat[i, ] <- sample(prodMat[i, ], ncol(prodMat))
				
			}
			
		}
			
	} else {

		# calculate unpermuted product
		prodMat <- value1 * value2
		
	}

	### calculate predictions
	#########################

	# initialize exponent
	S <- numeric(nrow(data))

	# for each datum
	for (i in 1:nrow(data)) {

		# for each featType of feature
		for (featType in c('linear', 'quadratic', 'product', 'forward hinge', 'reverse hinge', 'threshold')) {

			# if feature featType occurs in lambdas, evaluate
			if (featType %in% feats$featType) {

				thisValue1 <- value1[which(feats$featType==featType), i]
				thisValue2 <- value2[which(feats$featType==featType), i]
				if ('product' %in% feats$featType) thisProd <- prodMat[which(feats$featType==featType), i]
			
				thisLambda <- feats$lambda[which(feats$featType==featType)]
				thisCrit <- feats$crit[which(feats$featType==featType)]
				
				thisBeta1 <- feats$beta1[which(feats$featType==featType)]
				thisBeta2 <- feats$beta2[which(feats$featType==featType)]
				
				thisRange <- thisBeta2 - thisBeta1
				
				if (featType=='linear') S[i] <- S[i] + sum(thisLambda * (thisValue1 - thisBeta1) / thisRange, na.rm=T)
				if (featType=='quadratic') S[i] <- S[i] + sum(thisLambda * (thisValue1^2 - thisBeta1) / thisRange, na.rm=T)
				if (featType=='product') S[i] <- S[i] + sum(thisLambda * (thisProd - thisBeta1) / thisRange, na.rm=T)
				# if (featType=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=T)
				
				
				if (featType=='forward hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 >= thisCrit] * (thisValue1[thisValue1 >= thisCrit] - thisCrit[thisValue1 >= thisCrit]) / thisRange[thisValue1 >= thisCrit], na.rm=T)
				if (featType=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=T)
				# if (featType=='reverse hinge') S[i] <- S[i] + sum(thisLambda[thisValue1 < thisCrit] * (thisCrit[thisValue1 < thisCrit] - thisValue1[thisValue1 < thisCrit]) / thisRange[thisValue1 < thisCrit], na.rm=T)
				
				
				if (featType=='threshold') S[i] <- S[i] + sum(thisLambda[which(thisValue1 > thisCrit)], na.rm=T)

			}
			
		}

	}

	raw <- exp(S - meta$value[match('linearPredictorNormalizer', meta$item)]) / meta$value[match('densityNormalizer', meta$item)]

	if (type=='raw') {
		out <- raw
	} else if (type=='logistic') {
		out <- (raw * exp(meta$value[match('entropy', meta$item)])) / (1 + raw * exp(meta$value[match('entropy', meta$item)]))
	} else if (type=='cloglog') {
		entropy <- x[length(x)]
		entropy <- as.numeric(substr(entropy, regexpr(',', entropy) + 2, nchar(entropy)))
		out <- 1 - exp(-exp(entropy) * raw)
	}

	out

}

