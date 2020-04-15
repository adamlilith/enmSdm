# source('C:/Ecology/Drive/R/enmSdm/working/work_parallel.r')

scratchDir <- 'D:/ecology/!Scratch/_temp/'

	library(doParallel)

### model red-bellied lemurs
data(mad0)
data(lemurs)

# climate data
bios <- c(1, 5, 12, 15)
clim <- raster::getData('worldclim', var='bio', res=10)
clim <- raster::subset(clim, bios)
clim <- raster::crop(clim, mad0)

# occurrence data
occs <- lemurs[lemurs$species == 'Eulemur rubriventer', ]
occsEnv <- raster::extract(clim, occs[ , c('longitude', 'latitude')])
occsEnv <- as.data.frame(occsEnv) # need to do this for prediction later

# background sites
bg <- 2000 # too few cells to locate 10000 background points
bgSites <- dismo::randomPoints(clim, 2000)
bgEnv <- extract(clim, bgSites)

# collate
presBg <- rep(c(1, 0), c(nrow(occs), nrow(bgSites)))
env <- rbind(occsEnv, bgEnv)
env <- cbind(presBg, env)
env <- as.data.frame(env)

preds <- paste0('bio', bios)	

regMult = c(seq(0.5, 5, by = 0.5), 7.5, 10)
classes = 'default'	
	
	data <- env
	
	###########
	## setup ##
	###########

	if (!is.null(anyway)) dropOverparam <- anyway
	
	# response and predictors
	if (class(resp) %in% c('integer', 'numeric')) resp <- names(data)[resp]
	if (class(preds) %in% c('integer', 'numeric')) preds <- names(data)[preds]

	# get response and predictors
	presentBg <- data[ , resp]
	data <- data[ , preds, drop=FALSE]

	### get combinations of features to test for each regularization multiplier
	classesToTest <- if (classes == 'default') {
		c('l', 'p', 'q', 'h')
	} else {
		unlist(strsplit(classes, ''))
	}

	if (any('p' %in% classesToTest) & ncol(data) == 1) {
		product <- FALSE
		warning('Data has only one variable so forcing product features to FALSE.')
	}

	## collate all presences
	allPres <- data[presentBg == 1, , drop=FALSE]

	## collate all background sites
	allBg <- data[presentBg == 0, , drop=FALSE]

	# create df of 1/0 to indicate each combination of classes to test
	if (testClasses) {
		classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
		classGrid <- classGrid[-which(rowSums(classGrid) == 0), , drop=FALSE]
	} else {
		classGrid <- data.frame(matrix(rep(1, length(classesToTest)), nrow=1))
	}

	names(classGrid) <- classesToTest

	if (forceLinear & any(classGrid$l == 0)) classGrid <- classGrid[-which(classGrid$l == 0), , drop=FALSE]

	##########
	## MAIN ##
	##########

	if (verbose) omnibus::say('Evaluating models with regularization multiplier:', post=0)

	# remember all models and evaluation data
	models <- list()
	tuning <- data.frame()

	tuning <- data.frame()
	
	# by BETA
	for (thisRegMult in regMult) {

		# by FEATURE COMBINATION
		for (countParam in 1:nrow(classGrid)) {

			# classes
			classes <- paste0(
				ifelse('l' %in% names(classGrid) && classGrid$l[countParam] == 1, 'l', ''),
				ifelse('p' %in% names(classGrid) && classGrid$p[countParam] == 1, 'p', ''),
				ifelse('q' %in% names(classGrid) && classGrid$q[countParam] == 1, 'q', ''),
				ifelse('h' %in% names(classGrid) && classGrid$h[countParam] == 1, 'h', ''),
				ifelse('t' %in% names(classGrid) && classGrid$t[countParam] == 1, 't', '')
			)
			
			tuning <- rbind(
				tuning,
					data.frame(
					regMult=thisRegMult,
					numPres=sum(presentBg),
					classes=classes
				)
			)
			
		}
		
	}
	
	cores <- 4
	cl <- makeCluster(cores)
	registerDoParallel(cl)

	mcOptions <- list(preschedule=TRUE, set.seed=FALSE, silent=FALSE)
	
		
	.trainMaxEntWorker <- function(
		i,
		scratchDir,
		tuning,
		presentBg,
		allPres,
		allBg
	) {
			
			thisRegMult <- tuning$regMult[i]
			thisClasses <- tuning$classes[i]
			
			# create scratch directory
			thisScratchDir <- if (is.null(scratchDir)) {
				base::tempfile(pattern='/_maxentTempFiles/')
			} else {
				base::tempfile(pattern='/_maxentTempFiles/', tmpdir=scratchDir)
			}
			
			omnibus::dirCreate(thisScratchDir, '/plots')

			# get parameters
			params <- c(
				paste0('betamultiplier=', thisRegMult),
				paste0('linear=', ifelse(grepl('l', thisClasses), 'true', 'false')),
				paste0('product=', ifelse(grepl('p', thisClasses), 'true', 'false')),
				paste0('quadratic=', ifelse(grepl('q', thisClasses), 'true', 'false')),
				paste0('hinge=', ifelse(grepl('h', thisClasses), 'true', 'false')),
				paste0('threshold=', ifelse(grepl('t', thisClasses), 'true', 'false')),
				paste0('jackknife=', ifelse(jackknife, 'true', 'false'))
			)

			if (args != '') params <- c(params, args)

			# train model
			model <- dismo::maxent(
				x=data,
				p=as.vector(presentBg),
				path=thisScratchDir,
				args=params
			)

			## predict to training (and maybe test presences)
			predPres <- dismo::predict(
				object=model,
				x=allPres,
				na.rm=TRUE,
				args='outputformat=raw'
			)

			## predict to background
			predBg <- dismo::predict(
				object=model,
				x=allBg,
				na.rm=TRUE,
				args='outputformat=raw'
			)

			rawSum <- sum(c(predPres, predBg), na.rm=TRUE)

			## log likelihood
			ll <- sum(log(predPres / rawSum), na.rm=TRUE)

			## number of parameters
			numCoeff <- 0
		
			for (thisLambda in model@lambdas) { # for each line in lambda object

				# if not a meta-data line
				if (!grepl(thisLambda, pattern='linearPredictorNormalizer') & !grepl(thisLambda, pattern='densityNormalizer') & !grepl(thisLambda, pattern='numBackgroundPoints') & !grepl(thisLambda, pattern='entropy')) {

					split <- strsplit(thisLambda, ', ')
					paramValue <- as.numeric(split[[1]][2])
					if (paramValue !=0) numCoeff <- numCoeff + 1 # increment number of parameters

				}

			}

			# AICc
			AICc <- -2 * ll + 2 * numCoeff + (2 * numCoeff * (numCoeff + 1)) / (sum(presentBg) - numCoeff - 1)

			# remove temporary files... note that "species.lambda" file cannot be removed unless R is closed, so we'll just make it smaller to reduce disk space usage
			write.csv(NULL, paste0(thisScratchDir, '/species.lambdas'))
			if (file.exists(paste0(thisScratchDir, '/presences'))) write.csv(NULL, paste0(scratchDir, '/', i, '/presences'))
			if (file.exists(paste0(thisScratchDir, '/absences'))) write.csv(NULL, paste0(scratchDir, '/', i, '/absences'))
			unlink(paste0(thisScratchDir, '/plots'), recursive=TRUE, force=TRUE)
			unlink(thisScratchDir, recursive=TRUE, force=TRUE)

			out <- list(list(
				model=model,
				thisTuning=data.frame(
					regMult=thisRegMult,
					classes=thisClasses,
					numClasses=nchar(thisClasses),
					numCoeff=numCoeff,
					logLik=ll,
					AICc=AICc
				)
			))
		
		}
		
	
	work <- foreach(i=1:nrow(tuning), .options.multicore=mcOptions, .combine='c', .packages = c('omnibus', 'dismo', 'rJava')) %dopar% 
		.trainMaxEntWorker(
			i = i,
			scratchDir = scratchDir,
			tuning = tuning,
			presentBg = presentBg,
			allPres = allPres,
			allBg = allBg
		)
			
	stopCluster(cl)
			
	