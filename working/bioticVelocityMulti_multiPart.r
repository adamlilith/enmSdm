
tic()		
		if (exists('work')) rm(work)
		cores <- 1

		# divvy up time periods among cores
		repAtTimes <- repIndicesFrom <- list()
		repSize <- floor(length(atTimes) / cores)
		numReps <- floor(length(atTimes) / repSize)
		for (i in 1:numReps) {
			
			extra <- ifelse(i > 1, 1, 0)
			repAtTimes[[i]] <- atTimes[(1 + (i - 1) * repSize - extra):(i * repSize)]
			repIndicesFrom[[i]] <- (1 + (i - 1) * repSize - extra):(i * repSize)
		
		}
		
		# add times rounded out
		if (tail(repAtTimes[[length(repAtTimes)]], 1) < tail(atTimes, 1)) {
			repAtTimes[[numReps + 1]] <- atTimes[(i * repSize)]:tail(atTimes, 1)
			repIndicesFrom[[numReps + 1]] <- (i * repSize):length(atTimes)
		}

		# multi-core
		if (cores > 1) {
			`%makeWork%` <- foreach::`%dopar%`
			cl <- parallel::makeCluster(cores)
			doParallel::registerDoParallel(cl)
		} else {
			`%makeWork%` <- foreach::`%do%`
		}
		
		mcOptions <- list(preschedule=TRUE, set.seed=FALSE, silent=FALSE)
		
		export <- c('.euclid', '.cardinalDistance', '.interpolateLatFromMatrix', '.interpolateLongFromMatrix')
		
		work <- foreach::foreach(i=seq_along(repAtTimes), .options.multicore=mcOptions, .combine='rbind', .inorder=FALSE, .export=export, .packages = c('raster'), .verbose=TRUE) %makeWork%
			bioticVelocity(
				x = raster::subset(x, repIndicesFrom[[i]]),
				times = repAtTimes[[i]],
				atTimes = repAtTimes[[i]],
				longitude = longitude,
				latitude = latitude,
				elevation = elevation,
				metrics = metrics,
				quants = quants,
				onlyInSharedCells = onlyInSharedCells,
				cores = 1,
				warn = FALSE#,
				# ...
			)
				
		if (cores > 1) parallel::stopCluster(cl)

		work
toc()
