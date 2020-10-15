# source('E:/Ecology/Drive/R/enmSdm/working/bioticVelocityTest.r')

library(raster)
source('E:/Ecology/Drive/R/enmSdm/R/bioticVelocity.r')
# source('C:/Ecology/Drive/R/enmSdm/R/private_euclid.r')
# source('C:/Ecology/Drive/R/enmSdm/R/private_cardinalDistance.r')
# source('C:/Ecology/Drive/R/enmSdm/R/private_interpolateLatFromMatrix.r')
# source('C:/Ecology/Drive/R/enmSdm/R/private_interpolateLongFromMatrix.r')

# metrics <- c("centroid", "nsCentroid", "ewCentroid", "nCentroid", "sCentroid", "eCentroid", "wCentroid", "nsQuants", "ewQuants", "mean", "quants", "prevalence", 'similarity')

# ### movement in north-south directions
# mat <- matrix(0, nrow=5, ncol=5)
# mat1 <- mat2 <- mat
# mat1[3, 3] <- 1
# mat2[2, 3] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # note that nCentroidVelocity is NaN because just one cell is occupied
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs)) # north movement
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs)) # south movement


# ### movement in east-west directions
# mat <- matrix(0, nrow=5, ncol=5)
# mat1 <- mat2 <- mat
# mat1[3, 3] <- 1
# mat2[3, 2] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # movement east
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))

# # movement west
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))


# ### movement of portions of range northward/southward
# mat <- matrix(0, nrow=11, ncol=11)
# mat1 <- mat2 <- mat
# mat1[6, 5] <- 1 # bottom
# mat1[5, 5] <- 1 # center
# mat1[5, 4] <- 1 # west
# mat1[5, 6] <- 1 # east
# mat1[4, 5] <- 1 # north

# mat2[6, 5] <- 1 # bottom
# mat2[5, 5] <- 1 # center
# mat2[5, 4] <- 1 # west
# mat2[5, 6] <- 1 # east
# mat2[3, 5] <- 1 # north

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # northern section moves north
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))
# # northern section moves south
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))

# ### quantile velocities: north/south movement
# mat <- matrix(0, nrow=11, ncol=11)
# mat1 <- mat2 <- mat

# mat1[2:10, 6] <- 1
# mat2[1:9, 6] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # shift north
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))
# # shift south
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))


# ### quantile velocities: east/west movement
# mat <- matrix(0, nrow=11, ncol=11)
# mat1 <- mat2 <- mat

# mat1[6, 2:10] <- 1
# mat2[6, 3:11] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # eastward shift
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))
# # westward shift
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))

# ### big block test
# mat <- matrix(0, nrow=7, ncol=7)
# mat1 <- mat2 <- mat

# mat1[3:5, 3:5] <- 1
# mat2[1:3, 1:3] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# # shift northwest
# (bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))
# # shift southeast
# (bioticVelocity(mats[ , , 2:1], metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs))


# ### big block test: common frame
# mat <- matrix(0, nrow=7, ncol=7)
# mat1 <- mat2 <- mat

# mat1[3:5, 3:5] <- 1
# mat1[1, ] <- NA
# mat1[ , 1] <- NA
# mat2[1:3, 1:3] <- 1

# mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
# rownames(mats) <- paste0('lat', nrow(mats):1)
# colnames(mats) <- paste0('long', 1:ncol(mats))
# lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
# longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

# mat1rast <- raster(mat1)
# mat2rast <- raster(mat2)
# matsRast <- stack(mat1rast, mat2rast)
# plot(matsRast, col=c('gray', 'darkgreen'))

# v1 <- bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs)
# v2 <- bioticVelocity(mats, metrics=metrics, times=1:2,
	# latitude=lats, longitude=longs, onlyInSharedCells=TRUE)
# (rbind(v1, v2))


### elevational centroid velocity

# single cell moving
# goes up, holds, then down
mat <- matrix(0, nrow=7, ncol=7)
mat1 <- mat2 <- mat3 <- mat4 <- elevation <- mat

mat1[4, 4] <- 1
mat2[3, 4] <- 1
mat3[3, 4] <- 1
mat4[2, 4] <- 1

elevation[4, 4] <- 1
elevation[3, 4] <- 4
elevation[2, 4] <- 2

mats <- array(c(mat1, mat2, mat3, mat4), dim=c(nrow(mat1), ncol(mat1), 4))
rownames(mats) <- paste0('lat', nrow(mats):1)
colnames(mats) <- paste0('long', 1:ncol(mats))
lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

mat1rast <- raster(mat1)
mat2rast <- raster(mat2)
mat3rast <- raster(mat3)
mat4rast <- raster(mat4)

rownames(elevation) <- paste0('lat', nrow(mats):1)
colnames(elevation) <- paste0('long', 1:ncol(mats))

elevRast <- raster(elevation)
plot(stack(elevRast, mat1rast, mat2rast, mat3rast, mat4rast))

(elevVel <- bioticVelocity(mats, metrics='elevCentroid',
times=1:4, latitude=lats, longitude=longs, elevation=elevation))

### elevational centroid quantiles

# quantiles mostly move up, hold, then down
mat <- matrix(0, nrow=4, ncol=4)
mat1 <- elevation <- mat

mat1[] <- runif(16)
mat2[] <- sort(mat1)
mat3 <- mat2
mat4 <- matrix(rev(mat3), nrow=4)

elevation[] <- 1:16

mats <- array(c(mat1, mat2, mat3, mat4), dim=c(nrow(mat1), ncol(mat1), 4))
rownames(mats) <- paste0('lat', nrow(mats):1)
colnames(mats) <- paste0('long', 1:ncol(mats))
lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

mat1rast <- raster(mat1)
mat2rast <- raster(mat2)
mat3rast <- raster(mat3)
mat4rast <- raster(mat4)

names(mat1rast) <- 'time1'
names(mat2rast) <- 'time2'
names(mat3rast) <- 'time3'
names(mat4rast) <- 'time4'

rownames(elevation) <- paste0('lat', nrow(mats):1)
colnames(elevation) <- paste0('long', 1:ncol(mats))

elevRast <- raster(elevation)
names(elevRast) <- 'elev'

plot(stack(elevRast, mat1rast, mat2rast, mat3rast, mat4rast))

(elevVel <- bioticVelocity(mats, metrics='elevQuants',
times=1:4, latitude=lats, longitude=longs, elevation=elevation))


