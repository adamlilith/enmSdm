##
source('E:/Ecology/Drive/R/enmSdm/R/bioticVelocity.r')

library(raster)

### movement in north-south directions
mat <- matrix(0, nrow=5, ncol=5)
mat1 <- mat2 <- mat
mat1[3, 3] <- 1
mat2[2, 3] <- 1

mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
rownames(mats) <- paste0('lat', nrow(mats):1)
colnames(mats) <- paste0('long', 1:ncol(mats))
lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)


### movement in east-west directions
mat <- matrix(0, nrow=5, ncol=5)
mat1 <- mat2 <- mat
mat1[3, 3] <- 1
mat2[3, 2] <- 1

mats <- array(c(mat1, mat2), dim=c(nrow(mat1), ncol(mat1), 2))
rownames(mats) <- paste0('lat', nrow(mats):1)
colnames(mats) <- paste0('long', 1:ncol(mats))
lats <- matrix(nrow(mats):1, nrow=nrow(mats), ncol=ncol(mats))
longs <- matrix(1:nrow(mats), nrow=nrow(mats), ncol=ncol(mats), byrow=TRUE)

mat1rast <- raster(mat1)
mat2rast <- raster(mat2)
matsRast <- stack(mat1rast, mat2rast)
plot(matsRast, col=c('gray', 'darkgreen'))

# movement east
(bioticVelocity(mats, times=1:2,
latitude=lats, longitude=longs, cores=2))

