### 
library(enmSdm)
library(omnibus)
library(terra)
library(tictoc)

rasts <- rast(listFiles('E:/Ecology/Drive/Research/ENMs - Calibration Region for ENMs & SDMs/Climate/+003 Century')[1:6])
newdata <- as.data.frame(rasts)
newdata <- newdata[rep(1:nrow(newdata), each=3), ]

source('E:/Ecology/Drive/R/enmSdm/R/predictEnmSdm.r')

load('E:/Ecology/Drive/Research/ENMs - Calibration Region for ENMs & SDMs/03 ENDMs/CLIMATE Dynamic BG REGION Dynamic CALIB DIST NA DISPERSAL DIST 1.5 LAMBDA 4/MaxEnt Models for Species 001.rda')

load('E:/Ecology/Drive/Research/ENMs - Calibration Region for ENMs & SDMs/02 Background Data with Folds/CLIMATE Dynamic BG REGION Dynamic CALIB DIST NA DISPERSAL DIST 1.5 LAMBDA 4/Background Sites with Folds Species 001.rda')

load('E:/Ecology/Drive/Research/ENMs - Calibration Region for ENMs & SDMs/00 Simulations/Niches/Niche Species 001.rda')

pca <- niche$pca

model <- models$bioclimUberModel

tic()
p2 <- predictEnmSdm(
	model = model,
	newdata = newdata,
	maxentFun='dismo',
	cores = 4
)
toc()

tic()
p1 <- predictEnmSdm(
	model = model,
	newdata = newdata,
	maxentFun='dismo',
	cores = 1
)
toc()

# pc2 <- predictEnmSdm(
	# model = pca,
	# newdata = newdata,
	# cores = 4
# )

# pc1 <- predictEnmSdm(
	# model = pca,
	# newdata = newdata,
	# cores = 1
# )


