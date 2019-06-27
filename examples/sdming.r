set.seed(123)

### contrived example
n <- 10000
x1 <- seq(-1, 1, length.out=n) + rnorm(n)
x2 <- seq(10, 0, length.out=n) + rnorm(n)
x3 <- rnorm(n)
y <- 2 * x1 + x1^2 - 10 * x2 - x1 * x2

y <- statisfactory::probitAdj(y, 0)
y <- y^3

presAbs <- runif(n) < y

trainData <- data.frame(presAbs=presAbs, x1=x1, x2=x2, x3=x3)

# # source('C:/Ecology/Drive/R/enmSdm/R/trainGlm.r')
# model <- trainGlm(trainData, out=c('model', 'models', 'tuning'))

source('C:/Ecology/Drive/R/enmSdm/R/trainMaxEnt.r')
model <- trainMaxEnt(trainData, out=c('model', 'models', 'tuning'), verbose=TRUE, regMult=1:2, classes='lp', dropOverparam=FALSE)

source('C:/Ecology/Drive/R/enmSdm/R/trainByCrossValid.r')
folds <- dismo::kfold(trainData, 3)
out <- trainByCrossValid(data=trainData, folds=folds, verbose=1)




#' set.seed(123)
#' 
#' ### contrived example
#' n <- 10000
#' x1 <- seq(-1, 1, length.out=n) + rnorm(n)
#' x2 <- seq(10, 0, length.out=n) + rnorm(n)
#' x3 <- rnorm(n)
#' y <- 2 * x1 + x1^2 - 10 * x2 - x1 * x2
#' 
#' y <- statisfactory::probitAdj(y, 0)
#' y <- y^3
#' hist(y)
#' 
#' presAbs <- runif(n) < y
#' 
#' trainData <- data.frame(presAbs=presAbs, x1=x1, x2=x2, x3=x3)
#' 
#' model <- trainMaxEnt(trainData, regMult=1:2, classes='lqp')
#' summary(model)
#' 
#' out <- trainByCrossValid(data=trainData, verbose=1)
#' out$tuning
