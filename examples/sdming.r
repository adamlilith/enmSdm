set.seed(123)

### contrived example
n <- 10000
x1 <- seq(-1, 1, length.out=n) + rnorm(n)
x2 <- seq(10, 0, length.out=n) + rnorm(n)
x3 <- rnorm(n)
y <- 2 * x1 + x1^2 - 10 * x2 - x1 * x2

y <- statisfactory::probitAdj(y, 0)
y <- y^3
hist(y)

presAbs <- runif(n) < y

trainData <- data.frame(presAbs=presAbs, x1=x1, x2=x2, x3=x3)

model <- trainGlm(trainData)

source('C:/Ecology/Drive/R/enmSdm/working/trainByCrossValid.r')
out <- trainByCrossValid(data=trainData, verbose=2)




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
#' model <- trainGlm(trainData)
#' 
#' out <- trainByCrossValid(data=trainData, verbose=2)
