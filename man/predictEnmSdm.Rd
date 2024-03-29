% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predictEnmSdm.r
\name{predictEnmSdm}
\alias{predictEnmSdm}
\title{Generic predict function for LMs, GLMs, GAMs, RFs, BRTs, CRFs, Maxent, and more}
\usage{
predictEnmSdm(
  model,
  newdata,
  maxentFun = "dismo",
  cores = 1,
  nrows = nrow(newdata),
  ...
)
}
\arguments{
\item{model}{Object of class \code{lm}, \code{glm}, \code{gam}, \code{randomForest}, \code{MaxEnt}, \code{MaxNet}, \code{prcomp}, \code{kde}, \code{gbm}, and possibly others (worth a try!).}

\item{newdata}{Data frame or matrix with data to which to predict}

\item{maxentFun}{This argument is only used if the \code{model} object is a MaxEnt model; otherwise, it is ignored. I takes a value of either \code{'dismo'}, in which case a MaxEnt model is predicted using the default \code{predict} function from the \pkg{dismo} package, or \code{'enmSdm'} in which case the function \code{\link[enmSdm]{predictMaxEnt}} function from the \pkg{enmSdm} package (this package) is used.}

\item{cores}{Number of cores to use. The default is 1. If >1 and \code{newdata} is a raster, only 1 core is used (i.e., basically, \code{cores} is ignored if you're writing to a raster... sorry!)}

\item{nrows}{Number of rows of \code{newdata} to predict at a time (assuming \code{newdata} is a \code{data.frame} or \code{matrix}). The default value is to predict all rows at once, but for very large data frames/matrices this can lead to memory issues in some cases. By setting the number of rows, \code{newdata} can be divided into chunks, and predictions made to each chunk, which may ease memory limitations. This can be combined with multi-coring (which will increase memory requirements). In this case, all cores combined will get \code{nrows} of data. How many rows are too many? You will have to decide depending on your data and the output. For example, predicting the outcome of a GLM on data with 10E6 rows ma be fine, but predicting a PCA (with multiple axes) to the data data may require too much memory. You can use \code{\link[omnibus]{memUse}} to see to help figure this out.}

\item{...}{Arguments to pass to the algorithm-specific \code{predict} function.}
}
\value{
Numeric.
}
\description{
This is a generic predict function that automatically uses the model common arguments for predicting models of the following types: linear models, generalized linear models (GLMs), generalized additive models (GAMs), random forests, boosted regression trees (BRTs)/gradient boosting machines (GBMs), conditional random forests, Maxent, and more.
}
\seealso{
\code{\link[stats]{predict}} from the \pkg{stats} package, \code{\link[dismo]{predict}} from the \pkg{dismo} package, \code{\link[raster]{predict}} from the \pkg{raster} package, \code{\link[terra]{predict}} from the \pkg{terra} package
}
