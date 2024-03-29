% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contBoyce.r
\name{contBoyce}
\alias{contBoyce}
\title{Continuous Boyce Index (CBI) with weighting}
\usage{
contBoyce(
  pres,
  contrast,
  numBins = 101,
  binWidth = 0.1,
  presWeight = rep(1, length(pres)),
  contrastWeight = rep(1, length(contrast)),
  autoWindow = TRUE,
  method = "spearman",
  dropZeros = TRUE,
  graph = FALSE,
  na.rm = FALSE,
  bg = NULL,
  bgWeight = NULL,
  ...
)
}
\arguments{
\item{pres}{Numeric vector. Predicted values at presence sites.}

\item{contrast}{Numeric vector. Predicted values at background sites.}

\item{numBins}{Positive integer. Number of (overlapping) bins into which to divide predictions.}

\item{binWidth}{Positive numeric value < 1. Size of a bin. Each bin will be \code{binWidth * (max - min)}. If \code{autoWindow} is \code{FALSE} (the default) then \code{min} is 0 and \code{max} is 1. If \code{autoWindow} is \code{TRUE} then \code{min} and \code{max} are the maximum and minimum value of all predictions in the background and presence sets (i.e., not necessarily 0 and 1).}

\item{presWeight}{Numeric vector same length as \code{pres}. Relative weights of presence sites. The default is to assign each presence a weight of 1.}

\item{contrastWeight}{Numeric vector same length as \code{contrast}. Relative weights of background sites. The default is to assign each presence a weight of 1.}

\item{autoWindow}{Logical. If \code{FALSE} calculate bin boundaries starting at 0 and ending at 1 + epsilon (where epsilon is a very small number to assure inclusion of cases that equal 1 exactly). If \code{TRUE} (default) then calculate bin boundaries starting at minimum predicted value and ending at maximum predicted value.}

\item{method}{Character. Type of correlation to calculate. The default is \code{'spearman'}, the Spearman rank correlation coefficient used by Boyce et al. (2002) and Hirzel et al. (2006), which is the "traditional" CBI. In contrast, \code{'pearson'} or \code{'kendall'} can be used instead.  See \code{\link[stats]{cor}} for more details.}

\item{dropZeros}{Logical. If \code{TRUE} then drop all bins in which the frequency of presences is 0.}

\item{graph}{Logical. If \code{TRUE} then plot P vs E and P/E versus bin.}

\item{na.rm}{Logical. If \code{TRUE} then remove any presences and associated weights and background predictions and associated weights with \code{NA}s.}

\item{bg}{Same as \code{contrast}. Included for backwards compatibility. Ignored if \code{contrast} is not \code{NULL}.}

\item{bgWeight}{Same as \code{contrastWeight}. Included for backwards compatibility. Ignored if \code{contrastWeight} is not \code{NULL}.}

\item{...}{Other arguments (not used).}
}
\value{
Numeric value.
}
\description{
This function calculates the continuous Boyce index (CBI), a measure of model accuracy for presence-only test data. This version uses multiple, overlapping windows, in contrast to \code{link{contBoyce2x}}, which covers each point by at most two windows.
}
\details{
CBI is the Spearman rank correlation coefficient between the proportion of sites in each prediction class and the expected proportion of predictions in each prediction class based on the proportion of the landscape that is in that class.  The index ranges from -1 to 1. Values >0 indicate the model's output is positively correlated with the true probability of presence.  Values <0 indicate it is negatively correlated with the true probability of presence.
}
\examples{

set.seed(123)
pres <- sqrt(runif(100))
contrast <- runif(1000)
contBoyce(pres, contrast)
contBoyce2x(pres, contrast)
presWeight <- c(rep(1, 10), rep(0.5, 90))
contBoyce(pres, contrast, presWeight=presWeight)
contBoyce2x(pres, contrast, presWeight=presWeight)
\dontrun{
# compare stability of CBI calculated with ecospat.boyce() in ecospat package
library(ecospat)
set.seed(123)
results <- data.frame()
for (perform in c(1, 1.5, 2)) {
	for (i in 1:30) {

   pres <- runif(100)^(1 / perform)
   contrast <- runif(1000)

   cbi_enmSdm <- contBoyce(pres, contrast)
   cbi_ecospat <- ecospat.boyce(contrast, pres, PEplot=FALSE)$Spearman.cor

   results <- rbind(
     results,
     data.frame(
       performance = rep(perform, 2),
       method = c('enmSdm', 'ecospat'),
       cbi = c(cbi_enmSdm, cbi_ecospat)
     )
   )

	}

}

results$performance[results$performance == 1] <- 'poor'
results$performance[results$performance == 1.5] <- 'OK'
results$performance[results$performance == 2] <- 'good'

results$category <- paste0(results$method, '\n', results$performance)

par(mfrow=c(1, 2))
boxplot(cbi ~ category,
	data=results,
	ylab='CBI',
	main='CBI of poor, OK, and good models',
	border=c(rep('darkred', 3),
	rep('darkblue', 3))
)
plot(results$cbi,
	pch=rep(c(21, 22, 23, 24), each=2),
	contrast=ifelse(results$method == 'ecospat', 'darkred', 'cornflowerblue'),
	main='Pairs of CBIs',
	ylab='CBI'
)
legend('bottomright', fill=c('darkred', 'cornflowerblue'), legend=c('ecospat', 'enmSdm'))
}

}
\references{
Boyce, M.S., Vernier, P.R., Nielsen, S.E., and Schmiegelow, F.K.A.  2002.  Evaluating resource selection functions.  \emph{Ecological Modeling} 157:281-300. \doi{https://doi.org/10.1016/S0304-3800(02)00200-4}

Hirzel, A.H., Le Lay, G., Helfer, V., Randon, C., and Guisan, A.  2006.  Evaluating the ability of habitat suitability models to predict species presences.  \emph{Ecological Modeling} 199:142-152. \doi{10.1016/j.ecolmodel.2006.05.017}
}
\seealso{
\code{\link[stats]{cor}}, \code{link[enmSdm]{contBoyce2x}}
}
