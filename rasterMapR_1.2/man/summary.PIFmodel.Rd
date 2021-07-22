\name{summary.PIFmodel}
\alias{summary.PIFmodel}
\title{summarize PIFmodel object}
\description{A set of summaries of a \code{\link{PIFmodel}} \code{object} is obtained.}
\usage{\method{summary}{PIFmodel}(object, 
    correlation = FALSE, 
    symbolic.cor = FALSE, 
    ...)}
\arguments{
  \item{object}{an \code{object} of class \code{"lm"}, usually, a result of
a call to \code{\link{lm}}.}
  \item{correlation}{\code{logical}; if \code{TRUE}, the
\code{correlation} matrix of the estimated
parameters is returned and printed}
  \item{symbolic.cor}{\code{logical}. If \code{TRUE}, print
the correlations in a symbolic form (see
symnum) rather than as numbers.}
  \item{\dots}{further arguments passed to or from other methods.}
}



\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
          'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,4:6))
strips <- RasterIntersection(stack)
## thrs <- thresraster(strips[[2L]], strips[[1L]])
## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
summary(model)
}
