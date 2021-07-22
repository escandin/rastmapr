\name{plot.PIFmodel}
\alias{plot.PIFmodel}
\title{PIF-model plots}
\description{This function plots linear models of Pseudo-Invariant Features
(\code{PIF} models).}
\usage{\method{plot}{PIFmodel}(x, ...)}
\arguments{
  \item{x}{\code{list}. Set of \code{PIF} models such as that produced
by \code{\link{PIFmodel}}.  }
  \item{\dots}{Additional arguments in \code{\link{plot}} other than
\code{axes}, \code{col}, \code{ylab}, \code{xlab}, or
\code{main}}
}

\value{\code{list}. Set of \code{\link{lm}} models}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
    tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
## thrs <- thresraster(strips[[2L]], strips[[1L]])
## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
## plot(model)
}
