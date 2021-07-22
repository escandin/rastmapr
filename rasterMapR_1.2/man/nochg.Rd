\name{nochg}
\alias{nochg}
\title{Change and no change mask}
\description{This function computes a binary change/no-change mask processing a
threshold raster, see \code{\link{thresraster}}}
\usage{nochg(thresraster, degfree, 
    pvalue = 1e-04)}
\arguments{
  \item{thresraster}{\code{RasterLayer}, such as that produced by
\code{\link{thresraster}}.}
  \item{degfree}{\code{numeric}. Degrees of freedom.}
  \item{pvalue}{\code{numeric}. Probability threshold}
}
\details{This function is implemented by
\code{\link{PIFmodel}}. The mask is
computed comparing values in the threshold
raster with the chi-square statistic, see
\code{\link{qchisq}}.}
\value{\code{RasterLayer}. }

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
thrs <- thresraster(strips[[2L]], strips[[1L]])
noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
## plot(noch)
}
