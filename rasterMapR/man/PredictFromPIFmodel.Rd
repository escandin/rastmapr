\name{PredictFromPIFmodel}
\alias{PredictFromPIFmodel}
\title{DEPRECATED PIF-model prediction}
\description{Use \code{\link{CalibrateRaster}} instead. The function predicts
hyper-spectral data using PIF model}
\usage{PredictFromPIFmodel(stacks, 
    pif.mod)}
\arguments{
  \item{stacks}{\code{list}. Set of \code{RasterBrick} data such as that
produced by \code{\link{EEstackWithoutMeta}}.}
  \item{pif.mod}{\code{list}. PIF models such as these
produced by
\code{\link{PIFmodel}}.}
}

\value{\code{RasterBrick}. Predicted features.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
model <- PIFmodel(strips, pvalue = 3E-1, brackets = 7)
plot(model)
calib <- CalibrateRaster(model, strips)
}
