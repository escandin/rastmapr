\name{CalibrateRaster2}
\alias{CalibrateRaster2}
\title{Calibrate raster}
\description{This function implements parameters from \code{\link{PIFmodel}} to
calibrate multi- or hiper-spectral layers.}
\usage{CalibrateRaster2(pifs, 
    mlayer, round. = 0)}
\arguments{
  \item{pifs}{\code{List}. Set of linear models such as these produced
by \code{\link{PIFmodel}}.}
  \item{mlayer}{\code{RasterLayer}. Mult-layer object such as that
returned by \code{\link{EEstackWithoutMeta}}. If missing
then the set used to compute the \code{pifs} is
recycled.}
  \item{round.}{\code{numeric} or \code{NULL}. Integer indicating
the number of decimal places. If \code{NULL} then this
argument is ignored.}
}



\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
## \\donttest{
##     tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
##                   'LT050060661988072201T1-SC20181031160603.tar.gz')
##     tarPaths <- system.file(tarFiles, package = 'aRn')
##     stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
##     ## model <- PIFmodel(stack, pvalue = 3E-1, brackets = 7)
##     ## or
##     mlayer <- RasterIntersection(stack)
##     model <- PIFmodel(mlayer, pvalue = 3E-1, brackets = 7)
##     ## then
##     calib <- CalibrateRaster(model, stack)
## }
 
}
