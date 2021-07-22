\name{calibrationParameters}
\alias{calibrationParameters}
\title{Calibration parameters}
\description{This function generates calibration parameters comparing bands in
two overlaying \code{Raster*} objects: reference and \code{target}, while
using a change/no-change mask, see \code{\link{nochg}}.}
\usage{calibrationParameters(ref, 
    target, threshMask, 
    nbrackets = 30)}
\arguments{
  \item{ref}{\code{RasterBrick}. Reference raster.}
  \item{target}{\code{RasterBrick}. Target raster.}
  \item{threshMask}{\code{Raster*}. Change/no-change mask such as that
produced by \code{\link{nochg}}.}
  \item{nbrackets}{\code{numeric}. Number of brackets, see Section
of \code{Details}.}
}
\details{This function is implemented by
\code{\link{PIFmodel}}. The function
develops four procedures. First, values in
the \code{Raster*} objects are extracted
using the binary mask. Second, the selected
values are reduced to classes according to
a number of brackets in
\code{nbrackets}. Third, average classes of
values are modeled generating intercepts
and slopes for each of the bands. Fourth,
non-significant models are subtracted.}
\value{\code{RasterLayer}. }

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1:4))
strips <- RasterIntersection(stack)
## thrs <- thresraster(strips[[2L]], strips[[1L]])
## noch <- nochg(thrs, degfree = nlayers(strips[[2L]]) - 1, pvalue = 4E-1)
## calp <- calibrationParameters(strips[[2L]], strips[[1L]], noch, nbrackets = 8)
model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
calib <- CalibrateRaster(model, stack)
## merged <- merge(calib, stack[[2L]][[names(calib)]])
## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
}
