\name{commonExtent}
\alias{commonExtent}
\title{Common extent between raster objects}
\description{This function calculates the common extent between all the raster
objects within a list. The objects should all overlap at least
partially and have the same projection, spatial resolution and be
aligned. This can be written in a function. It can calculate
either the largest extent that includes all the listed elements
(union) other minimum common extent (intersection).}
\usage{commonExtent(inlist, 
    method = "intersection")}
\arguments{
  \item{inlist}{.}
  \item{method}{.}
}
\details{This function ...}
\value{... }

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
