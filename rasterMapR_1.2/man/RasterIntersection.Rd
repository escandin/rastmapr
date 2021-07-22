\name{RasterIntersection}
\alias{RasterIntersection}
\title{Raster Intersection}
\description{This function extracts intersections between two Earth-Explorer
Landsat data sets.}
\usage{RasterIntersection(rasterImp)}
\arguments{
  \item{rasterImp}{Either \code{list} of two rasters, or
\code{caracter} vector of paths to the data.}
}
\details{The function implements
\code{\link{intersect}} to derive a
list of two raster objects sharing
the same region.}
\value{\code{list} of rasters.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
## model <- PIFmodel(strips, pvalue = 3E-1, brackets = 7)
## plot(model)
## calib <- CalibrateRaster(model, strips)
}
