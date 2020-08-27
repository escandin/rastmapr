\name{sampleRaster}
\alias{sampleRaster}
\title{Sample Raster}
\description{This function processes Raster* objects to sample \code{n}
square polygons. Sampling without replacement is
implemented. Polygons having \link{NA} pixels are excluded.}
\usage{sampleRaster(rl, side = sqrt(10) * 
    1000, n = 30)}
\arguments{
  \item{rl}{\code{Raster*}. Raster Object with \code{UTM} Coordinate
Reference System.}
  \item{side}{\code{numeric}. Side of the
square. Default produces squares of 10
km2.}
  \item{n}{\code{numeric}. Number of samples.}
}

\value{\code{list}. Set of \code{n} fixed-area rasters and the polygon geometry.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
set.seed(2)
sr <- sampleRaster(stack[[1L]], side = 200, n = 2)
}
