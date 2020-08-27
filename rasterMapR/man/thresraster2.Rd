\name{thresraster2}
\alias{thresraster2}
\title{Threshold raster}
\description{This function produces a threshold raster by computing
standardized differences between two overlaying \code{Raster*}
objects: reference and target, according to a range of valid pixel
values.}
\usage{thresraster2(refstack, 
    tarstack, cca = FALSE, 
    propsamp = 1)}
\arguments{
  \item{refstack}{\code{RasterBrick}. Reference raster stack}
  \item{tarstack}{\code{RasterBrick}. Target raster stack}
  \item{cca}{\code{logical}. Whether to perform a canonical
transformation to input rasters and then use it to
produce SSSD.}
  \item{propsamp}{\code{numeric}. proportion of pixels to sample for
fitting the \code{cca}.}
}

\value{\code{RasterLayer}. }

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
thrs <- thresraster(strips[[2L]], strips[[1L]])
}
