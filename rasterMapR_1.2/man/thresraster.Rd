\name{thresraster}
\alias{thresraster}
\title{Threshold raster}
\description{This function produces a threshold raster by computing
standardized differences between two overlaying \code{Raster*}
objects: reference and \code{target}, according to a range of valid pixel
values.}
\usage{thresraster(ref, target, 
    minvalid = 0, maxvalid = Inf)}
\arguments{
  \item{ref}{\code{RasterBrick}. Reference raster}
  \item{target}{\code{RasterBrick}. Target raster}
  \item{minvalid}{\code{numeric}. Minimum valid value.}
  \item{maxvalid}{\code{numeric}. Maximum valid value.}
}
\details{This function is implemented by
\code{\link{PIFmodel}}. The standardized
differences in the threshold raster will be
contrasted with the chi-square statistic to
identify significant differences between
the overlaying rasters, see
\code{\link{nochg}}.}
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
