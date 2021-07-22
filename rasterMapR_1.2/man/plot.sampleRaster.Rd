\name{plot.sampleRaster}
\alias{plot.sampleRaster}
\title{sampleRaster plots}
\description{This function plots \code{\link{sampleRaster}} objects.}
\usage{\method{plot}{sampleRaster}(x, 
    y = 1, ...)}
\arguments{
  \item{x}{Raster* object.}
  \item{y}{If \code{x} is a RasterStack or RasterBrick: integer,
character (layer name(s)), see \code{plot} method in the
\code{raster} package.}
  \item{\dots}{Additional arguments in \code{\link{plot}}.}
}

\value{\code{plot}. Raster with sample units.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
set.seed(2)
sr <- sampleRaster(stack[[1L]], side = 200, n = 2)
## plot(sr)
}
