\name{stack2df}
\alias{stack2df}
\title{Extracts into data frame }
\description{This function extracts into a data frame the pixel values for
all bands from different classes defined in a spatial dataframe}
\usage{stack2df(inrast, invec, 
    classcolname = "class")}
\arguments{
  \item{inrast}{.}
  \item{invec}{.}
  \item{classcolname}{.}
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
## model <- PIFmodel(strips, pvalue = 4E-1, brackets = 8)
## calib <- CalibrateRaster(model, stack)
## merged <- merge(calib, stack[[2L]][[names(calib)]])
## plotRGB(merged, r = 3, g = 2, b = 1, stretch = 'lin')
}
