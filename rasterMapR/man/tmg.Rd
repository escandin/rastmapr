\name{tmg}
\alias{tmg}
\title{Thematic map generator}
\description{This function produces a thematic map based on either training
data or previously calibrated randomForest object.}
\usage{tmg(inrast, calibObject, 
    ntrees = 1000, classcolname = "class", 
    validObject, maptype = "class", 
    savefiles = TRUE, 
    outfile = "outfilename", 
    plot = TRUE, tmproduce = TRUE)}
\arguments{
  \item{inrast}{\code{list}. List of input rasters}
  \item{calibObject}{\code{logical}.}
  \item{ntrees}{\code{numeric}. "none" for no normalization,
"aRn" for std diff normalization or "cor" for
pearson correlation (pending)}
  \item{classcolname}{. Item in the list corresponding
to the raster to use as reference}
  \item{validObject}{\code{numeric}. Items to be normalized and
mosaiced. If mosaicitems=1, no mosaics are
built.}
  \item{maptype}{\code{character}}
  \item{savefiles}{.}
  \item{outfile}{\code{character}.}
  \item{plot}{\code{logical}.}
  \item{tmproduce}{\code{numeric}.}
}

\value{return(tm) # not sure what object to return since it depends on the inputs}

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
