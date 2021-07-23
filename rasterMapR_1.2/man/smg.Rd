\name{smg}
\alias{smg}
\title{Seamless mosaic generator}
\description{This function produces a cloud masked and radiometrically normalized }
\usage{smg(inlist, cca = "FALSE", 
    method = "none", 
    refitem = NA, mosaicitems = c(1, 
        3), QAbandname = "pixel_qa", 
    cloudbuff = NA, savelcloudmsk = TRUE, 
    indem = NA, normbands = seq(1, 
        7), pval.pif = 0.01, 
    pval.chg = 0.99, 
    sensor = "OLI", verbose = TRUE, 
    norm.ext = NULL)}
\arguments{
  \item{inlist}{\code{list}. List of input rasters}
  \item{cca}{\code{logical}.}
  \item{method}{\code{character}. "none" for no normalization,
"s3d" for std diff normalization or "cor" for
pearson correlation (pending)}
  \item{refitem}{. Item in the list corresponding to the
raster to use as reference}
  \item{mosaicitems}{\code{numeric}. Items to be normalized and
mosaiced. If mosaicitems=1, no mosaics are
built.}
  \item{QAbandname}{\code{character}}
  \item{cloudbuff}{.}
  \item{savelcloudmsk}{\code{logical}.}
  \item{indem}{.}
  \item{normbands}{\code{numeric}.}
  \item{pval.pif}{\code{numeric}.}
  \item{pval.chg}{\code{numeric}.}
  \item{sensor}{\code{character}}
  \item{verbose}{\code{logical}.}
  \item{norm.ext}{.}
}

\value{... }

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'rasterMapR')
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
