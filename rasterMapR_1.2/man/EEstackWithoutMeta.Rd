\name{EEstackWithoutMeta}
\alias{EEstackWithoutMeta}
\title{Stack EE avoiding metadata}
\description{This function stacks Landsat data shared by Earth Explorer while
avoiding embedded metadata.}
\usage{EEstackWithoutMeta(tarPath, 
    bands = NULL, sat.nm = "LE07", 
    drop.unmatched = TRUE, 
    meta.wd = FALSE)}
\arguments{
  \item{tarPath}{\code{List}. File paths or names of the compressed
(\code{.tar}) data containing multi- or hyper-spectral
radiometric data.}
  \item{bands}{\code{numeric} or \code{NULL}. Bands to be
processed after matching the band names (see
\code{sat.nm}). Default \code{NULL} preserves the
\code{bands}.}
  \item{sat.nm}{\code{character} or \code{NULL}. Regular
expression indicating a Landsat satellite used
to match names of the \code{bands}. Default
\code{'LE07'}. If \code{NULL} then the band
names are not matched.}
  \item{drop.unmatched}{\code{logical}. Drop the unmatched
\code{bands}.}
  \item{meta.wd}{\code{logical}. Copy metadata in the working
directory.}
}
\details{Compressed files of Multi-
or hyper-spectral layers
(\code{.tar}) are recursively
decompressed, stacked, and
renamed. The data is processed using
the names of the layers while
ignoring the metadata. However, a
portion of the names in the \code{bands}
must match each other. Names of the
\code{bands} should also contain the band
numbers. Only \code{bands} with similar
extents are stacked while ignoring
data with different
extents. \code{\link{equateEEnames}}
is implemented to match data from
different satellites}
\value{RasterBrick.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
}
