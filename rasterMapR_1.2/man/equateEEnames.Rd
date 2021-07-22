\name{equateEEnames}
\alias{equateEEnames}
\title{Equate EE names}
\description{This function equates names of bands in Earth-Explorer (EE) data
sets using conventions of a specific Landsat satellite.}
\usage{equateEEnames(dtls, sat.nm = "LE07", 
    change.dt = TRUE)}
\arguments{
  \item{dtls}{\code{list}. Set of two \code{RasterBrick} objects:
target and reference, sharing a common region such as
these produced by \code{\link{EEstackWithoutMeta}} or
\code{\link{RasterIntersection}}.}
  \item{sat.nm}{\code{character}. Regular expression
indicating a Landsat satellite used as reference
to match names of the bands, see
\code{Details}.}
  \item{change.dt}{\code{logical}. Print the renamed data. If
\code{FALSE} then only new names are printed.}
}
\details{Landsat satellite in
\code{sat.nm} is matched with columns in
\code{\link{LBands}}. }
\value{\code{list} of either characters or \code{RasterBrick} depending
on \code{change.dt}}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'aRn')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6), sat.nm = NULL)
eqstack <- equateEEnames(stack)
}
