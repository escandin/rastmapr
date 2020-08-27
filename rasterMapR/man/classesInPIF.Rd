\name{classesInPIF}
\alias{classesInPIF}
\title{Classes in PIF}
\description{This function is under testing and has not been implemented
yet. For computing \code{PIF} see \code{\link{PIFmodel}},
\code{\link{thresraster}}, \code{\link{nochg}}, and
\code{\link{calibrationParameters}}. The function calculates
classes into \code{brackets} of Pseudo-Invariant Features.}
\usage{classesInPIF(pifdt, brackets = 15)}
\arguments{
  \item{pifdt}{\code{data.frame}. Data set of PIF}
  \item{brackets}{\code{numeric}. Number of \code{brackets}.}
}
\details{Selected \code{PIF} are
scaled by computing the square of the
scaled differences between two rasters,
see \code{\link{scale}}. The scaled
\code{PIF} are then averaged using classes
of reflectance in the compared images
(\code{PIF} classes). The classes are
computed using \code{\link{quantile}}.}
\value{\code{data.frame}.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
classesInPIF(NULL)
}
