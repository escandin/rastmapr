\name{rclMatrix}
\alias{rclMatrix}
\title{reclassification matrix}
\description{A binary reclassification matrix is produced by specifying a
constant (threshold). The constant is used to split the interval
\code{c(-Inf, Inf)} into two intervals.}
\usage{rclMatrix(thr = 0, oneFirst = TRUE)}
\arguments{
  \item{thr}{\code{numeric}. Constant (threshold). Default \code{0}}
  \item{oneFirst}{\code{numeric}. Condition in the first row
becomes \code{1}. If \code{FALSE} then the second
column will do.}
}
\details{See \code{rcl} argument in
\code{\link{reclassify}}. This
function is implemented by
in-package \code{\link{PIFmodel}}
and \code{\link{RasterIntersection}}}
\value{\code{Matrix}.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
rclMatrix(4, FALSE)
}
