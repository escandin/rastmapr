\name{EEstackMeta}
\alias{EEstackMeta}
\title{Stack EE data}
\description{This function is under development. The function implements
\code{\link{stackMeta}} to stack Landsat hyper-spectal layers
using embedded metadata. Use \code{\link{EEstackWithoutMeta}}
instead.}
\usage{EEstackMeta(tarPath, 
    bands = NULL, ...)}
\arguments{
  \item{tarPath}{\code{character} or \code{List}. Name(s) of the
commpresed data set(s). These can include the file
paths.}
  \item{bands}{\code{numeric}. Bands to be used.}
  \item{\dots}{Other arguments in \code{\link{stackMeta}}}
}

\value{RasterBrick.}

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
EEstackMeta(NULL)
}
