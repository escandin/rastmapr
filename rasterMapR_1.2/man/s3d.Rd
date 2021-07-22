\name{s3d}
\alias{s3d}
\title{Iterated Sum of the Squared Standaradized}

\usage{s3d(strips, cca = FALSE, 
    distype = "gamma", 
    distsamp = 0.01, 
    pval.pif = 0.01, 
    pval.chg = 0.99, 
    minPIF = 5, thres = 0.01, 
    maxiter = 20, prefix = "", 
    norm.ext = NULL, 
    fitline = TRUE, writemasks = TRUE)}
\arguments{
  \item{strips}{\code{list}. Target and reference raster stacks
respectively.}
  \item{cca}{\code{logical}. Same argument as for
\code{PIFmodel2}.}
  \item{distype}{\code{character}. Same as noch2}
  \item{distsamp}{\code{numeric}. Same as noch2}
  \item{pval.pif}{\code{numeric}.}
  \item{pval.chg}{\code{numeric}.}
  \item{minPIF}{\code{numeric}. Minimum number of pixels that
should be selected to produce the normalization. The
function will abort if the number of pixels selected
is lower.}
  \item{thres}{\code{numeric}.}
  \item{maxiter}{\code{numeric}.}
  \item{prefix}{\code{character}.}
  \item{norm.ext}{ or \code{NULL}.}
  \item{fitline}{\code{logical}. Label pseudoinvariant targets
and use them to fit a regression.}
  \item{writemasks}{\code{logical}. Write in the hard drive the no
change masks produced in each iteration}
}
\details{This function ...}
\value{... }

\author{Victor Gutierrez-Velez <victorhugo@temple.edu>, Wilson Lara <wilarhen@temple.edu>, Jeronimo Rodriguez <jeronimo.rodriguez@temple.edu >}




\examples{
tarFiles <- c('LT050070651987081201T1-SC20181031175314.tar.gz',
              'LT050060661988072201T1-SC20181031160603.tar.gz')
tarPaths <- system.file(tarFiles, package = 'rasterMapR')
stack <- EEstackWithoutMeta(tarPaths, c(1,3:6))
strips <- RasterIntersection(stack)
sssd <- s3d(strips, distsamp = 0.5,pval.pif = 0.4, pval.chg = 0.5, maxiter = 1)
}
