\name{LBands}
\docType{data}
\alias{LBands}
\title{Landsat Band designations
}
\description{Band designations for the Landsat satellites
}
\usage{data(LBands)}
\format{
    A data frame with the following seven variables.
    \describe{
      \item{\code{patt}}{Regular expression matching
	initial four letters of names in Earth Explorer data}
      \item{\code{band}}{Band number}
      \item{\code{name}}{Band name}
      \item{\code{wave_mmm_from}}{Mimimum wavelength (micrometers)}
      \item{\code{wave_mmm_to}}{Maximum wavelenth (micrometers)}
      \item{\code{res_m_initial}}{Original resolution (meters)}
      \item{\code{res_m_final}}{Final resolution (meters)}
      }
}
\details{The band designations are compared with regular expressions in
  \code{patt} to match names in Earth Explorer data.
}
\source{
\href{https://www.usgs.gov/faqs/what-are-band-designations-landsat-satellites?qt-news_science_products=0#qt-news_science_products}{https://www.usgs.gov/faqs/what-are-band-designations-landsat-satellites?qt-news_science_products=0#qt-news_science_products}
}
\examples{
str(LBands)
}
\keyword{data sets}
