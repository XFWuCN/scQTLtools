% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{load_biclassify_info}
\alias{load_biclassify_info}
\alias{load_biclassify_info,eQTLObject-method}
\title{Access biclassification information from an eQTLObject.}
\usage{
load_biclassify_info(x)

\S4method{load_biclassify_info}{eQTLObject}(x)
}
\arguments{
\item{x}{An \code{eQTLObject}.}
}
\value{
A character or list containing biclassification information.

biclassify information of eQTLObject.
}
\description{
Access biclassification information from an eQTLObject.

Method to access eQTLObject biclassify information.
}
\examples{
  data(EQTL_obj)
  load_biclassify_info(EQTL_obj)
}
