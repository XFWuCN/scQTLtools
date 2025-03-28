% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{get_result_info}
\alias{get_result_info}
\title{Generic to access the result of identifying eQTLs from scRNA-seq data}
\usage{
get_result_info(x)
}
\arguments{
\item{x}{A eQTLObject object.}
}
\value{
A dataframe.
}
\description{
Generic to access the result of identifying eQTLs from scRNA-seq data
}
\examples{
  data(testEQTL)
  get_result_info(testEQTL)

}
