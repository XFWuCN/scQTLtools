% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{set_result_info}
\alias{set_result_info}
\alias{set_result_info,eQTLObject-method}
\title{Set eQTL results in an eQTLObject.}
\usage{
set_result_info(x, value)

\S4method{set_result_info}{eQTLObject}(x, value)
}
\arguments{
\item{x}{An \code{eQTLObject}.}

\item{value}{A data frame where each row corresponds to a SNP–gene pair.}
}
\value{
An updated \code{eQTLObject}.

An updated \code{eQTLObject}.
}
\description{
Set eQTL results in an eQTLObject.

Method to set the result of identifying eQTLs from scRNA-seq data.
}
\examples{
  data(EQTL_obj)
  result <- data.frame(0, nrow = 3, ncol = 3)
  set_result_info(EQTL_obj, result)

  data(EQTL_obj)
  result <- matrix(0, nrow = 3, ncol = 3)
  set_result_info(EQTL_obj, result)

}
