% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSNPsLoc-callQTL.R
\name{createSNPsLoc}
\alias{createSNPsLoc}
\title{Create SNP location dataframe.}
\usage{
createSNPsLoc(snpList, snp_mart = NULL, snpDataset = "hsapiens_snp")
}
\arguments{
\item{snpList}{A list of SNPs IDs.}

\item{snp_mart}{An object of class Mart representing the BioMart SNP
database to connect to. If provided, this should be a Mart object obtained
by calling \code{useEnsembl()}, which allows specifying a mirror in case of
connection issues. If \code{NULL}, the function will create and use a Mart
object pointing to the Ensembl SNP BioMart, using the specified
\code{snpDataset} and a default mirror.}

\item{snpDataset}{A character string specifying the SNP dataset to use from
Ensembl. Default is \code{hsapiens_snp} for human SNPs.}
}
\value{
A data frame containing the SNP genomic locations.
}
\description{
Create SNP location dataframe.
}
\examples{
snpList <- c('rs546', 'rs549', 'rs568', 'rs665', 'rs672')
library(biomaRt)
snp_mart <- useEnsembl(biomart = "snps",
                        dataset = "hsapiens_snp",
                        mirror = 'asia')
snp_loc <- createSNPsLoc(snpList = snpList,
                         snp_mart = snp_mart)
}
