% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iTMEcell.R
\name{iTMEcell}
\alias{iTMEcell}
\title{Identification of abnormal tumor microenvironment (TME) cells}
\usage{
iTMEcell(ExpData,clinical,nperm=1000)
}
\arguments{
\item{ExpData}{A gene expression profile of interest (rows are genes, columns are samples).}

\item{clinical}{A dataframe with three columns which are "sample" (sample id),"status" (survival status of samples, "0" represents live and "1" represents dead) and "time" (survival time of samples).}

\item{nperm}{Number of random permutations (default: 1000).}
}
\value{
A dataframe with seven columns those are cell names, marker source, marker size, marker genes, centrality (eigenvector centrality), P-value and FDR.
}
\description{
The function "iTMEcell" is used to calculate the eigenvector centrality of TME cells and identify abnormal TME cells.
}
\examples{
library(igraph)
#Obtain input data
GEP<-GetExampleSet('GEP')
clinicaldata<-GetExampleSet('clinicaldata')
#Run the function
iTMEcellresult<-iTMEcell(ExpData=GEP,clinical=clinicaldata,nperm=1000)
}
