\name{athaliana}
\alias{athaliana}
\docType{data}
\title{
DNA samples from Arabidopsis Thaliana
}
\description{
A list object containing 241 DNA samples from Arabidopsis Thaliana individuals associated with their population and bioclim variables
}
\usage{data("athaliana")}
\format{List with 3 objects
\describe{
  \item{genotype}{A matrix of 241 rows and 10 000 columns. It corresponds to 241 swedish plant accessions from The 1,001 Genomes database for Arabidopsis thaliana. The 10,000 filtered SNPs were randomly selected from a large set of variants. The matrix contains no missing genotypes.}
  \item{pop}{A vector of length 241. The i-th element corresponds to the population of the i-th individuals. The  individuals  were  clustered  in  groups  based  on  analysis  of  population structure accounting for geographic proximity.}
  \item{bio}{A matrix of 241 rows and 18 columns. Global climate and weather data corresponding to individual geographic coordinates were downloaded from the World-Clim  database  (https://worldclim.org).   Eighteen  bioclimatic  variables,  derived from the monthly temperature and rainfall values,  were considered as representing the current environmental matrix.}
}
}

\keyword{datasets}
