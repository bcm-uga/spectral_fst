\name{fmodel}
\alias{fmodel}
\docType{data}
\title{
DNA samples from fmodel simulation
}
\description{
A list object containing a matrix genotype of size n = 150 and L = 18967 loci and a vector pop of length n = 150 referring to population for each individual.
}
\usage{data("fmodel")}

\format{List with 2 objects
\describe{
  \item{genotype}{A matrix of 150 rows and 18 967 columns. It corresponds to fmodel simulations of 3 populations of 50 individuals with drift coefficient F1=0.3, F2=0.2, F1=0.1}
  \item{pop}{A vector of length 150. The i-th element corresponds to the population of the i-th individuals. }
}
}

\keyword{datasets}
