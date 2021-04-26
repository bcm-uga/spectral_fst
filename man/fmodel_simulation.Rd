\name{fmodel_simulation}
\alias{fmodel_simulation}
\title{
F model simulation
}
\description{
This functions aims at simulating haploid genotype using the F model. It allows you to simulate as many populations as you want with their associated drift coefficients.
}
\usage{
fmodel_simulation(n = 150, L = 20000, c_vector = c(1/3, 1/3, 1/3), F_vector = c(0.3, 0.2, 0.1))
}
\arguments{
  \item{n}{
  An integer. Specify the number of individuals you want to simulate
}
  \item{L}{
  An integer. Specify the number of loci you want to simulate. You can possibly obtain a matrix with a different number of loci as column with zero variance will be remove
}
  \item{c_vector}{
  A vector of length K with K the number of population you want to simulate. This vector represents the proportion of individuals you want in each population. Hence the entries of this vector must add up to one.
}
  \item{F_vector}{
  A vector of length K with K the number of population you want to simulate. Each entry represents the drift coefficient for each population. Hence, each entry must be comprise between zero and one
}
}

\value{
\item{genotype}{A matrix of n rows and Lmod columns. Lmod can be sligthly different than L as column with zero variance has been removed.}
\item{pop}{A vector of length n. Each entry corresponds to the population of an individual}
}


\examples{

## We simulate a  fmodel with defaults parameters

fmodel <- fmodel_simulation(n = 150, L = 20000, c_vector = c(1/3, 1/3, 1/3), F_vector = c(0.3, 0.2, 0.1))

X <- fmodel$genotype
pop <- fmodel$pop

# We run spectral fst on the simulated datas
spectral_fst <- compute_partition(X, pop)
print("Fst compute with between population matrix Zst")
spectral_fst$Fst
}

