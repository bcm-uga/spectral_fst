\name{haploidisation}
\alias{haploidisation}
\title{
Convert diploid genotype matrix into haploid.
}
\description{
This function takes in input a diploid genotype matrix of n lines (=individuals) and L columns (= loci) and returns a 2n * L equivalent haploid genotype matrix
}
\usage{
haploidisation(genotype, population_assignment)
}
\arguments{
  \item{genotype}{
  A genotype matrix of diploid individuals with n rows (=individuals) and L columns (=loci).
}
  \item{population_assignment}{
    A vector of length n, each element is the population of the corresponding individual in genotype matrix
}
}

\value{
\item{haploid_genotype}{A genotype matrix of 2*n rows and L columns. The corresponding haploid matrix}
\item{haploid_population_assignment}{A vector of size 2*n. The new vector of population assignment corresponding to the new haploid matrix}
}



\examples{
# Creation of a diploid matrix

diploid_matrix <- rbind(c(0,1,2),c(2,2,2))
diploid_pop <-  c(1,2)

hap_process <- haploidisation(diploid_matrix, diploid_pop)

# New haploid matrix and population assignment vector

haploid_matrix <- hap_process$haploid_matrix
haploid_pop <- hap_process$haploid_population_assignment
}

