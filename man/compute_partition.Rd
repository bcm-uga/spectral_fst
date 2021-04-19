\name{compute_partition}
\alias{compute_partition}
\title{Computing partition matrices and related quantities}
\usage{
compute_partition ( genotype, population_assignment)
}
\description{
This function aims at computing matrices that we define in the paper as within population matrix Zs and between population matrix Zst
}

\arguments{
\item{genotype}{
A genotype matrix of haploid or transformed haploid individuals with n rows and p columns.
}
\item{population_assignment}{
A vector of length n, each element is the population of the corresponding individual in genotype matrix
}
}

\value{
\item{Zs}{The within population matrix Zs}
\item{Zst}{The between population matrix Zst}
\item{Fst}{Squared norm of Zst}
\item{Fst_estimate}{the sum of the first (nb_pop − 1) eigenvalues of scaled PCA for genotype matrix}
\item{leadingeigenZs}{The leading eigenvalue of Zs}
\item{MPestimate}{the estimation of Marchenko Pastur}
\item{eigenZs}{eigenvalues of scaled Zs}
\item{eigenZst}{eigenvalues of scaled Zst}
\item{eigenZ}{eigenvalues of scaled Z}
}

\examples{
See tutorials directory
}
