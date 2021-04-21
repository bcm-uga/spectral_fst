\name{compute_partition}
\alias{compute_partition}
\title{Computing partition matrices and related quantities from spectral theory}
\usage{
## No adjustment
compute_partition ( genotype, population_assignment)

## With adjustment
compute_partition ( genotype, population_assignment, make_adjustment = TRUE, adjusting_variables)
}
\description{
This function aims at computing matrices that we define in the paper as within population matrix Zs and between population matrix Zst. From these matrices, it computes related quantities such as the mean Fst over loci, the approximation of Fst using Z matrix, random matrix theory (RMT) prediction.. Here, n corresponds to the number of individuals and L to the number of locus.
}

\arguments{
\item{genotype}{
A genotype matrix of n rows and L columns. The matrix can be composed of haploid, transformed haploid or diploid individuals.
}
\item{population_assignment}{
A vector of length n. Each element is the population of the corresponding individual in genotype matrix
}
\item{make_adjustment}{
A logical value that specifies whether you want to compute some correction on your genotype or not.
}

\item{adjusting_variables}{
An optional matrix of n rows and d columns where d is the number of adjusting variables. This matrix has to be specify when make_adjustment is TRUE.
}
}

\value{
\item{Zs}{A matrix of n rows and L columns in the case of haploid genotype, 2*n rows and L columns in the case of diploid genotype. It corresponds to the within population matrix Zs}
\item{Zst}{A matrix of n rows and L columns in the case of haploid genotype, 2*n rows and L columns in the case of diploid genotype. It corresponds to the between population matrix Zst}
\item{Fst}{A real value. It is computed using the squared norm of Zst and it is equal to the average value of Wright's Fst over all loci included in the genotype matrix}
\item{Fst_approximation}{A real value. The sum of the first (nb_pop − 1) eigenvalues of scaled PCA for genotype matrix Z}
\item{leadingeigenZs}{A real value. The leading eigenvalue of Zs}
\item{RMTprediction}{A real value. The prediction of the leading eigenvalue of the residual matrix using RMT}
\item{eigenZs}{A matrix of size n x 2 for haploid genotype, 2n x 2 for diploid. The projection of Zs on 2 first axis of PCA}
\item{eigenZst}{A matrix of size n x 2 for haploid genotype, 2n x 2 for diploid. The projection of Zst on 2 first axis of PCA}
\item{eigenZ}{A matrix of size n x 2 for haploid genotype, 2n x 2 for diploid. The projection of Z on 2 first axis of PCA}
\item{pcZs}{A prcomp object. It corresponds to the PCA of scaled Zs}
\item{pcZst}{A prcomp object. It corresponds to the PCA of scaled Zst}
\item{pcZ}{A prcomp object. It corresponds to the PCA of scaled Zst}
}

\examples{
See tutorials directory
}
