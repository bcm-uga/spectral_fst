\name{compute_partition}
\alias{compute_partition}
\title{Spectral estimates of inbreeding coefficients from multilocus genotypes}
\usage{

compute_partition ( genotype, population_assignment, Y=NULL )
}
\description{
This function computes a within-population (residual) matrix, Zs, and a between-population matrix, Zst, for a genotype matrix, Z, and for K predefined populations. It allows the genotype matrix to  be adjusted for specific covariates, Y, representing technical effects, such as coverage, or environmental variables, such as bioclimatic data. The function returns a spectral estimate of the average value of Fst across loci, an approximation of Fst by PCA eigenvalues, and the leading eigenvalue of the residual matrix.
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

# Loading datas from package
data('athaliana')

# We compute Fst after removing adaptive genetic variation associated with environmental variables

spectral_adjusted_athaliana <- compute_partition(athaliana$genotype, athaliana$pop, athaliana$bio)

print("Fst compute with between population matrix Zst")
spectral_adjusted_athaliana$Fst
print("Fst approximation using Z matrix")
spectral_adjusted_athaliana$Fst_approximation
}

\references{
François, O., Gain, C., (2021). A Spectral Theory for Wright’s Inbreeding Coefficients and Related Quantities. Under Review.
}
