\name{mean_wright_fst}
\alias{mean_wright_fst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute mean Fst using Wright's formula
}
\description{
This function aims at computing mean value of Fst over loci on an haploid genotype matrix
}
\usage{
mean_wright_fst(genotype, population_assignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
  A genotype matrix of haploid individuals with n rows (=individuals) and p columns (=loci).
}
  \item{population_assignment}{
  A vector of length n, each element is the population of the corresponding individual in genotype   matrix
}
}

\value{
The mean value of Fst
}


\examples{
# Loading datas from package
data('genotype_athaliana')
data('pop_athaliana')

wright_fst_athaliana <- mean_wright_fst(genotype_athaliana, pop_athaliana)

print("mean Fst compute with Wright's formula :")
wright_fst_athaliana
}


