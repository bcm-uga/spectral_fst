\name{spectralfst}
\alias{spectralfst}
\title{
Compute Wright's inbreeding coefficient using coefficient of determination
}
\description{
This function computes wright's Fst using coefficient of determination obtained from the regression of genotype on the population labels. This allows you to compute Fst on adjusted genotypic matrices. To do so, you need to specify a matrix Y of adjusting variables.
}
\usage{
spectralfst(genotype, population_labels, Y = NULL, adjustment_type = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
A matrix of genotypes with n rows and L columns, where n is number of individuals and L is the number of loci. Genotypes can be
#' haploid (0,1) or diploid (0,1,2). Missing values are not allowed.
}
  \item{population_labels}{
A vector of length n containing population labels for each sample in the genotype matrix. Population labels are integers
}
  \item{Y}{
A matrix of covariates with n rows and d columns, where d is the number of covariates. The values are used to compute an adjusted estimate of Fst. Missing values are not allowed. If \code{NULL}, then spectral estimates are unadjusted.
}
  \item{adjustment_type}{
A value 1 or 2. If adjustment_type = 1, then you will obtain Fst from matrix BX. If adjustment_type = 2, then you will obtain Fst from matrix W + \epsilon .
}

}

\value{
A vector of size L, with the corresponding value of Fst for each locus.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(spectralfst)

# We will use athaliana data
# It is composed of genotype, population labels, and bioclimatic
# variables that will allow us to test the adjustment parameter
data(athaliana)

# In this first part of the example, we will adjust on matrix W + E

fst_adjusted_w <- spectralfst(athaliana$genotype, athaliana$pop, athaliana$bio, 2)

# You now have the Fst at all locus

# In this second part, we will compute Fst on matrix BY where B are
# effect size and Y are adjusting variables

fst_adjusted_by <- spectralfst(athaliana$genotype, athaliana$pop, athaliana$bio, 1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
