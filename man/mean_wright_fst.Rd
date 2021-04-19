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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The mean value of Fst
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

