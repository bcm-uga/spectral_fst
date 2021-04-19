\name{lfmm2}
\alias{lfmm2}




\title{Fitting Latent Factor Mixed Models (Least squares algorithm)}
\description{Latent Factor Mixed Models (LFMMs) are factor regression models in which the response variable is a genotypic matrix, and the explanatory variables are environmental measures of ecological interest or trait values. The \code{\link{lfmm2}} function estimates latent factors based on an exact least-squares approach.
}

\usage{
lfmm2 (input, env, K, lambda)
}

\arguments{
\item{input}{A genotypic matrix
}
\item{env}{A matrix of environmental covariates
}
\item{K}{
    An integer corresponding to the number of latent factors.
}
\item{lambda}{
    A positive numeric value for a ridge regularization parameter. The default value is set to 1e-5.
}
}
\value{
\code{lfmm2} returns an object of class \code{lfmm2Class} that contains $K$ estimated latent factors \code{@U} and their loadings  \code{@V}.


}
\author{
    Olivier Francois
}


\references{
Caye K, Jumentier B, Lepeule J, Francois O. (2019). LFMM 2: fast and accurate inference of gene-environment associations in genome-wide studies. Molecular biology and evolution, 36(4), 852-860.
}
\examples{
### Example of analysis using lfmm2 ###


# Simulation with 10 target loci, with effect sizes ranging between -10 an 10
# n = 100 individuals and L = 1000 loci

X <- as.matrix(rnorm(100)) # causal environmental variable
B <- rep(0, 1000)
target <- sample(1:1000, 10) # target loci
B[target] <- runif(10, -10, +10) # effect sizes

# Creating hidden factors and loadings

U <- t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X))+ matrix(rnorm(300), ncol = 3)
V <- matrix(rnorm(3000), ncol = 3)

# Simulating a binarized matrix containing haploid genotypes
# Simulation performed with the generative LFMM

Y <- tcrossprod(as.matrix(X), B) + tcrossprod(U, V) + matrix(rnorm(100000, sd = .5), nrow = 100)
Y <- matrix(as.numeric(Y > 0), ncol = 1000)

######################################
# Fitting an LFMM with K = 3 factors #
######################################

mod2 <- lfmm2(input = Y, env = X, K = 3)

# Computing P-values and plotting their minus log10 values
# Target loci are highlighted

pv <- lfmm2.test(object = mod2, input = Y, env = X, linear = TRUE)
plot(-log10(pv$pvalues), col = "grey", cex = .4, pch = 19)
points(target, -log10(pv$pvalues[target]), col = "red")

}


\keyword{lfmm2}
\keyword{tutorial}
