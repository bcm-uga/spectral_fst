
#' @title fmodel_simulation
#
#' @description This function aims at simulating a fmodel
#
#' @param n number of individuals
#' @param L number of loci
#' @param c_vector A vector of length K (with K = number of populations). Each entry represents the proportion of a given population among the total population
#' @param F_vector A vector of length K (with K = number of populations). Each entry represent the drift coefficient of a given population
#'
#' @return fmodel_simulation
#' genotype : A matrix of n rows and L columns
#' pop : A vector of length n corresponding to population for each individual


fmodel_simulation <- function(n = 150, L = 20000, c_vector=c(1/3,1/3,1/3), F_vector=c(0.3,0.2,0.1)){
  # genotype matrices
  X <- matrix(NA, nrow = n, ncol = L)
  K <- length(c_vector)

  if (K != length(F_vector)) stop("c_vector and F_vector must have same length")
  if (sum(c_vector) != 1) stop("c_vector must sum to 1")

  n_vector <- c()
  for (i in seq(1,K-1)){
    n_vector <- c(n_vector, round(c_vector[i] * n))
  }

  n_vector <- c(n_vector, n - sum(n_vector))

  for (l in 1:L){
    p <- runif(1)
    for (i in seq(1,K)){
      if (i == 1) {
        current_index <- 1:n_vector[i]
      }else{
        current_index <- (sum(n_vector[1:(i-1)]) + 1):sum(n_vector[1:i])
      }
      Fi <- F_vector[i]
      p_pop <- rbeta(1, p*(1-Fi)/Fi, (1-p)*(1-Fi)/Fi)
      X[current_index,l] <- rbinom(n_vector[i], size=1, prob=p_pop)
    }
  }

  ## SNPs filtering
  boo.unique <- apply(X, 2, FUN = function(x) length(unique(x))) == 1
  X <- X[,!boo.unique]
  L <- ncol(X)

  # Population_assignment vector creation

  pop <- c()
  for (i in seq(1,K)){
    pop <- c(pop, rep(i, n_vector[i]))
  }

  return(list(genotype=X, pop=pop))
}
