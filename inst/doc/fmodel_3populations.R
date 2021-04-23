## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(spectralfst)

## ------------------------------------------------------------------------
## parameters F model
F1 = 0.3
F2 = 0.2
F3 = 0.1

L = 20000
n = 150
c = c(1/3,1/3,1/3)


# genotype matrices
X <- matrix(NA, nrow = n, ncol = L)
n1 <- round(n*c[1])
n2 <- round(n*c[2])
n3 <- n - n1 - n2 


for (l in 1:L){
  p <- runif(1)
  p1 <- rbeta(1, p*(1-F1)/F1, (1-p)*(1-F1)/F1)
  p2 <- rbeta(1, p*(1-F2)/F2, (1-p)*(1-F2)/F2)
  p3 <- rbeta(1, p*(1-F3)/F3, (1-p)*(1-F3)/F3)  
  X[1:n1,l] <- rbinom(n1, size = 1, prob = p1)
  X[(n1+1):(n1+n2),l] <- rbinom(n2, size = 1, prob = p2)
  X[(n1+n2+1):n,l] <- rbinom(n3, size = 1, prob = p3)  
}

## SNPs filtering
boo.unique <- apply(X, 2, FUN = function(x) length(unique(x))) == 1
X <- X[,!boo.unique]
L <- ncol(X)

# Population_assignment vector creation

pop <- c(rep(1,n1), rep(2,n2), rep(3,n3))

## ------------------------------------------------------------------------
dim(X)

## ------------------------------------------------------------------------
fst_wright <- mean_wright_fst(X, pop)
print("Wright's formula gives a mean fst of : ")
fst_wright

## ------------------------------------------------------------------------
spectral_fst <- compute_partition(X, pop)
print("Fst compute with between population matrix Zst")
spectral_fst$Fst

## ------------------------------------------------------------------------
print("Fst approximation using Z matrix")
spectral_fst$Fst_approximation

## ------------------------------------------------------------------------
print("The leading eigenvalue of Zs matrix is :")
spectral_fst$leadingeigenZs

## ------------------------------------------------------------------------
print("The RMT prediction of the leading eigenvalue of Zs is : ")
spectral_fst$RMTprediction

## ---- fig.show='hold'----------------------------------------------------
plot(spectral_fst$eigenZ[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZ, col = pop, main="PC Plot")

## ---- fig.show='hold'----------------------------------------------------
plot(spectral_fst$eigenZst[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZst, col = pop, main="PC Plot")

## ---- fig.show='hold'----------------------------------------------------
plot(spectral_fst$eigenZs[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZs, col = pop, main="PC Plot")

