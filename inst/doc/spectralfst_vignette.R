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

## ------------------------------------------------------------------------
data('genotype_athaliana')
data('env_var')
data('pop_athaliana')

## ------------------------------------------------------------------------
spectral_athaliana <- compute_partition(genotype_athaliana, pop_athaliana)
wright_fst_athaliana <- mean_wright_fst(genotype_athaliana, pop_athaliana)

print("mean Fst compute with Wright's formula :")
wright_fst_athaliana
print("Fst compute with between population matrix Zst")
spectral_athaliana$Fst
print("Fst approximation using Z matrix")
spectral_athaliana$Fst_approximation

## ------------------------------------------------------------------------
spectral_adjusted_athaliana <- compute_partition(genotype_athaliana, pop_athaliana, make_adjustment = T, adjusting_variables = env_var)

print("Fst compute with between population matrix Zst")
spectral_adjusted_athaliana$Fst
print("Fst approximation using Z matrix")
spectral_adjusted_athaliana$Fst_approximation

## ------------------------------------------------------------------------
m <- cbind(spectral_athaliana$eigenZ[1:6],spectral_adjusted_athaliana$eigenZ[1:6])
plot(seq(1,12)/2,as.numeric(t(m)),
     col = rep(c("darkblue","orange"), length = 12), 
     type = "h", lwd = 10, 
     ylab = "Proportion of variance",
     xlab = "PC axis",
    # main = "PCA",
     cex.axis = .9,
     las = 1)

legend(x = 3.5, y = 0.081, 
       legend = c("All SNPs", "Adjusted data"),
       col = c("darkblue","orange"),
       lty = c(1,1),
       pch = 19,
       cex = 1)

## ------------------------------------------------------------------------
m <- cbind(c(spectral_athaliana$eigenZst[1], spectral_athaliana$eigenZs[1:5]),
           c(spectral_adjusted_athaliana$eigenZst[1], spectral_adjusted_athaliana$eigenZs[1:5]))
plot(seq(1,12)/2,as.numeric(t(m)), 
     col = rep(c("lightblue","orange"), length = 12), 
     type = "h", lwd = 10, 
     ylab = "Proportion of variance",
     xlab = "PC axis",
    # main = "Structure Components",
     cex.axis = .9,
     las = 1)

legend(x = 3.5, y = 0.077, 
       legend = c("All SNPs", "Adjusted data"),
       col = c("lightblue","orange"),
       lty = c(1,1),
       pch = 19,
       cex = 1)

