## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(spectralfst)

## ------------------------------------------------------------------------
data('fmodel')
X <- fmodel$genotype
pop <- fmodel$pop

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
data('athaliana')

## ------------------------------------------------------------------------
spectral_athaliana <- compute_partition(athaliana$genotype, athaliana$pop)
wright_fst_athaliana <- mean_wright_fst(athaliana$genotype, athaliana$pop)

print("mean Fst compute with Wright's formula :")
wright_fst_athaliana
print("Fst compute with between population matrix Zst")
spectral_athaliana$Fst
print("Fst approximation using Z matrix")
spectral_athaliana$Fst_approximation

## ------------------------------------------------------------------------
spectral_adjusted_athaliana <- compute_partition(athaliana$genotype, athaliana$pop, Y =athaliana$bio)

print("Fst compute with between population matrix Zst")
spectral_adjusted_athaliana$Fst
print("Fst approximation using Z matrix")
spectral_adjusted_athaliana$Fst_approximation

## ---- fig.width = 7, fig.height = 7--------------------------------------
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

## ---- fig.width = 7, fig.height = 7--------------------------------------
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

