## ----setup, include = FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------------------------------------------------------
library(spectralfst)

## ---- include=FALSE-----------------------------------------------------------------------------------------------------------
obj <-  fmodel_simulation(n=150, L=20000, c_vector=c(1/3,1/3,1/3), F_vector=c(0.3,0.2,0.1))

## -----------------------------------------------------------------------------------------------------------------------------
data('fmodel')
X <- fmodel$genotype
pop <- fmodel$pop
dim(X)

## -----------------------------------------------------------------------------------------------------------------------------
# Mean fst from Wright's formula
fst_wright <- mean_wright_fst(X, pop)
fst_wright

## -----------------------------------------------------------------------------------------------------------------------------
#Compute Fst from the between-population matrix Zst
spectral_fst <- compute_partition(X, pop)
spectral_fst$Fst

## -----------------------------------------------------------------------------------------------------------------------------
#Leading eigenvalue of Zs/sqrt(n)
spectral_fst$leadingeigenZs

## -----------------------------------------------------------------------------------------------------------------------------
#RMT prediction of the leading eigenvalue of Zs/sqrt(n)
round(spectral_fst$RMTprediction, digits = 4)

## -----------------------------------------------------------------------------------------------------------------------------
#Smallest non-null eigenvalue of Zst/sqrt(n) (= second one)
spectral_fst$eigenZst[2]

## -----------------------------------------------------------------------------------------------------------------------------
# Fst approximation using PCA
spectral_fst$Fst_approximation

## ---- fig.show='hold'---------------------------------------------------------------------------------------------------------
plot(spectral_fst$eigenZ[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZ, col = pop, main="PC Plot")

## ---- fig.show='hold'---------------------------------------------------------------------------------------------------------
plot(spectral_fst$eigenZst[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
spectral_fst$pcZst[,2] <- - spectral_fst$pcZst[,2] #PC signs are arbitrary
plot(spectral_fst$pcZst, col = pop, main="PC Plot")

## ---- fig.show='hold'---------------------------------------------------------------------------------------------------------
plot(spectral_fst$eigenZs[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZs, col = pop, main="PC Plot")

## -----------------------------------------------------------------------------------------------------------------------------
pop <- c(rep(1,100), rep(2,50))

## -----------------------------------------------------------------------------------------------------------------------------
#Wright's Fst
fst_wright <- mean_wright_fst(X, pop)
fst_wright

## -----------------------------------------------------------------------------------------------------------------------------
#Squared norm of the between population matrix Zst
spectral_fst <- compute_partition(X, pop)
spectral_fst$Fst

## -----------------------------------------------------------------------------------------------------------------------------
#Leading eigenvalue of Zs/sqrt(n)
spectral_fst$leadingeigenZs

## -----------------------------------------------------------------------------------------------------------------------------
#Smallest non zero eigenvalue of Zst/sqrt(n) (first one = Fst)
spectral_fst$eigenZst[1]
spectral_fst$eigenZst[1] > spectral_fst$leadingeigenZs

## -----------------------------------------------------------------------------------------------------------------------------
#RMT prediction of the leading eigenvalue of Zs/sqrt(n)
(1 - 0.05412755)*(1/sqrt(148) + 1/sqrt(18967))^2
spectral_fst$RMTprediction

## -----------------------------------------------------------------------------------------------------------------------------
#Fst approximation using PCA of the full matrix
spectral_fst$Fst_approximation

## ---- fig.show='hold'---------------------------------------------------------------------------------------------------------
plot(spectral_fst$eigenZs[1:10] * 100, ylab="variance (%)", main="Scree plot", ylim=c(0,8))
plot(spectral_fst$pcZs, col = pop, main="PC Plot")

## -----------------------------------------------------------------------------------------------------------------------------
data('athaliana')

## -----------------------------------------------------------------------------------------------------------------------------
spectral_athaliana <- compute_partition(athaliana$genotype, athaliana$pop)
wright_fst_athaliana <- mean_wright_fst(athaliana$genotype, athaliana$pop)

#mean Fst from Wright's formula
wright_fst_athaliana

#Fst estimate from the between population matrix Zst
spectral_athaliana$Fst

## -----------------------------------------------------------------------------------------------------------------------------
#Smallest non zero eigenvalue of Zst
spectral_athaliana$eigenZst[1]
#Leading eigenvalue of Zs
spectral_athaliana$leadingeigenZs
spectral_athaliana$eigenZst[1] > spectral_athaliana$leadingeigenZs

## -----------------------------------------------------------------------------------------------------------------------------
#RMT prediction of the leading eigenvalue of Zs
spectral_athaliana$RMTprediction

#Fst approximation using PCA of the full matrix
spectral_athaliana$Fst_approximation

## -----------------------------------------------------------------------------------------------------------------------------
spectral_adjusted_athaliana <-compute_partition(athaliana$genotype, athaliana$pop, Y = athaliana$bio)

## -----------------------------------------------------------------------------------------------------------------------------
# Adjusted Fst 
spectral_adjusted_athaliana$Fst

## -----------------------------------------------------------------------------------------------------------------------------
(spectral_athaliana$Fst - spectral_adjusted_athaliana$Fst)/spectral_athaliana$Fst

## ---- fig.width = 7, fig.height = 7-------------------------------------------------------------------------------------------
m <- cbind(c(spectral_athaliana$eigenZst[1], spectral_athaliana$eigenZs[1:5]),
           c(spectral_adjusted_athaliana$eigenZst[1], spectral_adjusted_athaliana$eigenZs[1:5]))
plot(seq(1,12)/2,as.numeric(t(m)), 
     col = rep(c("darkblue","orange"), length = 12), 
     type = "h", lwd = 12, 
     ylab = "Proportion of variance",
     xlab = "Axis",
    # main = "Structure Components",
     cex.axis = 2,
    cex = 2,
     las = 1)

legend(x = 3.5, y = 0.077, 
       legend = c("All SNPs", "Adjusted data"),
       col = c("darkblue","orange"),
       lty = c(1,1),
       pch = 19,
       cex = 2)

## -----------------------------------------------------------------------------------------------------------------------------
#Smallest non zero eigenvalue of Zst
spectral_adjusted_athaliana$eigenZst[1]
#leading eigenvalue of Zs
spectral_adjusted_athaliana$leadingeigenZs
spectral_adjusted_athaliana$eigenZst[1] > spectral_adjusted_athaliana$leadingeigenZs

