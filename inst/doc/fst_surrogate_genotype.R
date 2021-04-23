## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=7
)

## ------------------------------------------------------------------------
library(spectralfst)

## ------------------------------------------------------------------------
data('pop_label')
data('athaliana_n241_L10k')
data('climate')

genotype_athaliana <- as.matrix(athaliana_n241_L10k)
env_var <- as.matrix(climate)
pop_athaliana <- pop_label[,1]

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

