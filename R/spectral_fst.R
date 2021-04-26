
#' @title compute_partition
#
#' @description This function computes a within-population (residual) matrix, Zs, and a between-population matrix, Zst, for a genotype matrix, Z, and for K predefined populations. It allows the genotype matrix to  be adjusted for specific covariates, Y, representing technical effects, such as coverage, or environmental variables, such as bioclimatic data. The function returns a spectral estimate of the average value of Fst across loci, an approximation of Fst by PCA eigenvalues, and the leading eigenvalue of the residual matrix.
#'
#' @param genotype A matrix of genotypes with n rows and L columns, where n is number of individuals and L is the number of loci. Genotypes can be
#' haploid (0,1) or diploid (0,1,2). Missing values are not allowed.
#' @param population_labels A vector of length n containing population labels for each sample in the genotype matrix. Population labels are integers
#' @param Y A matrix of covariates with n rows and d columns, where d is the number of covariates. The values are used to compute an adjusted estimate of Fst. Missing values are not allowed. If \code{NULL}, then spectral estimates are unadjusted.

#'
#' @return spectral_result
#' Zs : A matrix of n rows and L columns in the case of haploid genotype, 2*n rows and L columns in the case of diploid genotype. It corresponds to the within population matrix Zs
#' Zst : A matrix of n rows and L columns in the case of haploid genotype, 2*n rows and L columns in the case of diploid genotype. It corresponds to the between population matrix Zst
#' Fst : A real value. It is computed using the squared norm of Zst and it is equal to the average value of Wright's Fst over all loci included in the genotype matrix
#' Fst_approximation : A real value. The sum of the first (nb_pop − 1) eigenvalues of scaled PCA for genotype matrix Z
#' leadingeigenZs : A real value. The leading eigenvalue of Zs
#' RMTprediction : A real value. The prediction of the leading eigenvalue of the residual matrix using RMT
#' eigenZs : A vector of length n for haploid genotype, 2*n for diploid. It corresponds to the eigenvalues of scaled Zs
#' eigenZst : A vector of length n for haploid genotype, 2*n for diploid. It corresponds to the eigenvalues of scaled Zst
#' eigenZ : A vector of length n for haploid genotype, 2*n for diploid. It corresponds to the eigenvalues of scaled Z
#' pcZs : A matrix of size n \times 2 for haploid genotype, 2n \times 2 for diploid. The projection of Zs on 2 first axis of PCA
#' pcZst : A matrix of size n \times 2 for haploid genotype, 2n \times 2 for diploid. The projection of Zst on 2 first axis of PCA
#' pcZ : A matrix of size n \times 2 for haploid genotype, 2n \times 2 for diploid. The projection of Z on 2 first axis of PCA


compute_partition <- function(genotype, population_labels, Y=NULL){

  # Several check before running the functions
  # non null values

  # Check if it is a diploid matrix
  if (2 %in% genotype){
    haploid_object <- haploidisation(genotype, population_labels)
    genotype <- haploid_object$haploid_matrix
    population_labels <- haploid_object$haploid_population_labels
  }
  # unique_pop stores all distinct population identifier

  unique_pop <- sort(unique(population_labels))
  nb_pop <- length(unique_pop)

  # pop_count count numbers of individuals in each population
  # pop_distribution stores
  pop_count <- table(population_labels)
  pop_distribution <- pop_count / sum(pop_count)

  # store n and L values
  n <- dim(genotype)[1]
  L <- dim(genotype)[2]

  # If you want to adjust for variables (e.g coverage or environmental variables)

  if (! is.null(Y)){
    nb_var <- ncol(Y)
    lfmm_genotype <- lfmm2(genotype, Y, n-(nb_var + 1))
    mod.lm <- lm(genotype ~ ., data=data.frame(Y, lfmm_genotype@U))
    sm <- summary(mod.lm)
    effect.size <- sapply(sm, FUN = function(x) x$coeff[1:(nb_var + 1), 1])
    Y <- cbind(rep(1.0, n), Y)
    genotype <- genotype - Y %*% effect.size
  }


  # store mean and sd for each column
  colmean <- colMeans(genotype)
  colsd <- apply(genotype, 2, sd)
  # We get rid of potential 0 variance column
  genotype <- genotype[,colsd != 0]
  colmean <- colmean[colsd != 0]
  colsd <- colsd[colsd != 0]
  # Updata of L values in case of column deletion
  L <- dim(genotype)[2]



  #==============================#
  # We compute Zst and Zs matrix #
  #==============================#

  Zst <- genotype
  Zs <- genotype

  for (i in seq(1, nb_pop)){
    pop_i_index <- population_labels == unique_pop[i]

    # We build a matrix of n_i rows (where ni is the number of individuals of population i) called mean_matrix_pop_i
    # This matrix has n_i identical lines. Each line is of length L with lk is the mean allele frequency in population
    # i at locus k
    # mean_matrix has the same format but each line is of length L with lk is the mean allele frequency in individuals
    # as a whole at locus k
    # sd_matrix has the same format but each line is of length L with lk is sd in individuals
    # as a whole at locus k. This matrix aims at realising normalization

    n_i <- pop_count[i]
    colmean_pop_i <- colMeans(genotype[pop_i_index,])
    mean_matrix_pop_i <- do.call("rbind", replicate(n_i, colmean_pop_i, simplify=FALSE))
    mean_matrix <- do.call("rbind", replicate(n_i, colmean, simplify=FALSE))
    sd_matrix <- do.call("rbind", replicate(n_i, colsd, simplify=FALSE))

    #       Zs matrix         #

    Zs[pop_i_index,] <- (genotype[pop_i_index,] - mean_matrix_pop_i) * (1/sd_matrix)

    #       Zst matrix        #

    Zst[pop_i_index,] <- (mean_matrix - mean_matrix_pop_i) * (1/sd_matrix)
  }

  #==========================================#
  # We compute pca on the obtained matrices #
  #==========================================#

  pc_z <- prcomp(genotype, scale=T, center=T)
  pc_z_st <- prcomp(Zst, scale=F)
  pc_z_s <- prcomp(Zs, scale=F)


  Fst <- sum(pc_z_st$sdev^2 / L)
  Fst_approximation <- sum((pc_z$sdev^2 / L)[1:(nb_pop-1)])
  leadingeigenZs <- (pc_z_s$sdev^2 / L)[1]

  #==========================================#
  #        Random Matrix prediction          #
  #==========================================#

  RMTprediction <- (1 - Fst)*((1/sqrt(n-nb_pop)) + (1/sqrt(L)))^2


  return(list(Zs = Zs, Zst=Zst, Fst=Fst, Fst_approximation=Fst_approximation, leadingeigenZs=leadingeigenZs, RMTprediction=RMTprediction, eigenZs=pc_z_s$sdev^2/L, eigenZst=pc_z_st$sdev^2/L, eigenZ=pc_z$sdev^2/L, pcZs = pc_z_s$x[,1:2], pcZst = pc_z_st$x[,1:2], pcZ = pc_z$x[,1:2]))
}


#' @title mean_wright_fst
#
#' @description This function aims at computing mean value of Fst over loci on an haploid genotype matrix
#
#' @param genotype A genotype matrix of haploid individuals with n rows (=individuals) and p columns (=loci).
#' @param population_labels A vector of length n containing population labels for each sample in the genotype matrix. Population labels are integers
#
#' @return mean_wright_fst
#' mean_wright_fst : The mean value of Fst


mean_wright_fst <- function(genotype, population_labels){
  # unique_pop stores all distinct population identifier

  unique_pop <- sort(unique(population_labels))
  nb_pop <- length(unique_pop)

  # pop_count count numbers of individuals in each population
  # pop_distribution stores
  pop_count <- table(population_labels)
  pop_distribution <- pop_count / sum(pop_count)

  # store n and L values
  n <- dim(genotype)[1]
  L <- dim(genotype)[2]

  # Allele frequency
  P <- colMeans(genotype)
  Q <- 1 - P
  sumcipiqi <- rep(0, L)
  for (i in seq(1, nb_pop)){
    ci <- pop_distribution[i]
    pi <- colMeans(genotype[(population_labels == unique_pop[i]),])
    qi <- 1 - pi
    sumcipiqi <- sumcipiqi + ci * pi * qi
  }

  # Fst at all loci
  Fst_loci <- (P*Q - sumcipiqi) / (P*Q)


  return(mean(Fst_loci))
}



#' @title haploidisation
#
#' @description This function takes in input a diploid genotype matrix of n lines (=individuals)
#' and L columns (= loci) and returns a 2n * L equivalent haploid genotype matrix
#
#' @param genotype A genotype matrix of diploid individuals with n rows (=individuals) and p columns (=loci).
#' @param population_labels A vector of length n containing population labels for each sample in the genotype matrix. Population labels are integers
#'
#' @return
#' haploid_genotype : The correspondent haploid matrix
#' haploid_population_labels : The new vector of population assignment corresponding to the new haploid matrix

haploidisation <- function(genotype, population_labels){

  # store n and L values
  n <- dim(genotype)[1]
  L <- dim(genotype)[2]

  haploid_matrix <- matrix(rep(0,(2*n*L)), nrow=2*n, ncol=L)
  haploid_population_labels <- c()
  for (i in seq(1,n)){
    current_line <- genotype[i,]
    one_current_line <- current_line == 1
    two_current_line <- current_line == 2
    random_one <- rbinom(n=L, size=1, prob=0.5)

    new_line_1 <- two_current_line + one_current_line*random_one
    new_line_2 <- two_current_line + (one_current_line - one_current_line*random_one)

    haploid_matrix[(2*(i-1) + 1),] <- new_line_1
    haploid_matrix[2*i,] <- new_line_2
    current_pop <- population_labels[i]
    haploid_population_labels <- c(haploid_population_labels, rep(current_pop, 2))
  }

  return(list(haploid_matrix=haploid_matrix, haploid_population_labels=haploid_population_labels))
}
