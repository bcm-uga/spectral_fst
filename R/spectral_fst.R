
#' @title compute_partition
#
#' @description This function aims at computing matrices that we define in the paper as within population matrix Zs and between population matrix Zst
#
#' @param genotype A genotype matrix of haploid or transformed haploid individuals with n rows and p columns.
#' @param population_assignement A vector of length n, each element is the population of the corresponding individual in genotype matrix
#
#' @return spectral_result
#' Zs : The within population matrix Zs
#' Zst : The between population matrix Zst
#' Fst : Squared norm of Zst
#' Fst_estimate : the sum of the first (nb_pop − 1) eigenvalues of scaled PCA for genotype matrix
#' leadingeigenZs : The leading eigenvalue of Zs
#' MPestimate : the estimation of Marchenko Pastur
#' eigenZs : eigenvalues of scaled Zs
#' eigenZst : eigenvalues of scaled Zst
#' eigenZ : eigenvalues of scaled Z

compute_partition <- function(genotype, population_assignment){

  # unique_pop stores all distinct population identifier

  unique_pop <- sort(unique(population_assignment))
  nb_pop <- length(unique_pop)

  # pop_count count numbers of individuals in each population
  # pop_distribution stores
  pop_count <- table(population_assignment)
  pop_distribution <- pop_count / sum(pop_count)

  # store n and L values
  n <- dim(genotype)[1]
  L <- dim(genotype)[2]

  # store mean and sd for each column
  colmean <- colMeans(genotype)
  colsd <- apply(genotype, 2, sd)

  #==============================#
  # We compute Zst and Zs matrix #
  #==============================#

  Zst <- genotype
  Zs <- genotype

  for (i in seq(1, nb_pop)){
    pop_i_index <- population_assignment == unique_pop[i]

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
  Fst_estimate <- sum((pc_z$sdev^2 / L)[1:(nb_pop-1)])
  leadingeigenZs <- (pc_z_s$sdev^2 / L)[1]

  #==========================================#
  #         Marchenko Pastur estimate        #
  #==========================================#

  MPestimate <- (1 - Fst)*((1/sqrt(n-nb_pop)) + (1/sqrt(L)))^2


  return(list(Zs = Zs, Zst=Zst, Fst=Fst, Fst_estimate=Fst_estimate, leadingeigenZs=leadingeigenZs, MPestimate=MPestimate, eigenZs=pc_z_s$sdev^2/L, eigenZst=pc_z_st$sdev^2/L, eigenZ=pc_z$sdev^2/L))
}


#' @title mean_wright_fst
#
#' @description This function aims at computing mean value of Fst over loci on an haploid genotype matrix
#
#' @param genotype A genotype matrix of haploid individuals with n rows (=individuals) and p columns (=loci).
#' @param population_assignement A vector of length n, each element is the population of the corresponding individual in genotype matrix
#
#' @return mean_wright_fst
#' mean_wright_fst : The mean value of Fst


mean_wright_fst <- function(genotype, population_assignment){
  # unique_pop stores all distinct population identifier

  unique_pop <- sort(unique(population_assignment))
  nb_pop <- length(unique_pop)

  # pop_count count numbers of individuals in each population
  # pop_distribution stores
  pop_count <- table(population_assignment)
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
    pi <- colMeans(genotype[(population_assignment == unique_pop[i]),])
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
#' @param population_assignement A vector of length n, each element is the population of the corresponding individual in genotype matrix
#'
#' @return
#' haploid_genotype : The correspondent haploid matrix
#' haploid_population_assignment : The new vector of population assignment corresponding to the new haploid matrix

haploidisation <- function(genotype, population_assignment){

  # store n and L values
  n <- dim(genotype)[1]
  L <- dim(genotype)[2]

  haploid_matrix <- matrix(rep(0,(2*n*L)), nrow=2*n, ncol=L)
  haploid_population_assignment <- c()
  for (i in seq(1,n)){
    current_line <- genotype[i,]
    one_current_line <- current_line == 1
    two_current_line <- current_line == 2
    random_one <- rbinom(n=L, size=1, prob=0.5)

    new_line_1 <- two_current_line + one_current_line*random_one
    new_line_2 <- two_current_line + (one_current_line - one_current_line*random_one)

    haploid_matrix[(2*(i-1) + 1),] <- new_line_1
    haploid_matrix[2*i,] <- new_line_2
    current_pop <- population_assignment[i]
    haploid_population_assignment <- c(haploid_population_assignment, rep(current_pop, 2))
  }

  return(list(haploid_matrix=haploid_matrix, haploid_population_assignment=haploid_population_assignment))
}
