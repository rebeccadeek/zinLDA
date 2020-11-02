#' zinLDA Simulated Data
#'
#' \code{simulateZINLDA} generates sparse count data according to a zero-inflated latent Dirichlet allocation model.
#'
#' @param D number of samples.
#' @param V number of unique taxa.
#' @param N vector of length D containing the total number of sequencing readings per sample.
#' @param K number of latent subcommunities.
#' @param Alpha scalar symmetric hyperparameter of the Dirichlet prior on theta.
#' @param Pi scalar symmetric zero-inflation hyperparameter of the ZIGD prior on beta.
#' @param a scalar symmetric hyperparameter of the ZIGD prior on beta.
#' @param b scalar symmetric hyperparameter of the ZIGD prior on beta.
#'
#' @return A list containing the following elements:
#' \item{cohort}{D-length list of character vectors containing the taxa assigned
#' to each sequencing read in each microbial sample.}
#' \item{z}{D-length list of vectors containing the subcommunity assignments for
#' each sequencing read in each microbial sample.}
#' \item{sampleTaxaMatrix}{matrix of counts, analogous to an OTU matrix.}
#' \item{theta}{matrix of subcommunity probabilities per sample.}
#' \item{beta}{matrix of taxa probabilities per subcommunity.}
#' \item{delta}{matrix of structural zero indicators for each taxa and subcommunity.}
#'
#' @export
#'
#' @examples
#' N.d = rdu(n = 50, min = 10000, max = 15000)
#' simData = simulateZINLDA(D = 50, V = 100, N = N.d, K = 5, Alpha = .1,
#'                          Pi = 0.4, a = .05, b = 10)
simulateZINLDA <- function(D, V, N, K, Alpha, Pi, a, b) {
  # error messages
  if(class(Alpha)!="numeric" | length(Alpha)!=1 |
     class(Pi)   !="numeric" | length(Pi)!=1 |
     class(a)    !="numeric" | length(a)!=1 |
     class(b)    !="numeric" | length(b)!=1){
    stop("ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  } else if(class(D)!="numeric" | length(D)!=1 | !(D %% 1 ==0)){
    stop("ERROR: D must be a length 1 integer.")
  } else if(class(V)!="numeric" | length(V)!=1 | !(V %% 1 ==0)){
    stop("ERROR: V must be a length 1 integer.")
  } else if(class(K)!="numeric" | length(K)!=1 | !(K %% 1 ==0)){
    stop("ERROR: K must be a length 1 integer.")
  } else if(class(N)!="integer" | length(N)!=D){
    stop("ERROR: N must be a length D numeric vector.")
  }

  # data labels
  taxa <- paste0("Taxa ",seq(V))
  communities <- paste0("Community ", seq(K))
  samples <- paste0("Sample ", seq(D))

  # simulate delta
  Delta = do.call(cbind, lapply(1:V, function(i){stats::rbinom(n = K, size = 1, prob = Pi)}))
  rownames(Delta) <- communities
  colnames(Delta) <- taxa

  # resample delta for i with delta_ij == 1 for every j until atleast 1 delta_ij = 0
  deltaCheck = T
  while(deltaCheck){
    V.rs = which(colSums(Delta) == K)
    n.V.rs = length(V.rs)
    Delta[,V.rs] = do.call(cbind, lapply(1:n.V.rs, function(i){stats::rbinom(n = K, size = 1, prob = Pi)}))
    deltaCheck = length(which(colSums(Delta) == K))
  }

  # sample beta ~ ZIGD(pi, a, b)
  Beta <- ifelse(Delta == 1, 0, NA)

  for (j in 1:K) {
    Q_j_index <- which(Delta[j,] == 0)
    Q_j <- stats::rbeta(length(Q_j_index)-1, a, b)

    for (l in seq_along(Q_j_index[1:length(Q_j_index)-1])) { # sequence along first K non-zero Q
      Beta[j, Q_j_index[l]] <- Q_j[l]*prod(1 - Q_j[0:(l-1)])  # beta_{ij} = Q_{ij} prod_{l=1}^{i-1} (1-Q_{lj})
    }
    # K+1 beta
    Beta[j, Q_j_index[length(Q_j_index)]] <- 1 - sum(Beta[j,], na.rm = T)
  }

  # sample theta ~ Dirichlet(alpha)
  Theta <- MCMCpack::rdirichlet(n = D, alpha = rep(Alpha, K))
  rownames(Theta) <- samples
  colnames(Theta) <- communities

  # function to generate samples
  generateSample <- function(N, theta, beta){

    sample <- vector(length = N)
    z_d <- vector(length = N)

    for (n in 1:N) {
      # For each N in the sample
      z_n <- stats::rmultinom(1, 1, theta)                  # First choose a community z_n ~ multinomial(theta)
      w_n <- stats::rmultinom(1, 1, beta[which(z_n == 1),]) # Second choose a word ~ Pr(w_n | z_n, beta)

      sample[n] <- colnames(beta)[which(w_n == 1)]
      z_d[n] <- which(z_n == 1)
      names(z_d)[n] <- sample[n]
    }

    return(list("sample" = sample,
                "z" = z_d))
  }

  # generate cohort
  cohort <- vector(mode = "list", length = D)
  z <- vector(mode = "list", length = D)
  for (d in 1:D) {
    sample.d <- generateSample(N[d], Theta[d,], Beta)
    cohort[[d]] <- sample.d[["sample"]]
    z[[d]] <- sample.d[["z"]]
  }

  # convert sample into a frequency count of each unique word
  sampleTaxaFreq <- lapply(cohort, table)

  # sample-taxa frequency matrix
  sampleTaxaMatrix <- matrix(data = 0, nrow = D, ncol = V)
  rownames(sampleTaxaMatrix) <- samples
  colnames(sampleTaxaMatrix) <- taxa

  # assign non-zero count to taxa in sample
  for (d in 1:D) {
    sampleTaxaMatrix[d, names(sampleTaxaFreq[[d]])] <- sampleTaxaFreq[[d]]
  }

  return(list("cohort" = cohort,
              "z" = z,
              "sampleTaxaMatrix" = sampleTaxaMatrix,
              "theta" = Theta,
              "beta" = Beta,
              "delta" = Delta))
}
