#' Zero-Inflated Latent Dirichlet Allocation
#'
#' Estimate a zinLDA model using Gibbs Sampling.
#'
#' @param sampleTaxaMatrix object of class \code{matrix} with integer entries.
#' @param K number of latent subcommunities.
#' @param a positive scalar symmetric hyperparameter of the ZIGD prior on beta.
#' @param b positive scalar symmetric hyperparameter of the ZIGD prior on beta.
#' @param pi scalar symmetric zero-inflation probability of the ZIGD prior on beta.
#' @param alpha positive scalar symmetric hyperparameter of the Dirichlet prior on theta.
#' @param runs length of the MCMC chain.
#'
#' @return An object of class \code{"zinLDA_gibbs"} containing parameter estimates
#' of \code{theta}, \code{beta}, and \code{delta} from each of the MCMC runs.
#'
#' @export
zinLDA <- function(sampleTaxaMatrix, K, alpha, pi, a, b, runs = 500) {
  # error messages
  if(class(alpha)!="numeric" | length(alpha)!=1 |
     class(pi)   !="numeric" | length(pi)!=1 |
     class(a)    !="numeric" | length(a)!=1 |
     class(b)    !="numeric" | length(b)!=1){
    stop("ERROR: This package assumes hyperparameters are symmetric. Each must be a length 1 double.")
  } else if(class(K)!="numeric" | length(K)!=1 | !(K %% 1 ==0)){
    stop("ERROR: K must be a length 1 integer.")
  } else if(class(sampleTaxaMatrix)[1]!="matrix" & typeof(sampleTaxaMatrix)!="double"){
    stop("ERROR: sampleTaxaMatrix must be an numeric matrix.")
  }

  pb = pbapply::timerProgressBar()

  # remove any columns with colSum() == 0ALL zero
  ##########################################################################################################################
  termCount <- colSums(sampleTaxaMatrix)
  termCount0 <- which(termCount == 0)

  if (length(termCount0) > 0) {
    sampleTaxaMatrix <- sampleTaxaMatrix[,-termCount0]
  }
  ##########################################################################################################################

  # initialize z
  ##########################################################################################################################
  LDAsave <- topicmodels::LDA(sampleTaxaMatrix, k = K, method = "Gibbs")

  termDist <- topicmodels::posterior(LDAsave)$terms   # beta-hat
  topicDist <- topicmodels::posterior(LDAsave)$topics # theta-hat

  D <- nrow(sampleTaxaMatrix)
  V <- ncol(sampleTaxaMatrix)

  z <- vector(mode = "list", length = D)  # Empty list to hold z assignment

  for (d in 1:D) {
    nonzeroV <- which(sampleTaxaMatrix[d,] != 0)               # Find non-zero counts
    initZd <- vector(mode = "list", length = length(nonzeroV)) # Temporary list to hold z assignment

    for (v in seq_along(nonzeroV)) { # For every unique non-zero taxa

      init <- function(d, i, sampleTaxaMatrix, termDist){
        w <- colnames(sampleTaxaMatrix)[i]                                                    # taxa name
        zCount <- stats::rmultinom(n = 1, size = sampleTaxaMatrix[d, w], prob = termDist[,w]) # sample from LDA posterior
        zTemp <- rep(as.numeric(row.names(zCount)), times = zCount)                           # multinom counts into z-vect
        names(zTemp) <- rep(w, sampleTaxaMatrix[d, w])                                        # named vector with taxa
        return(zTemp)
      }

      initZd[[v]] <- init(d, nonzeroV[[v]], sampleTaxaMatrix, termDist)                # sample initial Z using LDA est.
    }

    initZd <- unlist(initZd)
    z[[d]] <- initZd
  }
  ##########################################################################################################################

  # nw, nd, ndsum
  ##########################################################################################################################
  nwFunction <- function(z) { # the number of times each taxa is assigned to each subcommunity / sample

    nwLapply <- function(x) { # Function for the lapply
      zeroTerms <- colnames(sampleTaxaMatrix)[which(!colnames(sampleTaxaMatrix) %in% unique(names(x)))] # 0 counts all topics
      zeroRows <- matrix(data = 0,                      # a matrix
                         nrow = length(zeroTerms),      # with terms that have zero counts
                         ncol = K,                      # for every topic
                         dimnames = list(c(zeroTerms),
                                         c(1:K)))

      nwTemp <- table(tibble::enframe(x)) # table of of vector z

      nwCondition <- ncol(nwTemp) < K

      if(nwCondition) {
        zeroTopic <- which(! rownames(termDist) %in% colnames(nwTemp)) # missing topics
        zeroColumns <- matrix(data = 0 ,                               # matrix
                              nrow = nrow(nwTemp),                     # of zero counts
                              ncol = length(zeroTopic),                # for missing columns/topics
                              dimnames = list(c(rownames(nwTemp)),
                                              c(zeroTopic)))

        nwTemp <- cbind(nwTemp, zeroColumns)
        nwTempTerms <- rownames(nwTemp)
        nwTemp <- matrix(nwTemp[,order(colnames(nwTemp))], # reorder
                         ncol = K,
                         dimnames = list(c(rownames(nwTemp)),
                                         c(1:K)))
      }

      nwTemp <- rbind(nwTemp, zeroRows)
      nwTemp <- nwTemp[order(factor(rownames(nwTemp), levels = colnames(sampleTaxaMatrix))),] # order rows by terms

      return(nwTemp)
    }

    nwList <- lapply(z, nwLapply) # nw matrix for each d
    nw <- Reduce("+", nwList)     # total nw matrix

    return(nw)
  }

  nw <- nwFunction(z)

  ndFunction <- function(z) { # the number of times each subcommunity is assigned to in each sample

    ndLapply <- function(x) {
      ndTemp <- table(x)      # number of times each subcommunity occurs in each sample

      ndCondition <- length(ndTemp) < K

      if(ndCondition) {
        zeroTopic <- which(! rownames(termDist) %in% names(ndTemp)) # missing subcommunities
        zeroColumns <- matrix(data = 0 ,                # matrix
                              nrow = 1,                 # of zero counts
                              ncol = length(zeroTopic), # for missing subcommunities
                              dimnames = list(c(""),
                                              c(zeroTopic)))

        nd <- c(ndTemp, zeroColumns)
        names(nd) <- c(names(ndTemp), colnames(zeroColumns))
        nd <- nd[order(names(nd))]

        ndTemp <- nd
      }
      return(ndTemp)
    }

    ndList <- lapply(z, ndLapply) # nd matrix for each d
    nd <- do.call(rbind, ndList)  # total nw matrix
    rownames(nd) <- rownames(sampleTaxaMatrix)

    return(nd)
  }

  nd <- ndFunction(z)

  ndsum <- rowSums(nd) # total number of subcommunity assignments/sample
  ##########################################################################################################################

  # sampling
  ##########################################################################################################################
  bSum = function(b, nw, w, V, k){
    bij = b + sum(nw[(w+1):V, k])
    return(bij)
  }
  bSumVect = Vectorize(bSum, vectorize.args = "w")

  # source Rcpp
  # Rcpp::sourceCpp(here::here("/code/01bi_20200402_GibbsZ.cpp"), embeddedR = FALSE)

  # adjust indexing for c++
  kVect <- 0:(K-1)
  Vm1 <- V-1

  # initialize empty lists
  deltaGibbs <- vector(mode = "list", length = runs)
  zGibbs <- vector(mode = "list", length = runs)
  betaGibbs <- vector(mode = "list", length = runs)
  thetaGibbs <- vector(mode = "list", length = runs)

  for (i in 1:runs) {
    # sample delta
    ######################################################################################################################
    sampleDelta <- function(nw) {
      deltaTemp <- matrix(data = 0,        # zero delta matrix
                          nrow = nrow(nw), # of dimension equal to nw
                          ncol = ncol(nw), # for each taxa/subcommunity
                          dimnames = list(c(rownames(nw)),
                                          c(colnames(nw))))

      sampleDeltaSc <- function(w, k) { # scalar sample delta function
        p1 <- pi / (pi + (1 - pi) * beta(a + nw[w, k], b + sum(nw[(w+1):V, k])) / beta(a, b)) # Pr(delta = 1)

        return(stats::rbinom(1, 1, p1)) # sample ber(p = p1)
      }

      sampleDeltaVect <- Vectorize(sampleDeltaSc, vectorize.args = "w")

      for (k in 1:K) {
        temp <- which(nw[,k] == 0)       # only sample for nw counts = 0
        last <- utils::tail(temp, n = 1) # NA produced for last word if p1 used
        if (length(temp) != 0){
          if (length(temp) == 1) {
            deltaTemp[last, k] <- stats::rbinom(1, 1, pi) # sample the last from ber(pi)
          }
          else{
            deltaTemp[temp[-length(temp)], k] <- sampleDeltaVect(w = temp[-length(temp)], k = k)
            deltaTemp[last, k] <- stats::rbinom(1, 1, pi) # sample the last from bern(pi)
          }
        }
      }
      return(deltaTemp)
    }

    delta <- sampleDelta(nw)
    ######################################################################################################################

    # sample Z
    ######################################################################################################################
    z <- lapply(z, function(x){x-1}) # -1 in topic indicators for C++ indexing
    gibbsZRcpp <- gibbSampleZRcpp(D, z, nw, nd, ndsum, kVect, delta, K, a, b, alpha, Vm1)

    z <- gibbsZRcpp$z
    z  <- lapply(z, function(x){x+1}) # +1 in topic indicators for R indexing
    nw <- gibbsZRcpp$nw
    nd <- gibbsZRcpp$nd
    ######################################################################################################################

    # Estimate Beta and Theta
    ######################################################################################################################
    # beta
    estimateBeta <- function(nw){
      betaTemp <- as.data.frame(matrix(data = 0,        # zero beta matrix
                                       nrow = nrow(nw), # of dimension equal to nw
                                       ncol = ncol(nw), # for each taxa/subcommunity
                                       dimnames = list(c(rownames(nw)),
                                                       c(colnames(nw)))))


      estimateBetaSc <- function(w, k) {
        k1 <- min(which(delta[,k] == 0))
        k1Condition <- w == k1
        vk <- max(which(delta[,k] == 0))
        vkCondition <- w == vk

        if(k1Condition) {
          pTemp = ((a + nw[w, k])/(a + nw[w, k] + b + sum(nw[(w+1):V, k])))
        }

        else if(vkCondition){
          w0 = which(delta[,k] == 0)
          aij = a + nw[w0[w0 < w], k]
          bij = bSumVect(b = b, nw = nw, w = w0[w0 < w], V = V, k = k)
          bijRatio = bij/(aij + bij)
          pTemp = prod(bijRatio)
        }

        else {
          w0 = which(delta[,k] == 0)
          aij = a + nw[w0[w0 < w], k]
          bij = bSumVect(b = b, nw = nw, w = w0[w0 < w], V = V, k = k)
          bijRatio = bij/(aij + bij)
          pTemp = ((a + nw[w, k])/(a + nw[w, k] + b + sum(nw[(w+1):V, k]))) *
            prod(bijRatio)
        }
        return(pTemp)
      }

      estimateBetaVect <- Vectorize(estimateBetaSc, vectorize.args = "w")

      for (k in 1:K) {
        temp <- which(delta[,k] == 0) # only for delta = 0 (non-structural zero, n > 0)
        betaTemp[temp, k] <- estimateBetaVect(w = temp, k = k)
      }
      return(t(betaTemp))
    }

    betaGibbs[[i]] <- estimateBeta(nw)

    # theta
    estimateTheta <- function(nd) {
      thetaTemp <- as.data.frame(matrix(data = 0,        # zero theta matrix
                                        nrow = nrow(nd), # of dimension equal to nd
                                        ncol = ncol(nd), # for each subcommunity/sample
                                        dimnames = list(c(rownames(nd)),
                                                        c(colnames(nd)))))

      estimateThetaSc <- function(d, k) {
        theta1 = ((nd[d, k] + alpha)/(ndsum[[d]] + K*alpha))
      }

      estimateThetaVect <- Vectorize(estimateThetaSc, vectorize.args = "d")

      for (k in 1:K) {
        thetaTemp[rownames(nd), k] <- estimateThetaVect(d = rownames(nd), k = k)
      }
      return(thetaTemp)
    }

    thetaGibbs[[i]] <- estimateTheta(nd)

    deltaGibbs[[i]] <- delta

    # run progress bar
    pbapply::setTimerProgressBar(pb, i/runs)
    ######################################################################################################################
  }
  ##########################################################################################################################

  # output
  out = list(delta = deltaGibbs,
             beta = betaGibbs,
             theta = thetaGibbs)

  # output class
  class(out) = "zinLDA_gibbs"

  return(out)
}
