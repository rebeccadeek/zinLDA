#' Posterior Estimates
#'
#' Calculate the posterior estimates of \code{delta}, \code{beta}, and \code{theta} for a fitted zinLDA model.
#'
#' @param obj object of class \code{"zinLDA_gibbs"}.
#' @param burnin number of burn-in runs of the MCMC.
#'
#' @return A list with the following elements:
#' \item{delta}{posterior indicator of structural zeros.}
#' \item{beta}{posterior probabilities of each taxon belonging to each subcommunity.}
#' \item{theta}{posterior probabilities of each subcommunity in each sample.}
#'
#' @export
posterior = function(obj, burnin=100){
  # error messages
  if(class(obj)!="zinLDA_gibbs"){
    stop("ERROR: obj must be a list of class 'zinLDA_gibbs'.")
  } else if(!(burnin %% 1 ==0)){
    stop("ERROR: burnin must be an interger.")
  } else if(burnin==0){
      warning("Warning: burnin==0. It is recommend to specify 0 < burnin.")
  } else if(burnin > length(obj$delta)){
    stop("ERROR: burnin must be smaller than runs.")
  }

  # 3D array of estimated theta from each run removing burn-in
  tensorMat.ls <- lapply(obj, function(x) {
    tens = array(unlist(x[-c(1:burnin)]), c(dim(x[[1]]), length(x[-c(1:burnin)])), dimnames = dimnames(x[[1]])) })

  posteriorEst.ls <- lapply(tensorMat.ls, function(x) {apply(x, 1:2, mean)})

  posteriorEst.ls$delta = ifelse(posteriorEst.ls$delta >= 0.5, 1, 0) # posterior pi to delta
  posteriorEst.ls$delta = t(posteriorEst.ls$delta)

  posteriorEst.ls$beta[which(posteriorEst.ls$delta == 1)] = 0        # beta = 0 if delta = 1
  posteriorEst.ls$beta = sweep(posteriorEst.ls$beta, 1, rowSums(posteriorEst.ls$beta), "/")

  return(posteriorEst.ls)
}
