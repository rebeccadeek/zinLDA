#' The Discrete Uniform Distribution
#'
#' \code{rdu} generates a random sample of size \code{n} from a discrete uniform
#' distribution on the interval \code{min} to \code{max}.
#' @param n number of observations.
#' @param min lower limit of the distribution.
#' @param max upper limit of the distribution.
#'
#' @return A length n numeric vector of realizations from a DU(\code{min}, \code{max}) distribution.
#'
#' @export
#'
#' @examples
#' N = rdu(n = 50, min = 10000, max = 15000)
rdu <- function(n, min, max) {
  # error messages
  if(class(n)!="numeric" | length(n)!=1){
    stop("ERROR: n must be a length 1 numeric integer.")
  }else if(class(min)!="numeric" | length(min)!=1){
    stop("ERROR: min must be a length 1 numeric integer.")
  }else if(class(max)!="numeric" | length(max)!=1 | max<min){
    stop("ERROR: max must be a length 1 numeric integer such that max > min.")
  }

  sample(min:max, n, replace = TRUE)
}
