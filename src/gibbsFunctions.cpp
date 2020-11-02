#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Rcpp function for rmultinomial(1, 1, p)
// [[Rcpp::export]]
IntegerVector rmultinomRcpp(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

// function to find the index of a string in a vector of strings
// [[Rcpp::export]]
int findStringLoc(StringVector inputVector, String id){
  int idLoc; // init storage for IDs
  int i;
  for(i = 0; i < inputVector.size(); i++)
    if(inputVector[i] == id) // check if input matches target
      idLoc = i;

  return(idLoc);
}

// function to calculate bij
// [[Rcpp::export]]
NumericVector bSumRcpp(double b_, NumericMatrix nw_, IntegerVector w_, int V_, int k_){
  int lenw = w_.size(); NumericVector bij(lenw);
  int i;
  for (i=0; i<lenw; i++){
    NumericVector nwK = nw_.column(k_);          // subset nw into column vector based on community k
    NumericVector nwKwv = nwK[seq(w_[i]+1, V_)]; // elements of nwK > w
    int nwKwvSum = sum(nwKwv);                   // sum elements of nwK > w
    bij(i) = b_ + nwKwvSum;
  }
  return(bij);
}

// function to calculate probabilities to sample Z in zinLDA Gibbs sampler
// [[Rcpp::export]]
NumericVector sampleZProbRcpp(int d, int w, IntegerVector kVect, arma::mat delta,
                              NumericMatrix nw, IntegerMatrix nd, IntegerVector ndsum,
                              int K, double a, double b, double alpha, int V){
  NumericVector p(K); // define global pTemp

  int i;
  for (i=0; i<K; i++){
    int k = kVect(i);
    arma::uvec w0 = arma::find(delta.col(k) == 0); // indices with delta == 0 for given k
    int k1 = min(arma::find(delta.col(k) == 0));   // find smallest index with delta == 0 for given k
    bool k1Condition = w == k1;                    // check is w = k1
    int vk = max(arma::find(delta.col(k) == 0));   // find largest index with delta == 0 for given k
    bool vkCondition = w == vk;                    // check is w == vk

    bool deltaCondition = delta(w, k) == 0;        // check is delta == 0

    NumericVector nwK = nw.column(k);   // subset nw into column vector based on community k

    if(deltaCondition) { // if delta == 0
      if(k1Condition) {  // if w = k1
        NumericVector nwKwv = nwK[seq(w+1, V)];                 // elements of nwK > w
        int nwKwvSum = sum(nwKwv);                              // sum elements of nwK > w
        p(i) = ((a + nw(w, k))/(a + nw(w, k) + b + nwKwvSum)) * // replace global pTemp w calc value
          ((nd(d, k) + alpha)/(ndsum(d) + K*alpha));
      }

      else if(vkCondition){ // if w == vk
        arma::uvec posW0 = w0.elem(arma::find(w0 < w));         // find = indices of w0<w, elem = elements of w0 of indices
        IntegerVector posW0iv = as<IntegerVector>(wrap(posW0)); // convert arma::uvec to SEXP to IntegerVector
        NumericVector nwKw0 = nwK[posW0iv];                     // extract rows (w0 < w)
        NumericVector aij = a + nwKw0;                          // extract rows (w0 < w) of nw for given k and add a
        NumericVector bij = bSumRcpp(b, nw, posW0iv, V, k);     // see bij function
        NumericVector bijRatio = bij/(aij + bij);
        double bijRatioProd = std::accumulate(bijRatio.begin(), bijRatio.end(), 1.0, std::multiplies<double>());
        p(i) = bijRatioProd * ((nd(d, k) + alpha)/(ndsum(d) + K*alpha)); // replace global pTemp w calculated value
      }

      else { // w != k1 and w != vk
        NumericVector nwKwv = nwK[seq(w+1, V)];                 // elements of nwK > w
        int nwKwvSum = sum(nwKwv);                              // sum elements of nwK > w
        arma::uvec posW0 = w0.elem(arma::find(w0 < w));         // find = indices of w0<w, elem = elements of w0 of indices
        IntegerVector posW0iv = as<IntegerVector>(wrap(posW0)); // convert arma::uvec to SEXP to IntegerVector
        NumericVector nwKw0 = nwK[posW0iv];                     // extract rows (w0 < w)
        NumericVector aij = a + nwKw0;                          // extract rows (w0 < w) of nw for given k and add a
        NumericVector bij = bSumRcpp(b, nw, posW0iv, V, k);     // see bij function
        NumericVector bijRatio = bij/(aij + bij);
        double bijRatioProd = std::accumulate(bijRatio.begin(), bijRatio.end(), 1.0, std::multiplies<double>());
        p(i) = ((a + nw(w, k))/(a + nw(w, k) + b + nwKwvSum)) * // replace global pTemp w calculated value
          bijRatioProd * ((nd(d, k) + alpha)/(ndsum(d) + K*alpha));
      }
    }
    else p(i) = 0; // if delta == 1 then p(z) == 0, replace global pTemp
  }
  p.names() = kVect; // + 1 for R indexing
  double sumP = sum(p);
  p = p/sumP;
  return(p);
}

// function to to sample Z in zinLDA Gibbs sampler
// [[Rcpp::export]]
List gibbSampleZRcpp(int D, List z, NumericMatrix nw, IntegerMatrix nd, IntegerVector ndsum,
                     IntegerVector kVect, arma::mat delta, int K, double a, double b, double alpha, int V){
  // deep copies
  List zDC = clone(z);
  NumericMatrix nwDC = clone(nw);
  IntegerMatrix ndDC = clone(nd);

  int d;
  for (d=0; d<D; d++){
    IntegerVector zd = zDC[d];                 // z for sample d
    int N = zd.size();                         // N for the dth sample

    int n;
    for (n=0; n<N; n++){
      int community = zd(n);                        // community for the nth taxon in sample d; -1 for indexing from R to C++
      StringVector taxaD = zd.names();              // taxa for sample d
      String taxon = taxaD(n);                      // nth taxon in sample d
      int w = findStringLoc(rownames(nwDC), taxon); // int value (row location) for nth taxon in sample d
      nwDC(w, community) -= 1;                      // - 1 bc gibbs sampling for -ith assignment
      ndDC(d, community) -= 1;                      // - 1 bc gibbs sampling for -ith assignment
      ndsum(d)     -= 1;                            // - 1 bc gibbs sampling for -ith assignment

      NumericVector p = sampleZProbRcpp(d, w, kVect, delta, nwDC, ndDC, ndsum, K, a, b, alpha, V);

      IntegerVector zSampled = rmultinomRcpp(p); // rmultinomial only taking n = 1 size = 1

      int scommunity = which_max(zSampled);

      zd(n) = scommunity;       // update the community for the nth taxon in the dth sample
      zDC[d] = zd;              // update z with the community for the nth taxon in the dth sample
      nwDC(w, scommunity) += 1; // +1 for new sampled community assignment
      ndDC(d, scommunity) += 1; // +1 for new sampled community assignment
      ndsum(d)      += 1;       // +1 for new sampled community assignment
    }
  }
  return(Rcpp::List::create(Rcpp::Named("z") = zDC,
                            Rcpp::Named("nw") = nwDC,
                            Rcpp::Named("nd") = ndDC));
}
