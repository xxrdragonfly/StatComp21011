#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param n one parameter of  bivariate density
//' @param a one parameter of  bivariate density
//' @param b one parameter of  bivariate density
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(100,10,2,3)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N,int n,int a,int b) {
  NumericMatrix X(N,2);
  X(0,0)=0.2;
  X(0,1)=0.3;
  for(int i = 1; i < N; i++) {
    double y = X(i-1, 1);
    X(i, 0) = rbinom(1,n,y)[0];
    double x = X(i, 0);
    X(i, 1) = rbeta(1,x+a,n-x+b)[0];
  }
  return X;
}

