#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::colvec soft_thr(arma::colvec input, double lambda, int pos, arma::rowvec sp_locs) {
  arma::colvec out; out = input;
  int N = out.n_rows;
  double val;
  for(int i=0; i<N; i++){
    val = lambda * sp_locs(i); // 0 when sp_locs is 0, lambda when sp_locs is 1
    if(pos==1){
      if( (input(i)-val) > 0){
        out(i) = input(i) - val;
      } else {
        out(i) = 0;
      }
    } else{
      int sign_in;  sign_in = (input(i) > 0) - (input(i) < 0);
      if( (sign_in*input(i)-val) > 0){
        out(i) = sign_in*(sign_in*input(i)-val);
      } else {
        out(i) = 0;
      }
    }
  }
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
  
/*** R
soft_thr(-1*matrix(1:10), 3, 0, c(rep(1,8),rep(0,2)))
*/
  