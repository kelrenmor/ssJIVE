#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
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

// [[Rcpp::export]]
List sfpca_bic(arma::mat x, int K, arma::mat Omegv, arma::rowvec sp_locsv, 
                    arma::rowvec lamvs, arma::rowvec alphavs) {
  
  x = x.t(); // want x to have obs in rows, as given they are in cols
  int n = x.n_rows;
  int p = x.n_cols;
  
  // FIXME add in option for observation smoothness/sparsity later (for now just leave out)
  
  // Assume no non-negativity constraints unless otherwise noted
  int posv = 0;
  
  // By default, set starting values to SVD initialization
  // and set maxit at 1000, with 5 iters to determine alpha/lambda params
  int startv = 0;
  int maxit = 1e3;
  int inner_maxit = 1e8;
  int iterS = 5;
  int inner_iter;
  arma::mat U(n,K);
  arma::mat V(p,K);
  arma::vec D(K);
  U.zeros(); V.zeros(); D.zeros();
  
  int rv = alphavs.n_cols; 
  int rlv = lamvs.n_cols; 
  alphavs = arma::sort(alphavs,"descend");
  lamvs = arma::sort(lamvs,"descend");
  double Lu = max(eig_sym(arma::eye(n,n)+n*arma::zeros(n,n))) + .01; // max(alphaus)*Omegu
  double Lv = max(eig_sym(arma::eye(p,p)+p*max(alphavs)*Omegv)) + .01;
  arma::rowvec optavs(K); arma::rowvec optlvs(K);
  arma::mat Xhat = x;
  
  // Initialize U/V/D/etc
  arma::mat U_tmp(n,K);
  arma::mat V_tmp(p,K);
  arma::vec D_tmp(K);
  arma::mat u(n,1); 
  arma::vec d_tmp(1);
  arma::mat v(p,1); 
  arma::mat bicv(rv,rlv);
  arma::mat oldu(n,1);
  arma::mat oldv(p,1);
  
  // Things saved along the way
  double optlv;
  double optav;
  arma::colvec tmpui(n);
  arma::colvec oldui(n);
  arma::colvec utild(n);
  arma::colvec tmpvi(p);
  arma::colvec oldvi(p);
  arma::colvec vtild(p);
  
  for (int k=0; k < K; k++) {
    bicv.zeros();
    // FIXME should have this depend on startv so you don't need to initialize with SVD
    // Initialize to regular SVD solution
    arma::svd(U_tmp,D_tmp,V_tmp,Xhat);
    u = U_tmp.col(0);
    v = V_tmp.col(0);
    d_tmp = D_tmp(0);
    
    double indo = 1.0; int iter = 1; double thr = 1e-6;
    while( (indo>thr) & (iter<maxit)  & (iter<iterS)){
      oldu = u; // save old u
      oldv = v; // save old v
      
      // FIXME to allow u smooth/sparse need to change this
      int ru = 1; int rlu = 1; // FIXME to allow u smooth/sparse need to change this
      arma::cube us(n,ru,rlu); us.zeros(); 
      for(int j=0; j<ru; j++){
        arma::mat Su(n,n,arma::fill::eye); // FIXME + alphaus(j)*n*Omegu;
        for(int i=0; i<rlu; i++){
          double indiu = 1.0;
          inner_iter = 0;
          while((indiu>thr)  & (inner_iter<inner_maxit)){
            tmpui = us.subcube(arma::span(),arma::span(j),arma::span(i));
            oldui = tmpui;
            utild = tmpui + (Xhat*v - Su*tmpui)/Lu;
            //us(:,j,i) = soft_thr(utild,lamus(i)/Lu,posu); // FIXME if add functionality
            us.subcube(arma::span(),arma::span(j),arma::span(i)) = utild;
            tmpui = utild;
            if(arma::norm(tmpui, 2) > 0) {
              us.subcube(arma::span(),arma::span(j),arma::span(i))  = tmpui/sqrt(as_scalar(tmpui.t()*Su*tmpui));
            } else{
              us.subcube(arma::span(),arma::span(j),arma::span(i)) = arma::zeros(n,1);
            }
            tmpui = us.subcube(arma::span(),arma::span(j),arma::span(i));
            indiu = arma::norm(tmpui - oldui, 2)/arma::norm(oldui, 2);
            inner_iter += 1;
            if( inner_iter==inner_maxit ){ printf("Warning innter_maxit reached during phase 1."); }
          }
          arma::uvec actu; actu = arma::find(tmpui != 0); // indices of non-zero elements
          double dfu = arma::trace( (arma::eye(actu.n_elem,actu.n_elem)).i() ); // FIXME + n*alphaus(j)*Omegu(actu,actu)));
          //bicu(j,i) = log( norm(Xhat*v - utild)^2/n ) + .5*log(n)*dfu/n;
        }
      }
      //iu = bicu==min(min(bicu(sum(abs(us))~=0)));
      //[indau,indlu] = ind2sub([ru rlu],find(iu==1));
      //optlu = lamus(indlu); optau = alphaus(indau);
      u = us.subcube(arma::span(),arma::span(0),arma::span(0)); // FIXME

      // Update right singluar vector (v)
      arma::cube vs(p,rv,rlv); vs.zeros(); 
      for(int j=0; j<rv; j++){
        arma::mat Sv(p,p,arma::fill::eye); Sv = Sv + alphavs(j)*p*Omegv;
        for(int i=0; i<rlv; i++){
          double indiv = 1.0;
          arma::colvec tmpvi(p); arma::colvec oldvi(p); arma::colvec vtild(p);
          inner_iter = 0;
          while((indiv>thr) & (inner_iter<inner_maxit)){
            tmpvi = vs.subcube(arma::span(),arma::span(j),arma::span(i));
            oldvi = tmpvi;
            vtild = tmpvi + (Xhat.t()*u - Sv*tmpvi)/Lv;
            //printf("%.2f ", accu(vtild));
            vs.subcube(arma::span(),arma::span(j),arma::span(i)) = soft_thr(vtild,lamvs(i)/Lv,posv,sp_locsv);
            //printf("%.2f ", accu(vs.subcube(arma::span(),arma::span(j),arma::span(i))));
            tmpvi = vs.subcube(arma::span(),arma::span(j),arma::span(i));
            if(arma::norm(tmpvi, 2) > 0) {
              vs.subcube(arma::span(),arma::span(j),arma::span(i))  = tmpvi/sqrt(as_scalar(tmpvi.t()*Sv*tmpvi));
            } else{
              vs.subcube(arma::span(),arma::span(j),arma::span(i)) = arma::zeros(p,1);
            }
            tmpvi = vs.subcube(arma::span(),arma::span(j),arma::span(i));
            double old_norm = arma::norm(oldvi, 2);
            indiv = arma::norm(tmpvi - oldvi, 2)/old_norm;
            inner_iter += 1;
            if( inner_iter==inner_maxit ){ printf("Warning innter_maxit reached during phase 1."); }
          }
          arma::uvec actv; actv = arma::find(tmpvi != 0); // indices of non-zero elements
          double dfv = arma::trace( arma::inv( arma::eye(actv.n_elem,actv.n_elem) + p*alphavs(j)*Omegv(actv,actv) ) );
          double tmp = arma::norm(Xhat.t()*u - tmpvi,2);
          bicv(j,i) = log( .5*tmp*tmp/p ) + .5*log(p)*dfv/p;
        }
      }
      double min_bicv = bicv.max();
      int indav = rv-1; // set to no smoothing (if equal, will keep this)
      int indlv = rlv-1; // set to no sparsity (if equal, will keep this)
      for(int j=0; j<rv; j++){
        for(int i=0; i<rlv; i++){
          //printf("%.2f ", bicv(j,i));
          if(accu(vs.subcube(arma::span(),arma::span(j),arma::span(i))) != 0){ // sum all elements
            if(bicv(j,i) < min_bicv){ // if you beat best bic so far, save indices
              indav = j;
              indlv = i;
              min_bicv = bicv(j,i);
              
            }
          }
        }
        //printf("\n ");
      }
      optlv = lamvs(indlv); optav = alphavs(indav);
      //printf("\n ");
      //printf("%.4f ", optlv); printf("%.4f\n", optav); printf("\n ");
      v = vs.subcube(arma::span(),arma::span(indav),arma::span(indlv));
      
      iter = iter + 1;
      indo = arma::norm(u - oldu)/arma::norm(oldu) + arma::norm(v- oldv)/arma::norm(oldv);
      
    } // while( (indo>thr) & (iter<maxit)  & (iter<iterS)){
    
    arma::mat Su(n,n,arma::fill::eye); // FIXME if add capability
    arma::mat Sv(p,p,arma::fill::eye); Sv = Sv + optav*p*Omegv;
    
    while((indo>thr) & (iter<maxit)){
      oldu = u; oldv = v;
      double indiu = 1.0;
      inner_iter = 1;
      while((indiu>thr) & (inner_iter<inner_maxit)){
        oldui = u;
        utild = u + (Xhat*v - Su*u)/Lu;
        //u = soft_thr(utild,optlu/Lu,posu);
        if(norm(u)>0){
          u  = u/sqrt((u.t()*Su*u).eval()(0));
        }else{
          u.fill(0.0);
        }            
        indiu = arma::norm(u - oldui, 2)/arma::norm(oldui, 2);
        inner_iter += 1;
        if( inner_iter==inner_maxit ){ printf("Warning innter_maxit reached during phase 2."); }
      }
      double indiv = 1.0;
      inner_iter = 1;
      while((indiv>thr) & (inner_iter<inner_maxit)){
        oldvi = v;
        vtild = v + (Xhat.t()*u - Sv*v)/Lv;
        v = soft_thr(vtild,optlv/Lv,posv,sp_locsv);
        if(arma::norm(v,2)>0){
          v  = v/sqrt((v.t()*Sv*v).eval()(0));
        }else{
          v.fill(0.0);
        }            
        indiv = norm(v - oldvi)/norm(oldvi);
        inner_iter += 1;
        if( inner_iter==inner_maxit ){ printf("Warning innter_maxit reached during phase 2."); }
      }
      indo = arma::norm(oldu - u, 2)/arma::norm(oldu,2) + arma::norm(oldv - v,2)/arma::norm(oldv,2);
      iter = iter + 1;
      if( iter==maxit ){ printf("Warning maxit reached during phase 2."); }
    }
    u = u/arma::norm(u,2);
    v = v/arma::norm(v,2);
    double dd = (u.t()*Xhat*v).eval()(0);
    Xhat = Xhat - dd*u*v.t();
    U.col(k) = u; V.col(k) = v; D(k) = dd;
    optavs(k) = optav; optlvs(k) = optlv;
    
    //printf("%d iter",iter);
    
  } // for (int k=0; k < K; k++) {
  
  List ret;
  ret["U"] = U;
  ret["D"] = D;
  ret["V"] = V;
  ret["optlvs"] = optlvs;
  ret["optavs"] = optavs;
  return ret;
}

/*** R
# Code for testing that my version of SFPCA and the original
# are getting to the same result!

# Load libraries
library(R.matlab)

# Read in Matlab file
setwd('/Users/Kelly/toxic/papers/dim_reduction/code/ssJIVE/')
data <- readMat('matlab_comparison/test_inputs_struct.mat')
data <- data[[1]][,,1]

# Quick-run version (fix Lambda and Alpha)
res <- sfpca_bic(t(data[["datBoth"]]), c(data[["K"]]), data[["OmegBoth"]], 
                 t(matrix(1,ncol(data[["datBoth"]]))), 
                 c(0.5,0.1), 0.5)
*/
