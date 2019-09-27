ssJIVE = function(data_list,r,rIndiv,SpSm=rep('N',length(data_list)),
                  lambda='BIC',alpha='BIC',scale='y', SpSm_indiv=SpSm,
                  ConvergeThresh=1e-6,MaxIter=1e3){
# [J,A] = JIVE(data_list,r,rIndiv,scale,ConvergeThresh,MaxIter)
# Input:
#   - data_list: List of data_list matrices with matching columns
#   (data_list[[i]] gives the i'th data_list matrix with column c holding observation c)
#   - r: Given rank for joint structure
#   - rIndiv: Given vector of ranks for individual structure'
#   - SpSm: Given vector specifying whether sparsity/smoothness is to be imposed.
#           "N"=none, "Sp"=sparsity, "Sm"=smoothness, "SpSm"=sparsity and smoothness
#   - lambda: penalty parameter (or vector of parameters) for sparsity (lambda = 'BIC' estimates via BIC)
#   - alpha: penalty parameter (or vector of parameters) for smoothness (lambda = 'BIC' estimates via BIC)
#   - scale: Should the data_listsets be centered and scaled by total variation? (default = 'y')
#   - ConvergeThresh:  Convergence threshold (default = 10^(-8))
#   - MaxIter: Maximum numeber of iterations (default = 1000)
# Output:
#   - J: Matrix of sparse joint structure, of dimension (d1+d2+...+dk, n)
#   - A: Matrix of sparse individual structure
#Note: Use ssvd on J and A to extract sparse components.

if(F){ # TO use during debugging
  data_list = list(X1=matrix(rnorm(100),ncol=20,nrow=5), X2=matrix(rnorm(200),ncol=20,nrow=10))
  r = 2
  rIndiv = c(1,1)
  SpSm = c("Sp","Sm")
  lambda = 'BIC'
  alpha = 'BIC'
  scale = 'y'
  ConvergeThresh = 1e-6
  MaxIter = 1e3
}
  
library(ssvd)
library(Rcpp)
library(RcppArmadillo)
  
if(!is.list(data_list)){
  print('Error: data_list should be of type list');
  print('For data_list matrices X1,X2,...,Xk, input data_list = list(X1,X2,...,Xk)');
}
ndata = length(data_list);

n <- d <- rep(NA,ndata)
for (i in 1:ndata){
  dim_tmp <- dim(data_list[[i]])
  d[i] <- dim_tmp[1]; n[i] <- dim_tmp[2];
}
n = unique(n);
if(length(n) != 1){
  print('Error: elements in data_list do not have same number of columns');
}

if(scale =='y'){
  for(i in 1:ndata){
    data_list[[i]] = scale(data_list[[i]], center=T, scale=F);
    data_list[[i]] = data_list[[i]]/norm(data_list[[i]], type='F');
  }
}

if(lambda == 'BIC'){
  lambda = (0.5)^(1:10)[c(FALSE, TRUE)]
} else{ lam_fix=T }
if(alpha == 'BIC'){
  alpha = (0.5)^(1:10)[c(FALSE, TRUE)]
} else{ alp_fix=T }
lambda_all = rep(list(lambda),ndata)
alpha_all = rep(list(alpha),ndata)

# Make row-bound concatenation of all data_list matrices
Tot = do.call(rbind, data_list) 
# Create empty matrices to hold joint (J) and individual (A) matrix components
J = matrix(0,nrow=sum(d),ncol=n);
A = matrix(0,nrow=sum(d),ncol=n);
X_est = 0;
# Get Omega matrix for where sparsity/smoothness should be imposed
# and specify sparsity_locs_all vector for vars where sparsity should be imposed.
Omeg_all = matrix(0,nrow=nrow(Tot),ncol=nrow(Tot));
sparsity_locs_all = rep(0,nrow(Tot))
for(i in 1:ndata){
  if(i!=1){ tmp <- sum(d[1:(i-1)]) }else{ tmp <- 0 }
  rows = (tmp+1):sum(d[1:i]);
  p_tmp = length(rows)
  if( (p_tmp>2) & (SpSm[i]%in%c("Sm","SpSm")) ){
    Omeg_all[rows, rows] <- getOmega(p_tmp)
  }
  if( SpSm[i]%in%c("Sp","SpSm") ){
    sparsity_locs_all[rows] <- rep(1,p_tmp)
  }
}

# Initialize empty matrices to hold results
UV1 = US = UV2 = list()

# Run JIVE algorithm 
mod_num <- 20 # point at which you stop allowing modifications
for(j in 1:MaxIter){
  print(paste("Iter",j))
  # Note that sfpca_bic takes in n x p matrix (i.e., obs in rows)
  if( sum(SpSm=="N")==length(SpSm) ){ # All regular (no sparsity or smoothness)
    tmp = svd(t(Tot), nu=r, nv=r)
    res = list(U=tmp$u, D=tmp$d[1:r], V=tmp$v)
  } else if(sum(SpSm=="Sp")==length(SpSm)){ # All sparse (no smoothness)
    # https://cran.r-project.org/web/packages/ssvd/ssvd.pdf
    tmp = ssvd(t(Tot), method="theory", r=r)
    res = list(U=tmp$u, D=tmp$d, V=tmp$v)
  }else{
    res = sfpca_bic(Tot, r, Omeg_all, sparsity_locs_all, lambda, alpha)
    if( min(c(res[["optavs"]]))==min(alpha) & j<mod_num & (!alp_fix) ){ alpha = c(alpha[-1], min(alpha)/2); print("Decreasing alp!") }
    if( max(c(res[["optavs"]]))==max(alpha) & j<mod_num & (!alp_fix) ){ alpha = c(max(alpha)*1.2, head(alpha, -1)); print("Increasing alp!") }
    if( min(c(res[["optlvs"]]))==min(lambda) & j<mod_num & (!lam_fix) ){ lambda = c(lambda[-1], min(lambda)/2); print("Decreasing lam!") }
    if( max(c(res[["optlvs"]]))==max(lambda) & j<mod_num & (!lam_fix) ){ lambda = c(max(lambda)*1.2, head(lambda, -1)); print("Increasing lam!") }
  }
  V1 = as.matrix(res[['V']]); V2 = as.matrix(res[['U']])
  D_tmp = c(res[['D']]); S = diag(length(D_tmp)); diag(S) = D_tmp
  J = V1 %*% S %*% t(V2);
  for(i in 1:ndata){
    if(i!=1){tmp <- sum(d[1:(i-1)])}else{tmp <- 0}
    rows = (tmp+1):sum(d[1:i]);
    U = data_list[[i]] - J[rows,];
    if( !(rIndiv[i]==0) ){
      if( SpSm_indiv[i]=="Sm" | SpSm_indiv[i]=="SpSm" ){
        res = sfpca_bic(U - U %*% V2 %*% t(V2), rIndiv[i], Omeg_all[rows,rows], 
                        sparsity_locs_all[rows], lambda, alpha)
        if( min(c(res[["optavs"]]))==min(alpha_all[[i]]) & j<mod_num & (!alp_fix) ){ alpha_all[[i]] = c(alpha_all[[i]][-1], min(alpha_all[[i]])/2); }
        if( max(c(res[["optavs"]]))==max(alpha_all[[i]]) & j<mod_num & (!alp_fix) ){ alpha_all[[i]] = c(max(alpha_all[[i]])*1.2, head(alpha_all[[i]], -1)); }
        if( min(c(res[["optlvs"]]))==min(lambda_all[[i]]) & j<mod_num & (!lam_fix) ){ lambda_all[[i]] = c(lambda_all[[i]][-1], min(lambda_all[[i]])/2); }
        if( max(c(res[["optlvs"]]))==max(lambda_all[[i]]) & j<mod_num & (!lam_fix) ){ lambda_all[[i]] = c(max(lambda_all[[i]])*1.2, head(lambda_all[[i]], -1)); }
      } else if(SpSm_indiv[i]=="Sp"){
        # https://cran.r-project.org/web/packages/ssvd/ssvd.pdf
        tmp = ssvd(t(U - U %*% V2 %*% t(V2)), method="theory", r=rIndiv[i])
        res = list(U=tmp$u, D=tmp$d, V=tmp$v)
      } else{ # Only type left is "None"
        tmp = svd(t(U - U %*% V2 %*% t(V2)), nu=rIndiv[i], nv=rIndiv[i])
        res = list(U=tmp$u, D=tmp$d[1:rIndiv[i]], V=tmp$v)
      }
      if( length(c(res[['D']]))==1 ){
        US[[i]] = as.matrix(c(res[['D']]))
      } else{
        US[[i]] = as.matrix(diag(c(res[['D']])))
      }
      UV1[[i]] = as.matrix(res[['V']]); UV2[[i]] = as.matrix(res[['U']])
    } else{ # Just fill in zero-matrices of correct dimension.
      UV1[[i]] = matrix(0, nrow=d[i], ncol=1); US[[i]] = matrix(0, nrow=1, ncol=1); UV2[[i]] = matrix(0, nrow=n, ncol=1)
    }
    A[rows,] = UV1[[i]] %*% US[[i]] %*% t(UV2[[i]]);
    Tot[rows,] = data_list[[i]] - A[rows,];
  }
  if(norm(X_est-J-A,'F')^2 < ConvergeThresh){
    break
  }
  print(paste("F norm diff:", norm(X_est-J-A,'F')^2))
  X_est = J+A;
  if(j==MaxIter){
    print('Warning: MaxIter iterations reached before convergence');
  }
}

res = list("UV1"=UV1, "US"=US, "UV2"=UV2, "V1"=V1, "S"=S, "V2"=V2)
#save(res, file="/Users/Kelly/jive_res.Rdata")
return(res)

}
