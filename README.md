# About ssJIVE

Joint and Individual Variation Explained (JIVE) is an exploratory method that separates joint and individual effects from multiple data sets. See paper by [Lock et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23745156) for details. Sparse and smooth JIVE (ssJIVE) allows for sparsity and/or smoothness to be imposed on the columns of V<sub>X</sub> and/or V<sub>Y</sub>. 

Suppose `X` and `Y` are data sets having `px` and `py` rows, respectively, and both having `n` columns. Let `dat_list=list(X,Y)`. Then running `ssJIVE(data_list=dat_list, r=2, rIndiv=c(4,3))` will run JIVE (no sparsity or smoothness imposed) assuming the shared latent space is 2-dimensional, and `X` and `Y` each have individual (unshared) latent spaces of 4 and 3 dimensions, respectively. Running `ssJIVE(data_list=dat_list, r=2, rIndiv=c(4,3), SpSm=c('Sp','Sm'))` will keep the latent dimensions but impose sparsity on the columns of V<sub>X</sub> and smoothness on the columns of V<sub>Y</sub>.

# Installing the package

After downloading and installing **R** (if necessary) and Rtools (if installing on a Windows machine) run the following from within **R**:
```
# Install the ssJIVE package
remotes::install_github("kelrenmor/ssJIVE")
```

