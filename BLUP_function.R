#' This is the script for BLUP function. It relies on the Conjugate Gradient method in 
#' file: "Conjugate_Gradient_Method.R"
#' 
#' version 1.05, on Oct 21, 2017. 

setwd("~/Documents/GitHub/Conjugate_Gradient_Method/")
source(file = "Conjugate_Gradient_Method.R")


BLUP <- function(geno=NULL, pheno = NULL, Vg = 0, Ve = 1, sig_thres = 1e-5, iter_num=50){
    if(is.null(geno) | is.null(pheno)){
        stop("please input the correct genotype/phenotype data frame.")
    }else if(nrow(geno) != length(pheno)){
        stop("the dimensions of the genotype and the phenotype data do not match.")
    }
    
    ## construct the A matrix
    lambda <- Ve/Vg*ncol(geno)
    geno <- as.matrix(geno)
    A <- t(geno)%*%geno + lambda*diag(x=1, nrow=ncol(geno))
    
    ## construct the b vector as in Ax = b
    b <- t(geno) %*% pheno
    
    x0 <- rep(0, ncol(geno))
    
    blup_vec <- ConjugateGradient(A, b, x0, sig_thres = sig_thres, iter_num = iter_num)
    blup_vec
}
