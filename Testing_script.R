#' This is a script to test the performance of CG-BLUP
#' 
#' by Longda Jiang, on Oct 20, 2017
#' 

setwd("~/Documents/GitHub/Conjugate_Gradient_Method/")
source(file = "Conjugate_Gradient_Method.R")
source(file = "REML_for_Blup.R")
source(file = "BLUP_function.R")

# load the geno and pheno data
geno <- read.table(file = "clean_geno_pheno/clean_genotype.geno", head=T, stringsAsFactors = F)
pheno <- read.table(file="clean_geno_pheno/simulated_phenotype.pheno", head=T, stringsAsFactors = F)

idx <- which(apply(geno, 2, sd, na.rm=T) == 0 )
geno <- geno[, -idx]

geno <- as.matrix(geno)

geno[!is.na(geno) & geno == 1] <- 2

# Calculation

## !!!    One thing to be noticed is that the mouse data contians only pure strains, whose genotype is always coded as 0/1, 
##        where 1 is for 2 minor alleles and 0 is for 2 minor alleles. 
##        Thus in the standardisation step, we will use 0-p/sqrt(p(1-p))

## standardise the genotype:

for( i in 1:ncol(geno)){
  p <- mean(geno[,i], na.rm=T)/2
  geno[is.na(geno[,i]), i] <- p
  
  geno[,i] <- (geno[,i]-p)/sd(geno[,i])
}


## construct the genetic correlation matrix:
GRM <- geno %*% t(geno)/ncol(geno)


## calculate the Vg and Vp
reml_result <- reml(grm = GRM, y=pheno$pheno, sigma_g = 0.99, sigma_e = 0.01, iter_num = 200)

# Vg = 0.4201491, Ve = 0.3969231


# BLUP: general formula: (XXt + Im*lambda)-1 * Xty = beta

Vg <- reml_result$Vcomp[1]
Ve <- reml_result$Vcomp[2]

beta_hat <- BLUP(geno = geno, pheno=pheno$pheno, Vg = Vg, Ve = Ve)  # the iteration number is 19,
system.time(BLUP(geno = geno, pheno=pheno$pheno, Vg = Vg, Ve = Ve))  
## >   user  system elapsed 
## >   15.085   0.154  15.257 


estimate <- geno %*% beta_hat

cor(estimate, pheno$pheno)
## The correlation r is 0.8256544, looks pretty good. 


######################################################################
##       But we should use corss validation to evaluate the performance.
######################################################################

# we will do a 5-k validation, with ~40 observations in each subset. 

# first, create 5 list of random idx
set.seed(10)
random_idx <- sample(nrow(pheno))
idx_idx <- cut(seq(1, nrow(pheno)), breaks=5, labels = F)

# Example:
# the training set
random_idx[which(idx_idx != 1)]
# test set
random_idx[which(idx_idx == 1)]


# cross-validation:

r_hat_vec <- vector()

for(i in 1:5){
    
    training_idx <- random_idx[which(idx_idx != i)]
    test_idx <- random_idx[which(idx_idx == i)]
    
    pheno_train <- pheno[training_idx, ]
    geno_train <- geno[training_idx, ]
    
    pheno_test <- pheno[test_idx, ]
    geno_test <- geno[test_idx, ]
    
    
    # REML
    GRM <- geno_train %*% t(geno_train)/ncol(geno_train)
    reml_result <- reml(grm = GRM, y=pheno_train$pheno, sigma_g = 0.95, sigma_e = 0.05, iter_num = 100)
    Vg <- reml_result$Vcomp[1]
    Ve <- reml_result$Vcomp[2]
    
    # BLUP
    beta_hat <- BLUP(geno = geno_train, pheno=pheno_train$pheno, Vg = Vg, Ve = Ve)  # the iteration number is 19,
    estimate <- geno_test %*% beta_hat
    r_hat <- cor( estimate, pheno_test$pheno) 
    
    # save the correlation r
    r_hat_vec <- c(r_hat_vec, r_hat)
}


r_hat_vec

# >  0.7039108 0.6393375 0.7450056 0.6928418 0.7785923  r_squ is ~0.5, looks legit.










