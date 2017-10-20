#' This is a script to test the performance of CG-BLUP
#' 
#' by Longda Jiang, on Oct 20, 2017
#' 

setwd("~/Documents/GitHub/Conjugate_Gradient_Method/")
source(file = "Conjugate_Gradient_Method.R")

# load the geno and pheno data
geno <- read.table(file = "clean_geno_pheno/clean_genotype.geno", head=T, stringsAsFactors = F)
pheno <- read.table(file="clean_geno_pheno/simulated_phenotype.pheno", head=T, stringsAsFactors = F)


# Actual calculation

## !!!    One thing to be noticed is that the mouse data contians only pure strains, whose genotype is always coded as 0/1, 
##        where 1 is for 2 minor alleles and 0 is for 2 minor alleles. 
##        Thus in the standardisation step, we will use 0-psqrt(2p(1-p))

## standardise the genotype:
scale_geno <- function(x){
  p <- mean(x, na.rm=T)/2
  
}

geno <- 







