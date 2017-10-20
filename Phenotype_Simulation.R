#' This is a test script for BLUP implementation of the CG method. 
#' The genotype data are from "Mouse HapMap Imputed Genotype Resource"
#'  URL: http://mouse.cs.ucla.edu/mousehapmap/emma.html
#' We will use the Chr19 genotype (~3k SNPs and 251 strains) and use it to
#'  simulate an additive phenotype. We will then run BLUP using CG on 200 strains (training set)
#'  and evaluate the prediction power using the BLUP predictor. 
#'  
#' by Longda Jiang (longda.jiang@uq.edu.au), on Oct 20 2017 .
#' 

setwd("~/Documents/GitHub/Conjugate_Gradient_Method")
raw_geno <- read.table(file="chr_19.snps", head=T, stringsAsFactors = F, na.strings = "N")  # the ACTG genotype with SNPid, strain names, etc
geno2 <- read.table(file="chr_19.emma", head=F, stringsAsFactors = F)   # the 0, 1, 2 format genotype 

## the raw genotype contains 3247 SNPs and 251 indi

h2 <- 0.6
nSNPs <- 50


## first we need to do some cleaning job for the raw genotype.

# 1. recode the genotype into 0, 1, 2
names(geno2) <- names(raw_geno)[-(1:4)]
row.names(geno2) <- raw_geno$SNP_ID

# 2. delete those SNPs with high missingness and strains with high missingness
cal_miss <- function(x){
  sum(is.na(x))/length(x)
}
# individual
idx <- which(apply(X = geno2, MARGIN = 2, FUN = cal_miss) < 0.1)   # 217 indi left
geno2_clean <- geno2[, idx]

# SNPs
idx <- which(apply(X = geno2_clean, MARGIN = 1, FUN = cal_miss) < 0.05)  # 2846 SNPs left
geno2_clean <- geno2_clean[idx,]


# 3. select the causal variants
geno2_clean <- t(geno2_clean)
geno2_clean <- as.data.frame(geno2_clean)

# 50 causal SNPs
set.seed(50)
causal_idx <- sample(1:2846, size=nSNPs)
causal_idx <- sort(causal_idx)

# we need to exclude those causal SNPs with low MAF since they are less representative
cal_maf <- function(x){
  sum(x, na.rm=T)/length(x)
}

eval_vector <- apply(X = geno2_clean[, causal_idx], MARGIN = 2, FUN = cal_maf)
idx <- which(eval_vector < 0.01 | eval_vector >0.99)   # 6 SNPs with low MAF need to be excluded
causal_idx <- causal_idx[-idx]

set.seed(60)
causal_idx_add <- sample(1:2846, size= nSNPs - lenght(idx))
causal_idx <- sort(union(causal_idx_add, causal_idx))

eval_vector <- apply(X = geno2_clean[, causal_idx], MARGIN = 2, FUN = cal_maf)
idx <- which(eval_vector < 0.01 | eval_vector >0.99)   # None

# thus "causal_idx" will be used as our calsual variants


# 4. simulate phenotype,  we will set the h2 to 0.6 with 50 causal snps.
set.seed(70)
effect <- rnorm(n = length(causal_idx))  # effect size for 50 SNPs

causal_geno <- geno2_clean[, causal_idx]
causal_geno <- as.matrix(causal_geno)
causal_geno[is.na(causal_geno)] <- 0

g <- causal_geno %*% effect    # the total genetic value for each individual
var_e <- var(g)/h2 - var(g)

set.seed(80)
env <- rnorm(n = nrow(geno2_clean), sd=sqrt(var_e))

pheno <- scale(g +  env)


# 5. we will now save this simulated pheno file and the clean genotype file. we will do the BLUP in my next script.

pheno <- cbind(row.names(geno2_clean), pheno, g)
pheno <- as.data.frame(pheno)
names(pheno) <- c("FID", "pheno", "genetic_value")
pheno$FID <- as.character(pheno$FID)
pheno$pheno <- as.numeric(as.character(pheno$pheno))
pheno$genetic_value <- as.numeric(as.character(pheno$genetic_value))

write.table(pheno, file="simulated_phenotype.pheno", quote=F, row.names = F, sep="\t")

write.table(geno2_clean, file="clean_genotype.geno", quote=F, row.names = F, sep="\t")





