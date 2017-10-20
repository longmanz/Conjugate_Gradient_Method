# Conjugate Gradient Method
This is a R script to try out the conjugate gradient (CG) method, and implemented it into genomic BLUP scheme.

## A short description of the CG method:
- To solve linear equation Ax = b often requires solving the inverse of A, which is usually computationally intensive.
- Yet, for a special kind of Ax = b, where A is a symmetric, positive-definite matrix, we can solve this equation by doing simple Matrix X vector iteration.  IT IS SUPER FAST!  
- Of all thid kind of methods, conjugate gradient method often requires the least iteration time. 
- For a thorough description of this method, please refer to   Shewchuk, Jonathan Richard. "An introduction to the conjugate gradient method without the agonizing pain." (1994).


## Implementation of CG into BLUP

### Genotype data
For evaluation, we downloaded a public available data from Mouse HapMap project (URL: http://mouse.cs.ucla.edu/mousehapmap/emma.html). Only Chr19 genotype was downloaded, with ~ 3000 SNPs and 251 strains. After QC (MAF <= 0.01, missingness of indi > 0.1 and missingness of SNPs > 0.05), we obtained 2864 SNPs and 217 strains. 


### Phenotype data
We used the above clean genotype data to generate a simulated additive phenotype, with h2 = 0.6 and num of causal variants/SNPs = 50. 


### todo



-P.S. The reason why I play around it is that BOLT-LMM applied this method in its MCMC REML algorithm, and it only requires ~5 iteration to solve the Mixed model equation  V-1 * y = x. 

