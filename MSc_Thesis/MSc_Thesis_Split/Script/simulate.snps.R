## Simulate independent SNPs to use as instruments in MR simulations
##
## author: athina spiliopoulou
## date: 6 January 2025
##==============================================================================

## set seed for random number generator (get same values when re-running)
set.seed(678)

## number of samples
n <- 10000

## number of SNPs (instruments) -- these are assumed to be independent
num.snps <- 25

## -- sample genotypes for the instruments --
## sample minor allele frequencies for the SNPs in the range [0.01, 0.5]
maf.snps <- runif(num.snps, min = 0.01, max = 0.5)

## sample genotypes for each individual and each SNP
## sample probability of having the minor allele uniformly in [0, 1], then
## threshold according to the MAF for that SNP to get the allele, then add up
## the two alleles for each SNP and each individual to get the genotype
allele.probs <- matrix(runif(2* n * num.snps, min = 0, max = 1), nrow = 2 * n, ncol = num.snps)
alleles <- allele.probs < matrix(rep(maf.snps, 2 * n), nrow = 2 * n, ncol = num.snps, byrow = TRUE)
genos <- alleles[1:n, ] + alleles[(n+1):(2*n), ]

## check the empirical MAF is roughly equal to the sampled value
cbind(colSums(genos) / (2*n), maf.snps)
