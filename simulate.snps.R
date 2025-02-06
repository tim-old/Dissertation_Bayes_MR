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
## sample allele frequencies for the SNPs in the range [0.01, 0.99]
af.snps <- runif(num.snps, min = 0.01, max = 0.99)

## sample genotypes for each individual and each SNP
## sample probability of having the minor allele uniformly in [0, 1], then
## threshold according to the allele frequency for that SNP to get the allele,
## then add up the two alleles for each SNP and each individual to get the
## genotype
allele.probs <- matrix(runif(2* n * num.snps, min = 0, max = 1), nrow = 2 * n, ncol = num.snps)
alleles <- allele.probs < matrix(rep(af.snps, 2 * n), nrow = 2 * n, ncol = num.snps, byrow = TRUE)
genos <- alleles[1:n, ] + alleles[(n+1):(2*n), ]

## check the empirical MAF is roughly equal to the sampled value
cbind(colSums(genos) / (2*n), af.snps)


## -- generate X --
## sample genetic effects on the exposure
gamma <- runif(num.snps, min = 0.03, max = 0.1)
## sample random noise for the exposure value for each individual
epsilon_x <- rnorm(n)
## generate X
X <- genos %*% gamma + epsilon_x


## -- generate Y (scenario 2 with no causal effect) --
beta_x <- 0 # null (no causal effect)
alpha <- runif(num.snps, min = -0.2, max = 0.2)
epsilon_y <- rnorm(n)

Y <- genos %*% alpha + beta_x * X + epsilon_y


## -- fit linear model to estimate the effects of the genetic instruments on the
##    exposure --
x.lm <- lm(X ~ genos)
betaXG <- coef(summary(x.lm))[2:(num.snps + 1), 1]
sebetaXG <- coef(summary(x.lm))[2:(num.snps + 1), 2]

## sanity check: betaXG should approximate gamma as we increase the sample size
plot(gamma, betaXG)


## -- fit linear model to estimate the effects of the genetic instruments on the
##    outcome --
y.lm <- lm(Y ~ genos)
betaYG <- coef(summary(y.lm))[2:(num.snps + 1), 1]
sebetaYG <- coef(summary(y.lm))[2:(num.snps + 1), 2]

## sanity check: in this scenario (no causal effect) betaYG should approximate
## alpha as we increase the sample size
plot(alpha, betaYG)


