source("TADA.R")

#################################################################
# Application of TADA-Denovo
#################################################################

# Model parameters: two categories of mutations - LoF and mis3 mutations ("probably damaging" by PolyPhen2)
mu.frac <- c(0.074, 0.32)
gamma.mean.dn <- c(20, 4.7)
beta.dn <- c(1,1)
l <- 100
pi0 <- 0.94 # the fraction of non-risk genes

# ASC (Autism Sequencing Consortium) data
# The file name contains the sample size information
# The only relevant counts are dn.LoF and dn.mis3
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
ntrio <- 2231  # number of trios  

# Running TADA-denovo
counts.dn <- as.array(cbind(data$dn.LoF, data$dn.mis3))
rs.dn <- TADA.denovo(counts.dn, ntrio, data$mut.rate, mu.frac, gamma.mean.dn, beta.dn)
data$BF.dn<- rs.dn$BF.total

# Estimating p-values of BFs (optional and slow)
rsp.dn <- TADAp.denovo(counts.dn, ntrio, data$mut.rate, mu.frac, gamma.mean.dn, beta.dn, l=l)
data$pval.TADA.dn <- rsp.dn$pval

# FDR estimation
data <- data[order(-data$BF.dn),]
data$qvalue.dn <- Bayesian.FDR(data$BF.dn, pi0)$FDR
write.csv(data, "data/TADA_denovo_results.csv", row.names=FALSE)

#################################################################
# Estimation of de novo parameters using Method of Moment approach
#################################################################

# ASC data 
# We start with LoF mutations 
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
ntrio <- 2231  # number of trios
m <- nrow(data)
mu <- (data$mut.rate) * 0.074
C <- sum(data$dn.LoF) # total number of de novo LoF mutations
M.obs <- sum(data$dn.LoF > 1) # number of multiple-hit genes

# MOM estimation: choose the value of k (the number of risk genes) so that the expected number of multiple-hit genes is closest to the observed number 
# The results are: #risk genes and the average RR
# Note: the value of beta usually has a small influence, so we fix its value. But one could try different values, e.g. from 0.1 to 10.
k.vec <- 100*6:14
beta <- 1.0   
M <- numeric(length(beta.vec))
for (k.idx in 1:length(k.vec)) {
  k <- k.vec[k.idx] 
  M[k.idx] <- denovo.MOM(ntrio, mu, C, beta, k)$M
}  
idx.k <- which.min(abs(M-M.obs))
k.est <- k.vec[idx.k]  # estimated number of risk genes
gamma.mean.est <- denovo.MOM(ntrio, mu, C, beta, k.est)$gamma.mean # estimated RR of de novo LoF mutations

# RR for mis3 mutations
mu.mis3 <- (data$mut.rate) * 0.32
C.mis3 <- sum(data$dn.mis3)
beta.mis3 <- 1.0
gamma.mean.mis3.est <- denovo.MOM(ntrio, mu.mis3, C.mis3, beta.mis3, k.est)$gamma.mean

#################################################################
# Simulation to assess the power of TADA.denovo
#################################################################

# Mutation rates (from real data) and parameters
# Note that we need two sets of parameters: one used for simulation, the other used by TADA for scoring genes (those with the suffix "est")
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
mu <- data$mut.rate
mu.frac <- c(0.074, 0.32)
pi <- 0.06
gamma.mean <- c(18, 5.4)
beta <- c(1, 0.5)
gamma.mean.est <- c(18, 5.4)
beta.est <- c(1, 0.5)

# Use simulation to assess the power of TADA.denovo
nr <- 10
rs <- numeric(nr)
N <- 1000*1:5
power.mean <- numeric(length(N))
power.sd <- numeric(length(N))
for (i in 1:length(N)) {
  rs <- replicate(nr, eval.TADA.denovo(N[i], mu, mu.frac, pi, gamma.mean, beta, gamma.mean.est, beta.est, FDR=0.1))
  power.mean[i] <- mean(rs)
  power.sd[i] <- sd(rs)
}
power.dn <- data.frame(mean=power.mean, sd=power.sd)
