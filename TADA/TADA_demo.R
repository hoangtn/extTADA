source("TADA.R")

#################################################################
# Application of TADA
#################################################################

# Model parameters: two categories of mutations - LoF and mis3 mutations ("probably damaging" by PolyPhen2)
mu.frac <- c(0.074, 0.32)
gamma.mean.dn <- c(20, 4.7)
beta.dn <- c(1,1)
gamma.mean.CC <- c(2.3, 1.00)
beta.CC <- c(4.0, 1000)
rho1 <- c(0.1, 0.5)
nu1 <- c(200, 100)
rho0 <- c(0.1, 0.5)
nu0 <- c(200, 100)
hyperpar <- as.array(rbind(gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0))
l <- 100
pi0 <- 0.94 # the fraction of non-risk genes

# ASC (Autism Sequencing Consortium) data
# The file name contains the sample size information
# The only relevant counts are dn.LoF and dn.mis3
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
ntrio <- 2231  # number of trios
ncase <- 1601  # number of cases
nctrl <- 5397  # number of controls
ntrans <- 1333 # number of subjects with transmission data
N <- list(dn=ntrio, ca=ntrans+ncase, cn=ntrans+nctrl)

# Running TADA
counts <- as.array(cbind(data$dn.LoF, data$case.LoF+data$trans.LoF, data$ctrl.LoF+data$ntrans.LoF, data$dn.mis3, data$case.mis3+data$trans.mis3, data$ctrl.mis3+data$ntrans.mis3))
rs <- TADA(counts, N, data$mut.rate, mu.frac, hyperpar)
data$BF <- rs$BF.total

# Estimating p-values of BFs (this is optional and slow)
rsp <- TADAp(counts, N, data$mut.rate, mu.frac, hyperpar, l=100)
data$pval.TADA <- rsp$pval

# FDR estimation
data <- data[order(-data$BF),]
data$qvalue <- Bayesian.FDR(data$BF, pi0)$FDR
write.csv(data, "data/TADA_results.csv", row.names=FALSE)

#################################################################
# Estimation of parameters of the prior distributions
#################################################################

# For estimation of RR of de novo LoF and mis3, see TADA_denovo_demo.R

# ASD data
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
ntrio <- 2231  
ncase <- 1601  
nctrl <- 5397  
ntrans <- 1333

# Estimation of RR of inherited variants using a set of known disease genes. The RR is approximately the enrichment of the LoF and mis3 counts in cases vs. in controls
ASDgenes.table <- read.csv("data/known_ASD_genes.csv", header=TRUE, as.is=TRUE)
ASDgenes <- data.frame(Gene=subset(ASDgenes.table, label==1)$gene)
data.ASDgenes <- merge(data, ASDgenes, by="Gene")
# Case-control
q.LoF.case <- sum(data.ASDgenes$case.LoF)/ncase
q.LoF.ctrl <- sum(data.ASDgenes$ctrl.LoF)/nctrl
q.LoF.case/q.LoF.ctrl
q.mis3.case <- sum(data.ASDgenes$case.mis3)/ncase
q.mis3.ctrl <- sum(data.ASDgenes$ctrl.mis3)/nctrl
q.mis3.case/q.mis3.ctrl

# Estimating the prior parameters of q: q ~ Gamma(q.mean*nu, nu). 
# Note: we assume that q|H_1 and q|H_0 have the same prior. 
# Only estimate the prior mean of q (q.mean), and choose nu to be a small number (e.g. 200)
q.LoF.mean <- mean(data$ctrl.LoF)/nctrl
q.mis3.mean <- mean(data$ctrl.mis3)/nctrl


#################################################################
# Simulation to assess the power of TADA
#################################################################

# To generate the simulation data, we need the prior parameters of q, the allele frequency. We estimate this using the negative binomial distribution (see the Section "Simulation to assess the power of TADA" in the Supplement of our paper)
# The prior distribution: q ~ Gamma(rho, nu) (assuming that q|H1 and q|H0 follow the same distribution, i.e. rho1=rho0=rho, nu1=nu0=nu)
library(MASS)
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
nctrl <- 5397
x <- data$ctrl.LoF
nbpar <- as.vector(fitdistr(x, "negative binomial", list(size=0.5, mu=0.5))$estimate)
mu <- nbpar[2]
theta <- nbpar[1]
rho <- theta
nu <- nctrl * theta/ mu

# Mutation rates (from real data) and parameters
# Note that we need two sets of parameters: one used for simulation, the other used by TADA for scoring genes (those with the suffix "est")
data <- read.csv("data/ASC_2231trios_1333trans_1601cases_5397controls.csv", header=TRUE, as.is=TRUE)
mu <- data$mut.rate
mu.frac <- c(0.074, 0.32)
pi <- 0.06
gamma.mean.dn <- c(18, 5.4)
beta.dn <- c(1, 0.5)
gamma.mean.CC <- c(2.3, 1.0)
beta.CC <- c(4.0, 1000)
rho1 <- c(0.66, 0.6)
nu1 <- c(1947, 123)
rho0 <- c(0.66, 0.64)
nu0 <- c(1947, 123)
rho1.est <- c(0.1, 0.5)
nu1.est <- c(200, 100)
rho0.est <- c(0.1, 0.5)
nu0.est <- c(200, 100)
hyperpar.est <- as.array(rbind(gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1.est, nu1.est, rho0.est, nu0.est))

# Use simulation to assess the power of TADA
nr <- 10
rs <- numeric(nr)
N <- 1000*1:5
power.mean <- numeric(length(N))
power.sd <- numeric(length(N))
for (i in 1:length(N)) {
  N.curr <- list(dn=N[i], ca=N[i], cn=N[i])
  rs <- replicate(nr, eval.TADA(N.curr, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, hyperpar.est, FDR=0.1, tradeoff=TRUE)$M1)
  power.mean[i] <- mean(rs)
  power.sd[i] <- sd(rs)
}
power <- data.frame(mean=power.mean, sd=power.sd)

