#################################################################
# Some useful functions
#################################################################

# Fisher's method of combining p values
# Input: a vector of p-values
combpval <- function(p) {
  k <- length(p)
  T <- -2 * sum(log(p))
  return (1 - pchisq(T,2*k))
}

# Genomic control
# Input: a vector of p-values, the quantile (0.5 or 0.75)
genom.ctrl <- function(p, quant) {
  chisq.obs <- qchisq(1-p, 1) # convert p-values to chi-square statistics (dof 1)
  chisq.obs.quant <- quantile(chisq.obs, probs=quant, names=FALSE) # 0.75 if first quantile
  chisq.exp.quant <- qchisq(quant, 1)
  return (chisq.obs.quant / chisq.exp.quant)
}

# Bayesian FDR control (PMID:19822692, Section2.3)
# BF: a sorted vector of BFs (in decreasing order)
# pi0: the prior probability that the null model is true
# alpha: the FDR target
# Return: the q-value of each BF, and the number of findings with q below alpha. 
Bayesian.FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}

# draw QQ plot
# Input: a vector of p-values
plotQQ <- function(p.obs) {
  obs <- -log10(p.obs)
  theo <- -log10(ppoints(length(obs)))
  qqplot(theo, obs, xlab=expression("Theoretical " * -log[10](p)), ylab=expression("Observed "*-log[10](p)))
  abline(0,1,col='red')
}

# Similar, but show QQ plot in the original scale (not log. scale)
plotQQ.unif <- function(p.obs) {
  obs <- (p.obs)
  theo <- (ppoints(length(obs)))
  qqplot(theo, obs, xlab="Theoretical p-values", ylab="Observed p-values")
  abline(0,1,col='red')
}

#################################################################
# Bayes Factor Computation for a Single Gene
#################################################################

# model evidence of de novo data: P(x_d|H_0) 
# Input: the count data x, the sample size N, the mutation rate mu 
evidence.null.dn <- function(x, N, mu) {
  return (dpois(x, 2*N*mu))
}

# model evidence of de novo data: P(x_d|H_1) 
# Input: the count data x, the sample size N, the mutation rate mu and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
evidence.alt.dn <- function(x, N, mu, gamma.mean, beta) {
  return (dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu)))
}


# model evidence of case-control data: P(x_1,x_0|H_0) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of q|H0: Gamma(rho0, nu0)
evidence.null.CC <- function(x, N, rho0, nu0) {
  marglik0.ctrl.log <- log(dnbinom(x$cn, rho0, nu0/(nu0+N$cn)))
  marglik0.case.log <- log(dnbinom(x$ca, rho0+x$cn, (nu0+N$cn)/(nu0+N$cn+N$ca)))
  marglik0.log <- marglik0.ctrl.log + marglik0.case.log
  
  return (list(ctrl=exp(marglik0.ctrl.log), case=exp(marglik0.case.log), total=exp(marglik0.log)))
}

# model evidence of case-control data: P(x_1,x_0|H_1) 
# Input: the count data x, the sample size N and the parameters
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
evidence.alt.CC <- function(x, N, gamma.mean, beta, rho1, nu1, q.lower=1e-8, q.upper=0.1, debug=FALSE) {
  integrand <- function(u) {
    q <- exp(u)
    return (dnbinom(x$ca, gamma.mean*beta, beta/(beta+N$ca*q)) * dgamma(q, rho1+x$cn, nu1+N$cn) * exp(u))
  }
  
  marglik1.ctrl <- dnbinom(x$cn, rho1, nu1/(nu1+N$cn))
  marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value
  
  marglik1 <- marglik1.ctrl * marglik1.case
  
  #   return (exp(marglik1.ctrl.log+marglik1.case.log))
  return (list(ctrl=marglik1.ctrl, case=marglik1.case, total=marglik1))
}

# Bayes factor of the case-control data
# BF.cn and BF.ca: contribution from control and case data, respectively
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor.CC <- function(x, N, gamma.mean, beta, rho1, nu1, rho0, nu0) {
  marglik0.CC <- evidence.null.CC(x, N, rho0, nu0)
  marglik1.CC <- evidence.alt.CC(x, N, gamma.mean, beta, rho1, nu1)
  
  BF.cn <- marglik1.CC$ctrl / marglik0.CC$ctrl
  BF.ca <- marglik1.CC$case / marglik0.CC$case
  BF <- BF.cn * BF.ca
  
  return (list(BF=BF, BF.cn=BF.cn, BF.ca=BF.ca))
}

# Bayes factor of de novo counts of a gene 
# x: the de novo count
# N: the sample size (number of families)
# mu: the mutation rate (of this type of mutational events)
# Prior distribution of RR: gamma ~ Gamma(gamma.mean*beta, beta)
bayes.factor.denovo <- function(x, N, mu, gamma.mean, beta) {
  marg.lik0 <- dpois(x, 2*N*mu)
  marg.lik1 <- dnbinom(x, gamma.mean*beta, beta/(beta+2*N*mu))
  BF <- marg.lik1/marg.lik0
  
  return (BF)
}

# Bayes factor of the gene combining de novo and case-control
# x: a list of (dn, ca, cn), counts in de novo, cases and controls
# N: a list of (dn, ca, cn), sample sizes
# hyperpar: (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)
# Prior distribution of RR in de novo: gamma.dn ~ Gamma(gamma.mean.dn*beta.dn, beta.dn)
# Prior distribution of RR in C/C data: gamma.CC ~ Gamma(gamma.mean.CC*beta.CC, beta.CC)
# Prior distribution of q|H1: Gamma(rho1, nu1)
# Prior distribution of q|H0: Gamma(rho0, nu0)
bayes.factor <- function(x, N, mu, hyperpar, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  BF.dn <- bayes.factor.denovo(x$dn, N$dn, mu, gamma.mean.dn, beta.dn)
  x.CC <- list(ca=x$ca, cn=x$cn)
  N.CC <- list(ca=N$ca, cn=N$cn)
  BF.CC <- bayes.factor.CC(x.CC, N.CC, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0)$BF
  BF <- BF.dn * BF.CC
  
  if (debug) {
    cat("BF.dn = ", BF.dn, "\n")
    cat("BF.CC = ", BF.CC, "\n")
  }
  
  return (BF)
}

#################################################################
# TADA-Denovo: analysis of de novo data 
#################################################################

# Genome-wide application of TADA for denovo data
# Input: counts, N, mu, gamma.mean, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  # Compute BFs
  for (i in 1:m) {
    for (j in 1:K)  BF[i,j] <- bayes.factor.denovo(counts[i,j], N, mu[i]*mu.frac[j], gamma.mean[j], beta[j])
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene: de novo only
# mu.gene: the mutation rates of a gene (K-dim. vector)
# N: the sample size
# l: the number of permutations
# gamma.mean, beta: RR of de novo mutations (vectors)
# Output: BF - l BFs from permutation; sample - permutate data
permute.gene.denovo <- function(mu.gene, N, l, gamma.mean, beta, dn.max=5) {
  K <- length(mu.gene)
  BF.gene <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  
  # permutation of l times
  count.gene <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N, mu.gene[j], gamma.mean[j], beta[j])
    }
    
    # permutation
    count.gene[,j] <- rpois(l, 2*N*mu.gene[j])
    for (i in 1:l) {
      x <- count.gene[i,j]
      cond.range.dn <- (x <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.gene[i,j] <- table.dn[x+1]
      } else {
        BF.gene[i,j] <- bayes.factor.denovo(x, N, mu.gene[j], gamma.mean[j], beta[j]) 
      }
    }
  }
  
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample=count.gene))
}

# Genome-wide application of TADA for denovo data: the difference with TADA.denovo is the report of p-values. 
# Input: counts, N, mu, mu.frac, gamma.mena, beta 
# counts: m x K matrix, where m is the number of gene, and K is the number of mutational categories. counts[i,j] is the number of de novo mutation in the j-th category of the i-th gene. 
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# gamma.mean, beta: the parameters of RR. Vector (one value for each category)
# l: the number of permutations to obtain the null distribution
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs. 
TADAp.denovo <- function(counts, N, mu, mu.frac, gamma.mean, beta, l=100, dn.max=5) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs
  rs <- TADA.denovo(counts, N, mu, mu.frac, gamma.mean, beta)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene.denovo(mu[i]*mu.frac, N, l, gamma.mean, beta, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# TADA: analysis of de novo and inherited data
#################################################################

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# denovo.only: whether using only de novo data (Boolean vector)
# pi.gene: for each gene, the estimated fractions of causal variants, one for each class of variants. These fractions will be used to set gene-specific RR (case-control)
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector)
TADA <- function(counts, N, mu, mu.frac, hyperpar, denovo.only=FALSE, pi.gene=1) {
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  
  if (length(denovo.only)==1) { denovo.only <- rep(denovo.only, m) }
  if (length(pi.gene)==1) { pi.gene <- array(1, dim=c(m,K))}
  
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
    
  # Compute BFs: BF[i,j] is the BF of the i-th gene in the j-th category
  for (i in 1:m) {
    if (denovo.only[i]==FALSE) {
      for (j in 1:K)  {
        # set hyperparameters
        hyperpar.gene <- hyperpar[,j]
        RR.product <- hyperpar.gene[3]*hyperpar.gene[4]
        hyperpar.gene[3] <- hyperpar.gene[3]*pi.gene[i,j] + (1-pi.gene[i,j])
        hyperpar.gene[4] <- RR.product/hyperpar.gene[3]
        
        # compute BF  
        start <- 3*(j-1)+1
        x <- list(dn=counts[i, start], ca=counts[i, start+1], cn=counts[i, start+2])
        BF[i,j] <- bayes.factor(x, N, mu[i]*mu.frac[j], hyperpar.gene)
      }
    } else {
      for (j in 1:K) {
        start <- 3*(j-1)+1
        x <- counts[i,start]
        BF[i,j] <- bayes.factor.denovo(x, N$dn, mu[i]*mu.frac[j], gamma.mean.dn[j], beta.dn[j])
      } 
    }
  }
  
  # Total BF per gene
  BF.total <- exp(rowSums(log(BF)))
  
  return (list(BF=BF, BF.total=BF.total))
}

# Compute permutation BFs of one gene
# mu.gene: the mutation rates of a gene (K-dim. vector), and the case-control counts (to be permuted): vectors (one value per category)
# N: sample size, three values for de novo, case and control
# l: number of permutations
# gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0: parameters. 
# table.CC: precomputed BFs. A 3-dim table, table.CC[i, j, k] stores the BF of (j-1, k-1) in the i-th category. 
# Output: BF - l BFs from permutation; sample.dn, sample.ca, sample.cn - permutate data
permute.gene <- function(mu.gene, count.ca, count.cn, N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=5) {
  K <- length(mu.gene)
  BF.dn <- array(1, dim=c(l,K))
  BF.CC <- array(1, dim=c(l,K))
  BF.gene.total <- numeric(l)
  ca.max <- dim(table.CC)[1] - 1
  cn.max <- dim(table.CC)[2] - 1
  
  # permutation l times for each category
  sample.dn <- array(0, dim=c(l,K))
  sample.ca <- array(0, dim=c(l,K))
  sample.cn <- array(0, dim=c(l,K))
  for (j in 1:K) {
    # pre-compute the de novo table for the j-th category
    table.dn <- numeric(dn.max+1)
    for (k in 0:dn.max) {
      table.dn[k+1] <- bayes.factor.denovo(k, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j])
    }
    
    # generate permutation data
    sample.dn[,j] <- rpois(l, 2*N$dn*mu.gene[j])
    sample.ca[,j] <- rhyper(l, count.ca[j] + count.cn[j], N$ca+N$cn-count.ca[j]-count.cn[j], N$ca)
    sample.cn[,j] <- count.ca[j] + count.cn[j] - sample.ca[,j]
    
    # compute the BFs
    for (i in 1:l) {
      x <- list(dn=sample.dn[i,j], ca=sample.ca[i,j], cn=sample.cn[i,j])
      cond.range.dn <- (x$dn <= dn.max)
      if (cond.range.dn==TRUE) {
        BF.dn[i,j] <- table.dn[x$dn+1]
      } else {
        BF.dn[i,j] <- bayes.factor.denovo(x$dn, N$dn, mu.gene[j], gamma.mean.dn[j], beta.dn[j]) 
      }
      
      cond.range.CC <- (x$ca <= ca.max) & (x$cn <= cn.max)
      if ( cond.range.CC==TRUE ) {
        BF.CC[i,j] <- table.CC[j, x$ca + 1, x$cn + 1]
      } else {
        BF.CC[i,j] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  BF.gene <- BF.dn * BF.CC
  BF.gene.total <- exp(rowSums(log(BF.gene)))
  
  return (list(BF=BF.gene.total, sample.dn=sample.dn, sample.ca=sample.ca, sample.cn=sample.cn))
}

# Genome-wide application of TADA 
# counts: m x 3K matrix, where m is the number of gene, and K is the number of mutational categories. Each category has three numbers: de novo, case and control. 
# N: sample size, three values for de novo, case and control
# mu: mutation rates of the genes
# mu.frac: the fraction of each category (vector)
# hyperpar: 8*K matrix, where each row is (gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0), and each column corresponds to one mutation type
# l: the number of permutations to obtain the null distribution
# dn.max, ca.max, cn.max: if counts are below these values, the BFs will be pre-computed. 
# Output: the statistics (BFs), an mxK matrix, and the total BFs (m-dim. vector). pval: the p-values of all genes. BF.null: the null distribution of BFs.
TADAp <- function(counts, N, mu, mu.frac, hyperpar, l=100, dn.max=5, ca.max=10, cn.max=10) {
  gamma.mean.dn <- hyperpar[1,]
  beta.dn <- hyperpar[2,]
  gamma.mean.CC <- hyperpar[3,]
  beta.CC <- hyperpar[4,]
  rho1 <- hyperpar[5,]
  nu1 <- hyperpar[6,]
  rho0 <- hyperpar[7,]
  nu0 <- hyperpar[8,]
  
  m <- length(mu)
  K <- length(mu.frac)
  BF <- array(1, dim=c(m,K))
  BF.total <- numeric(m)
  pval <- numeric(m)
  
  # Compute BFs of the observed data
  rs <- TADA(counts, N, mu, mu.frac, hyperpar)
  BF <- rs$BF
  BF.total <- rs$BF.total
  
  # Pre-compute the bayes-factors of the case-control data
  table.CC <- array(1, dim=c(K, (ca.max+1), (cn.max+1)))
  for (j in 1:K) {
    for (x1 in 0:ca.max) {
      for (x0 in 0:cn.max) {
        x <- list(ca=x1,cn=x0)
        table.CC[j, x1+1, x0+1] <- bayes.factor.CC(x, N, gamma.mean.CC[j], beta.CC[j], rho1[j], nu1[j], rho0[j], nu0[j])$BF
      }
    }
  }
  
  # Create the null distribution
  BF.null <- numeric(m*l)
  for (i in 1:m) {
#     print(i)
    BF.null[((i-1)*l+1):(i*l)] <- permute.gene(mu[i]*mu.frac, counts[i, 3*(1:K)-1], counts[i, 3*(1:K)], N, l, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, table.CC, dn.max=dn.max)$BF
  }
  
  # p-values of each gene
  BF.null.srt <- sort(BF.null, decreasing=TRUE)
  pval <- findInterval(-BF.total, -BF.null.srt)/length(BF.null.srt)
  pval[pval==0] <- 0.5/length(BF.null.srt)
  
  return (list(BF=BF, BF.total=BF.total, pval=pval, BF.null=BF.null))
}

#################################################################
# MOM estimation of hyperprior parameters from de novo data
#################################################################

# Prob. of having d or more de novo mutations under H1 
# Use simulation, but could also use analytic form 
multihit.prob <- function(N, mu, gamma.mean, beta, d=2, S=100) {
  p <- numeric(S)
  gamma <- rgamma(S, gamma.mean*beta, rate=beta)
  for (i in 1:S) {
    p[i] <- 1 - ppois(d-1, 2*N*mu*gamma[i])
  }
  return (mean(p))
}

# Estimate the number of multihit genes in a genome. 
# d: the parameter of the multiplicity test. 
# Returns: M0 - the number of multihit genes from non-risk genes; M1 - the number from risk genes. 
count.multihit <- function(N, mu, pi, gamma.mean, beta, d=c(2,3), S=2) {
  m <- length(mu)
  M0 <- numeric(length(d))
  M1 <- numeric(length(d))
  
  # M1: the number of causal genes having d or more de novo mutations
  p.alt <- array(0, dim=c(m, length(d)))  # p1[i,j] = P(X_i >= d_j|H1)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.alt[i,j] <- multihit.prob(N, mu[i], gamma.mean, beta, d=d[j], S=S)
    }
  }
  for (j in 1:length(d)) { 
    M1[j] <- m * pi  * mean(p.alt[,j]) 
  }
  
  # M0: the number of non-causal genes having d or more de novo mutations
  p.null <- array(0, dim=c(m, length(d)))  # p0[i,j] = P(X_i >= d_j|H0)
  for (i in 1:m) {
    for (j in 1:length(d)) {
      p.null[i,j] <- 1 - ppois(d[j] - 1, 2*N*mu[i])
    }
  }
  for (j in 1:length(d)) { 
    M0[j] <- m * (1-pi) * mean(p.null[,j]) 
  }
  
  result <- data.frame(d=d, M0=M0, M1=M1)
  return (result)
}

# Estimating relative risk and the number of multiple hits from de novo data
# Input: sample size (N), mutation rates of all genes (mu), observed number of de novo events (C), beta (parameter of the prior distribution of gamma), k (number of disease genes)
# Output: the average relative risk (gamma.mean), the expected number of multi-hit genes (M)
denovo.MOM <- function(N, mu, C, beta, k) {
  m <- length(mu) # number of genes
  
  # enrichment of de novo events
  nu <- C / (2 * N * sum(mu))
  
  # MOM estimator of gamma.mean
  gamma.mean <- (nu-1)*m/k +1
  
  # expected M (choose d = 2)
  rs <- count.multihit(N, mu, k/m, gamma.mean, beta, d=2)    
  M <- sum(rs$M1) + sum(rs$M0)
  
  return (list(gamma.mean=gamma.mean, M=M))
}

#################################################################
# Empirical Bayes estimation of hyperprior parameters
#################################################################

# Evalaute the marginal log-likelihood at given parameters
# pi: the fraction of causal genes
# counts: the count data (of one mutational category), a date frame
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the negative log-likelihood, and posterior, BF of each gene
marginal <- function(hyperpar, pi, counts, N, mu, prior.weight=0, denovo.only=FALSE, debug=FALSE) {
  gamma.mean.dn <- hyperpar[1]
  beta.dn <- hyperpar[2]
  gamma.mean.CC <- hyperpar[3]
  beta.CC <- hyperpar[4]
  rho1 <- hyperpar[5]
  nu1 <- hyperpar[6]
  rho0 <- hyperpar[7]
  nu0 <- hyperpar[8]
  
  n <- nrow(counts)
  prob <- numeric(n)
  posterior <- numeric(n)
  bf <- numeric(n)
  for (i in 1:n)  {
    if (debug) cat("i = ", i, "\tdn = ", counts[i,]$dn, "\tca = ", counts[i,]$ca, "\tcn = ", counts[i,]$cn, "\n")
    if (denovo.only==TRUE) {
      prob.M1 <- evidence.alt.dn(counts[i,"dn"], N$dn, mu[i], gamma.mean.dn, beta.dn) 
      prob.M0 <- evidence.null.dn(counts[i,"dn"], N$dn, mu[i])
    } else {
      prob.M1 <-  evidence.alt(counts[i,],N,mu[i],gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1)
      prob.M0 <-  evidence.null(counts[i,],N,mu[i],rho0, nu0) 
    }
    if (debug) cat("prob.M1 = ", prob.M1,"\n")
    if (debug) cat("Prob.M0 = ", prob.M0, "\n")
    prob[i] <- pi * prob.M1 + (1 - pi) * prob.M0
    posterior[i] <- pi*prob.M1 / (prob[i])
    bf[i] <- prob.M1 / prob.M0
  }
  marg.nll <- -sum(log(prob[!is.na(prob)]))
  
  # add the prior
  if (prior.weight > 0) {
    u <- prior.weight
    q1.mean <- rho1/nu1
    q0.mean <- rho0/nu0
    marg.nll <- marg.nll - log(dgamma(q1.mean, q0.mean*u, u)) 
  }
  
  result <- list(prob=prob, marginal=marg.nll, posterior=posterior, bf=bf)
  return (result)
}

# Convert the full set of parameters (pi, hyperpar) to a subset based on the option
# est.option: a Boolean vector (size 8), one for each hyperpar. If FALSE, the corresponding parameter will be fixed. 
# est.pi: Boolean. 
fullpar2subpar <- function(hyperpar, pi, est.option, est.pi) {
  subpar <- NULL
  for (i in 1:length(est.option)) {
    if (est.option[i]) subpar <- c(subpar, hyperpar[i])
  }
  
  if (est.pi==TRUE)  subpar <- c(subpar, pi)
  
  return (subpar)
}

# Convert the subset of parameters (subpar) to the full set based on the option and intial values
# hyperpar.init, pi.init: used to set parameters not to be estimated
subpar2fullpar <- function(subpar, hyperpar.init, pi.init, est.option, est.pi) {
  hyperpar <- numeric(length(hyperpar.init))
  curr <- 1
  for (i in 1:length(est.option)) {
    if (est.option[i]) { hyperpar[i] <- subpar[curr]; curr <- curr+1 }
    else hyperpar[i] <- hyperpar.init[i]
  }
  if (est.pi==TRUE)  pi <- subpar[curr]
  else pi <- pi.init
  
  return (list(hyperpar=hyperpar, pi=pi))
}

# Empirical Bayes estimation of hyperparameters 
# counts: the count data - the number of mutations in cases, controls and de novo data respectively. 
# N: sample size (dn, ca, cn, respectively). 
# mu: vector of mutation rates
# lower, upper: the search range of the parameters 
# est.option: a Boolean vector specifying whether to estimate each parameter. If FALSE, use the corresponding parameter in hyperpar.init
# est.pi: whether to estimate pi. If FALSE, use the given value of pi; if TRUE, use the given value as initial
# prior.weight: putting a prior so that q1 and q0 tend to be close (recommended value: 5000). Default 0. 
# Output: the hyperparameters and the NLL (negative log likelihood) at the parameters. 
empBayes <- function(counts, N, mu, hyperpar.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, est.option=rep(FALSE, 8), est.pi=FALSE, prior.weight=0, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(subpar.log) {
    subpar <- exp(subpar.log)
    allpar <- subpar2fullpar(subpar, hyperpar.init, pi.init, est.option, est.pi)
    hyperpar <- allpar$hyperpar
    pi <- allpar$pi
    result <- marginal(hyperpar,  pi, counts, N, mu, prior.weight=prior.weight)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(subpar, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # no parameter to estimate
  if ( sum(est.option==TRUE) == 0 & est.pi == FALSE ) {
    value <- marginal(hyperpar.init, pi.init, counts, N, mu, prior.weight=prior.weight)$marginal
    return (list(hyperpar=hyperpar.init, value=value))
  }
  
  # initialize the sub. parameters
  subpar.init <- fullpar2subpar(hyperpar.init, pi.init, est.option, est.pi)
  sublower <- fullpar2subpar(lower, lower.pi, est.option, est.pi)
  subupper <- fullpar2subpar(upper, upper.pi, est.option, est.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(subpar.init), marginal.loglike, method="L-BFGS-B", lower=log(sublower), upper=log(subupper))
  subpar.est.log <- like.optim$par
  value.max <- like.optim$value
  allpar.est <- subpar2fullpar(exp(subpar.est.log), hyperpar.init, pi.init, est.option, est.pi)
  hyperpar.est <- allpar.est$hyperpar
  pi.est <- allpar.est$pi
  
  if (est.pi==TRUE) { return (list(pi=pi.est, hyperpar=hyperpar.est, value=value.max)) }
  else { return (list(hyperpar=hyperpar.est, value=value.max)) }
}

# Empirical Bayes estimation of hyperparameters (gamma.mean, beta and pi) for de novo data
# counts: count data - the number of de novo mutations. 
# N - sample size (trios). 
# mu: vector of mutation rate
# lower, upper: the search range of the parameters 
empBayes.denovo <- function(counts, N, mu, gamma.mean.init, beta.init, pi.init, lower, upper, lower.pi=1e-10, upper.pi=1, debug=FALSE) {  
  # marginal likelihood function (parameters in log-scale)
  marginal.loglike <- function(par.log) {
    par <- exp(par.log)
    gamma.mean <- par[1]
    beta <- par[2]
    pi <- par[3]
    hyperpar <- c(gamma.mean, beta, 0, 0, 0, 0, 0, 0)
    m <- length(counts)
    counts.full <- data.frame(dn=counts, ca=numeric(m), cn=numeric(m))
    N.full <- list(dn=N, ca=0, cn=0)
    result <- marginal(hyperpar, pi, counts.full, N.full, mu, denovo.only=TRUE)$marginal  
    
    if (debug) {
      cat("Parameters = (", paste(par, collapse=", "), ")\tValue = ", result, "\n", sep="") 
    }
    
    return (result)
  }
  
  # initialize the parameters
  par.init <- c(gamma.mean.init, beta.init, pi.init)
  par.lower <- c(lower, lower.pi)
  par.upper <- c(upper, upper.pi)
  
  # maximization of marginal likelihood
  like.optim <- optim(log(par.init), marginal.loglike, method="L-BFGS-B", lower=log(par.lower), upper=log(par.upper))
  par.est <- exp(like.optim$par)
  value.max <- like.optim$value
  gamma.mean.est <- par.est[1]
  beta.est <- par.est[2]
  pi.est <- par.est[3]
  
  return (list(gamma.mean=gamma.mean.est, beta=beta.est, pi=pi.est, value=value.max)) 
}

#################################################################
# Functions for simulation
#################################################################

# Generate simulation data of a set of genes (multiple mutational categories): de novo mutations only
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean*beta, beta). Vectors (one per category) 
# Output: sample matrix (m by K), where m is the number of genes and K the number of variant categories. sample.info: more information of the samples, including the indicator (risk gene or not) and the RR. 
simulator.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  z <- rbinom(m, 1, pi)
  gamma <- array(1, dim=c(m,K))
  x <- array(0, dim=c(m,K))
  k <- sum(z==1)
  for (j in 1:K) {
    gamma[z==1, j] <- rgamma(k, gamma.mean[j]*beta[j], beta[j])
    x[,j] <- rpois(m, 2 * mu * mu.frac[j] * gamma[,j] * N)
  }
  
  sample.info <- cbind(mu, z, gamma, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Generate simulation data of a set of genes (multiple mutational categories)
# N: sample size (number of trios)
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: Relative risk of de novo mutation: gamma|M1 ~ Gamma(gamma.mean.dn*beta.dn, beta.dn). Vectors.
# gamma.mean.CC, beta.CC: Relative risk of inherited mutation (case/control): gamma.CC|M1 ~ Gamma(gamma.mean.CC*beta.CC, beta.CC). Vectors
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# tradeoff option: if TRUE, implement q-gamma tradeoff (i.e. higher gamma means lower q). Suppose, gamma_i is the RR, then q_i is proportional to mu_i / gamma_i, where the constant is determined from the mean of q, mu and gamma. 
# Output: sample matrix (m by 3K), where m is the number of genes and K the number of variant categories.  sample.info: more information of the samples, including the indicator (risk gene or not) and the RR.
simulator <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=FALSE) {
  m <- length(mu) # number of genes
  K <- length(mu.frac) # number of mutational categories
  
  # the tradeoff parameter (delta:=mu.mean/q.mean)
  delta <- mean(mu) * mu.frac / (rho0 / nu0)
  
  z <- rbinom(m, 1, pi)
  gamma.dn <- array(1, dim=c(m,K))
  gamma.CC <- array(1, dim=c(m,K))
  q <- array(0, dim=c(m,K))
  x <- array(0, dim=c(m,3*K))
  k <- sum(z==1)
  for (j in 1:K) {
    # sample de novo 
    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
    col <- 3*(j-1)+1
    x[,col] <- rpois(m, 2 * mu * mu.frac[j] * gamma.dn[,j] * N$dn)
    
    # sample case-control
    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
    if (tradeoff==FALSE) {
      q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
    } else {
      q[z==1, j] <- mu[z==1] * mu.frac[j] / (delta[j] * gamma.CC[z==1, j])
    }
    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca)
    x[,col+2] <- rpois(m, q[,j] * N$cn)
    
  }
  
  sample.info <- cbind(mu, z, gamma.dn, gamma.CC, q, x)
  
  return (list(sample=x, sample.info=sample.info))
}

# Evaluation of the power (FDR) of TADA.denovo
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean, beta: RR parameters for de novo mutations (vectors)
# gamma.mean.est, beta.est: the parameters used by TADA.denovo
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA.denovo <- function(N, mu, mu.frac, pi, gamma.mean, beta, gamma.mean.est, best.est, FDR=0.1) {
  # sample the simulation data
  sample <- simulator.denovo(N, mu, mu.frac, pi, gamma.mean, beta)$sample
  
  # run TADA.denovo
  sampleBF <- TADA.denovo(sample, N, mu, mu.frac, gamma.mean.est, beta.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  return (M1)
}

# Evaluation of the power (FDR) of TADA
# N: sample size, the number of families of the de novo study
# mu: mutation rate of each gene (a vector)
# mu.frac: for each type of mutation, its fraction of the gene's mutation rate
# pi: the fraction of risk genes 
# gamma.mean.dn, beta.dn: RR parameters for de novo mutations (vectors)
# gamma.mean.CC, beta.CC: RR parameters for inherited mutations (vectors)
# Frequency parameter of risk genes: q|M1 ~ Gamma(rho1, nu1)
# Frequency parameter of non-risk genes: q|M0 ~ Gamma(rho0, nu0)
# hyperpar.est: the parameters used by TADA
# FDR: the FDR level to be controled
# Output: the number of discoveries at the specified FDR level
eval.TADA <- function(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, hyperpar.est, FDR=0.1, tradeoff=FALSE) {
  # sample the simulation data
  sample <- simulator(N, mu, mu.frac, pi, gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=tradeoff)$sample
  
  # run TADA
  sampleBF <- TADA(sample, N, mu, mu.frac, hyperpar.est)$BF.total
  
  # Bayesian FDR control
  M1 <- Bayesian.FDR(sort(sampleBF, decreasing=TRUE), 1-pi, FDR)$ND
  
  # sample information (counts and BFs)
  sim.results <- data.frame(counts=sample, BF=sampleBF)
  sim.results <- sim.results[order(-sim.results$BF),]
  
  return (list(M1=M1, sample=sim.results))
}
