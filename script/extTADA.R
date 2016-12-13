##########There are three models written in this script: DNandCCextTADA (de novo + case control), DNextTADA (only de novo), CCextTADA (only case-control)
##########Users can change prior information as they wish
library("rstan")
#library("coda")
########################################
#######De novo + Case control
#August 1, 2016: Hoang Nguyen

DNandCCextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCdn; //Number of de novo classes
    int<lower=1> NCcc; //Number of case-control classes
    int Ndn[NCdn]; //Number of trios
    int Ncase[NCcc]; //Number of cases
    int Ncontrol[NCcc]; //Number of controls
    int Ntotal[NCcc]; //Number of cases + controls

    int dataDN[NN, NCdn]; //denovo data: Kdn classes
    real mutRate[NN, NCdn]; //mutation rates: Kdn classes
    int dataCCcase[NN, NCcc]; //case-control data: Kcc classes
    int dataCCtotal[NN, NCcc]; //case+control data: Kcc classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; // Ncase/(Ncase + Ncontrol)

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
    real<lower=0> hyper2GammaMeanDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
    real<lower=0> hyper2BetaDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
    real<lower=0> hyper2GammaMeanCC[NCcc]; //hyperGammaMeanCC ~ gamma(hyper2GammaMeanCC, hyper2BetaCC)
    real<lower=0> hyper2BetaCC[NCcc]; //hyperGammaMeanCC ~ gamma(hyper2GammaMeanCC, hyper2BetaCC)




    }

parameters {
    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks

//    real<lower=lowerBeta> betaDN[Kdn]; //Rate parameters for de novo relative risks
    //real<lower=lowerBeta> hyperBetaDN[Kdn]; //Hyper rate parameters for de novo relative risks
    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    //ordered[Kcc] hyperGammaMeanCC;
    real<lower=lowerGamma> gammaMeanCC[NCcc]; //parameters (in the sampling process) for de novo relative risks: gammaMeanCC ~ gamma(hyperGammaMeanCC*hyperBetaCC, hyperBetaCC)


}

transformed parameters {
    real hyperBetaCC[NCcc];
    real hyperBetaDN[NCdn];
    if (adjustHyperBeta != 0) {
      for (i1i in 1:NCcc){
            hyperBetaCC[i1i] = exp(betaPars[1]*hyperGammaMeanCC[i1i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }

      for (i2i in 1:NCdn){
            hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }
   }
    else {
        hyperBetaCC = hyperBetaCC0;
        hyperBetaDN = hyperBetaDN0;
        }
    }
 model {
     int newIndex;
     real ps[K];
     real sumDN[2];
     real sumCC[2];
     pi0 ~ beta(1, 5); //prior for the proportion of risk genes

     //Both CC + DN

     //Case control data: sample for hyper priors (NPcc populations and Kcc categories)
     for (ip in 1:NCcc){
         hyperGammaMeanCC[ip]  ~ gamma(hyper2GammaMeanCC[ip], hyper2BetaCC[ip]); //gamma(1, 0.05); //normal(15, 10); //gamma(1, 0.1);

     }

     for (ip in 1:NCcc){
         gammaMeanCC[ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
         }

  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
          hyperGammaMeanDN[ip]	~ gamma(hyper2GammaMeanDN[ip], hyper2BetaDN[ip]); //gamma(1, 0.02); //gamma(1, 0.05); //gamma(1, 0.1); //normal(15, 10);
    }

    for (ip in 1:NCdn){
          gammaMeanDN[ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0);
         ps[2] = log(pi0);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
    //For case-control data
         for (jj in 1:NCcc){
             ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
             ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[jj]*Ncase[jj]/(gammaMeanCC[jj]*Ncase[jj] + Ncontrol[jj])); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"


########################################
#######De novo only
DNextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCdn; //Number of de novo classes
    int Ndn[NCdn]; //Number of trios

    int dataDN[NN, NCdn]; //denovo data: Kdn classes
    real mutRate[NN, NCdn]; //mutation rates: Kdn classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaDN0[NCdn];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
    real<lower=0> hyper2GammaMeanDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
    real<lower=0> hyper2BetaDN[NCdn]; //hyperGammaMeanDN ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)

    }

parameters {
    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanDN[NCdn]; //Hyper parameters for de novo relative risks
    real<lower=lowerGamma> gammaMeanDN[NCdn]; //parameters (in the sampling process) for de novo relative risks

}

transformed parameters {
    real hyperBetaDN[NCdn];
    if (adjustHyperBeta != 0) {
      for (i2i in 1:NCdn){
            hyperBetaDN[i2i] = exp(betaPars[1]*hyperGammaMeanDN[i2i]^(betaPars[2]) + betaPars[3]);

       }
   }
    else {
        hyperBetaDN = hyperBetaDN0;
        }
    }
 model {
     int newIndex;
     real ps[K];
     real sumDN[2];
     pi0 ~ beta(1, 5); //prior for the proportion of risk genes

     //Both CC + DN


  //De novo data: sample for hyper priors (NPdn populations and Kdn categories)
    for (ip in 1:NCdn){
          hyperGammaMeanDN[ip]	~ gamma(hyper2GammaMeanDN[ip], hyper2BetaDN[ip]); //gamma(1, 0.1); //normal(15, 10);
    }

    for (ip in 1:NCdn){
          gammaMeanDN[ip] ~ gamma(hyperGammaMeanDN[ip]*hyperBetaDN[ip], hyperBetaDN[ip]);
           }

////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0);
         ps[2] = log(pi0);
   //For de novo data
         for (jj in 1:NCdn){
             ps[1] = ps[1] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]); //Add Null hypothesis
             ps[2] = ps[2] + poisson_lpmf(dataDN[ii, jj] | Ndn[jj]*2*mutRate[ii, jj]*gammaMeanDN[jj]); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"


################CC model
########################
CCextTADA <- "
data {
    int<lower=1> NN; //Number of genes
    int<lower=1> K; //Number of classes
    int<lower=1> NCcc; //Number of case-control classes
    int Ncase[NCcc]; //Number of cases
    int Ncontrol[NCcc]; //Number of controls
    int Ntotal[NCcc]; //Number of cases + controls

    int dataCCcase[NN, NCcc]; //case-control data: Kcc classes
    int dataCCtotal[NN, NCcc]; //case+control data: Kcc classes
    real betaPars[3]; //These parameters are used to adjust beta values as a non-linear function of gamma means

    real<lower=0> thetaH0[NCcc]; // Ncase/(Ncase + Ncontrol)

    //These below parameters should be default
    real<lower=0> upperPi0;
    real<lower=0> lowerPi0;
    real<lower=0> lowerHyperGamma; //Low limit for mean of relative risks
    real<lower=0> lowerGamma; //Low limit for relative risks
    real<lower=0> lowerBeta;
    real<lower=0> hyperBetaCC0[NCcc];
    int<lower=0> adjustHyperBeta; //Option to adjust betas or not; if this option is used (!=0) then beta is a non-linear function of gamma means
    real<lower=0> hyper2GammaMeanCC[NCcc]; //hyperGammaMeanCC ~ gamma(hyper2GammaMeanCC, hyper2BetaCC)
    real<lower=0> hyper2BetaCC[NCcc]; //hyperGammaMeanCC ~ gamma(hyper2GammaMeanCC, hyper2BetaCC)

    }

parameters {
    real<lower=lowerPi0,upper=upperPi0> pi0; //Proportion of risk genes
    real<lower=lowerHyperGamma> hyperGammaMeanCC[NCcc]; //Hyper parameters for case-control relative risks
    real<lower=lowerGamma> gammaMeanCC[NCcc]; //parameters (in the sampling process) for de novo relative risks: gammaMeanCC ~ gamma(hyperGammaMeanCC*hyperBetaCC, hyperBetaCC)


}

transformed parameters {
    real hyperBetaCC[NCcc];

    if (adjustHyperBeta != 0) {
      for (i1i in 1:NCcc){
            hyperBetaCC[i1i] = exp(betaPars[1]*hyperGammaMeanCC[i1i]^(betaPars[2]) + betaPars[3]); //,  hyperBetaMax);

       }

   }
    else {
        hyperBetaCC = hyperBetaCC0;
            }
    }
 model {
     int newIndex;
     real ps[K];
     real sumDN[2];
     real sumCC[2];
    pi0 ~ beta(1, 5); //prior for the proportion of risk genes

     //Both CC + DN

     //Case control data: sample for hyper priors (NPcc populations and Kcc categories)
     for (ip in 1:NCcc){
         hyperGammaMeanCC[ip]  ~ gamma(hyper2GammaMeanCC[ip], hyper2BetaCC[ip]); //gamma(1, 0.05); //normal(15, 10); //gamma(1, 0.1);

     }

     for (ip in 1:NCcc){
         gammaMeanCC[ip] ~ gamma(hyperGammaMeanCC[ip]*hyperBetaCC[ip], hyperBetaCC[ip]);
         }


////Main program
//Loop through data points
////
     for (ii in 1:NN){
         ps[1] = log1m(pi0);
         ps[2] = log(pi0);

    //For case-control data
         for (jj in 1:NCcc){
             ps[1] = ps[1] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], thetaH0[jj]); //Add Null hypothesis
             ps[2] = ps[2] + binomial_lpmf(dataCCcase[ii, jj] | dataCCtotal[ii, jj], gammaMeanCC[jj]*Ncase[jj]/(gammaMeanCC[jj]*Ncase[jj] + Ncontrol[jj])); //Add alternative hypothesis
             }
         target += log_sum_exp(ps);
         //increment_log_prob(log_sum_exp(ps));
         }
}
"

###Bayes Factor For Case-control
BayesFactorCC3 <- function(x.case, x.control, Nsample,
                           gamma.meanCC, betaCC, rhoCC, nuCC){
    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

                                        #    print(range(gAll))
    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
                                        #                print(bControl)
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
                                        #print(bCase)
            return(exp(bCase + bControl))
        }
        t1 <- evidence.alt.cc3()

        return(t1)
    })
                                        #    print(altCC)
    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
        bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, rho1 + x$cn, (N$cn + nu1)/(N$cn + nu1 + N$ca)))


        t1 <- exp(bCase + bControl)

        return(t1)
    })
                                        #    print(nullCC)

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(tempBF)
}

################


library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}








loc2plot.old <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}


#gives the HPD value for an observation px,py in a density constructed from x, y.
gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
#	d1 <- (x-px)^2+(y-py)^2
#	best <- d1 == min(d1)
#	lev <- mean(fitted(fit)[best])
	lev <- predict.locfit(fit,list(px,py))
	slev <- sort(fitted(fit))
	indic <- slev <= lev
	sum(indic)/length(x)
}

