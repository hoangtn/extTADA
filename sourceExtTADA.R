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

