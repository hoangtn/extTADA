calculateFDR <- function(pars,
                         caseData = NULL,
                         controlData = NULL,
                         dnData = NULL,
                         mutData = NULL,
                         geneName){
    outData <- data.frame(geneName)
    if (!is.null(dnData))
        outData <- cbind(outData, dnData)
    if (!is.null(mutData))
        outData <- cbind(outData, mutData)
    if (!is.null(caseData))
        outData <- cbind(outData, caseData)
    if (!is.null(controlData))
        outData <- cbind(outData, controlData)


    bfAll <- rep(1, dim(outData)[1])

    if ( length(pars$gammaMeanDN) == 0) {
        message("No parameters for de novo data; therefore, these categories are not calculated in this step.\n")
        }  else {
        bfDN <- matrix(1, nrow = dim(dnData)[1], ncol = dim(dnData)[2]) ##De novo bayes factors
        for (j2 in 1:dim(bfDN)[2]) {
            e.hyperGammaMeanDN <- pars$gammaMeanDN[j2]
            e.hyperBetaDN <- pars$betaDN[j2]
            e.bf <- bayes.factor.denovo(x =  dnData[, j2],
                                    N = pars$nfamily[j2],
                                    mu =  mutData[, j2],
                                    gamma.mean = e.hyperGammaMeanDN,
                                    beta = e.hyperBetaDN)
            bfDN[, j2] <- e.bf
        }
        bfAll <- bfAll*apply(bfDN, 1, prod)
    }

    if (length(pars$gammaMeanCC) == 0) {
        message("No parameters for case-control data;  therefore, these categories are not calculated in this step.\n")
    } else {

        bfCC <- matrix(1, ncol = dim(caseData)[2], nrow = dim(caseData)[1])
        for (cc3 in 1:dim(bfCC)[2]){
            e.hyperGammaMeanCC <- pars$gammaMeanCC[cc3]
            e.hyperBetaCC <- pars$betaCC[cc3]
            e.nu <- 200
            t.case <- caseData[, cc3]
            t.control <- controlData[, cc3]
            e.rho <- e.nu*mean(t.case + t.control)/(pars$ncase[cc3] + pars$ncontrol[cc3])
            e.bf <- BayesFactorCC3(Nsample = list(ca = pars$ncase[cc3], cn = pars$ncontrol[cc3]),
                                   x.case = t.case, x.control = t.control,
                                   gamma.meanCC = e.hyperGammaMeanCC, betaCC = e.hyperBetaCC,
                                   rhoCC = e.rho, nuCC = e.nu)
            bfCC[, cc3] <- e.bf
        }

        bfAll <- bfAll*apply(bfCC, 1, prod)
    }
    outData$BF <- bfAll
    outData <- outData[order(-outData$BF),]
    outData$qvalue <- Bayesian.FDR(outData$BF, 1 - pars$pi0)$FDR

    return(outData)
}

