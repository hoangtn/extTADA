x <- rnorm(100)
source("script/TADA.R")
source("script/extTADA.R")
fileN <- "DN.adjustBeta.1.3.nI20000.nThin.1.index.15_42_Dec_01_2016.lof_DD.damaging_DD.nGroupDN.2.nGroupCC.2.adjustRatio.1.RData"
fileN <- "DN.adjustBeta.1.3.nI20000.nThin.1.index.15_42_Dec_01_2016.lof_AST.damaging_AST.nGroupDN.2.nGroupCC.2.adjustRatio.1.RData"
fileN <- "GuipponiArrayNoFINcombined.DNandCC.adjustHyperBeta.0..adjustRatioMut.1.3.nIteration20000.nThin.20.index.10_46_24_November.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData"
#load(paste0("../NewTestR/TestR/HBproject/MainPaper/RDataNewDec1/", fileN))
load(paste0("../NewTestR/TestR/HBproject/MainPaper/RDataNew1/", fileN))



fileR <- dir("script", ".R$")
for (ii in fileR)
        source(paste0("script/", ii))
data <- read.table("data/data_mut_DD.csv", header = TRUE, as.is = TRUE)
head(data
     allDNData <- data[, paste0("dn_", c("damaging", "lof"), "_DD")]
     allMutData <- data[,paste0("mut_", c("damaging", "lof"))]
     head(data.frame(allMutData, allDNData))

extTADA <- function(modelName  ,
                         dataDN = NULL, mutRate = NULL, Ndn = NULL,
                         dataCCcase = NULL, dataCCcontrol = NULL, Ncase = NULL, Ncontrol = NULL,
                    nIteration = NULL, nIteration2 = NULL,
                    nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta = TRUE)
     {

         if (is.null(nIteration))
             stop("======\nNo input for the nIteration parameter\n=====")
         if ((is.null(dataDN) | is.null(mutRate)) & (is.null(dataCCcase) | is.null(dataCCcontrol)))
             stop("Need to have input data: only DN, only CC or DN + CC")
         if (is.null(nThin))
             nThin = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
         if (nCore != nChain)
             warning("nCore is different from nChain")
         if (!is.null(dataDN))
             Ngene <- dim(dataDN)[1]
         if (!is.null(dataCCcase))
             Ngene <- dim(dataCCcase)[1]
         message("\nThere are ", Ngene, " in this analysis")

         if (is.null(hyper2GammaMeanDN) & !is.null(dataDN))
             hyper2GammaMeanDN <- rep(1, dim(dataDN)[2])
         if (is.null(hyper2BetaDN) & !is.null(dataDN))
             hyper2BetaDN <- rep(0.025, dim(dataDN)[2])
        if (is.null(hyper2GammaMeanCC) & !is.null(dataCCcase))
             hyper2GammaMeanCC <- rep(1, dim(dataCCcase)[2])
         if (is.null(hyper2BetaCC) & !is.null(dataCCcase))
             hyper2BetaCC <- rep(0.2, dim(dataCCcase)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaDN0))
             hyperBetaDN0 <- rep(1, dim(dataDN)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaCC0))
             hyperBetaCC0 <- rep(1, dim(dataCCcase)[2])

         NCdn <- dim(dataDN)[2]
         NCcc <- dim(dataCCcase)[2]
         if (is.null(Ncase))              Ncase = 0
         if (is.null(Ncontrol)) Ncontrol = 1
         if (is.null(Ndn)) Ndn = 0
         if (is.null(hyperBetaDN0)) hyperBetaDN0 <- rep(1, length(dataDN[1, ]))
         if (is.null(hyperBetaCC0)) hyperBetaCC0 <- rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2GammaMeanDN)) hyper2GammaMeanDN = rep(1, length(dataDN[1, ]))
         if (is.null(hyper2GammaMeanCC)) hyper2GammaMeanCC = rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2BetaDN)) hyper2BetaDN = rep(0.05, length(dataDN[1, ]))
         if (is.null(hyper2BetaCC)) hyper2BetaCC = rep(0.2, length(dataCCcase[1, ]))



         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                          hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC)  )

         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
####Re-sample using new hyper betas
         if ( autoAdjustHyperBeta){
             if (!is.null(nIteration2))
                 nIteration = nIteration2

             if (is.null(nThin2))
                 nThin2 = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
             nThin = nThin2

         mcmcData <- as.data.frame(mcmcModel)
         cName <- colnames(mcmcData)
         hyperGammaMeanNameCC <- cName[grep("hyperGammaMeanCC", cName)]
         hyperGammaMeanNameDN <- cName[grep("hyperGammaMeanDN", cName)]

             hyperGammaMeanDN0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameDN]), 2, median))
             hyperGammaMeanCC0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameCC]), 2, median))

             hyperBetaDN0 <- exp(betaPars[1]*hyperGammaMeanDN0^(betaPars[2]) + betaPars[3])
             hyperBetaCC0 <- exp(betaPars[1]*hyperGammaMeanCC0^(betaPars[2]) + betaPars[3])

             adjustHyperBeta <- 0

         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                           hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC)  )

         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
         }



         return(mcmcModel)
          }


mcmcDD <- extTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData, Ndn = rep(4293, 2),
                  nIteration = 1000, nIteration2 = 2000)

load("../DataFile/Test4testPrior1077n2CC/TestArray.50008.gammaMeanDN.1.5.gammaMeanDenovo.5.gammaMeanDN2.1.5.gammaMeanDenovo2.5.pi0.0.02.rhocC.0.0002161572.rhoCC2.0.000123642505820651.lowerHyperGamma.1.casecontrolAndDenovo.1classes.RData")
mcmcCC <- extTADA(modelName = CCextTADA, dataCCcase = data.frame(y.case.lof), dataCCcontrol = data.frame(y.control.lof),
                  Ncase = 3157, Ncontrol = 4672,
                  nIteration = 1000, nIteration2 = 10000)
stan_trace(mcmcCC)

stan_trace(testIntegratedModel)

estimatePars <- function(pars, mcmcResult, nThin = NULL){
    mcmcDataFrame <- as.data.frame(mcmcResult)

    if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)
    if (length(pars[!is.element(pars, colnames(mcmcDataFrame))]) > 0)
        warning((pars[!is.element(pars, colnames(mcmcDataFrame))]), " is/are not in mcmc chains")
    pars <- pars[is.element(pars, colnames(mcmcDataFrame))]

    hpdList <- NULL
    for (iPar in pars) {
        xData <-  mcmcDataFrame[, iPar]
        if ((sd(xData) > 10^-3) & (length(grep("Beta", iPar)) == 0))
            hpdList <- rbind(hpdList, loc1stats(xData))
        else
            hpdList <- rbind(hpdList, rep(median(xData), 3))
    }
    rownames(hpdList) <- pars

return(hpdList)
}

pars0 <- estimatePars(pars = colnames(as.data.frame(testIntegratedModel)),
                     mcmcResult = testIntegratedModel)

plotParHeatmap <- function(pars, mcmcResult, nThin = NULL, color = "blue",
                           xLim = NULL, yLim = NULL,
                           mainLab = NULL, xLab = NULL, yLab = NULL,
                           maxk0 = 500, cprob = c(0.5, 0.05)){
    mcmcDataFrame <- as.data.frame(mcmcResult)
    if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)

    x <- mcmcDataFrame[, pars[1]]
    y <- mcmcDataFrame[, pars[2]]
    sc1<-sd(x)
    sc2<-sd(y)
    fit <- locfit(~x+y,scale=c(sc1,sc2), maxk = maxk0)
    lev <- sort(fitted(fit))[floor(cprob*length(x))]
    max.i <- fitted(fit)==max(fitted(fit))
    mode.x<-x[max.i][[1]]
    mode.y<-y[max.i][[1]]
    if (is.null(yLim))
        yLim <- c(0, max(y))
    if (is.null(xLim))
        xLim <- c(0, max(x))
    if (is.null(xLab))
        xLab <- pars[1]
    if (is.null(yLab))
        yLab <- pars[2]

    plot(mode.x,mode.y,
     #xaxt="n",
         main = mainLab, #diseaseName,
         xlim = xLim, ylim = yLim,
         pch=-1, xlab = xLab, ylab = yLab)

my.color.palette=colorRampPalette(c("white", color), space = "Lab")
plot(fit, type="image", m=300, add=TRUE,
     col=my.color.palette(25)[-c(2,3,5,6,8,10,12)])

plot(fit,add=TRUE,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""), labcex=0.75,vfont=c("sans serif","bold"))
points(mode.x,mode.y,pch=3)

}

plotParHeatmap(pars = c("pi0", "hyperGammaMeanDN[1]"), mcmcResult = testIntegratedModel)

#pars = list(gammaMeanDN, gammaMeanCC, betaDN, betaCC, pi0)
allPar = pars0[, 1]
pars = list(gammaMeanDN = allPar[grep("hyperGammaMeanDN", names(allPar))],
            betaDN = allPar[grep("hyperBetaDN", names(allPar))],
            gammaMeanCC = allPar[grep("hyperGammaMeanCC", names(allPar))],
            betaCC = allPar[grep("hyperBetaCC", names(allPar))],
            pi0 = allPar[grep("pi0", names(allPar))],
            nfamily = rep(ntrio, 3),
            ncase = allNcase,
            ncontrol = allNcontrol)
caseData <- allCaseData
controlData <- allControlData
dnData <- allDenovoData
mutData <- allMutationData
geneName <- data[, 1]
calculateFDR <- function(pars,
                         caseData,
                         controlData,
                         dnData,
                         mutData,
                         geneName){
    outData <- data.frame(geneName, dnData, mutData, caseData, controlData)
    bfAll <- rep(1, dim(outData)[1])

    if ( length(pars$gammaMeanDN) == 0) {
        message("No parameters for de novo data")
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
        message("No parameters for case-control data")
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

data22 <- calculateFDR(pars = pars,
                       caseData = caseData, controlData = controlData,
                       dnData = dnData, mutData = mutData,
                       geneName = geneName)
