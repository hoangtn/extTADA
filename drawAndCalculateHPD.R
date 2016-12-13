x <- rnorm(100)
source("script/TADA.R")
source("script/extTADA.R")
fileN <- "DN.adjustBeta.1.3.nI20000.nThin.1.index.15_42_Dec_01_2016.lof_DD.damaging_DD.nGroupDN.2.nGroupCC.2.adjustRatio.1.RData"
fileN <- "DN.adjustBeta.1.3.nI20000.nThin.1.index.15_42_Dec_01_2016.lof_AST.damaging_AST.nGroupDN.2.nGroupCC.2.adjustRatio.1.RData"
fileN <- "GuipponiArrayNoFINcombined.DNandCC.adjustHyperBeta.0..adjustRatioMut.1.3.nIteration20000.nThin.20.index.10_46_24_November.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData"
#load(paste0("../NewTestR/TestR/HBproject/MainPaper/RDataNewDec1/", fileN))
load(paste0("../NewTestR/TestR/HBproject/MainPaper/RDataNew1/", fileN))

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
