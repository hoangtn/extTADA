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
