estimatePars <- function(pars, mcmcResult, nThin = NULL){
    mcmcDataFrame <- as.data.frame(mcmcResult)
    pars <- pars[grep("hyper|pi", pars)]
    message("====\nOnly pi and hyper parameters are estimated in this step\n",
            "extTADA does not calculate HPDs for hyper betas, just their medians\n===\n")

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

    colnames(hpdList) <- c("Mode", "lCI", "uCI")

return(hpdList)
}
