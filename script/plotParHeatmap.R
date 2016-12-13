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

