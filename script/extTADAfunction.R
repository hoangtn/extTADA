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
         message("\nThere are ", Ngene, " genes in this analysis")

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
         if ((adjustHyperBeta == 0) & is.null( hyperBetaCC0) & !is.null(dataCCcase))
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

         message("\n=============FIRST TIME ===============\n")
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

             message("\n=============SECOND TIME ===============\n")
             message("\nThis process is using beta values estimated from the previous process\n")
         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
         }



         return(mcmcModel)
          }

