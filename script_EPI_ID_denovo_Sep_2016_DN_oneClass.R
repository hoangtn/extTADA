args <- commandArgs(TRUE)
nUpdate <- as.numeric(args[1])
nIteration <-   as.numeric(args[2])
inDex <- as.character(args[3])
nThin <- as.numeric(args[4])
nChain <- as.numeric(args[5])
annotationType <- as.character(args[6])
annotationType2 <- as.character(args[7])
outDir <- as.character(args[9])
lowerHyperGamma <- as.numeric(args[10])
annotationType3 <- as.character(args[11])
resultDir <- as.character(args[12])
nGroupDN <- as.numeric(args[13])
nGroupCC <- as.numeric(args[14])
nCore <-  nChain ###cores = chains
adjustHyperBeta <- as.numeric(args[15])
swapData <- as.numeric(args[16])

dn1 <- unlist(strsplit(annotationType, "_", fixed = TRUE))[1]
dn2 <- unlist(strsplit(annotationType2, "_", fixed = TRUE))[1]
dn3 <- unlist(strsplit(annotationType3, "_", fixed = TRUE))[1]

print(annotationType)
print(annotationType2)
print(annotationType3)


library("rstan")
library("coda")

#data <- read.table("/hpc/users/nguyet26/psychen/methods/extTADA/data/scz_data_Sep_2016_exac_noexac_oldDataBases_and_allSilentFCPK.txt",
 #                  header = TRUE)
data <- read.table("/hpc/users/nguyet26/psychen/methods/extTADA/Re_annotate/DenovoData/epi_id_dd_scz_data_Sep_2016_addNewID.txt",
                   header = TRUE, as.is = TRUE, sep = ",")
data2 <- read.table("~/psychen/methods/extTADA/scripts/TADA/data/ASC_2231trios_1333trans_1601cases_5397controls.csv", sep = ",", header = TRUE)

t1 <- pmatch(data2[, 1], data[, 1])
data3 <- data.frame(data2, index = t1)
data3 <- data3[!is.na(data3$index),]

data[, "dn_lof_AST_TADA"] <- rep(0, dim(data)[1])
data[, "dn_lof_AST_TADA"][data3[, "index"]] <- data3[, "dn.LoF"]

data[, "dn_damaging_AST_TADA"] <- rep(0, dim(data)[1])
data[, "dn_damaging_AST_TADA"][data3[, "index"]] <- data3[, "dn.mis3"]

data[, "dn_missense_AST_TADA"] <- rep(0, dim(data)[1])
data[, "dn_missense_AST_TADA"][data3[, "index"]] <- data3[, "dn.mis3"]

dim(data)
#data <- read.table("/hpc/users/nguyet26/psychen/methods/extTADA/data/scz_data_April_2016_exac_noexac_oldDataBases_and_allSilentFCPK.txt",
 #                  header = TRUE)

silentData <- read.table("~/psychen/resources/genic_mutation_w_annotation/result_mutation_gene_silent.txt",  header = FALSE, as.is = TRUE)
sumSilent <- sum(silentData[, 2])
sA <- pmatch(silentData[, 1], data[, 1])
length(sA); sA <- sA[!is.na(sA)]; length(sA)

sum(silentData[sA,][, 2])

allData <- read.table("~/psychen/methods/TADA_runs/scz.dat", header = TRUE, as.is = TRUE)
nTrioMenachem=617
apply(allData[, -1], 2, sum)*2*nTrioMenachem

sum(silentData[, 2])*2*nTrioMenachem #n = 207.0952
#nSilent = 227; nSilentOfFromer = 156
#wc -l /hpc/users/nguyet26/psychen/resources/denovos/data/dnenrich_inputs/SCZ.silent.mut

#data <- read.table("/hpc/users/nguyet26/psychen/methods/extTADA/data/scz_data_April_2016_exac_noexac_oldDataBases_and_allSilentFCPK.txt",
 #                  header = TRUE)

adjustRatioMut <- 156/(sum(silentData[, 2])*2*nTrioMenachem)


print(annotationType)
print(annotationType2)
print(annotationType3)
#.damaging_SingletonTransmittedDataAddCC_noexac
dim(data)
#data <- data[data$mut_lof > 0, ]
t1 <- as.numeric(data$mut_lof)
t1 <- t1[t1 > 0]
data[data$mut_lof == 0, ]$mut_lof <- min(t1)/10
#data <- data[data$mut_missense >0, ]
t2 <- as.numeric(data$mut_missense)
t2 <- t2[t2 > 0]
data[data$mut_missense == 0, ]$mut_missense <- min(t2)/10
dim(data)

#data <- data[data$mut_missense >0, ]
t2 <- as.numeric(data$mut_damaging)
t2 <- t2[t2 > 0]
data[data$mut_damaging == 0, ]$mut_damaging <- min(t2)/10

##Silent
#t3 <- as.numeric(data$mut_silentCFPK)
#t3 <- t3[t3 > 0]
#data[data$mut_silentCFPK == 0, ]$mut_silentCFPK <- min(t3)/10


dim(data)


ntrio <- 617 #Menachem
ntrio <- ntrio + 14 #Girard
ntrio <- ntrio + 105 #Gulsuner
ntrio <- ntrio + 57 #McCarthy
ntrio <- ntrio + 231 #Xu



ntrioID <- 100 #deLight
ntrioID <- ntrioID + 0 #Gilissen
ntrioID <- ntrioID + 41 #Hamdan
ntrioID <- ntrioID + 51 #Rauch
ntrioID <- ntrioID + 820 ##Add Lei: 461 boys + 359 girls



ntrioEPI <- 0 #Epi4k
ntrioEPI <- ntrioEPI + 356 #EuroEPINOMICS-RES

ntrioDD <- 4293 #http://biorxiv.org/content/biorxiv/early/2016/06/16/049056.full.pdf
ntrioSCZ = 1024
ntrioAST = 3985
ntrioAST_SSC <- 2508
ntrioAST_DeRubeis2014 <- 2270
ntrioAST_TADA <- 2231


if (annotationType == "lof_CHD_allN1575") {
    ntrio <- 1575
    adjustRatioMut <- sum(data$dn_silent_CHD_allN1575)/(sum(silentData[sA,][, 2])*2*ntrio) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}

if (annotationType == "lof_AST_allN13303") {
    ntrio <- 5122 #13303
    adjustRatioMut <- sum(data$dn_silent_AST_allN13303)/(sum(silentData[sA,][, 2])*2*ntrio) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}

if (annotationType == "lof_AST_TADA") {
    ntrio <- ntrioAST_TADA
    adjustRatioMut <- sum(data$dn_silent_AST_DeRubeis2014)/(sum(silentData[sA,][, 2])*2*ntrio) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}


if (annotationType == "lof_AST_DeRubeis2014") {
    ntrio <- ntrioAST_DeRubeis2014
    adjustRatioMut <- sum(data$dn_silent_AST_DeRubeis2014)/(sum(silentData[sA,][, 2])*2*ntrio) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}


if (annotationType == "lof_AST_SSC") {
    ntrio <- ntrioAST_SSC
    adjustRatioMut <- sum(data$dn_silent_AST_SSC)/(sum(silentData[sA,][, 2])*2*ntrio) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}


if (annotationType == "lof_AST") {
    ntrio <- ntrioAST
    adjustRatioMut <- sum(data$dn_silent_AST)/(sum(silentData[sA,][, 2])*2*ntrioAST) #1103/1204 #1133/(2*ntrioAST*sumSilent)
}
    
if (annotationType == "lof_EPI") {
    ntrio <- ntrioEPI
    adjustRatioMut <- 91/(2*ntrioEPI*sumSilent)
}
if (annotationType == "lof_ID") {
    ntrio <- ntrioID
    adjustRatioMut <- (51 + 213)/(2*ntrioID*sumSilent)
}
if (annotationType == "lof_EPIandID")
    ntrio <- ntrioEPI + ntrioID

if (annotationType == "lof_DD")
    ntrio <- ntrioDD

if (annotationType == "lof_SCZ") {
    ntrio <- ntrioSCZ
    adjustRatioMut <- 263/(2*ntrioSCZ*sumSilent)
}
N <- list(dn=ntrio)
data <- data[data$mut_lof != 0, ]
dim(data)

#counts <- cbind(sCountCC2[, 1], sCountCC2[, 2:6])


if (dn1 == "disruptive")
	dn1 <- "lof"

allDNData <- data[, paste0("dn_", c(annotationType2, annotationType))][, 2]
allMutData <- data[,paste0("mut_", c(dn2, dn1))][, 2]

mutationLoF <- paste0("mut_", dn1)
mutationMis3 <- paste0("mut_", dn2)
mutationCFPK <- paste0("mut_", dn3)
#dn1 <- paste(dn1, "_noexac", sep = "")
#dn2 <- paste(dn2, "_noexac", sep = "")


print(dn1)
head(data, 2)
dnLoF <- paste("dn_", annotationType, sep = "")
#mutationLoF <- paste("mut_", dn1, sep = "")


#length(protectLoF[protectLoF < 1])/length(y.case.lof)

yLoF <- data[, dnLoF]
mutLoF <- data[, mutationLoF]


dnMis3 <- paste("dn_",  annotationType2, sep = "")
#mutationMis3 <- paste("mut_", dn2, sep = "")

yMis3 <- data[, dnMis3]
mutMis3 <- data[, mutationMis3]

##Silent

dnCFPK <- paste("dn_",  annotationType3, sep = "")

yCFPK <- data[, dnCFPK]
mutCFPK <- data[, mutationCFPK]

message("Number of categories: ", nGroupDN, " DN and ", nGroupCC, " CC")

adjustRatioMut <- 0.8719 #1

allMutData <- adjustRatioMut*allMutData
mutRateArray <-  adjustRatioMut*cbind(mutCFPK, mutMis3, mutLoF)
dataDNArray <-  as.matrix(cbind(yCFPK, yMis3, yLoF))


lDN <- ifelse(nGroupDN == 3, 1,
              ifelse(nGroupDN == 2, 2, 3))
lCC <- ifelse(nGroupCC == 3, 1,
              ifelse(nGroupCC == 2, 2, 3))

swapColumn <- function(x, c1 = 2, c2 = 3){
    tempC <- x[, c1]
    x[, c1] <- x[, c2]
    x[, c2] <- tempC
    return(x)
}

if (swapData == 1){
    dataDNArray <- swapColumn(dataDNArray)
    mutRateArray<- swapColumn(mutRateArray)
    
}
    

adjustHyperBeta0 <- 0
adjustHyperBeta0 <- adjustHyperBeta

nGroupDN = 1

mixDataKclasses <- list(K = 2, Kcc = nGroupCC, Kdn = nGroupDN,  hyperBetaMax = 20,
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = length(yLoF),
                                                Ndn = N$dn, 
                                                dataDN = data.frame(allDNData),
                                                mutRate = data.frame(allMutData),
                                                upperPi0 = 0.5, lowerPi0 = 0.001,
                                                lowerHyperGamma = 1, lowerBeta = 1,
                                                hyperBetaDN0 = array(c(rep(10, nGroupDN -1), 1)),
                                                hyperBetaCC0 = array(rep(10, nGroupCC)),
                        betaPars = c(6.7771073, -1.7950864, -0.2168248))



nCore <- nChain
source("U_publishModel/extTADA.R")

testIntegratedModel <- stan(model_code = DNextTADA,
                        data = mixDataKclasses, iter = nIteration,
                                                chains = nChain, cores = nCore,
                                                thin = nThin)

save.image(paste0(resultDir, "/DNlof.adjustBeta.", adjustHyperBeta, ".",
                         nChain, ".nI", nIteration, ".nThin.",
        nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".adjustRatio.", round(adjustRatioMut, 2),
        ".RData"))


b1 <- as.data.frame(testIntegratedModel)
bMCMC <- mcmc(b1)
bHPD <-  HPDinterval(bMCMC)
medianHPD <- NULL


for (ii in 1:dim(bHPD)[1]){
    t1 <- b1[, ii]

    t2 <- t1[(t1>=bHPD[ii, 1]) & (t1<=bHPD[ii, 2])]

    medianHPD[ii] <- mean(t2)

     }

names(medianHPD) <- colnames(b1)

t3 <- c(medianHPD)

write.table(t(t3), paste0(resultDir, "/MCMC.", ntrio,
                         nChain, ".nIteration", nIteration, ".nThin.", 
        nThin, ".index.", inDex, ".adjustHyperBeta.", adjustHyperBeta,
        ".", annotationType, ".", annotationType2, ".nGroup.", nGroupDN,
        ".nGroupCC.", nGroupCC,
        ".notChangeDNgammaPrior.txt"),
         col.names = FALSE, quote= FALSE, row.names = FALSE)

save.image(paste0(resultDir, "/DNlof.adjustBeta.", adjustHyperBeta, ".",
                         nChain, ".nI", nIteration, ".nThin.",
        nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".adjustRatio.", round(adjustRatioMut, 2),
        ".RData"))

