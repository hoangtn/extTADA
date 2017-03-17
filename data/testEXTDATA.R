 fileR <- dir("../script", ".R$")
for (ii in fileR)
    source(paste0("../script/", ii))
data <- read.table("../data/data_mut_DD.csv", header = TRUE, as.is = TRUE)
head(data)

allDNData <- data[, paste0("dn_", c("damaging", "lof"), "_DD")]
allMutData <- data[,paste0("mut_", c("damaging", "lof"))]
head(data.frame(allMutData, allDNData))

mcmcDD <- extTADA(modelName = DNextTADA, #extTADA for only de novo data
                  dataDN = allDNData, mutRate = allMutData,
                  Ndn = array(rep(4293, 2)),
                  nIteration = 1000,
                  nIteration2 = 2000)
