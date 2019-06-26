library(heatmap3)
library(gplots)

#grid accuracy: 1e5

source("Knockoff.R")

regions <- read.csv("RealDataAnalysis/Data/BrainConnectivity/AAL_YAN.csv", head = TRUE, sep = ",") 
brain <-read.csv("RealDataAnalysis/Data/BrainConnectivity/BrainRegions.csv", sep = ";", head = FALSE)
# Removes commas in column 6
brain[,6] <- gsub(",", "", brain[,6], fixed = TRUE)
regnames <- brain[,3][regions$no]
p <- 116
FDRtarget <- 0.2

###NC, 10 patients
temp.nc <- matrix(0, p, p)
# set.seed(4)
# sam <- sample(28:37, 10)
set.seed(90)
for(i in 28:37)
{
  weights <- (1/2)^(i-27)
  FDRtarget <- 0.2*weights
  # index <- sam[i-27]
  dat.temp <- readRDS(file= "RealDataAnalysis/Data/BrainConnectivity/detrendedData")[ , i, ]
  edge.mat.temp <- GraphEstimation(dat.temp, FDRtarget, FALSE)$continuous.thresh
  temp.nc <- temp.nc + abs(edge.mat.temp)/10
}
temp.nc <- temp.nc/max(temp.nc)
diag(temp.nc) <- NA
rownames(temp.nc) <- regnames
colnames(temp.nc) <- regnames


ind <- which(regions$X42region == "Y")
oma.default <- c(0, 0, 0, 0)
par(oma = oma.default + c(5, 4, 3, 6.5))
heatmap.2(temp.nc[ind, ind], Rowv = FALSE, Colv = FALSE, 
  dendrogram = "none", trace = 'none', density.info = "none",
  cexRow = 0.47, srtRow = 330, offsetRow = 0.06, 
  cexCol = 0.47, srtCol = 60, offsetCol = 0.02,
  keysize = 0.2, key.xlab = "",key.title = "",
  key.par = list(mar=c(3.3, 1, 3.35, 1), cex = 0.5),
  lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(1, 4.5), lwid = c(1, 10, 1),
  col = gray((100:0)/100), na.color="tomato", margins=c(3.5, 0),
  add.expr=rect(c(0.6, 12.5, 20.5, 26.5),
                c(30.5, 22.5, 16.5, 0.6),
                c(12.5, 20.5, 26.5, 42.4),
                c(42.4, 30.5, 22.5, 16.5), 
                border = "tomato", lwd = 1.1))
mtext(text = "Frontal lobe", side = 2, line = 4.5,las = 1,cex = 0.5, padj = -10)
mtext(text = "Parietal lobe", side = 2, line = 4.5,las = 1,cex = 0.5, adj = 1)
mtext(text = "Occipital lobe", side = 2, line = 4.5,las = 1,cex =0.5, padj = 6)
mtext(text = "Temporal lobe", side = 2, line = 4.5,las = 1, cex = 0.5, padj = 16)
mtext(text = "NC", side = 1,line = -12.3, adj = 0.45)



