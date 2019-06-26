library(heatmap3)
library(gplots)
library(plotrix)
library(matrixcalc)
library(compositions)
source("Knockoff.R")

f <- "RealDataAnalysis/Data/AmericanGut/07-taxa/100nt/ag-cleaned_L2.txt"
data.raw <- read.table(f , sep = "\t", header = TRUE, fill = TRUE, comment.char="", quote = "")
bacteria.indx <- which(startsWith(colnames(data.raw), "k__"))
bacteria <- apply(data.raw[,bacteria.indx], 2, as.numeric)
names.temp <- colnames(bacteria)
for (i in 1:length(names.temp)) 
{
  name <- strsplit(names.temp[i], "_")[[1]]
  name <- name[length(name)]
  if(startsWith(name, "."))
  {
    colnames(bacteria)[i] <- substr(name, 2, nchar(name) - 1)
  }else{
    colnames(bacteria)[i] <- name
  }
}
dat <- as.matrix(as.data.frame(bacteria))
n <- dim(dat)[1]  # n=17348 number of observations
p <- dim(dat)[2]  # p=67 number of types of bacteria

# select the bacterias which has # non-zero observations > 0.5%*n
ratio <- 0.005
del <- c()
for(i in 1:p)
{
  if(sum(dat[,i] != 0) < n*ratio){
    del <- c(del,i)
  }
}
dat <- dat[,-c(del)] # n=17348  p=32
dat <- clr(dat) #CLT transformation
dat.smoke <- dat[which(data.raw$SMOKING_FREQUENCY=="Daily"
                       |data.raw$SMOKING_FREQUENCY=="Occasionally (1-2 times/week)"
                       |data.raw$SMOKING_FREQUENCY=="Rarely (a few times/month)"
                       |data.raw$SMOKING_FREQUENCY == "Regularly (3-5 times/week)" ), ]
dat.no.smoke <- dat[which(data.raw$SMOKING_FREQUENCY=="Never"), ]


###non-smoker
set.seed(3)
temp.results <- matrix(0, dim(dat)[2], dim(dat)[2])
for(i in 1:10)
{
  dat.temp <- dat.no.smoke[sample(dim(dat.no.smoke)[1], 1234), ]
  weights <- (1/2)^(i)
  FDRtarget <- 0.2*weights
  edge.mat.temp <- GraphEstimation(dat.temp, FDRtarget, FALSE)$continuous.thresh
  temp.results <- temp.results + abs(edge.mat.temp)
}
temp.results.nonsmoker <- temp.results/max(temp.results)


###smoker
set.seed(5)
temp.results1 <- GraphEstimation(dat.smoke, 0.2, FALSE)$continuous.thresh
temp.results1 <- abs(temp.results1)
temp.results.smoker <- temp.results1/max(temp.results1)


### making plots

# non-smoker
diag(temp.results.nonsmoker) <- NA
rownames(temp.results.nonsmoker) <- colnames(dat)
colnames(temp.results.nonsmoker) <- colnames(dat)
par(oma =  c(2, 0.1, 1.5, 5.8))
heatmap.2(abs(temp.results.nonsmoker),
  Rowv=FALSE, Colv=FALSE, 
  dendrogram = "none", trace=  'none', density.info = "none",
  cexRow = 0.55, srtRow = 330, offsetRow = 0.06, 
  cexCol = 0.55, srtCol = 60, offsetCol = 0.02,
  keysize = 0.2, key.xlab = "",key.title = "",
  key.par = list(mar=c(2.9,0.1,4.5,1.5), cex = 0.5),
  lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(1, 4.5), lwid = c(1, 10, 1),
  col = gray((100:0)/100), na.color = "tomato", margins = c(3.5, 0))
mtext(text = "Non-smoker", side = 1, line = -15.3, adj = 0.5)

dev.off()

# smoker
diag(temp.results.smoker) <- NA
rownames(temp.results.smoker) <- colnames(dat)
colnames(temp.results.smoker) <- colnames(dat)

par(oma =  c(2, 0.1, 1.5, 5.8))
heatmap.2(abs(temp.results.smoker),
  Rowv = FALSE, Colv = FALSE, 
  dendrogram = "none", trace = 'none', density.info = "none",
  cexRow = 0.55, srtRow = 330,offsetRow = 0.06, 
  cexCol = 0.55, srtCol = 60, offsetCol = 0.02,
  keysize = 0.2, key.xlab = "",key.title = "",
  key.par = list(mar = c(2.9,0.1,4.5,1.5), cex = 0.5),
  lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(1, 4.5), lwid = c(1, 10, 1),
  col = gray((100:0)/100), na.color = "tomato", margins = c(3.5, 0))
mtext(text="Smoker", side = 1, line = -15.3, adj = 0.5)

dev.off()

#making histograms and boxplot 

temp.smoker <- temp.results.smoker[upper.tri(temp.results.smoker)]
temp.non.smoker <- temp.results.nonsmoker[upper.tri(temp.results.nonsmoker)]
temp.list <- list(temp.non.smoker, temp.smoker)

par(mfrow = c(1, 2), mai = c(0.60, 0.65, 0.25, 0.2))

# histograms
l <- lapply(temp.list, hist, plot = FALSE)
mids <- unique(unlist(lapply(l, function(x)x$mids)))
counts <- lapply(l, function(x)x$counts[match(x = mids, table = x$mids, nomatch = NA)])
breaks <- unique(unlist(lapply(l, function(x)x$breaks)))
h <- do.call(rbind, counts)
barplot(h, beside = TRUE, col = c("lightskyblue", "indianred1"), xaxt = "n",cex.axis = 0.75)
atX <- axTicks(side = 1, axp = c(par("xaxp")[1:2], length(breaks)-1))
labels <- seq(min(breaks), max(breaks), length.out = 1 + par("xaxp")[3])
labels <- round(labels, digits = 1)
axis(side = 1, at = atX, labels = breaks, cex.axis = 0.75)
title(ylab = "Frequency", xlab = "Non-zero signal strengths", line = 2, cex.lab = 0.75)
legend(16, 330, c("Non-smoker", "Smoker"), 
       fill = c("lightskyblue", "indianred1"), cex = 0.6, bty = "n")

# boxplot
names(temp.list) <- c("non-smoker","smoker")
boxplot(temp.list, col = c("lightskyblue", "indianred1"), cex.axis = 0.75)
title(ylab = "Signal strengths", line = 2.2, cex.lab = 0.75)

dev.off()

