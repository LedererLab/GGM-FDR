library(heatmap3)
source("Knockoff.R")
dat.sites <- read.csv("RealDataAnalysis/Data/AtlanticAmphibians/ATLANTIC_AMPHIBIANS_sites.csv", header = TRUE)
dat.species <- read.csv("RealDataAnalysis/Data/AtlanticAmphibians/ATLANTIC_AMPHIBIANS_species.csv", header = TRUE)
dat.species <- dat.species[which(dat.species$order == "Anura"),]

temp <- dat.species[, c(1,8,9,10)]
#1 "id":id for study sites
#8 "valid_name" : valid name of species 
#9 "individuals": number of counts
#10 "endemic" 

dat.raw <- temp[complete.cases(temp), ] 
dat.raw.end <- dat.raw[which(dat.raw$endemic == 1), ]   
dat.raw.no.end <- dat.raw[which(dat.raw$endemic == 0), ] 
end.spec <- dat.raw.end$valid_name
nonend.spec <- dat.raw.no.end$valid_name


####endemic
sites <- dat.sites[which(dat.sites$record=="ab"), ]$id
num.obs <- length(sites) # n=346
#sort(table(dat.raw.end$valid_name),decreasing=TRUE)[1:30] 
num.spec <- 30
new.table <- data.frame(matrix(nrow = num.obs,ncol = num.spec))
rownames(new.table) <- as.vector(sites)
colnames(new.table) <- as.vector(as.data.frame(sort(table(dat.raw.end$valid_name), decreasing = TRUE)[1:num.spec])$Var1)
spec.names <- colnames(new.table)
for(i in 1:num.obs)
{
  for(j in 1:num.spec)
  {
    temp.row <- dat.raw[which(dat.raw$id == rownames(new.table)[i]
                              &dat.raw$valid_name == spec.names[j]),]
    if(nrow(temp.row) == 2){
      print(i)
      print(j)
      break
    }
    if(nrow(temp.row) == 0)
    {
      new.table[i,j] <- 0
    }else{
      new.table[i,j] <- temp.row$individuals
    }
  }
}

trans.data <- clr(new.table)
end.dat <- trans.data

set.seed(19)
temp.results.end <-  GraphEstimation(end.dat, 0.2, FALSE)$continuous.thresh
temp.results.end <- abs(temp.results.end)/max(abs(temp.results.end))

# making plot
diag(temp.results.end) <- NA
rownames(temp.results.end) <- colnames(end.dat)
colnames(temp.results.end) <- colnames(end.dat)
par(oma = c(3.5, 3, 5, 5.5))
heatmap.2(temp.results.end,
          Rowv = FALSE, Colv = FALSE, 
          dendrogram = "none",trace = 'none', density.info = "none",
          cexRow = 0.5, srtRow = 330, offsetRow = 0.06, 
          cexCol = 0.5, srtCol = 60, offsetCol = 0.02, 
          keysize = 0.2, key.xlab = "", key.title = "",
          key.par = list(mar = c(3.5,1,3,1), cex = 0.5),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(1, 4.5), lwid = c(1, 10, 1),
          col = gray((100:0)/100), na.color = "tomato", margins = c(3.5, 0))
mtext(text = "Endemic species", side = 1, line = -13, adj = 0.35)

dev.off()


####non-endemic

num.obs <- length(sites)
#sort(table(dat.raw.no.end$valid_name),decreasing=TRUE)[1:30] #>=1%
num.spec <- 30
new.table <- data.frame(matrix(nrow = num.obs, ncol = num.spec))
rownames(new.table) <- as.vector(sites)
colnames(new.table) <-as.vector(as.data.frame(sort(table(dat.raw.no.end$valid_name), decreasing = TRUE)[1:num.spec])$Var1)
spec.names<- colnames(new.table)

for(i in 1:num.obs)
{
  for(j in 1:num.spec)
  {
    temp.row <- dat.raw[which(dat.raw$id == rownames(new.table)[i]
                              &dat.raw$valid_name == spec.names[j]),]
    if(nrow(temp.row) == 2){
      print(i)
      print(j)
      break
    }
    if(nrow(temp.row) == 0)
    {
      new.table[i,j] <- 0
    }else{
      new.table[i,j] <- temp.row$individuals
    }
  }
}
trans.data <- clr(new.table)
noend.dat <- trans.data

set.seed(50)
results <- GraphEstimation(noend.dat, 0.2, FALSE)
temp.results <-  results$continuous.thresh
temp.results.nonend <- abs(temp.results)/max(abs(temp.results))

# making plot
diag(temp.results.nonend) <-NA
rownames(temp.results.nonend) <- colnames(noend.dat)
colnames(temp.results.nonend) <- colnames(noend.dat)
par(oma =  c(3.5, 3, 5, 5.5))
heatmap.2(temp.results.nonend,
          Rowv = FALSE, Colv = FALSE, 
          dendrogram = "none", trace = 'none', density.info = "none",
          cexRow = 0.5, srtRow = 330,offsetRow = 0.06, 
          cexCol = 0.5, srtCol = 60, offsetCol = 0.02,
          keysize = 0.2, key.xlab = "",key.title = "",
          key.par = list(mar = c(3.5,1,3,1), cex = 0.5),
          lmat = rbind(c(5, 4, 2), c(6, 1, 3)), lhei = c(1, 4.5), lwid=c(1, 10, 1),
          col = gray((100:0)/100),
          na.color = "tomato", margins = c(3.5, 0))
mtext(text = "Non-endemic species", side = 1,line = -13, adj = 0.35)

dev.off()

