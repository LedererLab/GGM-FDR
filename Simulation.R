library(base) # used to symmetrize matrices
library(MASS) # used to generate normally distributed observations
library(Matrix) # used to symmetrize matrices
library(spcov) # used to generate target correlation matrix
library(matrixcalc) # used to symmetrize matrices
library(compare)
library(glasso)
library(ppcor)
library(plyr)
library(graphics)
library(huge)

p <- 200
n <- 300
model <-  huge.generator( n = n, d = p, graph="band", vis = FALSE, verbose = FALSE, g = 4, v = 151, u = 11) #generate band graph
model$sparsity #sparsity level
Sigma <- model$sigma #true covariance matrix 
kappa(Sigma) #condition number
adj <- model$theta #true adjacency matrix
trueConnections <- which(adj>0) # find the true connections
trueNoneConnections <- which(adj==0)


rho.list <- c(seq(0.001, 0.1, 0.005), seq(0.105,1.25,0.005))
mu <- rep(0, p)
M <- 20
gl.fdp.mat <- gl.pwr.mat <- matrix(0, M, length(rho.list))
mb.or.fdp.mat <- mb.or.pwr.mat <- matrix(0, M, length(rho.list))
mb.and.fdp.mat <- mb.and.pwr.mat <- matrix(0, M, length(rho.list))
ct.fdp.mat <- ct.pwr.mat <- matrix(0, M, length(rho.list))
pt.fdp.mat <- pt.pwr.mat <- matrix(0, M, length(rho.list))
dat <- list()
for(k in 1:M){
  X <- mvrnorm(n, mu, Sigma)
  dat[[k]] <- X
  s <- var(X)
  pcorr <- pcor(X)$estimate
  gl.fdp <- mb.or.fdp <- mb.and.fdp <- ct.fdp <- pt.fdp <- NA
  gl.pwr <- mb.or.pwr <- mb.and.pwr <- ct.pwr <- pt.pwr <- NA
  
  for(i in 1:length(rho.list))
  {
    rho.temp <- rho.list[i]
    gl.adj <- huge(X, lambda = rho.temp, nlambda = 1, method = "glasso", verbose = FALSE)$path[[1]]
    mb.or.adj <- huge(X, lambda = rho.temp, nlambda = 1, method = "mb", sym = "or", verbose = FALSE)$path[[1]]
    mb.and.adj <- huge(X, lambda = rho.temp, nlambda = 1, method = "mb", sym = "and", verbose = FALSE)$path[[1]]
    ct.adj <- huge(X, lambda = rho.temp, nlambda = 1, method = "ct", verbose = FALSE)$path[[1]]
    pt.adj <- huge(pcorr, lambda = rho.temp, nlambda = 1, method = "ct", verbose = FALSE)$path[[1]]
    
    gl.FP <- gl.adj[trueNoneConnections] # false positives
    gl.err <- length(which(gl.FP == 1))
    gl.fdp[i] <- gl.err / max(1, length(which(gl.adj == 1)) )
    gl.pwr[i] <- (length(which(gl.adj == 1)) - gl.err) / max(1, length(which(adj == 1)))
    
    mb.or.FP <- mb.or.adj[trueNoneConnections] # false positives
    mb.or.err <- length(which(mb.or.FP == 1))
    mb.or.fdp[i] <- mb.or.err / max(1, length(which(mb.or.adj == 1)) )
    mb.or.pwr[i] <- (length(which(mb.or.adj == 1)) - mb.or.err) / max(1, length(which(adj == 1)))
    
    mb.and.FP <- mb.and.adj[trueNoneConnections]
    mb.and.err <- length(which(mb.and.FP == 1))
    mb.and.fdp[i] <- mb.and.err / max(1, length(which(mb.and.adj == 1)) )
    mb.and.pwr[i] <- (length(which(mb.and.adj == 1)) - mb.and.err) / max(1, length(which(adj == 1)))
    
    ct.FP <- ct.adj[trueNoneConnections]
    ct.err <- length(which(ct.FP == 1))
    ct.fdp[i] <- ct.err / max(1, length(which(ct.adj == 1)) )
    ct.pwr[i] <- (length(which(ct.adj == 1)) - ct.err) / max(1, length(which(adj == 1)))
    
    pt.FP <- pt.adj[trueNoneConnections]
    pt.err <- length(which(pt.FP == 1))
    pt.fdp[i] <- pt.err / max(1, length(which(pt.adj == 1)) )
    pt.pwr[i] <- (length(which(pt.adj == 1)) - pt.err) / max(1, length(which(adj == 1)))
  }
  
  gl.fdp.mat[k,] <- gl.fdp
  mb.or.fdp.mat[k,] <- mb.or.fdp
  mb.and.fdp.mat[k,] <- mb.and.fdp
  ct.fdp.mat[k,] <- ct.fdp
  pt.fdp.mat[k,] <- pt.fdp
  
  gl.pwr.mat[k, ] <- gl.pwr
  mb.or.pwr.mat[k, ] <- mb.or.pwr
  mb.and.pwr.mat[k, ] <- mb.and.pwr
  ct.pwr.mat[k, ] <- ct.pwr
  pt.pwr.mat[k,] <- pt.pwr
}

source("Knockoff.R")

set.seed(100)
FDRtarget.list <- seq(0.001, 1, 0.01)
ko.fdp.list <- ko.pwr.list <- c()
for(k in 1:length(FDRtarget.list)){
  FDRtarget <- FDRtarget.list[k]
  ko.fdp <- ko.pwr <- NA
  
  for(i in 1:M){
    X <- dat[[i]]
    ko.est <- GraphEstimation(X, FDRtarget, plus = FALSE)
    ko.adj <- ko.est$adjacency.matrix
    
    if(dim(ko.est$edge.set)[1]==0){
      ko.err<-0
    }else{
      ko.FP <- ko.adj[trueNoneConnections] # false positives
      ko.err <- length(which(ko.FP == 1))
    }
    ko.fdp[i] <- ko.err/ max(1, length(which(ko.adj == 1)) )
    ko.pwr[i] <- (length(which(ko.adj == 1)) - ko.err) / max(1, length(which(adj == 1)))
  }
  ko.fdp.list[k] <- mean(ko.fdp)
  ko.pwr.list[k] <- mean(ko.pwr)
}



### making plots

par(mfrow=c(1,2),mai = c(0.6, 0.6, 0.4, 0.2))

par(xpd = FALSE)
plot(NULL, NULL, xlim = c(0,1),cex = 0.5,xlab="", ylab = "",cex.axis = 0.8, type = "n", ylim = c(0,1), las = 1)
points(colMeans(gl.fdp.mat), colMeans(gl.pwr.mat), xlim = c(1,0),pch = 15, cex = 0.5, col = "darkorchid")
points(colMeans(mb.and.fdp.mat), colMeans(mb.and.pwr.mat), xlim = c(1,0), pch = 16, cex = 0.5, col = "orange")
points(colMeans(mb.or.fdp.mat), colMeans(mb.or.pwr.mat), xlim = c(1,0), pch = 17, cex = 0.5, col = "deepskyblue")
points(colMeans(ct.fdp.mat), colMeans(ct.pwr.mat), xlim=c(1,0), pch = 5,cex = 0.4, col = "deeppink")
points(colMeans(pt.fdp.mat), colMeans(pt.pwr.mat), xlim=c(1,0), pch = 6,cex = 0.4, col = "green3")
points(ko.fdp.list, ko.pwr.list, col = "firebrick1", pch = 8, cex = 0.4)
title(ylab = "Power", xlab = "FDR", line = 2)
grid(col = "darkgray")
title(main = "n=300, p=200", cex.main = 1.2, line = 0.6, font.main = 2)
legend(0.66, 0.5, legend = c("KO","GLASSO","MB(and)","MB(or)","CT","PT"),  col = c("firebrick1","darkorchid","orange","deepskyblue","deeppink","green3"),
       pch = c(8,seq(15,17),5,6), bty="n", y.intersp = 0.75, cex = 0.5)


plot(NULL, NULL, xlim = c(0,1), pch = 20, cex = 0.5, xlab = "",ylab = "", cex.axis = 0.8, type="n", ylim = c(0,1), las = 1)
points(FDRtarget.list, ko.fdp.list, xlim = c(0,1), pch = 20, cex = 0.5, xlab = "", ylab = "", col="firebrick1")
abline(a = 0, b = 1)
title(ylab = "Actual FDR",xlab = "Target FDR", cex.lab = 0.9, line = 2)
title(main = "n=300, p=200", cex.main = 1.2, line=0.6, font.main=2)

dev.off()