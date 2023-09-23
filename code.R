library(IMIFA)
library(fields)
library(mvtnorm)
library(msm)
library(SphericalCubature)
library(kernlab)
library(MASS)
library(huge)
library(mvtnorm)
library(MASS)
library(pracma)
library(expm)

source('Latentfactor.R')
source('multiRBFjointP2MHhigh.R')
source("otherdependent.R")


load("data.rda")
#Xred is the matrix of dimension 100 \times 207 with transformed predictor vectors for screened set of regions.
#     the columns with indices 1:9 + 9*(i-1) is associated with the region named, ROIname[i]
#     the nodal attributes are according to the "nodalname" variable
#datasca is the mean-centered and normalized data on cognitive scores of dimension 7 \times 100.
#Ytr is the training data matrix with 30% entries missing of dimension 7 \times 100.
#Yte is the test data matrix with missing entries at training locations of dimension 7 \times 100.
#ROInames have the screened ROI names.
#nodalname have the names of nodal atttributes.
#tasknames have the names of the cognitive tasks i.e. rownames of datasca.

runcode = T

if(runcode){
  #Fitting unsupervided Latent factor model
  fitPFA <- LFMPFA(datasca, d=100)
  #Fitting the supervised model with training data
  fit <- MultiRBFjointHigh(Ytr=Ytr, Xtr=Xred, Xte=Xred, cutpval=0.01, grpindex=grpindex, pindexprior=T, pindexGprior=F)
  
  #Fitting the supervised model with the complete data
  fitFsca <- MultiRBFjointHigh(Ytr=datasca, Xtr=Xred, Xte=NULL, cutpval=0.01, grpindex=grpindex, pindexprior=T, pindexGprior=F)
}

# if(!runcode){
#   load("outputs.rda")
# }


tasknames <- c("Oral Reading", "Picture Vocabulary", "Processing Speed", "Working Memory", "Episodic Memory", "Executive Function", "Flanker Task")

#Unsupervised case

lambdapsam <- array(0, dim=c(7,7, length(fitPFA$lambdals)))#array(unlist(lambdapp), dim = c(dim(lambdap), 1500))

for(i in 1:length(fitPFA$lambdals)){
  sigma2p <- fitPFA$sigma2ls[[i]]
  lambdapsam[,,i] <- fitPFA$lambdals[[i]] %*% diag(sigma2p)
}

lambdapp <- t(apply(lambdapsam, 1:2, mean))

lambdapOP <- sp_OP(lambdapsam, itermax = 10)
lambdapp  <- t(lambdapOP$Phi)

colnames(lambdapp) <- tasknames

out <- plotload(t(lambdapp[1:4, ])); Infolam <- out$info


################Figure 1##################
par(mfrow=c(1, 1), mar = c(4.15, 8, 2.15, 5))
image(Infolam$x, Infolam$y, lambdapp[1:4, ], col = out$col, breaks = out$brk, main = "", yaxt="n", xlab="Factors", ylab="")
axis(2, at=seq_along(tasknames),labels=as.character(tasknames[1:7]), las=2)
image.plot(Infolam$x, Infolam$y, lambdapp[1:4, ], col = out$col, breaks = out$brk, legend.only = T)



#####################Supervised case

G <- length(names)
grpindex <- list()
for(i in 1:G){
  grpindex[[i]] <- 1:9 + 9*(i-1)
}

lenlist <- list()
for(i in 1:G){
  lenlist[[i]] <- length(grpindex[[i]])
}

lambdapsam <- array(0, dim=c(7,7, length(fitFsca$lambdap)))#array(unlist(lambdapp), dim = c(dim(lambdap), 1500))

fitU <- fitFsca

for(i in 1:length(fitU$lambdap)){
  sigma2p <- fitU$sigma2ls[[i]]
  lambdapsam[,,i] <- fitU$lambdap[[i]] %*% diag(sigma2p)
}

lambdapp <- t(apply(lambdapsam, 1:2, mean))

colnames(lambdapp) <- tasknames

lambdapOP <- sp_OP(lambdapsam, itermax = 10)
lambdapp  <- t(lambdapOP$Phi)

out <- plotload(t(lambdapp[1:4, ])); Infolam <- out$info

################Figure 5##################
par(mfrow=c(1, 1), mar = c(4.15, 8, 2.15, 5))
image(Infolam$x, Infolam$y, (lambdapp[1:4, ]), col = out$col, breaks = out$brk, main = "", yaxt="n", xlab="Factors", ylab="")
axis(2, at=seq_along(tasknames),labels=as.character(tasknames[1:7]), las=2)
image.plot(Infolam$x, Infolam$y, (lambdapp[1:4, ]), col = out$col, breaks = out$brk, legend.only = T)


#################################Prediction MSE due to our proposed method in Table 2################
Ytp <- Reduce('+', fit$Ytep)/length(fit$Ytep)
signif(mean(rowMeans((Ytp-Yte)^2, na.rm=T)), digits = 2)
