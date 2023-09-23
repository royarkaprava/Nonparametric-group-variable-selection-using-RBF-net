#Xred is the matrix of dimension n \times p with transformed predictor vectors for screened set of regions.
#     the columns with indices 1:mg + mg*(i-1) is associated with the i-th region.
#    
#Xte is the test data matrix
#Ytr is the training data matrix with NA's at the missing entries of dimension v \times n or at the test places.
#The Ytep variable in the the output of this function will be of same dimention as Ytr with NA's replaced by predicted values.

source("Latentfactor.R") #It is needed if runFA = T

MultiRBFjointHigh <- function(Ytr, Xtr, K=NULL, Xte = NULL, corel=T,lasso=F, runFA=F, noupR=1000, noupRL=1000, grpindex, indz=NULL, fitL=NULL, updateLow=F, Shr = T, cutpval=0.01, pindexprior=F,pindexGprior=F, pindex=NA, pindexG=NA, Total_itr = 10000, burn = 5000){
  
  set.seed(1) 
  
  HMC = function (U, grad_U, epsilon, L = 30, current_q, k, arc, st=F)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, k) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, k)
      if(st){
        D <- 1
        if (j!=L) p = p - epsilon * grad_U(q, k)-D*p+sqrt(2*D)*rnorm(length(q),0,1)}
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, k) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, k)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, k)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  HMCm = function (U, grad_U, epsilon, L = 30, current_q, k, j, arc)
  {
    q = current_q
    p = rnorm(length(q),0,1)  # independent standard normal variates
    current_p = p
    
    # Make a half step for momentum at the beginning
    
    p = p - epsilon * grad_U(q, k, j) / 2
    
    # Alternate full steps for position and momentum
    
    for (j in 1:L)
    {
      # Make a full step for the position
      
      q = q + epsilon * p
      q = (q>tl)*tl + (q < tl1) * (tl1) + q * (q <= tl) * (q >= tl1)
      
      # Make a full step for the momentum, except at end of trajectory
      
      if (j!=L) p = p - epsilon * grad_U(q, k, j)
    }
    
    # Make a half step for momentum at the end.
    
    p = p - epsilon * grad_U(q, k, j) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    
    p = -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    
    current_U = U(current_q, k, j)
    current_K = sum(current_p^2) / 2
    proposed_U = U(q, k, j)
    proposed_K = sum(p^2) / 2
    
    # Accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    
    R <- current_U-proposed_U+current_K-proposed_K
    if(is.nan(R)){R = -Inf}
    if (log(runif(1)) < R)
    {
      up = q  # accept
      arc <- arc + 1
    }
    else
    {
      up = current_q
    }
    return(list(up = up, arc = arc))  # reject
  }
  
  UL <- function(x, sigm){
    tempmat <- diag(p)
    tempmat[lower.tri(tempmat, diag = T)] <- x
    temp <- tempmat %*% diag(sigm) %*% t(tempmat)
    
    return(temp)
  }
  
  Umu <- function(muc, k, j){
    Yred <- eta[j, ]
    compo <- compoBig[,,j]
    coef  <- coefBig[, j]
    if(K > 1){
      if(K==2){tvec <- compo[, -k] * coef[-k]}
      if(K>2){
        tvec <- compo[, -k]%*%coef[-k]
      }
      Yred <- eta[j, ] - tvec 
    }
    tempmat <- sigmaM1
    tempmat <- tempmat[Ones, Ones]
    
    if(length(Ones) > 1 ){
      Qtemp <- tX[Ones, ]-matrix(muc[Ones], length(Ones), n)
      Par <- tempmat %*% (Qtemp);
      part1 <- colSums((Qtemp) * Par)
    }
    if(length(Ones) == 1){
      temp <- tempmat
      Qtemp <- array((X[, Ones]-muc[Ones])^2*temp)#-matrix(mu[k, Ones], length(Ones), n)
      #Par <- Qtemp
      part1 <- sum(Qtemp)
    }
    if(sum(zindex)==0){
      part1 <- 0
    }
    #part1 <- apply(X-matrix(muc, n, p, byrow=T), 1, FUN = function(x){wcrossprod(x, x, )})
    
    ret <- Yred - coef[k] * exp(- part1)
    ret <- sum(ret ^ 2) / sigma1[j]^2
    
    return(ret / 2 )#+ sum((muc[Ones]-mumean)^2) / (2*slab1^2))
  }
  
  grad_Umu <- function(muc, k, j, compoc=compo){
    Yred <- eta[j, ]
    compo <- compoBig[,,j]
    coef  <- coefBig[, j]
    if(K > 1){
      if(K==2){tvec <- compo[, -k] * coef[-k]}
      if(K>2){
        tvec <- compo[, -k]%*%coef[-k]
      }
      Yred <- eta[j, ] - tvec 
    }
    tempmat <- sigmaM
    #tempmat <- tempmat[Ones,]
    
    Par <- tempmat %*% (tX-matrix(muc, d, n));
    part1 <- colSums((tX-matrix(muc, d, n)) * Par)
    
    ret <- Yred - coef[k] * exp(- part1)
    retder <- 2*array((Yred - coef[k] * exp(- part1)) * coef[k] * exp(- part1))
    
    part1der <- -2*tempmat%*% t(X-matrix(muc, n, d, byrow=T))
    
    retf <- part1der %*% retder
    
    temp <- rep(0, p)
    temp[Ones] <- retf
    
    return(temp/(sigma1[j]^2))# + (muc-mumean) / slab1^2)
  }
  
  
  UsigD <- function(x, k=1){
    sum <- 0
    tMat <- choL
    temp1 <- tcrossprod(tMat[Ones, ]*x[Ones])
   
    for(j in 1:p){
      Y <- eta[j, ]
      #compo <- compoBig[,,j]
      coef  <- coefBig[, j]
      k = 1
      if(length(Ones) > 1 ){
        Qtemp <- tX[Ones, ]-matrix(mu[k, Ones, j], length(Ones), n)
        Par <- temp1 %*% (Qtemp);
        compo <- colSums((Qtemp) * Par)
      }
      if(length(Ones) == 1){
        Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp1)#-matrix(mu[k, Ones], length(Ones), n)
        #Par <- Qtemp
        compo <- sum(Qtemp)
      }
      if(sum(zindex)==0){
        compo <- 0
      }
      
      part1 <- coef[k] * exp(-compo)
      for(k in 2:K){
        tMat <- choL
        
        if(length(Ones) > 1 ){
          Qtemp <- tX[Ones, ]-matrix(mu[k, Ones, j], length(Ones), n)
          Par <- temp1 %*% (Qtemp);
          compo <- colSums((Qtemp) * Par)
        }
        if(length(Ones) == 1){
          Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp1)#-matrix(mu[k, Ones], length(Ones), n)
          #Par <- Qtemp
          compo <- sum(Qtemp)
        }
        if(sum(zindex)==0){
          compo <- 0
        }
        part1 <- part1 + coef[k] * exp( - compo)
      }
      
      ret <- Y - part1
      ret <- sum(ret ^ 2) / sigma1[j]^2
      sum <- sum + ret
    }
    baad <- 0
    
    if(length(Ones)<d){baad <- (itr>noupS1)*sum(x[-Ones]^2/(2*betavar[-Ones]^2))} 
    return(sum / 2 + sum(x^2/(2*betavar^2))-baad)
  }
  
  
  grad_UsigD <- function(x, k=1){
    k = 1
    dersum <- 0
    tMat <- choL
    temp1 <- tcrossprod(tMat[Ones, ]*x[Ones])
    temp2 <- matrix(0,d,d)
    if(length(Ones)>1){temp2[Ones, Ones] <- tcrossprod(tMat[Ones, ]*x[Ones], tMat[Ones, ])}
    #temp3 <- tcrossprod(tMat*x*zindex, tMat)
    for(j in 1:p){
      if(length(Ones) > 1 ){
        Qtemp <- t(X[, Ones])-matrix(mu[k, Ones, j], length(Ones), n)
        Par <- temp1 %*% (Qtemp);
        compo <- colSums((Qtemp) * Par)
      }
      if(length(Ones) == 1){
        Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp1)#-matrix(mu[k, Ones], length(Ones), n)
        #Par <- Qtemp
        compo <- sum(Qtemp)
      }
      if(sum(zindex)==0){
        compo <- 0
      }
      
      part1 <- coefBig[k, j] * exp( - compo)
      matvec   <- (X-matrix(mu[k, ,j], n, d, byrow=T)) %*% temp2
      matvec   <- (X-matrix(mu[k, ,j], n, d, byrow=T)) * matvec
      sum      <-  - 2*(matvec)*coefBig[k, j] * exp( - compo)
      #retder <-  coef[k] * exp(- compo)
      for(k in 2:K){
        tMat <- choL
        
        if(length(Ones) > 1 ){
          Qtemp <- tX[Ones, ]-matrix(mu[k, Ones, j], length(Ones), n)
          Par <- temp1 %*% (Qtemp);
          compo <- colSums((Qtemp) * Par)
        }
        if(length(Ones) == 1){
          Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp1)#-matrix(mu[k, Ones], length(Ones), n)
          #Par <- Qtemp
          compo <- sum(Qtemp)
        }
        if(sum(zindex)==0){
          compo <- 0
        }
        
        part1 <- part1 + coefBig[k, j] * exp( - compo)
        matvec   <- (X-matrix(mu[k, ,j], n, d, byrow=T)) %*% temp2
        matvec   <- (X-matrix(mu[k, ,j], n, d, byrow=T)) * matvec
        sum      <- sum - 2*(matvec)*coefBig[k, j] * exp( - compo)
      }
      
      ret <- eta[j,] - part1
      retder <- ret / sigma1[j]^2
      
      sumr <- sum * retder#matrix(retder, n, p)
      
      der <- - colSums(sumr)
      dersum <- dersum + der
    }
    
    
    ret <- dersum + x / betavar^2
    
    #if(itr <= noupMM){ret[-indz] <- 0}
    
    return(ret)
  }
  
  UZ <- function(zindexc){
    sum <- 0
    Ones <- which(zindexc==1)
    tempmat <- choL
    if(length(Ones)>0){temp <- tcrossprod(tempmat[Ones, ]*sigD[Ones])}
    for(j in 1:p){
      mean <- 0
      for(k in 1:K){
        #ind <- which((tX[Ones, ]-mu[k, Ones, j])^2)
        #temp <- diag(sigD[k, Ones]) %*% tempmat[Ones, ] %*%  t(tempmat[Ones, ]) %*% diag(sigD[k,Ones])
        if(length(Ones) > 1 ){
          Qtemp <- tX[Ones, ]-matrix(mu[k, Ones, j], length(Ones), n)
          Par <- temp %*% (Qtemp);
          part1 <- colSums((Qtemp) * Par)
        }
        if(length(Ones) == 1){
          Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp)#-matrix(mu[k, Ones], length(Ones), n)
          #Par <- Qtemp
          part1 <- sum(Qtemp)
        }
        
        if(sum(zindexc)==0){
          part1 <- 0
        }
        mean <- mean + coefBig[k,j] * exp(- part1)
      }
      
      ret <- eta[j, ] - mean 
      ret <- sum(ret ^ 2) / sigma1[j]^2
      sum <- sum + ret
    }
    return(sum / 2)# + sum(x^2)/200)
  }
  
  UsigL <- function(x){
    sum <- 0
    Ones <- which(zindex==1)
    tempmat <- diag(d)
    tempmat[lower.tri(tempmat, diag = T)] <- x
    if(length(Ones)>0){temp <- tcrossprod(tempmat[Ones, ]*sigD[Ones])}
    
    for(j in 1:p){
      mean <- 0
      for(k in 1:K){
        #temp <- diag(sigD[k, Ones]) %*% tempmat[Ones, ] %*%  t(tempmat[Ones, ]) %*% diag(sigD[k,Ones])
        if(length(Ones) > 1 ){
          Qtemp <- tX[Ones, ]-matrix(mu[k, Ones, j], length(Ones), n)
          Par <- temp %*% (Qtemp);
          part1 <- colSums((Qtemp) * Par)
        }
        if(length(Ones) == 1){
          Qtemp <- array((X[, Ones]-mu[k, Ones, j])^2*temp)#-matrix(mu[k, Ones], length(Ones), n)
          #Par <- Qtemp
          part1 <- sum(Qtemp)
        }
        
        if(sum(zindex)==0){
          part1 <- 0
        }
        mean <- mean + coefBig[k,j] * exp(- part1)
      }
      
      ret <- eta[j, ] - mean 
      ret <- sum(ret ^ 2) / sigma1[j]^2
      sum <- sum + ret
    }
    return(sum / 2)# + sum(x^2)/200)
  }
  
  
  
  polrec <- function(angl, r =1){
    len  <- length(angl)
    temp <- exp(cumsum(log(sin(angl[1:(len-1)]))))
    temp <- c(1, temp)
    
    temp[1:len] <- temp[1:len] * cos(angl)
    temp <- c(temp, prod(sin(angl)))
    return(r * temp)
  }
  
  choL <- list()
  d <- ncol(Xtr)
  p <- nrow(Ytr)
  sigmaM <- list()
  sigmaM1 <- list()
  polthl <- list()
  
  
  if(is.null(Xte)){Yte <- NULL}
  
  Y  <- Ytr#c(Ytr, Ytest)
  #X2 <- rbind(Xtr, Xte)
  #X2 <- X2 - matrix(colMeans(X2), nrow(X2), p, byrow = T)
  
  X <- Xtr#X2[1:nrow(Xtr), ]
  tX <- t(X)
  #Xte <- X2[-c(1:nrow(Xtr)), ]
  #out  <- kmeans(X, centers = K)
  
  n <- nrow(X)#nrow(Xtr)
  nte <- nrow(Xte)
  
  slab  <- 10
  
  indzl <- list()
  
  Y <- Ytr
  if(!is.null(Xte)){Y <- t(apply(Ytr, 1, function(x){nax <- which(is.na(x)==T); x[nax] <- rep(mean(x[-nax]), length(nax)); return(x)}))}
  
  if(is.null(indz)){
    indzM <- matrix(0, p, d)
    for(j in 1:p){
      ntr <- which(is.na(Ytr[j, ])==T)
      if(length(ntr)>0){
        m3 <- huge.npn(cbind(Ytr[j, -ntr], Xtr[-ntr, ]), npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
        
        if(corel){
          result <- NULL#mvnTest::HZ.test(m3[,1:2])#apply(m3, 2, FUN=function(x){cor.test(m3[,1],x, method=c("pearson"))$p.value})#
          for(i in 1:ncol(Xtr)){
            result[i] <- cor.test(m3[,1], m3[,1+i], method=c("pearson"))$p.value
          }
          
          result[which(is.na(result))] <- Inf
          
          indzM[j, which(abs(result)<=0.01)] <- 1
        }
        if(lasso){
          lambda <- cv.glmnet(m3[,-1], m3[, 1])
          result <- glmnet(m3[,-1], m3[,1], lambda = lambda$lambda.min) 
          indzM[j, which(result$beta!=0)] <- 1
        }
      }
      
      if(length(ntr)==0){
        m3 <- huge.npn(cbind(Ytr[j, ], Xtr[, ]), npn.func = "shrinkage", npn.thresh = NULL, verbose = TRUE)
        
        if(corel){
          result <- NULL#mvnTest::HZ.test(m3[,1:2])#apply(m3, 2, FUN=function(x){cor.test(m3[,1],x, method=c("pearson"))$p.value})#
          for(i in 1:ncol(Xtr)){
            result[i] <- cor.test(m3[,1], m3[,1+i], method=c("pearson"))$p.value
          }
          
          result[which(is.na(result))] <- Inf
          
          indzM[j, which(abs(result)<=cutpval)] <- 1
        }
        if(lasso){
          lambda <- cv.glmnet(m3[,-1], m3[, 1])
          result <- glmnet(m3[,-1], m3[,1], lambda = lambda$lambda.min) 
          indzM[j, which(result$beta!=0)] <- 1
        }
      }
    }
    
    indz <- which(colMeans(indzM)!=0)
    if(is.null(K)){
    Kvec <- rep(0, p)
    for(j in 1:p){
      ntr <- which(is.na(Ytr[j, ])==T)
      rvmm <- NULL
      rvmerror <- NA
      if(length(ntr)>0){
        out <- try(rvmm <- rvm(Xtr[-ntr, indz], Ytr[j, -ntr],kernel="rbfdot", type="regression",list(sigma=0.1)), silent = T)} #kpar="automatic", iterations  =100000)#
      if(length(ntr)==0){
        out <- try(rvmm <- rvm(Xtr[, indz], Ytr[j, ],kernel="rbfdot", type="regression",list(sigma=0.1)), silent = T)} #kpar="automatic", iterations  =100000)#
      multi <- 1
      while(class(out)=="try-error"){
        multi <- multi+1
        if(length(ntr)>0){
          out <- try(rvmm <- rvm(Xtr[-ntr, indz], Ytr[j, -ntr],kernel="rbfdot", type="regression",list(sigma=0.1*multi)), silent = T)} #kpar="automatic", iterations  =100000)#
        if(length(ntr)==0){
          out <- try(rvmm <- rvm(Xtr[, indz], Ytr[j, ],kernel="rbfdot", type="regression",list(sigma=0.1*multi)), silent = T)} #kpar="automatic", iterations  =100000)#
        if(multi>10){break}
      }
      if(!is.null(rvmm)){
        Kvec[j] <- out@nRV}
    }
    
    K <- max(Kvec) 
    }
  }
  
  if(K==0){K=4}
  
  compoBig  <- array(0, dim=c(K, d, p))
  mu   <- array(0, dim=c(K, d, p))#out$centers
  
  if(is.null(fitL)){
    if(updateLow){
      if(!is.null(Xte)){
        fitL <- MultiRBFjointLow(Ytr=Ytr, Xtr=Xtr[, indz], K = K, Xte = Xte[, indz], Yte=Yte, Total_itr = 5000, burn = 2500)
      }
      if(is.null(Xte)){
        fitL <- MultiRBFjointLow(Ytr, Xtr[, indz], K = K, Xte = NULL, Yte=NULL, Total_itr = 5000, burn = 2500) # noupR = noupRL, noupM = 1000, noupS = 20000,
      } 
    }
  }
  print(K)
  if(!updateLow){
    vec1    <- rep(1, length(indz))
    choLp   <- diag(length(indz))
    Y <- t(apply(Ytr, 1, function(x){nax <- which(is.na(x)==T); x[nax] <- rep(mean(x[-nax]), length(nax)); return(x)}))
    
    for(j in 1:p){
      out     <- kmeans(cbind(Y[j,], Xtr), centers = K)
      temp    <- out$centers
      mu[,,j] <- temp[, -1]#matrix(rnorm(p*K),K, p)#
    }
    
    compoBig  <- array(0, dim=c(n, K, p))
    coefBig <- matrix(rnorm(K*p),K, p)
    
    vec <- rep(0, d)#1/sqrt(p)
    vec[indz] <- vec1#rep(1, length(indz))#rnorm(length(indz), sd = 1)#
    
    sigD <- vec#matrix(vec, K, d, byrow = T)
    
    V     <- diag(d)
    polth <- diag(d)
    
    index1 <- which(diag(d)==1)
    index2 <- index1 + 1
    index2 <- index2[-d]
    
    indplth <- which(lower.tri(polth)==T)
    
    index3 <- setdiff(indplth, union(index1, index2))
    index3len <- length(index3)
    
    
    V     <- diag(d)
    polth <- diag(d)
    
    temp <- matrix(0, length(indz), length(indz))
    temp <- choLp
    
    temp <- diag(1/sqrt(rowSums(temp^2))) %*% temp
    
    V[indz, indz] <- temp
    
    for(l in 2:d){
      angl <- rect2polar(V[l,l:1])
      angl <- angl$phi
      angl[which(angl<0)] <- angl[which(angl<0)]  + 2*pi
      polth[l, 1:(l-1)]  <- angl
      
      if(l>2) {V[l, l:1] <- polrec(polth[l, 1:(l-1)])} 
      if(l==2) {V[l, l:1] <- c(cos(polth[l, 1:(l-1)]), sin(polth[l, 1:(l-1)]))} 
    }
    
    polthl <- polth
    choL   <-  V
    
    tempmat <- tcrossprod(choL*sigD) #diag(sigD[i, ]) %*% choL[[i]] %*% t(choL[[i]]) %*% diag(sigD[i, ])
    sigmaM  <- tempmat
    
    FX <- Ytr   
    for(j in 1:p){ 
      for(i in 1:K){
        Par              <- tempmat %*% (tX-matrix(mu[i, ,j], d, n));
        compoBig[, i, j] <- exp(-colSums((tX-matrix(mu[i, , j], d, n)) * Par)) 
      }
      fit <- lm(Y[j, ] ~ compoBig[, , j]-1)
      coef <- fit$coefficients
      coef[is.nan(coef)] <- 0
      coef[is.na(coef)] <- 0
      coefBig[, j] <- coef
      FX[j, ] <- as.matrix(compoBig[, , j]) %*% array(coefBig[, j])
    }
    
    
    
    #tcov <- solve(cov(Xtr[, indz]))[1,]#cov(Xtr[, indz])[1,]
    
    
    #for(i in 1:K)
    
    Yr <- Y - matrix(rowMeans(Y), p, n)
    
    lambda <- Yr %*% ginv(crossprod(Yr)) %*% t(Yr) #matrix(rnorm(p*r, 0, sd = 1), p, r)#
    sigma1 <- rep(1, p)#apply(Y - lambda %*% eta, 1, sd) #
    sigma2 <- rep(1, p)#apply(eta - gamma %*% Y, 1, sd) #rep(1, p)#
   
    if(runFA){
      fitPFA <- LFMPFA(Yr, d=100, showupdate = F)
      sigma2 <- Reduce('+', fitPFA$sigma2ls)/length(fitPFA$lambdals)
      sigma1 <- Reduce('+', fitPFA$sigma1ls)/length(fitPFA$lambdals)
      lambda <- Reduce('+', fitPFA$lambdals)/length(fitPFA$lambdals) #%*% diag(sigma2p)
    } 
    
    r <- ncol(lambda)
    
    epsilon2 <- ginv(crossprod(lambda)) %*% (crossprod(lambda, Yr))
    eta <- Y - lambda %*% epsilon2
  }
  
  YFr <- Y-FX
  
  nu <- 1
  a1 <- 2.1
  a2 <- 3.1
  eta.var2 <- rep(0, r)
  var <- eta.var2
  delta <- rnorm(r)
  eta.var2[r:1] <- rep(1, r)#cumsum(exp(delta[r:1])) / sum(exp(delta))
  #set.seed(seed)
  phi <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)
  
  psi <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
  tau <- exp(cumsum(log(psi)))
  
  sigmaM1    <- sigmaM
  sigDls     <- list()
  cholls     <- list()
  muls       <- list()
  coefls     <- list()
  sigls      <- list()
  Ytels      <- list()
  indls      <- list()
  vecl       <- list()
  Y_trainhat <- list()
  gindex_p   <- list()
  sigma1_p   <- list()
  sigma2_p   <- list()
  lambda_p   <- list()
  
  itr    <- 0
  
  
  zindiindex <- list()
  lenlist <- list()
  G <- length(grpindex)
  sum <- 0
  for(i in 1:G){
    lenlist[[i]] <- length(grpindex[[i]])
    zindiindex[[i]] <- rep(1, lenlist[[i]])
    #grpindex[[i]]  <- 1:lenlist[[i]] + sum
  }
  
  gindex  <- rep(1, G)
  bgindex <- unlist(lapply(1:length(gindex), FUN = function(i){rep(gindex[i], lenlist[[i]])}))
  
  zindex <- bgindex*unlist(zindiindex)#indz#sample(0:1, p, replace = T, prob = c(1-p1, p1))
  #zindex[indz] <- 1
  Ones   <- which(zindex==1)
  
  arsigD <- rep(0, 1)#rep(0)
  arch   <- 0#rep(0, K)
  arm    <- matrix(0, K, p)
  armf   <- 0
  
  sdsigD <- rep(1e-5, 1)#10^(-2*co1)#
  sdch   <- 1e-3#rep(1, K)
  sdm    <- matrix(1e-4, K, p)
  sdmf   <- rep(1e-4, K)
  
  noupS1 <- Total_itr
  noupS <- Total_itr
  noupS2 <- 00
  noupD <- 2000
  noupM3 <- Total_itr#3000
  noupM <- 500 #2000
  noupM2 <- 500 #1500
  noupM1 <- 500 # 1000
  
  # noupS1 <- 000
  # noupS <- 00
  # noupS2 <- 00
  # noupD <- 2000
  # noupM3 <- Total_itr#3000
  # noupM <- 500
  # noupM2 <- 500
  # noupM1 <- 500
  noupR  <- 2300
  
  burn1 <- Total_itr#burn
  
  coefmean <- NULL
  coefsd <- NULL
  
  for(j in 1:p){
    coefmean[j] <- (max(Y[j, ]) + min(Y[j, ])) /(2*K)
    coefsd[j]  <- (max(Y[j, ]) - min(Y[j, ])) /(4*sqrt(K)) 
  }
  
  len <- 20
  len1 <- 1
  
  #if(is.na(p1)){p1 <- 1/p} #length(indz)
  
  if(is.na(pindex)){pindex <- rep(0.5, G)}##1-p1 rep(0.5, G)#1-1/unlist(lenlist)
  if(is.na(pindexG)){pindexG <- 1-1/G} #0.9#
  
  blen <- 1
  #len  <- 1
  
  sigma <- 1
  tl    <- 1#30*p1 #.3
  tl1   <- 0#-30*p1#0.3
  
  betavar <- rep(slab, d)
  burn1 <- Total_itr#burn
  pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)
  
  while(itr < Total_itr){
    itr <- itr + 1
    
    varg   <- lambda%*%diag(sigma2^2)%*%t(lambda)+diag(sigma1^2)
    
    for(j in 1:n){
      naindex <- which(is.na(Ytr[, j])==T)
      if(length(naindex)>0){
        temp    <- YFr[, j]
        vargn   <- matrix(varg[naindex, naindex], nrow = length(naindex) )
        vargg   <- matrix(varg[- naindex, - naindex], nrow = p - length(naindex))
        covng   <- matrix(varg[naindex, - naindex], nrow = length(naindex))
        conMean <- covng %*% solve(vargg) %*% temp[-naindex]
        temp[naindex] <- conMean#array(rmvnorm(1, conMean, sigma = convar))
        YFr[, j]        <- temp 
      }
    }
    
    for(j in 1:p){
      if(!is.null(Xte)){
        indNA <- which(is.na(Ytr[j,]))
        tXtej <- tX[, indNA]
        nte  <- length(indNA)
        compoBig1 <- array(0, dim=c(nte,K, p))
        for(k in 1:K){
          Par        <- sigmaM1 %*% (tXtej-matrix(mu[k, , j], d, nte));
          compoBig1[, k, j] <- exp(-colSums((tXtej-matrix(mu[k, , j], d, nte)) * Par)) 
          #compo1[, k] <- exp(-apply(Xte-matrix(mu[k, ], nrow(Xte), p, byrow=T), 1, FUN = function(x){wcrossprod(x, x, sigmaM1[[k]])}))
        }
        Y[j, indNA] <- YFr[j, indNA] + as.matrix(compoBig1[,,j]) %*% coefBig[, j]#rnorm(length(testind), compo[testind, ] %*% array(coef), sigma) #compo[testind, ] %*% array(coef) #
        #Y <- c(Ytr, Ytest)
      }
    }
    
    Yhatred <- Y - FX - lambda %*% epsilon2
    
    for(i in 1:p){
      al       <- 0.1 + n /2
      be       <- 0.1 + sum((Yhatred[i, ])^2)/2
      sigma1[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    for(i in 1:r){
      al       <- 100 + n /2
      be       <- 0.1 + sum((epsilon2[i, ])^2)/2
      sigma2[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    var.pm <- ginv(crossprod(lambda / sigma1) + diag(1 / sigma2^2))
    var.pm <- (var.pm + t(var.pm)) / 2
    
    YFr <- Y - FX
    mean.etai <- t(lambda/sigma1^2) %*% (YFr)
    mean.etai <- var.pm %*% mean.etai
    
    temp <- apply(mean.etai, 2, FUN=function(x){mvtnorm::rmvnorm(1, x, var.pm)})
    epsilon2  <- temp
    
    
    for(i in 1:r){
      mean.lami <- rowSums((YFr - lambda[, -i] %*% epsilon2[ - i, ])*matrix(epsilon2[i, ], p, n, byrow=T)/sigma1^2)
      var.lami  <- 1/(sum(epsilon2[i, ]^2)/sigma1^2 + phi[, i]*tau[i]) #phi[, i]*tau[i]) #rep(1/100, p)))
      #var.lami  <- ( var.lami + t(var.lami) ) / 2
      mean.lami <- var.lami * mean.lami
      
      lambda[, i] <- rnorm(p, mean.lami, sqrt(var.lami))
      
      phi[, i]    <- rgamma(p, (nu + 1) / 2, ((lambda[, i]^2) * tau[i] + nu) / 2)
    }
    
    eta <- Y - lambda %*% epsilon2
    
    temp1 <- colSums((lambda ^ 2) * phi)
    
    tauprime1 <- tau / (psi[1])
    
    psi[1] <- rgamma(1, a1 + p*r/2, 1 + sum(temp1 * tauprime1) / 2)
    tau    <- tauprime1 * psi[1]
    
    for(i in 2:r){
      tauprime1 <- tau / (psi[i])
      psi[i] <- rgamma(1, a2 + p*(r - i + 1) / 2, 1 + sum(temp1[i:r] * tauprime1[i:r]) / 2)
      tau    <- tauprime1 * psi[i]
    }
    
    eta.var2[r:1] <- sigma2[r:1]^2#cumsum(exp(delta[r:1])) / sum(exp(delta))
    
    Ones <- which(zindex==1)
    #Sparse Update vec
    if(itr > noupS1){
      temp        <- sigD
      temp[Ones]  <- sigD[Ones] + rnorm(length(Ones), 0, sdsigD)
      if(length(Ones)>1){temp[-Ones] <- rnorm(d-length(Ones), sd = slab)}
      D           <- UsigD(sigD) - UsigD(temp)
      
      vec <- sigD
      if(is.nan(D)){D = -Inf}
      if(is.na(D)){D = -Inf}
      if(D > log(runif(1))){
        arsigD <- arsigD + 1
        vec    <- temp
      } 
    }
    
    if((itr <= noupS)*(itr>noupS2)){
      temp   <- HMC(UsigD, grad_UsigD, sdsigD, L = 1, sigD, 1, arsigD)
      vec    <- temp$up
      arsigD <- temp$arc
    }
    
    sigD <- vec#matrix(vec, K, p, byrow = T)
    
    #if(itr > noupS1){
    for(g in 1:G){
      al   <- 0.1 + lenlist[[g]] / 2
      bl   <- 0.1 + sum(sigD[grpindex[[g]]] ^ 2) / 2
      slab <- sqrt(1/rgamma(1, al, bl)) 
      
      betavar[grpindex[[g]]] <- rep(slab, lenlist[[g]])
    }
    #}
    
    if(itr>noupD){
      if(itr %% len == 0){
        currentU <- UZ(zindex)
        for(g in 1:G){
          gindexc    <- gindex
          gindexc[g] <- 1-gindex[g]
          bgindexc <- unlist(lapply(1:length(gindexc), FUN = function(i){rep(gindexc[i], lenlist[[i]])}))
          
          zindexc <- bgindexc * unlist(zindiindex) 
          
          if(gindex[g] == 1){ratio <- - currentU + UZ(zindexc)}
          if(gindex[g] == 0){ratio <- - UZ(zindexc) + currentU}
          
          pratio <- 0
          if(pindexG < 1){pratio <- exp(log(1 - pindexG) + ratio)}
          prob <-  pratio / (pratio + pindexG)
          if(pratio==Inf){prob=1}
          
          gindex[g] <- rbinom(blen, 1, prob)
          if(gindex[g]==gindexc[g]) {currentU <- UZ(zindexc);bgindex<-bgindexc; zindex <- zindexc} 
          if(gindex[g] == 0){zindiindex[[g]] <- rbinom(lenlist[[g]], size = 1, prob = 1-pindex[g])}
          if(gindex[g] == 1){
            currentU <- UZ(zindex)
            newsamZ  <- grpindex[[g]]#sample(1:p, p)
            for(i in newsamZ){
              # indup   <- (1:blen) + (i-1)*blen
              zindexc    <- zindex
              zindexc[i] <- 1- zindex[i]
              
              if(zindex[i] == 1){ratio <- - currentU + UZ(zindexc)}
              if(zindex[i] == 0){ratio <- - UZ(zindexc) + currentU}
              
              pratio <- 0
              #ratio <- - UZ(zindexc) + UZ(zindexc1)
              if(pindex[g] < 1){pratio <- exp(log(1 - pindex[g]) + ratio)}
              #prob <- (1 - pindex)* exp(ratio)  / ((1 - pindex)* exp(ratio) + pindex)
              prob <-  pratio / (pratio + pindex[g])
              if(pratio==Inf){prob=1}
              
              zindex[i] <- rbinom(blen, 1, prob)
              if(zindex[i]==zindexc[i]) {currentU <- UZ(zindexc)} 
            }
            zindiindex[[g]] <- zindex[grpindex[[g]]] 
            if(pindexprior){pindex[g] <- rbeta(1, 1 + lenlist[[g]] - sum(zindiindex[[g]]), 1 + sum(zindiindex[[g]]))} # prob of spike}
          }
          zindex <- bgindex*unlist(zindiindex)
        }
        if(pindexGprior){pindexG <- rbeta(1, 1 + G - sum(gindex), 1 + sum(gindex))}
      }
    }
    
    Ones <- which(zindex==1)
    
    
    if(itr > noupR){
      polth  <- polthl
      polthc <- polth
      #polthc[indplth] <- polth[indplth] + rnorm(length(indplth), sd = sdch[k])#runif(length(indplth), - sdch[k], sdch[k])
      
      Vc <- choL
      
      l=2
      temp <- runif(1, 0, 2*pi)
      polthc[l, 1:(l-1)] <- temp
      if(zindex[l]==1){
        temp <- msm::rtnorm(1, polth[l, 1:(l-1)], sd = sdch, 0, 2*pi)
        polthc[l, 1:(l-1)] <- temp
      }
      
      Vc[l, l:1] <- c(cos(temp), sin(temp)) 
      
      for(l in 3:d){
        if(zindex[l]==1){
          temp <- polth[l, 1:(l-1)]
          temp[1:(l-2)] <- msm::rtnorm(l-2, polth[l, 1:(l-2)], sd = sdch, 0, pi) #polthc[l, 1:(l-1)]
          temp[l-1]     <- msm::rtnorm(1, polth[l, (l-1)], sd = sdch, 0, 2*pi)#(temp[l-1] %% (2 * pi))
          polthc[l, 1:(l-1)] <- temp
          Vc[l, l:1] <- polrec(temp) 
        }
        if(zindex[l]==0){
          temp <- polth[l, 1:(l-1)]
          temp[1:(l-2)] <- runif(l-2, 0, pi) #polthc[l, 1:(l-1)]
          temp[l-1]     <- runif(1, 0, 2*pi)#(temp[l-1] %% (2 * pi))
          polthc[l, 1:(l-1)] <- temp
          Vc[l, l:1] <- polrec(temp)
        }
      }
      
      mat <- choL
      
      temp <- Vc[lower.tri(Vc, diag = T)]
      D    <- UsigL(mat[lower.tri(mat, diag = T)]) - UsigL(temp)
      if(is.nan(D)){D=-Inf}
      if(is.na(D)){D=-Inf}
      
      if(D > log(runif(1))){
        arch                           <- arch + 1
        mat[lower.tri(mat, diag = T)]  <- temp
        choL                           <- Vc
        polthl                         <- polthc
      }
    }
    
    for(j in 1:p){
      for(k in 1:K){
        if((itr > noupM)*(itr <= noupM3)){
          if(itr %% len1 == 0){
            # temp    <- HMC(Umu, grad_Umu, sdm[k], L = 5, mu[k, ], k, arm[k])
            # mu[k, ] <- temp$up
            # arm[k]  <- temp$arc
            
            temp       <- mu[k, , j]
            temp[Ones] <- msm:: rtnorm(length(Ones), mu[k, Ones, j], sdm[k, j], tl1, tl) #mu[k, ] + rnorm(p, 0, sdm[k])
            if(length(Ones)>1){temp[-Ones] <- runif(d-length(Ones))}#msm:: rtnorm(p-length(Ones), mumean, sd = slab1, tl1, tl)}
            D    <- Umu(mu[k, ,j], k, j) - Umu(temp, k, j)
            
            if(is.nan(D)){D = -Inf}
            if(is.na(D)){D = -Inf}
            if(D > log(runif(1))){
              arm[k, j]  <- arm[k, j] + 1
              mu[k, , j] <- temp
            }
          }
        }
        
        if((itr <= noupM2)*(itr > noupM1)){
          temp      <- HMCm(Umu, grad_Umu, sdm[k, j], L = 1, mu[k, ,j], k, j, arm[k, j])
          mu[k, ,j] <- temp$up
          arm[k,j]  <- temp$arc
        }
        
        
        
        tempmat  <- tcrossprod(choL*sigD)
        tempmat1 <- tcrossprod(choL*sigD*zindex)
        sigmaM   <- tempmat
        sigmaM1  <- tempmat1
        
        if(itr <= noupD){
          Par        <- tempmat %*% (tX-matrix(mu[k, , j], d, n));
          compoBig[, k, j] <- exp(-colSums((tX-matrix(mu[k, , j], d, n)) * Par))}
        
        #compo[, k] <- exp(-apply(X-matrix(mu[k, ], n, p, byrow=T), 1, FUN = function(x){wcrossprod(x, x, tempmat)}))}
        if(itr > noupD){
          Par        <- tempmat1 %*% (tX-matrix(mu[k, , j], d, n));
          compoBig[, k, j] <- exp(-colSums((tX-matrix(mu[k, , j], d, n)) * Par)) 
        }
      }
      
      for(k in 1:K){
        if(K > 1){
          if(K==2){tvec <- compoBig[, -k, j] * coefBig[-k, j]}
          if(K>2){
            tvec <- compoBig[, -k, j]%*%coefBig[-k, j]
          }
          Yred <- eta[j, ] - tvec
        }
        mean <- sum(Yred*compoBig[,k,j]) / sigma^2 + coefmean[j] / coefsd[j]^2
        if(itr > burn1){var  <- 1/(sum(compoBig[,k, j]^2)/sigma^2+coefsd[j]^(-2)/sigma^2)}
        if(itr <= burn1){var  <- 1/(sum(compoBig[,k, j]^2)/sigma^2+coefsd[j]^(-2))}
        
        coefBig[k, j] <- rnorm(1, mean*var, sd = sqrt(var))
      }
      
      FX[j, ] <- as.matrix(compoBig[, , j]) %*% array(coefBig[, j])
    }
    
    
    #print(pindex)
    if(itr %% 100 == 0){
      #image(t(lambda))
      #if(!is.null(Xte)){print(mean(rowMeans((Y-Yte)^2, na.rm=T)))}
      
      if(itr > noupR){
        ar <- arch/ (itr-noupR)
        if(ar<.20){sdch <- sdch * (.5)}
        if(ar>.40){sdch <- sdch * (2)}
        if(sdch>15){sdch <- 15}
      }
      
      for(j in 1:K){
        if(itr > noupM){
          ar <- arm[j]/ (itr-noupM)
          if(ar<.10){sdm[j] <- sdm[j] * (.5)}
          if(ar>.30){sdm[j] <- sdm[j] * (2)}
          if(sdm[j]>15){sdm[j] <- 15}
        }
        
        if((itr <= noupM2)*(itr > noupM1)){
          ar <- arm[j]/ (itr-noupM1)
          if(ar<.45){sdm[j] <- sdm[j] * (.1)}
          if(ar>.70){sdm[j] <- sdm[j] * (10)}
        }
      }
      if(itr<=noupS){
        if(itr > noupS2){
          for(i in 1:1){
            ar <- arsigD[i]/ (itr-noupS2)
            #ar <- ar / p
            if(ar<.45){sdsigD[i] <- sdsigD[i] * (.1)}
            if(ar>.70){sdsigD[i] <- sdsigD[i] * (10)}  
          }
        }
      }
      if(itr>noupS1){
        ar <- arsigD/ (itr-noupS1)
        if(ar<.10){sdsigD <- sdsigD * (.1)}
        if(ar>.30){sdsigD <- sdsigD * (10)} 
      }
      cat(mean(arm)/(itr-noupM), "acceptance rate for muK")
      cat(mean(arsigD)/(itr-noupS2), "acceptance rate for sigK")
      cat(mean(arch)/(itr-noupR), "acceptance rate for chol")
    }
    
    if(itr==noupM){
      arm <- matrix(0, K, p)
      sdm <- matrix(1e-4, K, p)
    }
    if(itr==noupS1){
      arsigD = 0
      sdsigD <- 1e-4
      # 
    }
    
    if(itr > burn){
      sigma1_p[[itr-burn]] <- sigma1
      sigma2_p[[itr-burn]] <- sqrt(eta.var2)
      lambda_p[[itr-burn]] <- lambda
      sigDls[[itr - burn]] <- sigD
      #cholls[[itr - burn]] <- choL
      #muls[[itr - burn]]   <- mu
      coefls[[itr - burn]] <- coef
      sigls[[itr - burn]]  <- sigma
      if(!is.null(Xte)){
        Ytels[[itr - burn]] <- Y
      }
      if(Shr){indls[[itr - burn]] <- unlist(zindiindex);gindex_p[[itr - burn]] <- gindex}
    }
    setTxtProgressBar(pb, itr)
  }
  close(pb)
  #out <- list(sigDp = sigDls, cholp = cholls, mup = muls, coefp = coefls, sigp = sigls, Ytep = Ytels, indlp = indls, Y_trainhatp = Y_trainhat)
  out <- list(sigp = sigls,sigDlp=sigDls, coefp = coefls, Ytep = Ytels, indlp = indls, gindexp = gindex_p, lambdap = lambda_p, sigma1ls =  sigma1_p, sigma2ls =  sigma2_p)
  
  return(out)
}