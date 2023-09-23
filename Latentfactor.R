LFMPFA <- function(Y, d=10,showupdate=T, Total_itr = 10000, burn = 5000){

  set.seed(1)
  
  n <- ncol(Y)
  p  <- nrow(Y)
  r  <- p
  
  
  # parameters initialization
  
  r = p
  nu <- 1
  a1 <- 2.1
  a2 <- 3.1
  eta.var2 <- rep(0, r)
  var <- eta.var2
  delta <- rnorm(r)
  eta.var2[r:1] <- rep(1, r)#cumsum(exp(delta[r:1])) / sum(exp(delta))
  #set.seed(seed)
  phi <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)
  
  psi0 <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
  tau0 <- exp(cumsum(log(psi0)))
  
  psi <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
  tau <- exp(cumsum(log(psi)))
  
  phi2 <- matrix(rgamma(p * r, nu / 2, nu / 2), p, r)
  
  psi2 <- c(rgamma(1,5,1), rgamma(r-1, 5, 1))
  tau2 <- exp(cumsum(log(psi2)))
  
  #lambda <- Plambda0%*%matrix(rnorm(p*r, 0, sd = 1), p, r)#* RO #1 / sqrt(phi * matrix(tau, p, r, byrow = T)
  
  lambda <- Y %*% ginv(crossprod(Y)) %*% t(Y)#matrix(rnorm(p*r, 0, sd = 1), p, r)#Y %*% ginv(crossprod(Y)) %*% t(Y)
  eta <- ginv(crossprod(lambda)) %*% (crossprod(lambda, Y))
  lambdaginv <- ginv(crossprod(lambda)) %*% t(lambda)
  
  gamma <- lambdaginv #matrix(0, r, p) #ginv(tcrossprod(Y)) %*% (tcrossprod(Y, eta))
  
  sigma1 <- apply(Y - lambda %*% eta, 1, sd) #rep(1, p)#apply(Y - lambda %*% eta, 1, sd) #
  
  sigma2 <- apply(eta - gamma %*% Y, 1, sd) #rep(1, p)#apply(eta - gamma %*% Y, 1, sd) #
  
  Yhat  <- lambda %*% eta
  for(i in 1:p){
    al       <- 0.1 + n /2
    be       <- 0.1 + sum((Y[i, ]-Yhat[i, ])^2)/2
    sigma1[i] <- sqrt(1/rgamma(1, al, be))
  }
  
  etahat  <- gamma %*% Y
  
  for(i in 1:r){
    al       <- 0.1 + n /2
    be       <- 0.1 + sum((eta[i, ]-etahat[i, ])^2)/2
    sigma2[i] <- sqrt(1/rgamma(1, al, be))
  }
  
  
  eta_p   <- list()
  lambda_p <- list()
  gamma_p  <- list()
  sigma1_p <- list()
  sigma2_p <- list()
  itr <- 0
  
  R <- 100
  incre <- 4
  secparam <- 10*(1:r)
  ard1 <- 0
  ard2 <- 0
  arl <- 0
  arg <- 0
  L2 <- 10
  epsilon2 <- t(mvtnorm::rmvnorm(n, sigma = diag(sigma2)))
  Q <- diag(p)
  po <- 1
  
  Q_p <- list()
  if(showupdate){pb <- txtProgressBar(min = itr, max = Total_itr, style = 3)}
  
  while (itr < Total_itr) {
    itr <- itr + 1
    
    Yhatred <- Q %*% Y - lambda %*% epsilon2
    
    for(i in 1:p){
      al       <- 0.1 + n /2
      be       <- 0.1 + sum((Yhatred[i, ])^2)/2
      sigma1[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    for(i in 1:r){
      al       <- d + n /2
      be       <- 0.1 + sum((epsilon2[i, ])^2)/2
      sigma2[i] <- sqrt(1/rgamma(1, al, be))
    }
    
    var.pm <- ginv(crossprod(lambda / sigma1) + diag(1 / sigma2^2))
    var.pm <- (var.pm + t(var.pm)) / 2
    
    mean.etai <- t(lambda/sigma1^2) %*% (Q %*% Y)
    mean.etai <- var.pm %*% mean.etai
    
    temp <- apply(mean.etai, 2, FUN=function(x){mvtnorm::rmvnorm(1, x, var.pm)})
    epsilon2  <- temp
    
    sigma1_p[[itr]] <- sigma1
    
    eta <- epsilon2
    
    for(i in 1:r){
      mean.lami <- rowSums((Q%*%Y - lambda[, -i] %*% eta[ - i, ])*matrix(eta[i, ], p, n, byrow=T)/sigma1^2)
      var.lami  <- 1/(sum(eta[i, ]^2)/sigma1^2 + phi[, i]*tau[i]) #phi[, i]*tau[i]) #rep(1/100, p)))
      #var.lami  <- ( var.lami + t(var.lami) ) / 2
      mean.lami <- var.lami * mean.lami
      
      lambda[, i] <- rnorm(p, mean.lami, sqrt(var.lami))
      
      phi[, i]    <- rgamma(p, (nu + 1) / 2, ((lambda[, i]^2) * tau[i] + nu) / 2)
    }
    
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
    
    sigma2_p[[itr]] <- sqrt(eta.var2)
    
    
    
    #eta_p[[itr]]  <- eta
    lambda_p[[itr]] <- lambda
    #gamma_p[[itr]] <- gamma
    
    
    #print(itr)
    
    if(itr %% R==0){
      u <- runif(1)
      if(u < exp(-1 - itr * 5 * 10^(-4) )){
        temp    <- colMeans(abs(lambda))
        #print(temp)
        c <- (which(temp < 1e-4))
        if(r-length(c)<3){c <- NULL}
        if(r < 3) {c <- NULL}
        if(length(c)> 0){
          r <- r - length(c)
          lambda <- lambda[, -c]
          phi <- phi[, -c]
          tau   <- tau[-c]
          psi   <- psi[-c]
          phi2 <- phi2[, -c]
          tau2   <- tau2[-c]
          psi2   <- psi2[-c]
          tau0   <- tau0[-c]
          psi0   <- psi0[-c]
          eta  <- eta[-c, ]
          epsilon2  <- epsilon2[-c, ]
          eta.var2 <- eta.var2[ -c]
          sigma2 <- eta.var2
          gamma    <- gamma[-c, ]
          delta    <- delta[ - c]
        }
      }
      R     <- R + incre
      #incre <- 2*incre
    }
    #eta.var2 <- sigma2^2
    #image(t(lambda))
    if(showupdate){setTxtProgressBar(pb, itr)}
  }
  if(showupdate){close(pb)}
  out <- list(sigma1ls =  sigma1_p[(burn+1):Total_itr], sigma2ls =  sigma2_p[(burn+1):Total_itr], lambdals =  lambda_p[(burn+1):Total_itr])
  return(out)
}

# fitPFA <- LFMPFA(datasca, d=100)
# 
# lambdapsam <- array(0, dim=c(7,7, length(fitPFA$lambdals)))#array(unlist(lambdapp), dim = c(dim(lambdap), 1500))
# 
# for(i in 1:length(fitPFA$lambdals)){
#   sigma2p <- fitPFA$sigma2ls[[i]]
#   lambdapsam[,,i] <- fitPFA$lambdals[[i]] %*% diag(sigma2p)
# }
# 
# lambdapp <- t(apply(lambdapsam, 1:2, mean))
