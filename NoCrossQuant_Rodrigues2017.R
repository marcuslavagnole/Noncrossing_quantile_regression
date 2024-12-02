require(bayesQR); require(coda); require (plyr); require(grDevices); require(gptk); 
require(plot3D); require(Rearrangement); require(geometry);
require(VGAM); require(quadprog); require(BSquare); require(quantreg);

###################################################################
################### Main Function #################################
###################################################################
## Perform Two-stage Quantile Regression

## Arguments:
# yi: vector of responses
# xi: transpose of the matrix of predictors (do not include the intercept)
# gridtau: vector of quantiles to be estimated
# nmcmc: number of MCMC draws (for first stage BayesQR)
# burnmcm: number of effective samples discarded as burn-in (for first stage BayesQR)
# keep: thinning parameter (for first stage BayesQR)
# prior: an S3 object of class "prior" (for first stage bayesQR)
# sigma2.p: Variance hyperparameter of the GP kernel (See Eq. 6, notation from the paper: sigma^2_k)
# hzero: True if bandwith=0 is fixed (i.e, perform only Standard Bayesian QR) and False to perform the GPQR
# linear: True if linear regression and False otherwise

## Value
# band: bandwidth estimate
# quant.final.mean: matrix of final quantile estimates (for every tau and xi) 
# quant.final.sd: matrix of final quantile standard error (for every tau and xi)
# CI.final: confidence interval for final quantile (for every tau and xi)
# beta.final: matrix of final beta estimates (for every tau and covariate), if Linear=T
# beta.final.sd: matrix of final beta standard error (for every tau and covariate), if Linear=T
# CI.final.beta: confidence interval for final beta estimates (for every tau and covariate), if Linear=T
# initial.quant: matrix of quantile estimates from the first stage (for every tau and xi) 
# initial.quant.sd: matrix of quantile standard error from the first stage (for every tau and xi)
# CI.quant.yu: confidence interval for quantiles from the first stage (for every tau and xi)
# initial.beta: matrix of beta estimates from the first stage (for every tau and covariate), if Linear=T
# initial.beta.sd: matrix of beta standard error from the first stage (for every tau and covariate), if Linear=T
# CI.beta.yu: confidence interval for beta from the first stage (for every tau and covariate), if Linear=T

gpqr <- function(yi, xi, gridtau, nmcmc, burnmcmc, keep, prior, sigma2.p=100, hzero=F, linear=T){
  xi.data <- as.data.frame(t(xi))
  x.matrix <- rbind(1, xi)
  xg.matrix <- x.matrix
  gridp <- gridtau
  sample.size <- nmcmc/keep-burnmcmc
  lgridtau <- length(gridtau)
  lgridp <- lgridtau
  lbeta <- dim(x.matrix)[1]
  # Standard bayesian quantile regression
  out.qr <- bayesQR(t(yi)~., quantile=gridtau, ndraw=nmcmc, keep=keep, prior=prior, data=xi.data)
  initial.beta <- matrix(rapply(out.qr, function.mean, classes='matrix'), lgridp, lbeta, byrow=T)
  count.na <- which(is.na(initial.beta[,1]))
  while (any(is.na(initial.beta))){
    for(i in count.na){
      tau.na <- gridtau[i]
      out.na <- bayesQR(t(yi)~., quantile=tau.na, ndraw=nmcmc, keep=keep, prior=prior, data=xi.data)
      out.qr[[i]]$sigmadraw <- out.na[[1]]$sigmadraw
      out.qr[[i]]$betadraw <- out.na[[1]]$betadraw
    }
    initial.beta <- matrix(rapply(out.qr, function.mean, classes='matrix'), lgridp, lbeta, byrow=T)
    count.na <- which(is.na(initial.beta[,1]))
  }
  # Calculating quantile candidates (induced quantile posterior mean)
  initial.beta.sd <- matrix(rapply(out.qr, function.sd, classes='matrix'), lgridp, lbeta, byrow=T)
  initial.quant <- initial.beta %*%xg.matrix
  initial.sigma <- rep(0,lgridp)
  cov.beta <- array(0,dim=c(lgridp, lbeta, lbeta))
  var.sigma <- rep(0,lgridp)
  cov.beta.sigma <- array(0,dim=c(lgridp, lbeta))
  initial.quant.cov <- array(0,dim=c(lgridp, n, n))
  for(i in 1:lgridp){
    initial.sigma[i] <- mean((out.qr[[i]]$sigmadraw)[-(1:burnmcmc)])
    beta.sigma <- cbind(out.qr[[i]]$betadraw, out.qr[[i]]$sigmadraw)[-(1:burnmcmc),]
    cov.all <- cov(beta.sigma)
    cov.beta[i,,] <- cov.all[1:lbeta,1:lbeta]
    var.sigma[i] <- cov.all[lbeta+1,lbeta+1]
    cov.beta.sigma[i,] <- cov.all[lbeta+1,1:lbeta]
    initial.quant.cov[i,,] <- t(xg.matrix)%*%cov.beta[i,,]%*%xg.matrix
  }
  initial.quant.sd <- t(sqrt(apply(initial.quant.cov, 1, diag)))
  quant.cand <- array(0, dim=c(lgridtau, lgridp, lgridx))
  quant.cand.sd <-  array(0, dim=c(lgridtau, lgridp, lgridx))
  for (itau in 1:lgridtau){
    for (ip in 1:lgridp){
      if(gridp[ip]>=gridtau[itau]){
        quant.cand[itau,ip,] <- as.vector(initial.beta[ip,]%*%xg.matrix + initial.sigma[ip]/(1-gridp[ip])*log(gridtau[itau]/gridp[ip]))
        c <- 1/(1-gridp[ip])*log(gridtau[itau]/gridp[ip])
      }else{
        quant.cand[itau,ip,] <- as.vector(initial.beta[ip,]%*%xg.matrix - initial.sigma[ip]/(gridp[ip])*log((1-gridtau[itau])/(1-gridp[ip])))
        c <- -1/(gridp[ip])*log((1-gridtau[itau])/(1-gridp[ip]))
      }
      quant.cand.sd[itau, ip, ] <-  sqrt(diag(initial.quant.cov[ip,,]) + c^2*var.sigma[ip]+ 2*c*t(xg.matrix)%*%cov.beta.sigma[ip,])
    } 
  }
  CI.quant.yu <- array(0,c(lgridp,lgridx,2))
  CI.beta.yu <- array(0,c(lgridp,lbeta,2))
  for (ip in 1:lgridp){ 
    quantdraw2 <- (out.qr[[ip]]$betadraw)[-c(1:burnmcmc),]%*%xg.matrix
    CI.beta.yu[ip,,] <- t(apply(out.qr[[ip]]$betadraw[-c(1:burnmcmc),], 2, quantile, c(0.025, 0.975))) # initial CI for beta
    CI.quant.yu[ip,,] <- t(apply(quantdraw2, 2, quantile, c(0.025, 0.975))) } # initial CI for the quantiles 
  # Gaussian process regression adjustment
  res <- find.l(gridtau=gridtau, gridx=gridx, gridp=gridp, quant.cand=quant.cand,
                cov.or.var= quant.cand.sd^2/sample.size, sigma2.p, by.band, ess=sample.size, hzero, linear)
  band <- res[[1]]
  quant.final.mean <- res[[2]][ , ,1] # Q(tau,x)
  quant.final.sd <- sqrt(res[[2]][ , ,2]) # Q(tau,x) Variance of the predicted quantile
  beta.final <- array(0, dim=c(lgridtau, lbeta))
  beta0.final <- rep(0,lgridtau)
  ix=1
  for(itau in 1:lgridtau){
    beta.final[itau, ] <- sapply(1:lbeta, function(i) t(initial.beta[,i])%*%(res[[3]][itau,ix,]),simplify="array")
    c1 <- t(t(initial.sigma[gridp>=gridtau[itau]])*(log(gridtau[itau]/gridp[gridp>=gridtau[itau]])/(1-gridp[gridp>=gridtau[itau]])))
    c2 <- t(t(initial.sigma[gridp<gridtau[itau]])*(log((1-gridtau[itau])/(1-gridp[gridp<gridtau[itau]]))/gridp[gridp<gridtau[itau]]))
    beta0.final[itau] <- beta.final[itau,1] + t(rbind(-c2,c1))%*%(res[[3]][itau,ix,])
  }
  beta.final[,1] <- beta0.final
  beta.final.sd  <- initial.beta.sd
  CI.final <- array(0, dim=c(lgridtau, lgridx, 2))
  CI.final.beta <- array(0, dim=c(lgridtau, lbeta, 2))
  CI.final[,,1] <- quant.final.mean - qnorm(0.975)*quant.final.sd
  CI.final[,,2] <- quant.final.mean + qnorm(0.975)*quant.final.sd
  CI.final.beta[,,1] <- beta.final - qnorm(0.975)*beta.final.sd
  CI.final.beta[,,2] <- beta.final + qnorm(0.975)*beta.final.sd
  
  list(band=band, quant.final.mean=quant.final.mean, quant.final.sd=quant.final.sd, CI.final=CI.final,
       beta.final=beta.final, beta.final.sd=beta.final.sd, CI.final.beta=CI.final.beta,
       initial.quant=initial.quant, initial.quant.sd=initial.quant.sd, CI.quant.yu = CI.quant.yu,
       initial.beta=initial.beta, initial.beta.sd=initial.beta.sd, CI.beta.yu = CI.beta.yu)
}


################################################################
############ Auxiliary Functions ################################
################################################################

# Functions: function.mean, function.sd, rmse
function.mean <- function(x) apply(x[-(1:burnmcmc),], 2, mean)
function.sd <- function(x) apply(x[-(1:burnmcmc),], 2, sd)
rmse <- function(true,est) {sqrt(apply((true-est)^2,1,mean))}


# Call 'GPreg' to find the smallest bandwith such that the quantiles do not cross
find.l <- function(gridtau, gridx, gridp, quant.cand, cov.or.var, sigma2.p, by.band, ess, hzero=F, linear=T){
  dectest <- -1
  bandw <- 0
  lgridtau <- length(gridtau)
  lgridp <- lgridtau
  gridxt <- as.matrix(t(gridx))
  cov.meanx <- apply(cov.or.var, c(1,2), mean)
  if(linear==F) {pos.chull <-1:length(gridxt)
  }else{
    if(hzero==T) {pos.chull <-1:nrow(gridxt)
    }else{
      if(ncol(gridxt)==1){
        pos.chull <- c(which.min(gridxt),which.max(gridxt))
      }else{ 
        hull <- convhulln(gridxt)
        pos.chull <- sort(unique(c(hull)))
      }
    }
  }
  lgridx2 <- length(pos.chull)
  it <- 1
  while (sum(dectest<0)>0 | it==0) {
    if (it==0){
      pos.chull <- 1:nrow(gridxt)
      lgridx2 <- length(pos.chull)
      bandw <- bandw - by.band
    }
    q.tau.givenx <- array(0,c(lgridtau,lgridx2,3))
    w.gp <- array(0,c(lgridtau,lgridx2,lgridp))
    
    if(linear==T){
      for (ix in 1:lgridx2){  
        for (itau in 1:lgridtau){   
          if(nrow(cov.or.var)==lgridp^2){
            sigma2.n <- matrix(cov.or.var[,pos.chull[ix],itau], ncol=lgridp)
          }else{
            sigma2.n <- cov.meanx[itau,] 
            sigma.tau <- cov.or.var[itau,itau,pos.chull[ix]]*ess
          }
          res.GP <- GPreg(xdata=gridp, ydata=quant.cand[itau,,pos.chull[ix]], sigma2.n=sigma2.n, x.star=gridtau[itau], sigma2.p, lb=bandw, sigma.tau=sigma.tau)
          q.tau.givenx[itau,ix,] <- res.GP[[1]]
          w.gp[itau,ix,] <- res.GP[[2]]
        } 
      }
    }else{
      for (ix in 1:lgridx2){  
        for (itau in 1:lgridtau){   
          if(nrow(cov.or.var)==lgridp^2){
            sigma2.n <- matrix(cov.or.var[,pos.chull[ix],itau], ncol=lgridp)
          } else {
            sigma2.n <- cov.or.var[itau,,pos.chull[ix]]
            sigma.tau <- cov.or.var[itau,itau,pos.chull[ix]]*ess
          }
          res.GP <- GPreg(xdata=gridp, ydata=quant.cand[itau,,pos.chull[ix]], sigma2.n=sigma2.n, x.star=gridtau[itau], sigma2.p, lb=bandw, sigma.tau=sigma.tau)
          q.tau.givenx[itau,ix,] <- res.GP[[1]]
          w.gp[itau,ix,] <- res.GP[[2]]
        } 
      }
    }  
    if(bandw>1) by.band <- 0.5
    if(bandw>10) by.band <- 1
    if(bandw>30) by.band <- 10
    if(bandw>200) by.band <- 50
    if(bandw>1000) by.band <- 200
    bandw <- bandw + by.band
    if(hzero==F){
      dectest<- q.tau.givenx[2:lgridtau,,1]-q.tau.givenx[1:(lgridtau-1),,1]
      if(sum(dectest<0)==0 & it==1 & linear==T){it<-0
      }else it<-1
    }else {dectest <- 1; it <- 1}
  }
  bandw <- bandw-by.band
  list(l=bandw, q.tau.givenx=q.tau.givenx, w.gp=w.gp)
}


# Second stage of the GPR: calculate the adjusted posterior mean for quantile tau given x
GPreg <- function(xdata, ydata, sigma2.n, x.star, sigma2.p, lb, sigma.tau){
  options = gpOptions()
  options$kern$comp = list("rbf")
  yTrain = ydata
  xTrain = xdata
  kern = kernCreate(xTrain, "rbf")
  kern$variance = sigma2.p
  lb <- ifelse(lb==0,10^(-10),lb)
  kern$inverseWidth = 1/lb^2
  xTest = x.star
  pos = which.min(abs(xTrain-xTest))
  Kx = kernCompute(kern, xTest, xTrain)
  Ktrain = kernCompute(kern, xTrain)
  if(is.matrix(sigma2.n)){ 
    noiseVar = sigma2.n
  } else {
    noiseVar = diag(dim(Ktrain)[1])*sigma2.n
  }
  invKtrain = (.jitCholInv(Ktrain + noiseVar , silent = TRUE))[[1]] 
  weights = Kx %*% invKtrain
  yPred =  weights %*% yTrain
  meanVar = kernDiagCompute(kern, xTest) - rowSums(Kx %*% invKtrain * Kx)
  yPredVar = meanVar + sigma.tau
  list(y.gp=c(yPred, yPredVar, meanVar), w.gp=weights)
}


# Matrix inversion for GPreg
.jitCholInv <- 
  function (M, Num = 10, silent = FALSE) 
  {
    jitter <- 0
    jitter1 <- abs(mean(diag(M))) * 1e-06
    eyeM <- diag(1, nrow = length(M[, 1]), ncol = length(M[1, 
                                                           ]))
    for (i in 1:Num) {
      try(stop(""), TRUE)
      Ch <- try(chol(M + jitter * eyeM), silent = TRUE)
      nPos <- grep("not positive definite", geterrmessage())
      if (length(nPos) != 0) {
        jitter1 <- jitter1 * 10
        jitter <- jitter1
        if (!silent) {
          warnmsg <- paste("Matrix is not positive definite, adding", 
                           signif(jitter, digits = 4), "jitter!")
          warning(warnmsg)
        }
      }
      else break
    }
    invCh <- try(solve(Ch, eyeM), silent = TRUE)
    if (class(invCh) == "try-error") {
      return(NaN)
    }
    else {
      invM <- invCh %*% t(invCh)
      if (jitter == 0) {
        ans <- list(invM = invM, jitter = jitter, chol = Ch)
      }
      else ans <- list(invM = invM, jitM = M + jitter * eyeM, 
                       jitter = jitter, chol = Ch)
      return(ans)
    }
  }

