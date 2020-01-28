# Function converts from monthly to quarterly.
## Data must start from month 1 for this function to work 
getMonths <- function(z){return(as.integer(12*(time(z)+0.01-trunc(time(z)))+1L))}
sumvals <- function(x){
  m <- getMonths(x)
  print(m[1])
  if(!(m[1] %in% c(1L,4L,7L,10L))){
    return(sumvals(window(x,time(x)[2])))
  } else {
    return(aggregate(x,nfrequency=4))
  }
}

## Block Bootstrap function

##Multiple horizon
Mult.boot <- function(y,x,size,beg,freq,horizon) {
  N <- length(y)
  num_blocks <- (N/size)
  holder <- matrix(NA, nrow=N, ncol=(size+size))
  q <- seq(from = 1, to = ((N-size)+1), by = size)
  for (i in q){
    holder[i,1:size] <- y[i:(i+(size-1))]
    holder[i,(size+1):(size+size)] <- x[i:(i+(size-1))]
  }
  holder <- na.omit(holder)
  k <- 1:num_blocks
  newdata_y <- matrix(NA, nrow=N, ncol=1)
  newdata_x <- matrix(NA, nrow=N, ncol=1)
  for (i in q){
    w <- sample(k, 1)
    newdata_y[i:(i+(size-1)),1] <- holder[w,1:size]
    newdata_x[i:(i+(size-1)),1] <- holder[w,(size+1):(size+size)]
  }
  ts_y <- ts(newdata_y, start=beg, frequency=freq)
  ts_x <- ts(newdata_x, start=beg, frequency=freq)
  block.reg <- dynlm(ts_y ~ L(ts_x, 0:horizon))
  block.irf <- matrix(NA, nrow=(horizon+1), ncol=1)
  for (i in 1:(horizon+1)){
    block.irf[i,1] <- (sum(coef(block.reg)[2:(1+i)]))
  }
  return(block.irf)
}

## Bootstrap function
bootfn <- function(fit,shocks) {                                            
  f <-  fitted(fit)
  res <- residuals(fit)
  z <- sample(c(-1,1), length(res), replace=TRUE)
  simdata <- f + z*res
  return(cumsum(coef(dynlm(simdata ~ L(shocks,0:10)))))
}

intfit <- function(r){
  sfit <- dynlm(r ~ L(shocks[,"Supply"],0:10))
  dfit <- dynlm(r ~ L(shocks[,"Demand"],0:10))
  ofit <- dynlm(r ~ L(shocks[,"Oil"],0:10))
  jfit <- dynlm(r ~ L(shocks[,"Jet"],0:10))
  replicate(1000, {list(supply=bootfn(sfit,shocks[,"Supply"]),demand=bootfn(dfit,shocks[,"Demand"]), oil=bootfn(ofit,shocks[,"Oil"]), jet=bootfn(jfit,shocks[,"Jet"]))}, simplify=FALSE)
}

#### Functions for Historical Decompositions
VARhd <- function(Estimation){
  
  ## make X and Y
  nlag    <- Estimation$p   # number of lags
  DATA    <- Estimation$y   # data
  QQ      <- VARmakexy(DATA,nlag,1)
  
  
  ## Retrieve and initialize variables 
  invA    <- t(chol(as.matrix(summary(Estimation)$covres)))   # inverse of the A matrix
  Fcomp   <- companionmatrix(Estimation)                      # Companion matrix
  
  #det     <- c_case                                           # constant and/or trends
  F1      <- t(QQ$Ft)                                         # make comparable to notes
  eps     <- ginv(invA) %*% t(residuals(Estimation))          # structural errors 
  nvar    <- Estimation$K                                     # number of endogenous variables
  nvarXeq <- nvar * nlag                                      # number of lagged endogenous per equation
  nvar_ex <- 0                                                # number of exogenous (excluding constant and trend)
  Y       <- QQ$Y                                             # left-hand side
  #X       <- QQ$X[,(1+det):(nvarXeq+det)]                    # right-hand side (no exogenous)
  nobs    <- nrow(Y)                                          # number of observations
  
  
  ## Compute historical decompositions
  
  # Contribution of each shock
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(nlag-1)*nvar))
  HDshock_big <- array(0, dim=c(nlag*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  # for each variable
    eps_big <- matrix(0,nvar,(nobs+1)) # matrix of shocks conformable with companion
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
    
  } 
  
  HD.shock <- array(0, dim=c((nobs+nlag),nvar,nvar))   # [nobs x shock x var]
  
  for (i in 1:nvar){
    
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,nlag), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD.shock)
  
}


VARmakexy <- function(DATA,lags,c_case){
  
  nobs <- nrow(DATA)
  
  #Y matrix 
  Y <- DATA[(lags+1):nrow(DATA),]
  Y <- DATA[-c(1:lags),]
  
  #X-matrix 
  if (c_case==0){
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    } 
  } else if(c_case==1){ #constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    X <- cbind(matrix(1,(nobs-lags),1), X) 
  } else if(c_case==2){ # time trend and constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    trend <- c(1:nrow(X))
    X <-cbind(matrix(1,(nobs-lags),1), t(trend))
  }
  A <- (t(X) %*% as.matrix(X)) 
  B <- (as.matrix(t(X)) %*% as.matrix(Y))
  
  Ft <- ginv(A) %*% B
  
  retu <- list(X=X,Y=Y, Ft=Ft)
  return(retu)
}

companionmatrix <- function (x) 
{
  if (!(class(x) == "varest")) {
    stop("\nPlease provide an object of class 'varest', generated by 'VAR()'.\n")
  }
  K <- x$K
  p <- x$p
  A <- unlist(Acoef(x))
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  return(companion)
}

GIRF <- function(res, hor=20, shk=1, shVar = 1, replic=100) {
  #res...estimation result from TVAR
  #hor...impulse response horizon
  #shk...shock size (in standard deviations)
  #shVar...shocked Variable
  #replic...no. bootstrap replications
  
  
  threshV <- which(res$model.specific$transCombin == 1) #which variable is used as threshold variable
  datS <- res$model[,1:res$k] #get dataset
  
  avgDiff <- array(dim = c(hor, ncol(datS), (nrow(datS) - res$lag)), NA) #results for each history will be saved here
  #simulate over all histories (regime dependence is dealt with in the next step)
  for (i in 1 : (nrow(datS) - res$lag)) {
    resultDiff <- lapply(1:replic, simTVAR, res, datS, hor, shk, threshV, shVar, i) #call simTVAR (second function)
    avgDiff[, , i] <- Reduce("+", resultDiff) / replic #mean of particular history
  }
  
  
  #Regimes
  regimes <- regime(res)[(res$lag + 1) : length(regime(res))]
  
  #GIRF for regimes (final result)
  girf <- array(dim = c(hor, ncol(datS), max(regimes)), NA)
  for (i in 1 : max(regimes)) {
    selectReg <- avgDiff[, , which(regimes == i)]
    girf[, , i] <- apply(selectReg, MARGIN=c(1, 2), sum) / dim(selectReg)[3]
  }
  
  return(girf)
}  


#k=1; results=res;dat=datS;horizon=hor;sh=shk;shvar=shVar;history=1;
simTVAR <- function(k, results, dat, horizon, sh, threshV, shvar, history) {
  #k...for apply function (index)
  #dat...data
  #horizon...impulse response horizon
  #sh...shock size (in standard deviations)
  #threshV...which variable is used as threshold Variable
  #shVar...shocked Variable
  #history...history of variable
  
  
  #sample bootstrap residuals
  resid <- matrix(NA, nrow=horizon, ncol=results$k)
  for (i in 1:results$k) resid[,i] <- sample(results$residuals[,i], size=horizon, replace=TRUE)
  #resid <- sample(results$residuals, size=horizon * nrow(results$coeffmat), replace=TRUE)
  #dim(resid) <- c(horizon,nrow(results$coeffmat))
  
  #bootstrap residuals for shocked series (same residuals with additional shock at the beginning)
  shock <- sqrt(var(results$residuals[, shvar])) * sh #shock at t
  shock <- c(rep(0, shvar-1), shock, rep(0,nrow(results$coeffmat)-shvar))
  resid_delta <- rbind(resid[1,] + shock, resid[-1,])
  
  #simulation without addtional shock --> innov = resid
  simul <- TVAR.sim(B = results$coeffmat, 
                    Thresh = results$model.specific$Thresh, 
                    nthres = results$model.specific$nthresh,
                    n = horizon, lag = results$lag, include = results$include,
                    thDelay = results$model.specific$thDelay, mTh = threshV,
                    starting = dat[history : (history + results$lag - 1), ], innov = resid)
  
  #with resid_delta
  simul_delta <- TVAR.sim(B = results$coeffmat, 
                          Thresh = results$model.specific$Thresh, 
                          nthres = results$model.specific$nthresh,
                          n = horizon, lag = results$lag, include = results$include,
                          thDelay = results$model.specific$thDelay, mTh = threshV,
                          starting = dat[history : (history + results$lag - 1), ], innov = resid_delta)
  
  diff <- simul_delta - simul
  
  return(diff)  
}

