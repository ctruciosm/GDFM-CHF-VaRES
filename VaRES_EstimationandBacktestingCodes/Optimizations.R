### Optimizations
#JSF = function(VaR, ES, r, alpha){
#  (ifelse(r<= VaR,1,0) - alpha)*G1(VaR) 
#  - ifelse(r<=VaR,1,0)*G1(r)
#  + G2(ES)*(ES - VaR + ifelse(r<=VaR,1,0)*(VaR - r)/alpha) 
#  - G3(ES) + a(r)
#}

# Scoring Functions


QL = function(VaR, r, alpha){
  D = dim(VaR)
  if (is.null(D)){
    k = 1
    val_ =  (alpha - ifelse(r<=VaR,1,0))*(r-VaR)
  } else {
    k = D[2]
    n = D[1]
    val_ = matrix(0,ncol = k, nrow = n)
    for (i in 1:k){
      val_[,i] = (alpha - ifelse(r<=VaR[,i],1,0))*(r-VaR[,i])
    }
  }
  return(val_)
}

AL = function(VaR, ES, r, alpha){
  D = dim(VaR)
  if (is.null(D)){
    k = 1
    val_ = ((-1/ES)*(ES - VaR + ifelse(r<=VaR,1,0)*(VaR - r)/alpha) - (-log(-ES)) + (1-log(1-alpha)))
  } else {
    k = D[2]
    n = D[1]
    val_ = matrix(0,ncol = k, nrow = n)
    for (i in 1:k){
    val_[,i] = ((-1/ES[,i])*(ES[,i] - VaR[,i] + ifelse(r<=VaR[,i],1,0)*(VaR[,i] - r)/alpha) - (-log(-ES[,i])) + (1-log(1-alpha)))
    }
  }
  return(val_)
}

FZ0 = function(VaR, ES, r, alpha){
  D = dim(VaR)
  if (is.null(D)){
    k = 1
    val_ = ((-1/ES)*(ES - VaR + ifelse(r<=VaR,1,0)*(VaR - r)/alpha) - (-log(-ES)))
  } else {
    k = D[2]
    n = D[1]
    val_ = matrix(0,ncol = k, nrow = n)
    for (i in 1:k){
      val_[,i] = ((-1/ES[,i])*(ES[,i] - VaR[,i] + ifelse(r<=VaR[,i],1,0)*(VaR[,i] - r)/alpha) - (-log(-ES[,i])))
    }
  }
  return(val_)
}

NZ = function(VaR, ES, r, alpha){
  D = dim(VaR)
  if (is.null(D)){
    k = 1
    val_ = ((1/(2*sqrt(-ES)))*(ES - VaR + ifelse(r<=VaR,1,0)*(VaR - r)/alpha) + sqrt(-ES))
  } else {
    k = D[2]
    n = D[1]
    val_ = matrix(0,ncol = k, nrow = n)
    for (i in 1:k){
      val_[,i] = ((1/(2*sqrt(-ES[,i])))*(ES[,i] - VaR[,i] + ifelse(r<=VaR[,i],1,0)*(VaR[,i] - r)/alpha) + sqrt(-ES[,i]))
    }
  }
  return(val_)
}

FZG = function(VaR, ES, r, alpha){
  D = dim(VaR)
  if (is.null(D)){
    k = 1
    val_ = ((ifelse(r<= VaR,1,0) - alpha)*VaR - ifelse(r<=VaR,1,0)*r+ (exp(ES)/(1+exp(ES)))*(ES - VaR + ifelse(r<=VaR,1,0)*(VaR - r)/alpha) - (log(1+exp(ES))) + log(2))
  } else {
    k = D[2]
    n = D[1]
    val_ = matrix(0,ncol = k, nrow = n)
    for (i in 1:k){
      val_[,i] = ((ifelse(r<= VaR[,i],1,0) - alpha)*VaR[,i] - ifelse(r<=VaR[,i],1,0)*r+ (exp(ES[,i])/(1+exp(ES[,i])))*(ES[,i] - VaR[,i] + ifelse(r<=VaR[,i],1,0)*(VaR[,i] - r)/alpha) - (log(1+exp(ES[,i]))) + log(2))
    }
  }
  return(val_)
}


#########################################################
####### Relative Score Combining 
#########################################################
RSC = function(lambda, VaR, ES, r, alpha, S){
  if(length(r) == 1){
    M = length(VaR) 
    VaR = matrix(VaR, ncol = M, nrow = 1)
    ES = matrix(ES, ncol = M, nrow = 1)
  } else {
    if(is.data.frame(VaR)) VaR = as.matrix(VaR)
    if(is.data.frame(ES)) ES = as.matrix(ES)
    M = dim(VaR)[2]
  }

  omega =  matrix(0,ncol = 1, nrow = M)
  for (i in 1:M){
    omega[i,1] = exp(-lambda * sum(S(VaR[,i], ES[,i], r, alpha)))/sum(exp(-lambda * apply(S(VaR, ES, r, alpha),2,sum) ))
  }

  if(!is.nan(sum(sum(omega)))) {
    VaR_c = VaR%*%omega
    ES_c = ES%*%omega
    return(sum(S(VaR_c, ES_c, r, alpha)))
  } else{
    return(Inf)
  }
}

RSC_Eval = function(lambda, VaR, ES, r, alpha, S){
  if(length(r) == 1){
    M = length(VaR) 
    VaR = matrix(VaR, ncol = M, nrow = 1)
    ES = matrix(ES, ncol = M, nrow = 1)
  } else {
    if(is.data.frame(VaR)) VaR = as.matrix(VaR)
    if(is.data.frame(ES)) ES = as.matrix(ES)
    M = dim(VaR)[2]
  }
  
  omega =  matrix(0,ncol = 1, nrow = M)
  for (i in 1:M){
    omega[i,1] = exp(-lambda * sum(S(VaR[,i], ES[,i], r, alpha)))/sum( exp(-lambda * apply(S(VaR, ES, r, alpha),2,sum) ))
  }
  
return(omega)
}

RSC_grid = function(parini, VaR, ES, r, alpha, S){
  lambda = seq(from = 0.000001, to = 100, length.out = 10^4)
  val = c()
  for (i in 1:(length(lambda))){
    val[i] = RSC(lambda[i], VaR, ES, r, alpha, S)
  }
  lini = lambda[which(val == min(val))[1]]
  
  if (RSC(parini, VaR, ES, r, alpha, S) < RSC(lini, VaR, ES, r, alpha, S)) lini = parini
  
  return(c(lambda[which(val == Inf)[1]-1],lini))
}

RSC_opt = function(parini, VaR, ES, r, alpha, S){
  parini = RSC_grid(parini, VaR, ES, r, alpha, S)
  lambda = suppressWarnings(optim(par = parini[2], fn = RSC, method = "L-BFGS-B", VaR = VaR, ES = ES, r = r, alpha = alpha, S = S, lower = 0.000001, upper = parini[1])$par)
  return(lambda)
}

#########################################################
####### Minimum Score Combining 
#########################################################
equal <- function(omega, VaR, ES, r, alpha, S) {
  M = length(omega)
  N = M/2
  return(c(sum(omega[1:N]), sum(omega[(N+1):M])))
}

MSC = function(omega, VaR, ES, r, alpha, S){
  N = dim(VaR)[2]
  VaR_c = VaR%*%omega[1:N]
  ES_c = VaR_c + (ES-VaR)%*%omega[(N+1):(2*N)] 
  return(sum(S(VaR_c, ES_c, r, alpha)))
}

MSC_grid = function(parini, VaR, ES, r, alpha, S){
  N = dim(VaR)[2]
  omega = parini
  VaR_c = VaR%*%parini[1:N]
  ES_c = VaR_c + (ES-VaR)%*%parini[(N+1):(2*N)] 
  SF = sum(S(VaR_c, ES_c, r, alpha))
  
  for(i in 1:10^5){
    omega[1:N] = runif(N)
    omega[1:(N-1)] = omega[1:(N-1)]/sum(omega[1:N])
    omega[N] = 1-sum(omega[1:(N-1)])
    
    omega[(N+1):(2*N)]  = runif(N)
    omega[(N+1):(2*N-1)]  = omega[(N+1):(2*N-1)]/sum(omega[(N+1):(2*N)])
    omega[2*N] = 1-sum(omega[(N+1):(2*N-1)])
    
    VaR_c = VaR%*%omega[1:N]
    ES_c = VaR_c + (ES-VaR)%*%omega[(N+1):(2*N)] 
    
    if(sum(S(VaR_c, ES_c, r, alpha)) < SF){
      parini = omega
      SF = sum(S(VaR_c, ES_c, r, alpha))
    }
  }
  return(parini)
}

MSC_opt = function(parini, VaR, ES, r, alpha, S){
  N = dim(VaR)[2]
  parini = MSC_grid(parini, VaR, ES, r, alpha, S)
  param = solnp(pars = parini, fun = MSC, eqfun = equal, eqB = c(1,1), LB = rep(0,2*N), VaR = VaR, ES = ES, r = r, alpha = alpha, S = S)$pars
  return(param)
}




