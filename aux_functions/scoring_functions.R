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

