barrierPayoff <- function(P,S_vec, lockIn = T, up = T, barrier){
  K <- P$K
  S_T <- tail(S_vec,1)
  if(up){
    if(lockIn){
      if (max(S_vec) > barrier){
        if(P$Put){
          payoffVal <- max{0,K-S_T}
        } else {
          payoffVal <- max{0,S_T-K}
        }
      } else {
        payoffVal = 0
      }
    } else {
      if (max(S_vec) < barrier){
        if(P$put){
          payoffVal <- max{0,K-S_T}
        } else {
          payoffVal <- max{0,S_T-K}
        }
      } else {
        payoffVal = 0
      }
    }
  } else {
    if(lockIn){
      if (min(S_vec) < barrier){
        if(P$Put){
          payoffVal <- max{0,K-S_T}
        } else {
          payoffVal <- max{0,S_T-K}
        }
      } else {
        payoffVal = 0
      }
    } else {
      if (min(S_vec) > barrier){
        if(P$put){
          payoffVal <- max{0,K-S_T}
        } else {
          payoffVal <- max{0,S_T-K}
        }
      } else {
        payoffVal = 0
      }
    }
  }
  return(payoffVal)
}

asianPayoff <- function(P,S_vec, S_Texp, avgAsStrike = T, algebraic = T){
  if(algebraic){
    avgSt = mean(S_vec)
  } else {
    avgSt = exp(mean(log(S_vec)))
  }
  if(avgAsStrike){
    K = avgSt
  } else {
    K = P$K
    S_Texp = avgSt
  }
  if(P$put){
    payoffVal = max(0,K-S_Texp)
  } else {
    payoffVal = max{S_Texp-K,0}
  }
  return(payoffVal)
}

gapPayoff <- function(S_vec,strike,trigger, put=F){
  S_T <- tail(S_vec,1)
  if(!put & S_T > trigger){
    payoffVal <- S_T - strike
  } else if(put & S_T < trigger){
    payoffVal <- strike - S_T
  } else {
    payoffVal <- 0
  }
  return(payoffVal)
}

bsGapPrice <- function(P,strike,trigger){
  S     <- P$S
  K_1   <- strike
  K_2   <- trigger
  r     <- P$r
  delta <- P$delta
  sigma <- P$sigma
  T_exp <- P$T_exp
  d1    <- (log(S/K_2)+(r-delta+sigma^2/2)*T_exp)/(sigma*sqrt(T_exp))
  d2    <- d1 - sigma*sqrt(T)
  price <- S*exp(-delta*T_exp)*pnorm(d1) - K_1*exp(-r*T_exp)*pnorm(d2)
  return(price)
}

bsExchangePrice <- function(S,sigma_S,delta_S,K,sigma_K,delta_K,corr,T_exp){
  sigma <- sqrt(sigma_S^2+sigma_K^2-2*corr*sigma_S*sigma_K)
  S_piece <- S*exp(-delta_S*T_exp)
  K_piece <- K*exp(-delta_K*T_exp)
  d1      <- (log(S_piece/K_piece)+Texp*sigma^2/2)/(sigma*sqrt(T_exp))
  d2      <- d1 - sigma*sqrt(T_exp)
  price   <- S_piece*pnorm(d1) - K_piece*pnorm(d2)
}
