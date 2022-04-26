black_scholes <- function(P){
  #Black Scholes with or without dividend yeild
  # P is the parameters of the model
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!P$put){
    price <- S*exp(-delta*T_exp) * pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
  } else {
    price <- -S*exp(-delta*T_exp) * pnorm(-d1) + K*exp(-r*T_exp)*pnorm(-d2)
  }
  return(price)
}

blackscholesPriceRunSim <- function(P,time_v){
  # P is output of parameterize with the defined parameters for the binomial model
  # time_v is a vector of desired time values of S
  ##ex: time_v <- c(0,25,30,70) would result in the prices of S at time 0, 25, 30, 70
  ## time_v should always include initial time 0
  steps <- length(time_v)-1
  diff <- diff(time_v)
  r <- P$r[1]
  delta <- P$delta[1]
  sigma <- P$sigma[1]
  S_v <- c(P$S)
  for(i in 1:steps){
    S_a <- tail(S_v,1)
    mean = (r-delta-sigma^2/2)*diff[i]
    variance = sigma^2*diff[i]
    sd = sqrt(variance)
    S_b <- S_a*rlnorm(1,mean,sd)
    S_v <- c(S_v,S_b)
  }
  return(S_v)
}

barrierPayoff <- function(P,S_vec, lockIn = T, up = T, barrier){
  K <- P$K
  S_T <- tail(S_vec,1)
  if(up){
    if(lockIn){
      if (max(S_vec) > barrier){
        if(P$put){
          payoffVal <- max(0,K-S_T)
        } else {
          payoffVal <- max(0,S_T-K)
        }
      } else {
        payoffVal = 0
      }
    } else {
      if (max(S_vec) < barrier){
        if(P$put){
          payoffVal <- max(0,K-S_T)
        } else {
          payoffVal <- max(0,S_T-K)
        }
      } else {
        payoffVal = 0
      }
    }
  } else {
    if(lockIn){
      if (min(S_vec) < barrier){
        if(P$put){
          payoffVal <- max(0,K-S_T)
        } else {
          payoffVal <- max(0,S_T-K)
        }
      } else {
        payoffVal = 0
      }
    } else {
      if (min(S_vec) > barrier){
        if(P$put){
          payoffVal <- max(0,K-S_T)
        } else {
          payoffVal <- max(0,S_T-K)
        }
      } else {
        payoffVal = 0
      }
    }
  }
  return(payoffVal)
}



# McDonald 14.3 Table

P <- list(S=0.9,sigma=0.1,r=0.06,delta=0.03,T_exp=0.5, put = T)
N <- 100000
Strikes = c(0.8,0.9,1.0)

payoffs = c()
for(i in 1:N){
  P$K      = 0.8
  pricerun = blackscholesPriceRunSim(P,(0:126)/126)
  payoff   = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.8)
  payoffs = c(payoffs,payoff)
}
price = exp(-0.06*0.5)*mean(payoffs)

for(i in 1:3){
  P$K = Strikes[i]
  standard_prices[i] = black_scholes(P)
  DI_B0.80_prices[i] = exp(-P$r*P$T_exp)*mean(DI_B0.80_payoffs[i,])
  DI_B0.85_prices[i] = exp(-P$r*P$T_exp)*mean(DI_B0.85_payoffs[i,])
  UO_B0.95_prices[i] = exp(-P$r*P$T_exp)*mean(UO_B0.95_payoffs[i,])
  UO_B1.00_prices[i] = exp(-P$r*P$T_exp)*mean(UO_B1.00_payoffs[i,])
  UO_B1.05_prices[i] = exp(-P$r*P$T_exp)*mean(UO_B1.05_payoffs[i,])
}

table14.3 = cbind(Strikes,standard_prices,DI_B0.80_prices,DI_B0.85_prices,UO_B0.95_prices,UO_B1.00_prices,UO_B1.05_prices)
table14.3
