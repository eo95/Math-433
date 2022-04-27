P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=25,K=46,sigma=0.29,choice=1,D_CF=c(22/252,.09),recombine=F)
binomial_pricing(P,call_payoff)[[2]][1]

P <- bs_parameterize(S=46.01,r=0.0021,T_exp=29/252,K=46,sigma=0.29,D_CF=c(.09,22/252))
black_scholes_div(P)

# Building Table 14.3

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

payoffs = list("aA"=c(),"bB"=c(),"cC"=c(),"dD"=c(),"eE"=c(),"fF"=c(),"gG"=c(),"hH"=c(),"iI"=c(),"jJ"=c(),"kK"=c(),"lL"=c(),"mM"=c(),"nN"=c(),"oO"=c())
prices  = rep(0,15)

for(i in 1:N){
  pricerun = blackscholesPriceRunSim(P,(0:126)/256)
  P$K = 0.8
  payoffaA = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payoffbB = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffcC = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffdD = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffeE = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  P$K = 0.9
  payofffF = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payoffgG = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffhH = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffiI = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffjJ = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  P$K = 1.0
  payoffkK = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payofflL = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffmM = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffnN = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffoO = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  
  payoffs$aA = c(payoffs$aA,payoffaA)
  payoffs$bB = c(payoffs$bB,payoffbB)
  payoffs$cC = c(payoffs$cC,payoffcC)
  payoffs$dD = c(payoffs$dD,payoffdD)
  payoffs$eE = c(payoffs$eE,payoffeE)
  payoffs$fF = c(payoffs$fF,payofffF)
  payoffs$gG = c(payoffs$gG,payoffgG)
  payoffs$hH = c(payoffs$hH,payoffhH)
  payoffs$iI = c(payoffs$iI,payoffiI)
  payoffs$jJ = c(payoffs$jJ,payoffjJ)
  payoffs$kK = c(payoffs$kK,payoffkK)
  payoffs$lL = c(payoffs$lL,payofflL)
  payoffs$mM = c(payoffs$mM,payoffmM)
  payoffs$nN = c(payoffs$nN,payoffnN)
  payoffs$oO = c(payoffs$oO,payoffoO)
}
for(i in 1:15){
  prices[i] = exp(-0.06*0.5)*mean(payoffs[[i]])
}
table = matrix(prices,3,5,byrow=T)
table



# Black Scholes Eqn vs Black Scholes Price Run Sim
N = 10000
P.Call <- list(S=277.5,K=280,r=.004452,T_exp=52/252,put=F,sigma=0.24,delta=0)
payoffs <- c()
for(i in 1:N){
  priceRun <- blackscholesPriceRunSim(P.Call,c(0,52/252))
  S_T <- tail(priceRun,1)
  K   <- P.Call$K
  payoff = max(0,S_T-K)
  payoffs = c(payoffs,payoff)
}
price <- exp(-52/252*.004452)*mean(payoffs)
price
black_scholes(P.Call)
