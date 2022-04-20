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
    payoffVal = max(S_Texp-K,0)
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

asianAveragePriceOptions_bs <- function(P,N="Infinity"){
  # Uses Equations in Appendix 14.A to price an asian option.
  # Asian option calculated with Geometric Average as S_T, and at N evenly spaced time intervals in [0,T_exp]
  r     = P$r
  delta = P$delta
  sigma = P$sigma
  if(N == "Infinity"){
    deltaNew = (r+delta+sigma^2/6)/2
    sigmaNew = sigma * sqrt(1/3)
  } else {
    deltaNew = 0.5*(r*(N-1)/N+(delta+sigma^2/2)*(N+1)/N-sigma^2/N^2*(N+1)*(2*N+1)/(6))
    sigmaNew = sigma/N*sqrt((N+1)*(2*N+1)/(6))
  }
  P$delta = deltaNew
  P$sigma = sigmaNew
  return(black_scholes(P))
}

asianAverageStrikeOptions_bs <- function(P,N="Infinity"){
  # Uses Equations in Appendix 14.A to price an asian option.
  # Asian option calculated with Geometric Average as K, and at N evenly spaced time intervals in [0,T_exp]
  if(N=="Infinity"){N = 10^8}
  r        = P$r
  delta    = P$delta
  sigma    = P$sigma
  rNew     = 0.5*(r*(N-1)/N+(delta+sigma^2/2)*(N+1)/N-sigma^2/N^2*(N+1)*(2*N+1)/(6))
  corr     = 0.5*sqrt(6*(N+1)/(2*N+1))
  piece    = (N+1)*(2*N+1)/(6)/(N^2)
  sigmaNew = sigma*sqrt(1+piece-2*corr*sqrt(piece))
  strikeNew = P$S
  P$r     = rNew
  P$sigma = sigmaNew
  P$K     = strikeNew
  return(black_scholes(P))
}

# Uses two functions above along with black_scholes() to create table 14.1 in McDonald

P.Call <- list(S=40,K=40,r=.08,sigma=0.3,delta=0,T_exp=1,put=F)
P.Put <- list(S=40,K=40,r=.08,sigma=0.3,delta=0,T_exp=1,put=T)
theNs <- c(1,2,3,5,10,50,1000)
AveragePriceCallPrice  <- c()
AveragePricePutPrice   <- c()
AverageStrikeCallPrice <- c()
AverageStrikePutPrice  <- c()
for(i in theNs){
  AvgPriceCallPrem   <-  asianAveragePriceOptions_bs(P.Call,N=i)
  AvgPricePutPrem    <-  asianAveragePriceOptions_bs(P.Put, N=i)
  AvgStrikeCallPrem  <- asianAverageStrikeOptions_bs(P.Call,N=i)
  AvgStrikePutPrem   <- asianAverageStrikeOptions_bs(P.Put, N=i)
  AveragePriceCallPrice <- c(AveragePriceCallPrice,AvgPriceCallPrem)
  AveragePricePutPrice   <- c(AveragePricePutPrice,AvgPricePutPrem)
  AverageStrikeCallPrice <- c(AverageStrikeCallPrice,AvgStrikeCallPrem)
  AverageStrikePutPrice  <- c(AverageStrikePutPrice,AvgStrikePutPrem)
}
AvgPriceCallPrem   <-  asianAveragePriceOptions_bs(P.Call,N="Infinity")
AvgPricePutPrem    <-  asianAveragePriceOptions_bs(P.Put, N="Infinity")
AvgStrikeCallPrem  <- asianAverageStrikeOptions_bs(P.Call,N="Infinity")
AvgStrikePutPrem   <- asianAverageStrikeOptions_bs(P.Put, N="Infinity")
AveragePriceCallPrice <- c(AveragePriceCallPrice,AvgPriceCallPrem)
AveragePricePutPrice   <- c(AveragePricePutPrice,AvgPricePutPrem)
AverageStrikeCallPrice <- c(AverageStrikeCallPrice,AvgStrikeCallPrem)
AverageStrikePutPrice  <- c(AverageStrikePutPrice,AvgStrikePutPrem)

theNs <- c(theNs,"Infinity")
table14.1 <- data.frame(theNs,AveragePriceCallPrice,AveragePricePutPrice,AverageStrikeCallPrice,AverageStrikePutPrice)
colnames(table14.1) <- c("N","AvgPriceCall","AvgPricePut","AvgStrikeCall","AvgStrikePut")
table14.1

# Using the BS assumptions to simulate the asian option price. Demonstrating the function and the simulator are working.

AvgPriceCallPrem   <- 0
AvgPricePutPrem    <- 0
AvgStrikeCallPrem  <- 0
AvgStrikePutPrem   <- 0
for(i in 1:10000){
  S_vec <- blackscholesPriceRunSim(P.Call,time_v = (0:1000)/1000)
  AvgPriceCallPrem   <- AvgPriceCallPrem  + exp(-.08)*asianPayoff(P.Call,S_vec,S_Texp=tail(S_vec,1),avgAsStrike = F,algebraic=F)
  AvgPricePutPrem    <- AvgPricePutPrem   + exp(-.08)*asianPayoff(P.Put,S_vec,S_Texp=tail(S_vec,1),avgAsStrike = F,algebraic=F)
  AvgStrikeCallPrem  <- AvgStrikeCallPrem + exp(-.08)*asianPayoff(P.Call,S_vec,S_Texp=tail(S_vec,1),avgAsStrike = T,algebraic=F)
  AvgStrikePutPrem   <- AvgStrikePutPrem  + exp(-.08)*asianPayoff(P.Put,S_vec,S_Texp=tail(S_vec,1),avgAsStrike = T,algebraic=F)
}
AveragePriceCallPrice  <- c(AveragePriceCallPrice, AvgPriceCallPrem/10000)
AveragePricePutPrice   <- c(AveragePricePutPrice,  AvgPricePutPrem/10000)
AverageStrikeCallPrice <- c(AverageStrikeCallPrice,AvgStrikeCallPrem/10000)
AverageStrikePutPrice  <- c(AverageStrikePutPrice, AvgStrikePutPrem/10000)

theNs <- c(theNs,"BS Sim")
table14.1WithSim <- data.frame(theNs,AveragePriceCallPrice,AveragePricePutPrice,AverageStrikeCallPrice,AverageStrikePutPrice)
colnames(table14.1) <- c("N","AvgPriceCall","AvgPricePut","AvgStrikeCall","AvgStrikePut")
table14.1WithSim

