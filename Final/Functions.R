
####                                        ####
####                                        ####
# Binomial Model
####                                        ####
####                                        ####

one_step_claim <- function(S,Su,Sd,Cu,Cd,r,h,delta=0,D=0) {
  # The function solves for C, A, B as a one time step binomial tree
  # delta is the continuous dividend rate
  # D is the discrete dividend at time "t+h"
  # Outputs a numeric vector of C, A, B
  u <- Su/S
  d <- Sd/S
  A <- exp(-delta*h)*(Cu-Cd)/(Su-Sd)
  B <- exp(-r*h)*(u*Cd-d*Cu)/(u-d) - A*D*exp(-r*h)
  C <- A*S + B
  return(c(C,A,B))
}

generate_S_m <- function(S,n,u,d) {
  # Generates the matrix of S values for each node of the binomial tree using S = S[1,1], n, u, d
  # Output is S_m: an n+1 by n+1 numeric matrix.
  
  S_m <- matrix(0,n+1,n+1)
  S_m[1,1] <- S
  for(i in 1:n){
    S_a <- S_m[1:i,i]
    S_b <- c(S_a[1]*u,S_a*d)
    S_m[1:(i+1),i+1] <- S_b
  }
  return(S_m)
}

generate_S_v <- function(S,n,u,d,D_v){
  S_v <- vector("numeric",2^(n+1) - 1)
  S_v[1] <- S
  for(i in 1:n){
    D <- D_v[i]
    for (j in 1:2^(i-1)){
      k <- 2^(i-1) + j - 1
      S_v[2*k]     <- u[i]*(S_v[k] - D)
      S_v[2*k + 1] <- d[i]*(S_v[k] - D)
    }
  }
  return(S_v)
}

generate_ud <- function(choice, h, r, delta, sigma, K, S, n, u, d){
  if(choice == 0){
    u <- u
    d <- d
  } else if(choice==1){
    # McDonald Eqn (10.9)
    u <- exp((r-delta)*h+sigma*sqrt(h))
    d <- exp((r-delta)*h-sigma*sqrt(h))
  } else if(choice==2){
    # Cox-Ross-Rubinstein (CRR)
    u <- exp((r-sigma^2/2)*h+sigma*sqrt(h))
    d <- exp((r-sigma^2/2)*h-sigma*sqrt(h))
  } else if(choice==3){
    # Jarrow-Rudd (JR)
    u <- exp(sigma*sqrt(h))
    d <- exp(-sigma*sqrt(h))
  } else if(choice==4){
    # Modified Cox-Ross-Rubinstein (MCRR)
    K_n <- log(K/S)/n
    V_n <- sigma*sqrt(h)
    u   <- exp(K_n+V_n)
    d   <- exp(K_n-V_n)
  } else if(choice==5){
    #Trigeorgis Tree (TRG)
    a = sqrt(h*sigma^2 + (h*r)^2)
    u = exp(a)
    d = exp(-a)
  }
  return(list(u = u, d = d))
}

generate_D_v <- function(vec,n,h,r_v){
  # Decompose vector to time and value components
  output <- rep(0,n)
  len <- length(vec)
  CFs    <- matrix(0,nrow=len/2,ncol=2)
  for (i in 0:(len/2-1)){
    for (j in 1:2){
      CFs[i+1,j] <- vec[i*2+j]
    }
  }
  for (i in 1:(len/2)){
    # time of dividend
    a = CFs[i,1]
    # previous interval number (starting with 0)
    b = ceiling(CFs[i,1]/h) - 1
    # discount value to previous interval
    r = r_v[b+1]
    interest = exp(-r*(a - b*h))
    # update output vector
    end_div = interest*CFs[i,2]
    output[b+1] = output[b+1] + end_div
  }
  return(output)
}

discretize_r_t <- function(r_t,n,h){
  areas <- rep(0,n)
  for(i in 1:n){
    areas[i] = integrate(r_t,(i-1)*h,i*h)$value
  }
  return(areas/h)
}

solve_binomial_pricing <- function(S_v,C_v,A_v,B_v,r,h,n,eur,payoff,K,delta,D_v){
  # Solves the values of C, A, and B for each node of the binomial tree.
  # Requires the matrix S_m, C_m, A_m, and B_m. As well as r, h, n, delta.
  # If not a European option, requires the payoff function, K.
  # Output is a list of numeric matrices: S_m, C_m, A_m, and B_m
  
  for(i in seq(n,1,-1)){
    for(j in seq(1,2^(i-1))){
      k <- 2^(i-1) + j -1
      values   <- one_step_claim(S_v[k],S_v[2*k],S_v[2*k + 1],C_v[2*k],C_v[2*k+1],r[i],h,delta[i],D_v[i])
      C_v[k] <- values[1]
      A_v[k] <- values[2]
      B_v[k] <- values[3]
    }
    if(!eur){
      a <- 2^(i-1)
      b <- 2^(i)-1
      C_1 <- C_v[a:b]
      S_i <- S_v[a:b]
      C_2 <- payoff(S_i,K)
      C_i <- ifelse(C_1>C_2,C_1,C_2)
      C_v[a:b] <- C_i
    }
  }
  return(list(S_v,C_v,A_v,B_v))
}

solve_binomial_pricing_recombine <- function(S_m,C_m,A_m,B_m,r,h,n,eur,payoff,K,delta,D_v){
  # Solves the values of C, A, and B for each node of the binomial tree.
  # Requires the matrix S_m, C_m, A_m, and B_m. As well as r, h, n, delta.
  # If not a European option, requires the payoff function, K.
  # Output is a list of numeric matrices: S_m, C_m, A_m, and B_m
  
  for(i in seq(n,1,-1)){
    for(j in seq(1,i)){
      values   <- one_step_claim(S_m[j,i],S_m[j,i+1],S_m[j+1,i+1],C_m[j,i+1],C_m[j+1,i+1],r,h,delta)
      C_m[j,i] <- values[1]
      A_m[j,i] <- values[2]
      B_m[j,i] <- values[3]
    }
    if(!eur){
      C_1 <- C_m[1:i,i]
      S_i <- S_m[1:i,i]
      C_2 <- payoff(S_i,K)
      C_i <- ifelse(C_1>C_2,C_1,C_2)
      C_m[1:i,i] <- C_i
    }
  }
  collab <- paste("t = ",0:n,sep="")
  rowlab <- c()
  for (i in 1:n){
    rowlab <- c(paste(rowlab[1],"u",sep=""), paste(rowlab,"d",sep=""))
  }
  rownames(S_m) <- rowlab
  colnames(S_m) <- collab
  rownames(C_m) <- rowlab
  colnames(C_m) <- collab
  rownames(A_m) <- rowlab
  colnames(A_m) <- collab
  rownames(B_m) <- rowlab
  colnames(B_m) <- collab
  return(list("S_m" = S_m,"C_m" = C_m,"A_m" = A_m,"B_m" = B_m))
}

binomial_pricing <- function(P,payoff){
  # Function that takes a list of formatted parameters and a payoff function, and implements the binomial model.
  #   Formatted parameters P derived from parametrize()
  #
  #   payoff in the form function(S,K) where S is an asset price and
  #   K is the strike price though the claim need not actually depend on K
  # Output is a list of matrices for the node values of S, C, A, B.
  recombine = P$recombine
  if(recombine){
    u <- P$u[1]
    d <- P$d[1]
    S_m <- generate_S_m(P$S,P$n,u,d)
    C_m <- matrix(0,P$n+1,P$n+1)
    D_m <- matrix(0,P$n+1,P$n+1)
    B_m <- matrix(0,P$n+1,P$n+1)
    S_T <- S_m[,P$n+1]
    C_m[,P$n+1] <- payoff(S_T,P$K)
    solution <- solve_binomial_pricing_recombine(S_m,C_m,D_m,B_m,P$r_v[1],P$h,P$n,P$eur,payoff=payoff,P$K,P$delta,P$D_v)
  } else {
    S_v <- generate_S_v(P$S,P$n,P$u,P$d,P$D_v)
    n <- P$n
    C_v <- vector("numeric",2^(n+1)-1)
    A_v <- vector("numeric",2^(n+1)-1)
    B_v <- vector("numeric",2^(n+1)-1)
    a <- 2^(n)
    b <- 2^(n+1)-1
    S_T <- S_v[a:b]
    C_v[a:b] <- payoff(S_T,P$K)
    solution <- solve_binomial_pricing(S_v,C_v,A_v,B_v,P$r_v,P$h,P$n,P$eur,payoff=payoff,P$K,P$delta,P$D_v)
  }
  
  return(solution)
}

parameterize <- function(S,r,T_exp=0,n=0,h=0,K,sigma=0.1,delta=0,choice=0,u=10^8,d=10^8,eur=T,D_CF='None',recombine = T){
  # S: Time 0 value of the asset priced in desired currency
  # r: Risk-free rate with respect to desired currency in decimal form
  #   may be given in the form of an integer, a vector, or a function
  # T_exp: Time at which claim will expire in years
  # n: number of steps desired
  # h: time between steps
  # K: strike value of claim in desired currency
  #   though your claim need not have a strike in general this is necessary for
  #   MCRR pricing
  # sigma: volatility evaluated annually in decimal form
  #   may be given in the form of an integer or a vector
  # delta: dividend rate evaluated annually in decimal form
  #   may be given in the form of an integer or a vector
  # choice: desired pricing algorithm, see generate_ud()
  # u: rate of increase per step (not necessary is choice != 0) 
  # d:rate of decrease per step (not necessary is choice != 0) 
  # eur: American or European claim (T -> European, F -> American)
  # D_CF: vector of Dividends within claim expiration period
  #   vector in form (time, dividend value, t_2, dv_2, ...)
  #   where time is in years and value is in desired currency
  # recombine: specifies whether the algorithm should be performed using a recombining
  #    tree or a non-recombining tree. Where there are no discrete dividends
  #    or variables changing over time a recombining model is ideal
  
  
  # Find T_exp, n, and h
  if(T_exp==0){T_exp=n*h}
  else if(n==0){n=T_exp/h}
  else {h=T_exp/n}
  
  # Handles r_v:
  if(recombine){
    r_v = r
  }else if (!is.function(r)){
    r_v <- rep(r,n)
  } else {
    r_v <- discretize_r_t(r,n,h)
  }
  
  #vectorize delta
  if(!recombine){
    delta <- rep(delta,n)
  }
  
  
  #vectorize sigma
  if(!recombine){
    sigma <- rep(sigma,n)
  }
  
  # Creates D_v
  if(!is.numeric(D_CF)){
    D_v <- rep(0,n)
  } else {
    D_v <- generate_D_v(D_CF,n,h,r_v)
  }
  
  # Creates u and d
  ud <- generate_ud(choice,h,r_v,delta,sigma,K,S,n,u,d)
  u <- ud$u
  d <- ud$d
  
  P <- list("S"=S,"r_v"=r_v,"T_exp"=T_exp,"n"=n,"h"=h,"K"=K,"sigma"=sigma,"delta"=delta,"D_v"=D_v,"choice"=choice,"u"=u,"d"=d,"eur"=eur,"recombine"=recombine)
  return(P)
}

call_payoff            <- function(S,K){
  output <- ifelse(S>K,S-K,0)
  return(output)
}

put_payoff             <- function(S,K){
  output <- ifelse(K>S,K-S,0)
  return(output)
}



####                                        ####
####                                        ####
# Black Scholes Functions
####                                        ####
####                                        ####

bs_parameterize <- function(S,r,T_exp,K,sigma,delta=0,eur=T,put=F,D_CF=0){
  # S: Time 0 value of the asset priced in desired currency
  # r: Risk-free rate with respect to desired currency in decimal form
  #   may be given in the form of an integer, a vector, or a function
  # T_exp: Time at which claim will expire in years
  # K: strike value of claim in desired currency
  #   though your claim need not have a strike in general this is necessary for
  #   MCRR pricing
  # sigma: volatility evaluated annually in decimal form
  #   may be given in the form of an integer or a vector
  # delta: dividend rate evaluated annually in decimal form
  #   may be given in the form of an integer or a vector
  # choice: desired pricing algorithm, see generate_ud()
  # eur: American or European claim (T -> European, F -> American)
  # D_CF: vector of cash flows ordered (time, dividend value, t_2, dv_2, ...)
  if(length(D_CF) == 0){
    D_m <- 0
  } else {
    D_m <- matrix(D_CF,ncol=2,nrow=length(D_CF)/2,byrow=T)
  }
  P <- list("S"=S,"r"=r,"T_exp"=T_exp,"K"=K,"sigma"=sigma,"delta"=delta,"eur"=eur,"put"=put,"D_CF"=D_m)
  return(P)
}

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


###Black scholes on other assets
#For use in cases with differing prepaid forward value like discrete dividends and currencies
#This requires a lot of case work so we can't work this is with a general P list
#So we require youre prepaid formula and the extra parameters you will need
#we need F_S and F_K to be of the form function(S, P, extra_P)
#where both return the prepaid forward for cash(K) and the asset(S)
prepaid_black_scholes <- function(P, F_S, F_K = K_default, extra_P = "None"){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  F0_s = F_S(S,P,extra_P)
  F0_k = F_K(K,P,extra_P)
  d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!P$put){
    price <- F0_s*pnorm(d1) - F0_k*pnorm(d2)
    return(price)
  }
  else{
    price <- -F0_s* pnorm(-d1) + F0_k*pnorm(-d2)
    return(price)
  }
}

K_default <- function(K, P, extra_P){
  r = P$r
  t = P$T_exp
  return(K*exp(-r*t))
}

dividend_F <- function(S, P, extra_P){
  # Decompose vector to time and value components
  D <- 0
  CFs = P$D_CF
  len = length(CFs[,1])
  r   = P$r
  for (i in 1:(len/2)){
    # time of dividend
    a = CFs[i,1]
    # discount value to 0
    interest = exp(-r*a)
    # update D
    end_div = interest*CFs[i,2]
    D = D + end_div
  }
  return(S - D)
}

BlackScholes <- function(S, K, sigma, r, T_exp, delta = 0, put=F){
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!put){
    price <- S * pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
    return(price)
  }
  else {
    price <- -S * pnorm(-d1) + K*exp(-r*T_exp)*pnorm(-d2)
    return(price)
  }
}

black_scholes_discrete_div <- function(P){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  F0_s = F_S(S,P,extra_P)
  F0_k = F_K(K,P,extra_P)
  if (type == 0 ){
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
  }
  delta = 0
  D_amounts = P$D_CF[,1]
  D_times = P$D_CF[,2]
  PV_D = sum(D_amounts*exp(-r*D_times))
  F0_s = S-PV_D
  F0_k = K*exp(-r*T_exp)
  d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!P$put){
    price <- F0_s*pnorm(d1) - F0_k*pnorm(d2)
  } else {
    price <- -F0_s* pnorm(-d1) + F0_k*pnorm(-d2)
  }
  return(price)
}

####                                        ####
####                                        ####
# Volatility Analysis
####                                        ####
####                                        ####


nonannual_volatility <- function(prices){
  returns_d <- daily_returns(prices)
  sigma_d <- sqrt(sum(returns_d^2)/(length(returns_d)-1))
  return(sigma_d)
}

daily_returns <- function(daily_prices){
  n <- length(daily_prices)
  a <- daily_prices[-n]
  b <- daily_prices[-1]
  daily_returns <- (b-a)/b
  return(daily_returns)
}

ann_volatility <- function(prices,numperyear=252){
  sigma_d <- nonannual_volatility(prices)
  sigma <- sqrt(numperyear)*sigma_d
  return(sigma)
}

impliedVolatility <- function(S,K,r,T_exp,delta=0,put=F,D_CF = c(0,0),claimVal,steps = 30){
  # Numerically approximates the implied volatility from the given parameters using a version of the midpoint method
  sigmas <- c(0,1)
  for(i in 1:steps){
    sigma_mid <- mean(sigmas)
    err_mid <- BlackScholes(S,K,sigma_mid,r,T_exp,delta,put)-claimVal
    if(err_mid<0){
      sigmas <- c(sigma_mid,sigmas[2])
    } else {
      sigmas <- c(sigmas[1],sigma_mid)
    }
    
  }
  return(mean(sigmas))
}


####                                        ####
####                                        ####
# Exotic Simulation
####                                        ####
####                                        ####

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
        if(P$Put){
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
        if(P$Put){
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