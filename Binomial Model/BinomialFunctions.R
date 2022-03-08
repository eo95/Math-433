# Functions Used in the Binomial Pricing Model

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

generate_ud <- function(choice, h, r, delta, sigma, mu, K, S, n, u, d){
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
    #Equal Probabilities Tree (EQP)
    a = mu*h
    b = 0.5*sqrt(4*h*sigma^2 - 3*(h*mu)^2)
    u = exp(0.5*a + b)
    d = exp(b - 1.5*a)
  } else if(choice==6){
    #Trigeorgis Tree (TRG)
    a = sqrt(h*sigma^2 + (h*mu)^2)
    u = exp(a)
    d = exp(-a)
  }
  return(list(u = u, d = d))
}

generate_D_v <- function(vec,n,h,r){
  output <- rep(0,n)
  len <- length(vec)
  CFs    <- matrix(0,nrow=len/2,ncol=2)
  for (i in 0:(len/2-1)){
    for (j in 1:2){
      CFs[i+1,j] <- vec[i*2+j]
    }
  }
  for (i in 1:(len/2)){
    a = CFs[i,1]
    b = ceiling(CFs[i,1]/h)
    if(is.numeric(r)){
      interest = exp(r*(b*h-a))
    } else{
      interest = exp(integrate(r,a,b*h)$value)
    }
    end_div = interest*CFs[i,2]
    output[b] = output[b] + end_div
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

binomial_pricing <- function(P,payoff){
  # Function that takes a list of formatted parameters and a payoff function, and implements the binomial model.
  # Output is a list of matrices for the node values of S, C, A, B.
  ## The compiler seems to think the dollar symbol works like quotations. Copy/Paste into RStudio has no issues with this function.
  S_v <- generate_S_v(P$S,P$n,P$u,P$d, P$D_v)
  C_v <- vector("numeric",2^(n+1)-1)
  A_v <- vector("numeric",2^(n+1)-1)
  B_v <- vector("numeric",2^(n+1)-1)
  a <- 2^(n)
  b <- 2^(n+1)-1
  S_T <- S_v[a:b]
  C_v[a:b] <- payoff(S_T,P$K)
  solution <- solve_binomial_pricing(S_v,C_v,A_v,B_v,P$r_v,P$h,P$n,P$eur,payoff=payoff,P$K,P$delta,P$D_v)
  return(solution)
}

parameterize <- function(S,r,T_exp=0,n=0,h=0,K,sigma=0.1,delta=0,mu=0,choice=0,u=10^8,d=10^8,eur=T,D_CF='None'){
  # Find T_exp, n, and h
  if(T_exp==0){T_exp=n*h}
  else if(n==0){n=T_exp/h}
  else {h=T_exp/n}
  
  # Creates r_v
  if(is.numeric(r)){
    r_v <- rep(r,n)
  } else {
    r_v <- discretize_r_t(r,n,h)
  }
  
  #vectorize delta
  delta <- rep(delta,n)
  
  # Creates D_v
  if(D_CF[1] == 'None'){
    D_v <- rep(0,n)
  } else {
    D_v <- generate_D_v(D_CF,n,h,r)
  }
  
  # Creates u and d
  ud <- generate_ud(choice,h,r_v,delta,sigma,mu,K,S,n,u,d)
  u <- ud$u
  d <- ud$d
  
  P <- list("S"=S,"r_v"=r_v,"T_exp"=T_exp,"n"=n,"h"=h,"K"=K,"sigma"=sigma,"delta"=delta,"D_v"=D_v,"choice"=choice,"u"=u,"d"=d,"eur"=eur)
  return(P)
}

# VOLATILITY FUNCTIONS
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
