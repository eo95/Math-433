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
    
    
    ### Removed because we can't use it
    
    #Equal Probabilities Tree (EQP)
    #a = r*h
    #b = 0.5*sqrt(4*h*sigma^2 - 3*(h*r)^2)
    #u = exp(0.5*a + b)
    #d = exp(b - 1.5*a)
  #} else if(choice==6){
    
    
    
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
    r = r_v[b]
    interest = exp(-r*(a - b*h))
    # update output vector
    end_div = interest*CFs[i,2]
    output[b+1] = output[b] + end_div
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
  return(list(S_m,C_m,A_m,B_m))
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
  if(!is.function(r)){
    r_v <- rep(r,n)
  } else {
    r_v <- discretize_r_t(r,n,h)
  }
  if(recombine){
    r_v = r
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
  if(D_CF == 'None'){
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

# Generic Payoffs

call_payoff            <- function(S,K){
  output <- ifelse(S>K,S-K,0)
  return(output)
}

put_payoff             <- function(S,K){
  output <- ifelse(K>S,K-S,0)
  return(output)
}


# Exotic Options

# Exotic Option {S_t} Generators

binomialPriceRunSim <- function(P,step_v){
  # Simulates a run of S values at desired steps. FUNCTION IS IN STEPS not in TIME.
  # P is output of parameterize with the defined parameters for the binomial model
  # step_v is a vector of desired step values of S on a single price run
  ##ex: step_v <- c(0,25,30,70) would result in the prices of S at step 0, 25, 30, 70
  ## step_v should always include initial step 0
  ## if desire an entire run at every step, step_v should be 0:n
  steps <- length(step_v)-1
  diff <- diff(step_v)
  r <- P$r_v[1]
  u <- P$u[1]
  d <- P$d[1]
  delta <- P$delta[1]
  q <- (exp((r-delta)*P$h)-d)/(u-d)
  S_v <- c(P$S)
  for(i in 1:steps){
    numUps <- rbinom(1,diff[i],q)
    S_a <- tail(S_v,1)
    scale <- ((u)^numUps)*((d)^(diff[i]-numUps))
    S_b <- S_a*scale
    S_v <- c(S_v,S_b)
  }
  return(S_v)
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
  
  
# Barrier Option
barrier_pricing <- function(S_v, n, payoff, K, barrier, out=T){
  l     <- 2^n - 1
  up    <- S_v[1] < barrier
  cross <- F
  C_v   <- vector("numeric",2^(n+1)-1)
  for (i in 0:l){
    path <- i
    k <- 1
    j <- 1
    while(k <= n){
      #path chooses which direction we go
      bin   <- path%%2
      path  <- (path - bin)/2
      j     <- 2*j + 1
      #are we in a different position relative to the barrier
      pos   <- S_v[j] < barrier
      pos   <- pos != up
      #shift cross to true if we are, keep it as true from then on
      cross <- cross || pos
      #iterate k so we dont go beyond vector limits
      k     <- k + 1
    }
    #we xor because we only want to use payoff if it crosses and its not an out
    #or it doesnt cross and it is an out
    if(xor(cross, out)){
      C_v[j] = payoff(S_v[j], K)
    } else{
      C_v[j] = 0
    }
  }
}
  
  #Asian Option
  #we require that n = c(per-1) where per is the number of compoundings of the averaging
  #per period (assuming evenly distributed) and c is a random constant
asian_pricing <- function(S_v, n, payoff, K=0, per, geo=F, underlying=T){
  l     <- 2^n - 1
  C_v   <- vector("numeric",2^(n+1)-1)
  c     <- n/(per-1)
  for (i in 0:l){
    path <- i
    k    <- 0
    j    <- 1
    sum  <- S_v[j]
    while(k <= n){
      #path chooses which direction we go
      bin   <- path%%2
      path  <- (path - bin)/2
      j     <- 2*j + 1
      k     <- k + 1
      #are we in an averaging period
      if (k%%c = per){
        sum <- sum + S_v[j]
      }
    }
    avg <- sum/per
    #if the avg is the underlying use it as that, otherwise use as k
    if(underlying){
      C_v[j] <- payoff(avg, K)
    } else{
      C_v[j] <- payoff(S_v[j],avg)
    }
  }
}
