
# Functions

one_step_claim <- function(S,Su,Sd,Cu,Cd,r,h,delta) {
  # The function solves for C, A, B as a one time step binomial tree
  # Outputs a numeric vector of C, A, B
  
  u <- Su/S
  d <- Sd/S
  A <- exp(-delta*h)*(Cu-Cd)/(Su-Sd)
  B <- exp(-r*h)*(u*Cd-d*Cu)/(u-d)
  C <- B + A*S
  return(c(C,A,B))
}

one_step_claim_discrete <- function(S,Su,Sd,Cu,Cd,r,h,D) {
  # The function solves for C, A, B as a one time step binomial tree
  # with discrete dividends
  # Outputs a numeric vector of C, A, B
  
  u <- Su/S
  d <- Sd/S
  A <- (Cu-Cd)/(Su-Sd)
  B <- exp(-r*h)*(Su*Cd-Sd*Cu)/(Su-Sd) - A*D*exp(-r*h)
  C <- B + A*S
  return(c(C,A,B))
}

generate_S_m <- function(S,n,u,d) {
  # Generates the matrix of S values for each node of the binomial tree using S = S[1,1], n, u, d
  # Output is S_m: an n+1 by n+1 numeric matrix.
  
  S_m <- matrix(0,n+1,n+1)
  S_m[1,1] <- S
  for(i in 1:n){
    S_a <- S_m[1:i,i]
    S_b <- c(S_a[1]*u,S_a * d)
    S_m[1:(i+1),i+1] <- S_b
  }
  return(S_m)
}
generate_S_m_new <- function(S,n,u,d){
  S_m <- c(S)
  for (i in 1:n){
    a <- (i-1)*i/2 + 1
    b <- (i+1)*i/2
    S_a <- S_m[a:b]
    S_b <- c(S_a[1]*u,S_a*d)
    S_m <- append(S_m, S_b)
  }
  return(S_m)
}
generate_S_m_nonconstant <- function(S,n,u,d){
  S_m <- vector("numeric",2^(n+1) - 1)
  S_m[1] <- S
  for(i in 1:n){
    for (j in 1:2^(i-1)){
      k <- 2^(i-1) + j - 1
      S_m[2*k]     <- u[i]*S_m[k]
      S_m[2*k + 1] <- d[i]*S_m[k]
    }
  }
  return(S_m)
}
generate_ud <- function(choice, h, r, delta, sigma, mu, u, d){
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
solve_binomial_pricing <- function(S_m,C_m,A_m,B_m,r,h,n,eur,payoff,K,delta){
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
solve_binomial_pricing_nonconstant <- function(S_m,C_m,A_m,B_m,r,h,n,eur,payoff,K,delta){
  # Solves the values of C, A, and B for each node of the binomial tree.
  # Requires the matrix S_m, C_m, A_m, and B_m. As well as r, h, n, delta.
  # If not a European option, requires the payoff function, K.
  # Output is a list of numeric matrices: S_m, C_m, A_m, and B_m
  
  for(i in seq(n,1,-1)){
    for(j in seq(1,2^(i-1))){
      k <- 2^(i-1) + j -1
      values   <- one_step_claim(S_m[k],S_m[2*k],S_m[2*k + 1],C_m[2*k],C_m[2*k+1],r[i],h,delta[i])
      C_m[k] <- values[1]
      A_m[k] <- values[2]
      B_m[k] <- values[3]
    }
    if(!eur){
      a <- 2^i
      b <- b
      C_1 <- C_m[a:b]
      S_i <- S_m[a:b]
      C_2 <- payoff(S_i,K)
      C_i <- ifelse(C_1>C_2,C_1,C_2)
      C_m[a:b] <- C_i
    }
  }
  return(list(S_m,C_m,A_m,B_m))
}
solve_binomial_pricing_new <- function(S_m,C_m,A_m,B_m,r,h,n,eur,payoff,K,delta){
  # Solves the values of C, A, and B for each node of the binomial tree.
  # Requires the matrix S_m, C_m, A_m, and B_m. As well as r, h, n, delta.
  # If not a European option, requires the payoff function, K.
  # Output is a list of numeric matrices: S_m, C_m, A_m, and B_m
  
  for(i in seq(n,1,-1)){
    for(j in seq(1,i)){
      k <- (i-1)*i/2 + j 
      values   <- one_step_claim(S_m[k],S_m[k+i],S_m[k+i+1],C_m[k+i],C_m[k+i+1],r[i],h,delta[i])
      C_m[k] <- values[1]
      A_m[k] <- values[2]
      B_m[k] <- values[3]
    }
    if(!eur){
      a <- (i-1)*i/2 + 1
      b <- b
      C_1 <- C_m[a:b]
      S_i <- S_m[a:b]
      C_2 <- payoff(S_i,K)
      C_i <- ifelse(C_1>C_2,C_1,C_2)
      C_m[a:b] <- C_i
    }
  }
  return(list(S_m,C_m,A_m,B_m))
}

binomial_pricing <- function(P,payoff){
  # Function that takes a list of formatted parameters and a payoff function, and implements the binomial model.
  # Output is a list of matrices for the node values of S, C, A, B.
  ## The compiler seems to think the dollar symbol works like quotations. Copy/Paste into RStudio has no issues with this function.
  h <- P$T_exp/P$n
  n <- P$n
  if (typeof(P$r) == "closure"){
    discretize_r_t(r,n,h)
  }
  ud <- generate_ud(P$choice, h, P$r, P$delta, P$sigma, P$mu, P$u, P$d)
  u <- ud$u
  d <- ud$d
  if(length(P$r) == 1){
    r <- rep(P$r, n)
  } else{
    r <- P$r
  }
  if(length(P$delta) == 1){
    delta <- rep(P$delta, n)
  } else{
    delta <- P$delta
  }
  if(length(u) == 1){
    S_m <- generate_S_m_new(P$S,n,u,d)
    C_m <- vector("numeric",(n+1)*(n+2)/2)
    A_m <- vector("numeric",(n+1)*(n+2)/2)
    B_m <- vector("numeric",(n+1)*(n+2)/2)
    a <- (n+1)*n/2+1
    b <- (n+1)*(n+2)/2
    S_T <- S_m[a:b]
    C_m[a:b] <- payoff(S_T,P$K)
    solution <- solve_binomial_pricing_new(S_m,C_m,A_m,B_m,r,h,n,P$eur,payoff=payoff,P$K,delta)
  } else{
    S_m <- generate_S_m_nonconstant(P$S,n,u,d)
    C_m <- vector("numeric",2^(n+1)-1)
    A_m <- vector("numeric",2^(n+1)-1)
    B_m <- vector("numeric",2^(n+1)-1)
    a <- 2^(n)
    b <- 2^(n+1)-1
    S_T <- S_m[a:b]
    C_m[a:b] <- payoff(S_T,P$K)
    solution <- solve_binomial_pricing_nonconstant(S_m,C_m,A_m,B_m,r,h,n,P$eur,payoff=payoff,P$K,delta)
  }
  return(solution)
}

generate_P <- function(S,K,r,T_exp,n,sigma=0.1,delta=0,mu=0,choice=0,u=10^8,d=10^8,eur=T){
  # Generate P takes the given values and choices for the problem and outputs a formatted list for binomial_pricing to use
  # S:      The time 0 price of the underlying asset.
  # K:      The strike price of the option contract.
  # r:      The risk free rate of return. a(t)=exp(r*t) with t in years.
  # T_exp:  The time until expiration of the contract from time 0 in years.
  # n:      The number of steps for the binomial model to take. For notation purposes, h = T_exp/n
  # sigma:  The volatility of the underlying asset.
  # delta:  The continuous dividend rate for an index.
  # mu:     The drift
  # choice: An integer between 0 and 4. Choice determines how u and d are calculated.
  # choice = 0: Given parameters and will not be calculated
  # choice = 1: Determined by Equation (10.9) in McDonald
  # choice = 2: Determined by the CRR equations.
  # choice = 3: Determined by the JR equations.
  # choice = 4: Determined by the MCRR equations.
  # u:     u = Su/S. Enter as a parameter if no calculation choice is made.
  # d:     d = Sd/S. Enter as a parameter if no calculation choice is made.
  # eur:   A logical argument defaulted as TRUE. If the option is American, set eur to be FALSE.
  
  P <- list("S"=S,"K"=K,"r"=r,"T_exp"=T_exp,"n"=n,"sigma"=sigma,"delta"=delta,"mu"=mu, "choice"=choice,"u"=u,"d"=d,"eur"=eur)
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

#Represent a generic payoff with a vector
#Vector form is as follows
# Entry 1   -- payoff of first inflection point
# Entry i   -- Slope before inflection point
# Entry i+1 -- Inflection asset value
# Entry n   -- final slope
#Ethan
generic_payoff        <- function(S,payoff){
  n <- length(payoff)
  output <- payoff[1]
  for (i in seq(3,n-1,by=2)){
    x <- payoff[i]
    slope <- payoff[i-1]
    if (S <= x){
      output <- output - (x - S)*slope
      return(output)
    } else{
      dx <- ifelse(i+2<n, payoff[i+2], S) - x
      slope <- payoff[i+1]
      output <- output + dx*slope
    }
  }
  return(output)
}
#Converts the more easily determined payoff form
#Beginning Slope
#Inflection points in (x,y) order
#Ending Slope
#into the form needed above
#Ethan
payoff_form         <- function(vec){
  n <- length(vec)
  if (n < 3){
    return(c(0,0,0,0))
  }
  payoff <- c(vec[3], vec[1])
  if (n < 6){
    payoff <- append(payoff, c(vec[n-2], vec[n]))
    return(payoff)
  }
  for (i in seq(2, n-4, by = 2)){
    slope <- (vec[i+3] - vec[i+1])/(vec[i+2] - vec[i])
    payoff <- append(payoff, c(vec[i], slope))
  }
  payoff <- append(payoff, c(vec[n-2], vec[n]))
  return(payoff)
}
discretize_r_t <- function(r_t,n,h){
  areas <- rep(0,n)
  for(i in 1:n){
    areas[i] = integrate(r_t,(i-1)*h,i*h)$value
  }
  return(areas/h)
}


# Testing Functions with Exercises, Examples and Problems:

# Pole Exercises: ----

# Problem 2
# (a)
P_CRR.256  <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=256)
P_CRR.512  <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=512)
P_JR.256   <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=256)
P_JR.512   <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=512)
P_MCRR.256 <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=256)
P_MCRR.512 <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=512)

C_1 <- binomial_pricing(P_CRR.256,payoff=put_payoff)[[2]][1]
C_2 <- binomial_pricing(P_CRR.512,payoff=put_payoff)[[2]][1]
C_3 <- binomial_pricing(P_JR.256,payoff=put_payoff)[[2]][1]
C_4 <- binomial_pricing(P_JR.512,payoff=put_payoff)[[2]][1]
C_5 <- binomial_pricing(P_MCRR.256,payoff=put_payoff)[[2]][1]
C_6 <- binomial_pricing(P_MCRR.512,payoff=put_payoff)[[2]][1]

sol_2a <- matrix(c(C_1,C_2,C_3,C_4,C_5,C_6),3,2,byrow=T)
rownames(sol_2a) <- c("CRR","JR","MCRR")
colnames(sol_2a) <- c("N=256","N=512")

# (b)
P_CRR.256  <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=256)
P_CRR.512  <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=512)
P_JR.256   <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=256)
P_JR.512   <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=512)
P_MCRR.256 <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=256)
P_MCRR.512 <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=512)

C_1 <- binomial_pricing(P_CRR.256,payoff=call_payoff)[[2]][1]
C_2 <- binomial_pricing(P_CRR.512,payoff=call_payoff)[[2]][1]
C_3 <- binomial_pricing(P_JR.256,payoff=call_payoff)[[2]][1]
C_4 <- binomial_pricing(P_JR.512,payoff=call_payoff)[[2]][1]
C_5 <- binomial_pricing(P_MCRR.256,payoff=call_payoff)[[2]][1]
C_6 <- binomial_pricing(P_MCRR.512,payoff=call_payoff)[[2]][1]

sol_2b <- matrix(c(C_1,C_2,C_3,C_4,C_5,C_6),3,2,byrow=T)
rownames(sol_2b) <- c("CRR","JR","MCRR")
colnames(sol_2b) <- c("N=256","N=512")
sol_2a
sol_2b

# McDonald Examples: ----

# Figure 10.4

P_fig10.4 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=2,n=2,choice=1)
binomial_pricing(P_fig10.4,payoff=call_payoff)

# Figure 10.5

P_fig10.5 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.5,payoff=call_payoff)

# Figure 10.6

P_fig10.6 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.6,payoff=put_payoff)

# Figure 10.7

P_fig10.7 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1,eur=F)
binomial_pricing(P_fig10.7,payoff=put_payoff)

# Figure 10.8

P_fig10.8 <- generate_P(S=110,K=100,sigma=0.3,r=0.05,T_exp=1,n=3,choice=1,eur=F,delta=0.035)
binomial_pricing(P_fig10.8,payoff=call_payoff)

# McDonald Problems: ----

# Problem 10.1

P_10.1 <- generate_P(S=100,K=105,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)

#(a)
binomial_pricing(P_10.1,payoff=call_payoff)

#(b)
binomial_pricing(P_10.1,payoff=put_payoff)

# Problem 10.2
P_10.2 <- generate_P(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.2,payoff=call_payoff)[[2]][1]

# Problem 10.3

P_10.3 <- generate_P(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.3,payoff=put_payoff)[[2]][1]

# Problem 10.6

P_10.6 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.6,payoff=call_payoff)

# Problem 10.8

P_10.8 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.8,payoff=put_payoff)

# Problem 10.10

P_10.10 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,eur=F,n=2)
binomial_pricing(P_10.10,payoff=call_payoff)

