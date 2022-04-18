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
  # D_CF: vector of cash flows ordered (amount_1,time_1,amount_2,time_2,...)
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


#From book example McDonald 12.1 Table infinity
#Let S =$41 K=$40, sigma = 0.3, r = 8%, T = 0.25 (3 months), and delta = 0 (no dividend yield)
# Book Soln = 6.961
P <- bs_parameterize(S=41,K=40,sigma=0.3,r=0.08,T=1)
black_scholes(P)

#Call Premium McDonald Example 12.1
#BlackScholes(41, 40, 0.3, 0.08, 0.25, 0, 0)
#Should equal 3.399 from textbook
P <- bs_parameterize(S=41,K=40,sigma=0.3,r=0.08,T=.25)
black_scholes(P)

#Put Premium  McDonald Example 12.2
#BlackScholes(41, 40, 0.3, 0.08, 0.25, 0, 0)
#Should equal 1.607
P <- bs_parameterize(S=41,K=40,sigma=0.3,r=0.08,T=.25,put=T)
black_scholes(P)


###Black scholes on other assets
#For use in cases with differing prepaid forward value like discrete dividends and currencies
#This requires a lot of case work so we can't work this is with a general P list
#So we require youre prepaid formula and the extra parameters you will need
#we need F_S and F_K to be of the form function(S, P, extra_P)
#where both return the prepaid forward for cash(K) and the asset(S)
black_scholes_discrete_div <- function(P){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
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

# McDonald Example 12.3
P <- bs_parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=0.25,D_CF=c(3,1/12))
black_scholes_discrete_div(P)

# THE GREEKS
gimmeGreeks <- function(P){
  S     = P$S
  K     = P$K
  sigma = P$sigma
  r     = P$r
  T_exp = P$T_exp
  delta = P$delta
  put   = P$put
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if(!put){
    greekDelta = exp(-delta*T_exp) * qnorm(d1)
    greekGamma = exp(-delta*T_exp) * dnorm(d1) / (S*sigma*sqrt(T_exp))
    greekTheta = delta*S*exp(-delta*T_exp)*qnorm(d1) - r*K*exp(-r*T_exp)*qnorm(d2) - K*exp(-r*T_exp)*dnorm(d2)*sigma/(2*sqrt(T_exp)) 
    greekVega  = S*exp(-delta*T_exp)*dnorm(d1)*sqrt(T_exp)
    greekRho   = T_exp*K*exp(-r*T_exp)*qnorm(d2)
    greekPsi   = -T_exp*S*exp(-delta*T_exp)*qnorm(d1)
  } else {
    greekDelta = -exp(-delta*T_exp) * qnorm(-d1)
    greekGamma = exp(-delta*T_exp) * dnorm(d1) / (S*sigma*sqrt(T_exp))
    greekTheta = delta*S*exp(-delta*T_exp)*qnorm(d1) - r*K*exp(-r*T_exp)*qnorm(d2) - K*exp(-r*T_exp)*dnorm(d2)*sigma/(2*sqrt(T_exp)) + r*K*exp(-r*T_exp) - delta*S*exp(-delta*T_exp)
    greekVega  = S*exp(-delta*T_exp)*dnorm(d1)*sqrt(T_exp)
    greekRho   = -T_exp*K*exp(-r*T_exp)*qnorm(-d2)
    greekPsi   = T_exp*S*exp(-delta*T_exp)*qnorm(-d1)
  }
  theGreeks = list("delta" = greekDelta,"gamma"=greekGamma,"theta"=greekTheta,"vega"=greekVega,"rho"=greekRho,"psi"=greekPsi)
    return(theGreeks)
}
