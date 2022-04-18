bs_parameterize <- function(S,r,T_exp=0,n=0,h=0,K,sigma,delta=0,eur=T){
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

  P <- list("S"=S,"r_v"=r,"T_exp"=T_exp,"K"=K,"sigma"=sigma,"delta"=delta,"eur"=eur)
  return(P)
}

#Black Scholes with or without dividend yeild
#S = Stock Price
#K = Strike Price
#r = risk free rate
#sigma = volatility
#T = expiration
#type ~ Put or Call option
#delta = dividend yield
black_scholes <- function(P, type){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r_v[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  if (type == 0 ){
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- S*exp(-delta*T_exp) * pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- -S*exp(-delta*T_exp) * pnorm(-d1) + K*exp(-r*T_exp)*pnorm(-d2)
    return(price)
  }
}

#Comparing Binomial option prices. Comparing to figure 10.3. 
#Suppose Stock price~ S= $41, Strike price~ K=$40, volatility~sigma= 0.3
#risk free rate~ r=0.08, time to expiration~T_exp=1, dividend yield~delta = 0
BlackScholes(S = 41, K = 40, sigma = 0.3, r = 0.08, T = 1, delta = 0, type=0)
#Output = 6.960999
#Mcdonald Answer =6.961

#Example 12.1 McDonald. Let S=$41, K=$40, sigma=0.3, r=0.08, T = 0.25, delta =0
#Computing the Black-Scholes call price
BlackScholes(S = 41, K = 40, sigma = 0.3, r = 0.08, T = 0.25, delta = 0, type=0)
#Output = 3.399078
#McDonald Answer = 3.399

#Computing the Blac-Scholes put price keeping the same parameters
BlackScholes(S = 41, K = 40, sigma = 0.3, r = 0.08, T = 0.25, delta = 0, type = 1)
# Output = 1.607025
# MacDonald Answer = 1.607


###Black scholes on other assets
#For use in cases with differing prepaid forward value like discrete dividends and currencies
#This requires a lot of case work so we can't work this is with a general P list
#So we require youre prepaid formula and the extra parameters you will need
    #we need F_S and F_K to be of the form function(S, P, extra_P)
    #where both return the prepaid forward for cash(K) and the asset(S)
prepaid_black_scholes <- function(P, F_S, F_K, extra_P, type){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r_v[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  F0_s = F_S(S,P,extra_P)
  F0_k = F_K(S,P,extra_P)
  if (type == 0 ){
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- F0_s*pnorm(d1) - F0_k*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- -F0_s* pnorm(-d1) + F0_k*pnorm(-d2)
    return(price)
  }
}



###Black Scholes on other options
#This is in cases where we are using exotic options rather than a put or call
#Because these can, for the most part, not be found explicitly we need to 
    #simulate the black scholes model to find the needed values for pricing
