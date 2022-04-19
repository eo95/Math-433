#Black Scholes with or without dividend yield
#S = Stock Price
#K = Strike Price
#r = risk free rate
#sigma = volatility
#T = expiration
#type ~ Put or Call option
#delta = dividend yield
#type = 0 ~ Pricing a Call Option
#type = 1 ~ Pricing a Put Option
BlackScholes <- function(S, K, sigma, r, T_exp, delta, type){
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (type == 0 ){
    price <- S*exp(-delta*T_exp) * pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
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
