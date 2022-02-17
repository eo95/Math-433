#Black Scholes with Discrete Dividends
BlackScholesD <- function(So, K, sigma, r, T, type, delta, dividends,t){
  if (type == 0 ){
    S = So - dividends*exp(-r*t)
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- S * pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
    S = So - dividends*exp(-r*t)
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- -S * pnorm(-d1) + K*exp(-r*T)*pnorm(-d2)
    return(price)
  }
}
