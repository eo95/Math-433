#Black Scholes for Prepaid Forward
BlackScholesF <- function(S, K, sigma, r, T, type, delta){
  if (type == 0 ){
    S = So - dividends*exp(-r*t)
    d1 <- (log(S*exp(-delta*T)/K*exp(-r*T)) + (sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- S * pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
    S = So - dividends*exp(-r*t)
    d1 <- (log(S*exp(-delta*T)/K*exp(-r*T)) + (r + sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- -S * pnorm(-d1) + K*exp(-r*T)*pnorm(-d2)
    return(price)
  }
}