#Black Scholes with or without dividend yeild
BlackScholes <- function(S, K, sigma, r, T, delta, type){
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

#0 is call
#1 is put
#This formula does not account for dividends
