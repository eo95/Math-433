#Black Scholes with or without dividend yeild
BlackScholes <- function(S, K, sigma, r, T, delta, type){
  if (type == 0 ){
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- S * pnorm(d1) - K*exp(-r*T)*pnorm(d2)
    return(price)
  }
  if (type== 1 ){
    d1 <- (log(S/K) + (r - delta + sigma^2/2)*T) / (sigma*sqrt(T))
    d2 <- d1 - sigma*sqrt(T)
    price <- -S * pnorm(-d1) + K*exp(-r*T)*pnorm(-d2)
    return(price)
  }
}

#0 is call
#1 is put
#This formula does not account for dividends

#From book example
#Let S =$41 K=$40, sigma = 0.3, r = 8%, T = 0.25 (3 months), and delta = 0 (no dividend yield)

#Call Premium
#BlackScholes(41, 40, 0.3, 0.08, 0.25, 0, 0)
#Should equal 3.399 from textbook

#Put Premium
#BlackScholes(41, 40, 0.3, 0.08, 0.25, 0, 0)
#Should equal 1.607
