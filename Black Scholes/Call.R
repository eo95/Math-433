Call <- function(S, K, sigma, r, T_exp, delta) {
  d1 <- (log(S/K) + (r + sigma^2/2)*T) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  S *exp(-delta*T_exp)* pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
  
}


