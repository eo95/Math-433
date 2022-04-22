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
  # D_CF: vector of cash flows ordered (time, dividend value, t_2, dv_2, ...)
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
prepaid_black_scholes <- function(P, F_S, F_K = K_default, extra_P = "None"){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  F0_s = F_S(S,P,extra_P)
  F0_k = F_K(K,P,extra_P)
  if (!P$put){
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- F0_s*pnorm(d1) - F0_k*pnorm(d2)
    return(price)
  }
  else{
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
    price <- -F0_s* pnorm(-d1) + F0_k*pnorm(-d2)
    return(price)
  }
}

K_default <- function(K, P, extra_P){
  r = P$r
  t = P$T_exp
  return(K*exp(-r*t))
}

dividend_F <- function(S, P, extra_P){
  # Decompose vector to time and value components
  D <- 0
  CFs = P$D_CF
  len = length(CFs[,1])
  r   = P$r
  for (i in 1:(len/2)){
    # time of dividend
    a = CFs[i,1]
    # discount value to 0
    interest = exp(-r*a)
    # update D
    end_div = interest*CFs[i,2]
    D = D + end_div
  }
  return(S - D)
}


black_scholes_discrete_div <- function(P){
  S = P$S
  K = P$K
  sigma = P$sigma[1]
  r = P$r[1]
  T_exp = P$T_exp
  delta = P$delta[1]
  F0_s = F_S(S,P,extra_P)
  F0_k = F_K(K,P,extra_P)
  if (type == 0 ){
    d1 <- (log(F0_s/F0_k) + (sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
    d2 <- d1 - sigma*sqrt(T_exp)
  }
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

# McDonald Example 12.4 - Currency

# For currency, the normal black_scholes() function works. with:
## In Parameters: S_0 = x_0 and delta = r_f (the foreign risk free rate)

P.Call <- list(S=1.25,sigma=.1,K=1.2,r=.01,T_exp=1,delta=0.03,put=F)
P.Put <-list(S=1.25,sigma=.1,K=1.2,r=.01,T_exp=1,delta=0.03,put=T)

black_scholes(P.Call)
black_scholes(P.Put)

###Black Scholes on other options
#This is in cases where we are using exotic options rather than a put or call
#Because these can, for the most part, not be found explicitly we need to 
    #simulate the black scholes model to find the needed values for pricing
blackscholesPriceRunSim <- function(P,time_v){
  # P is output of parameterize with the defined parameters for the binomial model
  # time_v is a vector of desired time values of S
  ##ex: time_v <- c(0,25,30,70) would result in the prices of S at time 0, 25, 30, 70
  ## time_v should always include initial time 0
  steps <- length(time_v)-1
  diff <- diff(time_v)
  r <- P$r_v[1]
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
# THE GREEKS
theGreeks <- function(P){
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
    greekDelta = exp(-delta*T_exp) * pnorm(d1)
    greekGamma = exp(-delta*T_exp) * dnorm(d1) / (S*sigma*sqrt(T_exp))
    greekTheta = (delta*S*exp(-delta*T_exp)*pnorm(d1) - r*K*exp(-r*T_exp)*pnorm(d2) - K*exp(-r*T_exp)*dnorm(d2)*sigma/(2*sqrt(T_exp)))/365
    greekVega  = (S*exp(-delta*T_exp)*dnorm(d1)*sqrt(T_exp))/100
    greekRho   = (T_exp*K*exp(-r*T_exp)*pnorm(d2))/100
    greekPsi   = (-T_exp*S*exp(-delta*T_exp)*pnorm(d1))/100
  } else {
    greekDelta = -exp(-delta*T_exp) * pnorm(-d1)
    greekGamma = exp(-delta*T_exp) * dnorm(d1) / (S*sigma*sqrt(T_exp))
    greekTheta = (delta*S*exp(-delta*T_exp)*pnorm(d1) - r*K*exp(-r*T_exp)*pnorm(d2) - K*exp(-r*T_exp)*dnorm(d2)*sigma/(2*sqrt(T_exp)) + r*K*exp(-r*T_exp) - delta*S*exp(-delta*T_exp))/365
    greekVega  = (S*exp(-delta*T_exp)*dnorm(d1)*sqrt(T_exp))/100
    greekRho   = (-T_exp*K*exp(-r*T_exp)*pnorm(-d2))/100
    greekPsi   = (T_exp*S*exp(-delta*T_exp)*pnorm(-d1))/100
  }
  theGreeks = list("delta" = greekDelta,"gamma"=greekGamma,"theta"=greekTheta,"vega"=greekVega,"rho"=greekRho,"psi"=greekPsi)
  return(theGreeks)
}

# Verification of Greeks - McDonald Table 12.2 Replication

P.40Call <- list(S=40,sigma=0.3,r=0.08,T_exp=91/365,K=40,put=F,delta=0)
P.40Call.Price <- black_scholes(P.40Call)
P.40Call.Greeks <- theGreeks(P.40Call)
P.40Call.Greeks
P.45Call <- list(S=40,sigma=0.3,r=0.08,T_exp=91/365,K=45,put=F,delta=0)
P.45Call.Price <- black_scholes(P.45Call)
P.45Call.Greeks <- theGreeks(P.45Call)
P.45Call.Greeks


names    <- c("Weights","Price","Delta","Gamma","Vega","Theta","Rho")
Results40  <- round(c(1,P.40Call.Price,P.40Call.Greeks$delta,P.40Call.Greeks$gamma,P.40Call.Greeks$vega,P.40Call.Greeks$theta,P.40Call.Greeks$rho),5)
Results45  <- round(c(-1,P.45Call.Price,P.45Call.Greeks$delta,P.45Call.Greeks$gamma,P.45Call.Greeks$vega,P.45Call.Greeks$theta,P.45Call.Greeks$rho),5)
Combined   <- Results40-Results45
Combined[1] <- '--'
Table12.2   <- cbind(Results40,Results45,Combined)
rownames(Table12.2) <- names
Table12.2


# Plots of the greeks with Novartis Parameters
NOVN.table <- read.csv(file.choose()) # CHOOSE NOVN DATA TABLE
NOVN.prices <- NOVN.table$Close
NOVN.sigma  <- ann_volatility(NOVN.prices)
P.Call.NOVN.1mo  <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=F,T_exp=1/12)
P.Call.NOVN.3mo  <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=F,T_exp=3/12)
P.Call.NOVN.1yr  <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=F,T_exp=1)
P.Put.NOVN.1mo   <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=T,T_exp=1/12)
P.Put.NOVN.3mo   <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=T,T_exp=3/12)
P.Put.NOVN.1yr   <- list(S=77.64,K=78,r=-0.0075,delta=0.0363,sigma=NOVN.sigma,put=T,T_exp=1)
StockPrices <- seq(70,90,.1)

# Delta Call Graph
Delta.Call.1mo <- c()
Delta.Call.3mo <- c()
Delta.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  DeltaNew  <- theGreeks(P.Call.NOVN.1mo)$delta
  Delta.Call.1mo <- c(Delta.Call.1mo,DeltaNew)
  P.Call.NOVN.3mo$S <- i
  DeltaNew  <- theGreeks(P.Call.NOVN.3mo)$delta
  Delta.Call.3mo <- c(Delta.Call.3mo,DeltaNew)
  P.Call.NOVN.1yr$S <- i
  DeltaNew  <- theGreeks(P.Call.NOVN.1yr)$delta
  Delta.Call.1yr <- c(Delta.Call.1yr,DeltaNew)
}

# Delta Put Graph
Delta.Put.1mo <- c()
Delta.Put.3mo <- c()
Delta.Put.1yr <- c()
for(i in StockPrices){
  P.Put.NOVN.1mo$S <- i
  DeltaNew  <- theGreeks(P.Put.NOVN.1mo)$delta
  Delta.Put.1mo <- c(Delta.Put.1mo,DeltaNew)
  P.Put.NOVN.3mo$S <- i
  DeltaNew  <- theGreeks(P.Put.NOVN.3mo)$delta
  Delta.Put.3mo <- c(Delta.Put.3mo,DeltaNew)
  P.Put.NOVN.1yr$S <- i
  DeltaNew  <- theGreeks(P.Put.NOVN.1yr)$delta
  Delta.Put.1yr <- c(Delta.Put.1yr,DeltaNew)
}
plot(StockPrices,Delta.Call.1mo,type='l',
     main = "Call and Put Deltas v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Delta",
     ylim = c(-1,1))
lines(StockPrices,Delta.Call.3mo,type='l',col='blue')
lines(StockPrices,Delta.Call.1yr,type='l',col='darkgreen')
lines(StockPrices,Delta.Put.1mo,type='l',col='darkred')
lines(StockPrices,Delta.Put.3mo,type='l',col='purple')
lines(StockPrices,Delta.Put.1yr,type='l',col='darkorange')
abline(0,0)
legend(82,0.25,legend=c("Call with 1 Month Expiry","Call: T = 3 Month","Call with 1 Year Expiry","Put with 1 Month Expiry","Put with 3 Month Expiry","Put with 1 Year Expiry"),cex=0.6,col=c('black','blue','darkgreen','darkred','purple','darkorange'),lty=c(1,1,1,1,1,1))


# Gamma Put Graph
Gamma.Call.1mo <- c()
Gamma.Call.3mo <- c()
Gamma.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  GammaNew  <- theGreeks(P.Call.NOVN.1mo)$gamma
  Gamma.Call.1mo <- c(Gamma.Call.1mo,GammaNew)
  P.Call.NOVN.3mo$S <- i
  GammaNew  <- theGreeks(P.Call.NOVN.3mo)$gamma
  Gamma.Call.3mo <- c(Gamma.Call.3mo,GammaNew)
  P.Call.NOVN.1yr$S <- i
  GammaNew  <- theGreeks(P.Call.NOVN.1yr)$gamma
  Gamma.Call.1yr <- c(Gamma.Call.1yr,GammaNew)
}
plot(StockPrices,Gamma.Call.1mo,type='l',
     main = "Call and Call Gammas v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Gamma")
lines(StockPrices,Gamma.Call.3mo,type='l',col='blue')
lines(StockPrices,Gamma.Call.1yr,type='l',col='darkgreen')
legend(83,0.11,legend=c("Call: T = 1 Month","Call: T = 3 Month","Call: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

# Vega Call Graph
Vega.Call.1mo <- c()
Vega.Call.3mo <- c()
Vega.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  VegaNew  <- theGreeks(P.Call.NOVN.1mo)$vega
  Vega.Call.1mo <- c(Vega.Call.1mo,VegaNew)
  P.Call.NOVN.3mo$S <- i
  VegaNew  <- theGreeks(P.Call.NOVN.3mo)$vega
  Vega.Call.3mo <- c(Vega.Call.3mo,VegaNew)
  P.Call.NOVN.1yr$S <- i
  VegaNew  <- theGreeks(P.Call.NOVN.1yr)$vega
  Vega.Call.1yr <- c(Vega.Call.1yr,VegaNew)
}
plot(StockPrices,Vega.Call.1mo,type='l',
     main = "Call and Put Vegas v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Vega",
     ylim = c(0,.4))
lines(StockPrices,Vega.Call.3mo,type='l',col='blue')
lines(StockPrices,Vega.Call.1yr,type='l',col='darkgreen')
legend(83,0.4,legend=c("Call: T = 1 Month","Call: T = 3 Month","Call: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

# Theta Call Graph
Theta.Call.1mo <- c()
Theta.Call.3mo <- c()
Theta.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  ThetaNew  <- theGreeks(P.Call.NOVN.1mo)$theta
  Theta.Call.1mo <- c(Theta.Call.1mo,ThetaNew)
  P.Call.NOVN.3mo$S <- i
  ThetaNew  <- theGreeks(P.Call.NOVN.3mo)$theta
  Theta.Call.3mo <- c(Theta.Call.3mo,ThetaNew)
  P.Call.NOVN.1yr$S <- i
  ThetaNew  <- theGreeks(P.Call.NOVN.1yr)$theta
  Theta.Call.1yr <- c(Theta.Call.1yr,ThetaNew)
}
plot(StockPrices,Theta.Call.1mo,type='l',
     main = "Call Thetas v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Call Theta")
lines(StockPrices,Theta.Call.3mo,type='l',col='blue')
lines(StockPrices,Theta.Call.1yr,type='l',col='darkgreen')
abline(0,0)
legend(83,-.010,legend=c("Call: T = 1 Month","Call: T = 3 Month","Call: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

# Theta Put Graph
Theta.Put.1mo <- c()
Theta.Put.3mo <- c()
Theta.Put.1yr <- c()
for(i in StockPrices){
  P.Put.NOVN.1mo$S <- i
  ThetaNew  <- theGreeks(P.Put.NOVN.1mo)$theta
  Theta.Put.1mo <- c(Theta.Put.1mo,ThetaNew)
  P.Put.NOVN.3mo$S <- i
  ThetaNew  <- theGreeks(P.Put.NOVN.3mo)$theta
  Theta.Put.3mo <- c(Theta.Put.3mo,ThetaNew)
  P.Put.NOVN.1yr$S <- i
  ThetaNew  <- theGreeks(P.Put.NOVN.1yr)$theta
  Theta.Put.1yr <- c(Theta.Put.1yr,ThetaNew)
}
plot(StockPrices,Theta.Put.1mo,type='l',
     main = "Put Thetas v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Put Theta")
lines(StockPrices,Theta.Put.3mo,type='l',col='blue')
lines(StockPrices,Theta.Put.1yr,type='l',col='darkgreen')
abline(0,0)
legend(83,-.010,legend=c("Put: T = 1 Month","Put: T = 3 Month","Put: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

# Rho Call Graph
Rho.Call.1mo <- c()
Rho.Call.3mo <- c()
Rho.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  RhoNew  <- theGreeks(P.Call.NOVN.1mo)$rho
  Rho.Call.1mo <- c(Rho.Call.1mo,RhoNew)
  P.Call.NOVN.3mo$S <- i
  RhoNew  <- theGreeks(P.Call.NOVN.3mo)$rho
  Rho.Call.3mo <- c(Rho.Call.3mo,RhoNew)
  P.Call.NOVN.1yr$S <- i
  RhoNew  <- theGreeks(P.Call.NOVN.1yr)$rho
  Rho.Call.1yr <- c(Rho.Call.1yr,RhoNew)
}
plot(StockPrices,Rho.Call.1mo,type='l',
     main = "Call Rhos v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Call Rho",
     ylim = c(0,0.6))
lines(StockPrices,Rho.Call.3mo,type='l',col='blue')
lines(StockPrices,Rho.Call.1yr,type='l',col='darkgreen')
legend(70,0.5,legend=c("Call: T = 1 Month","Call: T = 3 Month","Call: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

# Psi Call Graph
Psi.Call.1mo <- c()
Psi.Call.3mo <- c()
Psi.Call.1yr <- c()
for(i in StockPrices){
  P.Call.NOVN.1mo$S <- i
  PsiNew  <- theGreeks(P.Call.NOVN.1mo)$psi
  Psi.Call.1mo <- c(Psi.Call.1mo,PsiNew)
  P.Call.NOVN.3mo$S <- i
  PsiNew  <- theGreeks(P.Call.NOVN.3mo)$psi
  Psi.Call.3mo <- c(Psi.Call.3mo,PsiNew)
  P.Call.NOVN.1yr$S <- i
  PsiNew  <- theGreeks(P.Call.NOVN.1yr)$psi
  Psi.Call.1yr <- c(Psi.Call.1yr,PsiNew)
}
plot(StockPrices,Psi.Call.1mo,type='l',
     main = "Call Psis v NOVN Prices",
     xlab = "NOVN (CHF)",
     ylab = "Call Psi",
     ylim = c(-0.8,0))
lines(StockPrices,Psi.Call.3mo,type='l',col='blue')
lines(StockPrices,Psi.Call.1yr,type='l',col='darkgreen')
legend(70,-0.55,legend=c("Call: T = 1 Month","Call: T = 3 Month","Call: T = 1 Year"),cex=0.6,col=c('black','blue','darkgreen'),lty=c(1,1,1))

