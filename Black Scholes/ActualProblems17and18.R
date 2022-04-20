BlackScholes <- function(S, K, sigma, r, T_exp, delta = 0, put=F){
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!put){
    price <- S * exp(-delta*T_exp) *pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
    return(price)
  }
  else {
    price <- -S * exp(-delta*T_exp) * pnorm(-d1) + K*exp(-r*T_exp)*pnorm(-d2)
    return(price)
  }
}

impliedVolatility <- function(S,K,r,T_exp,delta=0,put=F,claimVal,steps = 30){
  # Numerically approximates the implied volatility from the given parameters using a version of the midpoint method
  sigmas <- c(0,1)
  print(sigmas)
  print(mean(sigmas))
  for(i in 1:steps){
    sigma_mid <- mean(sigmas)
    print(sigma_mid)
    err_mid <- BlackScholes(S,K,sigma_mid,r,T_exp,delta,put)-claimVal
    if(err_mid<0){
      sigmas <- c(sigma_mid,sigmas[2])
    } else {
      sigmas <- c(sigmas[1],sigma_mid)
    }
    
  }
  return(mean(sigmas))
}

# Problem 17 ----

strikes = c(90,95,99,100,105,110)
pricesTrues  = c(3.433,5.647,5.476,5.829,9.661,14.073)
impVols = c()
modelPrices = c()

for(i in 1:length(strikes)){
  K = strikes[i]
  priceTrue = pricesTrues[i]
  impVol = impliedVolatility(S=100, K,       r=0.05,T_exp=200/252,put=T,claimVal=priceTrue)
  modelPrice <- BlackScholes(S=100, K,impVol,r=0.05,T_exp=200/252, delta = 0,put=T)
  impVols <- c(impVols,impVol)
  modelPrices <- c(modelPrices,modelPrice)
}

table = cbind(strikes,pricesTrues,impVols,modelPrices)

table

plot(strikes,impVols,
     main = "Implied Volatilities v Strikes of Eur. Put Options",
     xlab = "Strike Prices",
     ylab = "Implied Volatility",
     type = "o")

# Problem 18 ----

truePrices  = matrix(c(7.0,3.7,1.6,8.3,5.2,2.9,10.5,7.5,5.1),3,3)
impVols     = matrix(0,3,3)
modelPrices = matrix(0,3,3)
OutputList   = list("TruePrices" = truePrices,"ImpliedVolatilities"=impVols,"ModelPrices"=modelPrices)
strikes     = c(45,50,55)
maturities  = c(.25,.5,1)

for(i in 1:3){
  rownames(OutputList[[i]]) <- strikes
  colnames(OutputList[[i]]) <- maturities
}

for(i in 1:3){
  for(j in 1:3){
    S = 50
    K     = strikes[i]
    T_exp = maturities[j]
    r = 0.05
    priceTrue = OutputList$TruePrices[i,j]
    impVol    = impliedVolatility(S,K,r,T_exp,claimVal=priceTrue)
    modelPrice <- BlackScholes(S, K,impVol,r,T_exp)
    OutputList$ImpliedVolatilities[i,j] <- impVol
    OutputList$ModelPrices[i,j] <- modelPrice
  }
}
OutputList
xData  <- strikes
yData1 <- OutputList$ImpliedVolatilities[,1]
yData2 <- OutputList$ImpliedVolatilities[,2]
yData3 <- OutputList$ImpliedVolatilities[,3]
plot(xData,yData1,xlim = c(45,55),
     ylim = c(0.3,0.4),
     main = "Implied Volatilities v Strikes of Eur. Put Options",
     xlab = "Strike Prices",
     ylab = "Implied Volatility",
     type = "l")
lines(xData,yData2, col='green',lty=2)
lines(xData,yData3, col='red',lty=3)
legend(51,.40,cex = 0.8,legend = c("T = 3 Mo.","T=6 Mo.","T = 12 Mo."),col=c("black","green","red"),lty=1:3)

# Plots for Monotonic Increasing Nature of Option Prices v Volatility

sigmas <- seq(0,1,10^-5)
calls <- BlackScholes(100,110,sigmas,.05,3/12)

plot(sigmas,puts,type='l',
     main = "Call Value v Volatility",
     sub  = "Parameters: S=100 K=110 r=.05 T=3/12",
     xlab = "X: Volatility",
     ylab = "Y: Claim Value of a Call Option")

puts <- BlackScholes(100,90,sigmas,.05,3/12,put=T)

plot(sigmas,puts,type='l',
     main = "Put Value v Volatility",
     sub  = "Parameters: S=100 K=90 r=.05 T=3/12",
     xlab = "X: Volatility",
     ylab = "Y: Claim Value of a Put Option")

# Plot to show HW 3 Q 3 was correct solution
S <- 100
K <- seq(70,130,.001)
sigma <- .3
r <- .05
T_exp <- 1
delta = 0.01

callPrices <- BlackScholes(S,K,sigma,r,T_exp,delta)
putPrices  <- BlackScholes(S,K,sigma,r,T_exp,delta, put = T)
straddlePrices <- callPrices + putPrices
K_min <- S*exp(T*(r-delta-sigma^2/2))

plot(K,straddlePrices,type="l",
     xlim = c(90,110),
     ylim = c(23,25),
     main = "Price of Straddle over Strike",
     sub  = "Parameters: S=100,sigma=0.3,r=0.05,T_exp=1,delta=0.01",
     xlab = "Strike Price",
     ylab = "Value of the Straddle")
xpoints <- c(K_min,K_min)
ypoints <- c(0,100)
lines(xpoints,ypoints,col="blue",lty = 2,type = 'l')
legend(100,25,legend = c("Calculated Strike"),lty = 2)

