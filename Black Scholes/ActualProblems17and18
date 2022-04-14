BlackScholes <- function(S, K, sigma, r, T_exp, delta = 0, put=F){
  d1 <- (log(S/K) + (r - delta + sigma^2/2)*T_exp) / (sigma*sqrt(T_exp))
  d2 <- d1 - sigma*sqrt(T_exp)
  if (!put){
    price <- S * pnorm(d1) - K*exp(-r*T_exp)*pnorm(d2)
    return(price)
  }
  else {
    price <- -S * pnorm(-d1) + K*exp(-r*T_exp)*pnorm(-d2)
    return(price)
  }
}

impliedVolatility <- function(S,K,r,T_exp,delta=0,put=F,claimVal,steps = 30){
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
legend(51,.40,legend = c("3 Months","6 Months","12 Months"),col=c("black","green","red"),lty=1:3)
