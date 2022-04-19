# Volatility Analysis

# Adobe

adobe.table <- read.csv(file.choose()) # Import the ADOBE DATA FILE
head(adobe.table)
adobe.prices <- adobe.table$Close

adobe.sigma <- ann_volatility(adobe.prices)

timeLengths <- seq(252/2,756,252/2)

adobe.sigmas <- c()

for(i in seq(30,756,1)){
  sigma <- ann_volatility(tail(adobe.prices,i))
  adobe.sigmas <- c(adobe.sigmas,sigma)
}

adobe.P <- list(S=452.13,K=450,r=0.0033,T_exp=54/252)
adobe.impVol     <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=54/252,claimVal=mean(c(36.15,37.65)))

plot(x=seq(30/252,756/252,1/252),y=adobe.sigmas,type = 'l',
     ylim = c(0.25,.55),
     main = 'Adobe Historical and Implied Volatilities',
     xlab = 'Past Time Period in Years',
     ylab = 'Volatility',
     )
lines(x=c(0,10),y=c(adobe.impVol,adobe.impVol),type='l',col='green',lty=2)
legend(1.9,.55,legend = c("Historical Volatility","Implied Volatility"),col=c("black","green"),lty=1:2)


# Comparing Implied Vol of Call Options at Same Expiry different Strikes

adobe.options.strikes <- c(440,445,450,455,460)
adobe.options.bids    <- c(41.55,38.75,36.15,33.55,31.35)
adobe.options.asks    <- c(43.00,40.35,37.65,35.15,32.65)
adobe.options.prices  <- (adobe.options.bids+adobe.options.asks)/2
adobe.options.impVol  <- c()
adobe.options.BSprices <- c()
for(i in 1:5){
  sigma <- impliedVolatility(S=452.13,K=adobe.options.strikes[i],r=0.0033,T_exp=54/252,claimVal=adobe.options.prices[i])
  adobe.options.impVol <- c(adobe.options.impVol,sigma)
  BSprice <- BlackScholes(S=452.13,K=adobe.options.strikes[i],r=0.0033,T_exp=54/252,sigma = sigma)
  adobe.options.BSprices <- c(adobe.options.BSprices,BSprice)
}



tableOfAdobeCallsDifferentK <- cbind(adobe.options.strikes,adobe.options.prices,adobe.options.impVol,adobe.options.BSprices)
tableOfAdobeCallsDifferentK

# Comparing Implied Vol of Put Options at Same Expiry different Strikes

adobe.options.strikes <- c(440,445,450,455,460)
adobe.options.bids    <- c(30.60,32.70,35.10,37.50,40.00)
adobe.options.asks    <- c(31.75,33.95,36.25,38.70,41.30)
adobe.options.prices  <- (adobe.options.bids+adobe.options.asks)/2
adobe.options.impVol  <- c()
for(i in 1:5){
  sigma <- impliedVolatility(S=452.13,K=adobe.options.strikes[i],r=0.0033,T_exp=54/252,claimVal=adobe.options.prices[i],put=T)
  adobe.options.impVol <- c(adobe.options.impVol,sigma)
  
}

tableOfAdobePutsDifferentK <- cbind(adobe.options.strikes,adobe.options.prices,adobe.options.impVol)

tableOfAdobePutsDifferentK


# Comparing Implied Vol of Call Options at Different Maturities

adobe.options.T_exp   <- c(10,29,54,73,91)/252
adobe.options.bids    <- c(16.40,28.45,36.15,42.00,45.80)
adobe.options.asks    <- c(17.35,29.45,37.65,43.15,47.05)
adobe.options.prices  <- (adobe.options.bids+adobe.options.asks)/2
adobe.options.impVol  <- c()
for(i in 1:5){
  sigma <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=adobe.options.T_exp[i],claimVal=adobe.options.prices[i])
  adobe.options.impVol <- c(adobe.options.impVol,sigma)
  
}

tableOfAdobePutsDifferentK <- cbind(adobe.options.T_exp,adobe.options.prices,adobe.options.impVol)

tableOfAdobePutsDifferentK

# Comparing Implied Vol of Put Options at Different Maturities

adobe.options.T_exp   <- c(10,29,54,73,91)/252
adobe.options.bids    <- c(16.05,28.05,35.10,40.30,43.65)
adobe.options.asks    <- c(16.75,28.70,36.25,41.50,45.10)
adobe.options.prices  <- (adobe.options.bids+adobe.options.asks)/2
adobe.options.impVol  <- c()
for(i in 1:5){
  sigma <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=adobe.options.T_exp[i],claimVal=adobe.options.prices[i],put=T)
  adobe.options.impVol <- c(adobe.options.impVol,sigma)
  
}

tableOfAdobePutsDifferentK <- cbind(adobe.options.T_exp,adobe.options.prices,adobe.options.impVol)

tableOfAdobePutsDifferentK


# TESCO Making a Smile

tesco.options.strikes <- c(265,270,275,280,285)
tesco.options.bids    <- c(11.75,8.00,5.00,3.00,1.50)
tesco.options.asks    <- c(16.75,12.50,9.50,6.75,4.75)
tesco.options.prices  <- (tesco.options.bids+tesco.options.asks)/2
tesco.options.impVol  <- c()
tesco.options.BSprices <- c()
for(i in 1:5){
  sigma <- impliedVolatility(S=277.5,K=tesco.options.strikes[i],r=0.004452,T_exp=10/252,claimVal=tesco.options.prices[i])
  tesco.options.impVol <- c(tesco.options.impVol,sigma)
  BSprice <- BlackScholes(S=277.5,K=tesco.options.strikes[i],r=0.004452,T_exp=10/252,sigma = sigma)
  tesco.options.BSprices <- c(tesco.options.BSprices,BSprice)
}



tableOftescoCallsDifferentK <- cbind(tesco.options.strikes,tesco.options.prices,tesco.options.impVol,tesco.options.BSprices)
tableOftescoCallsDifferentK

plot(tesco.options.strikes,tesco.options.impVol, type = 'o',
     main ="Tesco Implied Volatility vs Strike",
     xlab = "Strike Price",
     ylab = "Implied Volatility",
     sub  ="\nParameters: S=277.5,r=.004452,T_exp=10/252")

# Beautiful Daily Price and Returns Plots for NYT Data

NYT.table  <- read.csv(file.choose()) # OPEN NYT DATA FILE
NYT.prices <- NYT.table$Close
plot(ts(NYT.prices),
     main="Daily Prices of NYT",
     xlab="Time (in days)",
     ylab="Price of NYT ($)")

NYT.returns <- daily_returns(NYT.prices)
hist(NYT.returns,
     main = "Histogram of NYT Daily Returns",
     xlab = "NYT Daily Returns")
ystd = c(1,1)
NYT.sigma <- ann_volatility(NYT.prices,numperyear=1)
plot((NYT.returns),type='l',
     main = "NYT Daily Returns vs Time",
     xlab = "Time (in days)",
     ylab = "Daily Returns")
abline(0,0)
lines(x=c(-10^6,10^6),y=ystd*NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*2*NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*-NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*-2*NYT.sigma,col= 'red',lty=2)
