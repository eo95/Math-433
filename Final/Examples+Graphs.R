###############################################
###############################################
# Normality
###############################################
###############################################

NYT.table  <- read.csv(file.choose()) # OPEN NYT DATA FILE
NYT.prices <- NYT.table$Close
plot(ts(NYT.prices),
     main="Daily Prices of NYT",
     xlab="Time (in days)",
     sub ="t=0 is Mar. 8, 2018",
     ylab="Price of NYT ($)")


NYT.returns <- daily_returns(NYT.prices)
hist(NYT.returns,
     main = "Histogram of NYT Daily Returns",
     xlab = "NYT Daily Returns")
ystd = c(1,1)
NYT.sigma <- ann_volatility(NYT.prices,numperyear=1)
plot((NYT.returns),type='l',
     main = "NYT Daily Returns vs Time",
     xlab = "Time (in days)\n",
     ylab = "Daily Returns",
     sub  = expression(paste("red lines at ±",sigma," and ±2",sigma)))
abline(0,0)
lines(x=c(-10^6,10^6),y=ystd*NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*2*NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*-NYT.sigma,col = 'red',lty=2)
lines(x=c(-10^6,10^6),y=ystd*-2*NYT.sigma,col= 'red',lty=2)


###############################################
###############################################
# Binomial
###############################################
###############################################

# Time to Compile
T_exp = 1
n = 1
S = 50
K = 0
r = 0.05
choice = 0
u = 1.5
d = 0.75
x <- c()
y_1 <- c()
for (i in 1:23){
  n <- i
  x <- append(x,n)
  P <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  avg <- 0
  for (j in 1:3){
    start_time <- Sys.time()
    binomial_pricing(P,call_payoff)
    end_time   <- Sys.time()
    avg        <- avg + as.numeric(end_time-start_time)
  }
  avg <- avg/5
  y_1   <- append(y_1, avg)
}
y_2 <- c()
for (i in 1:23){
  n <- i
  P <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n,recombine=F)
  for (j in 1:3){
    start_time <- Sys.time()
    binomial_pricing(P,call_payoff)
    end_time   <- Sys.time()
    avg        <- avg + as.numeric(end_time-start_time)
  }
  avg <- avg/5
  y_2   <- append(y_2, avg)
}
plot(x,y_2, col= "blue", type="l", xlab="Number of Steps", ylab="Time (seconds)", main="Time to Compile vs. Number of Steps")
lines(x,y_1, col="green")


# Convergence to BS (our example)
x <- 11:400
y <- c()
for (i in 11:400){
  P <- parameterize(S = 107, T_exp = 0.5, n = i, r = 0.035,K = 100, sigma = 0.12, choice = 2, delta=0.3, recombine=T)
  y <- append(y, binomial_pricing(P, call_payoff)[[2]][1])
}
plot(x,y,xlab="Steps",ylab="Price (dollars)", main="Binomial Price with Increasing Steps", ylim=c(0.7,1.1))
P <- bs_parameterize(S = 107, T_exp = 0.5, r = 0.035, K = 100, sigma = 0.12, delta=0.3, put=F)
bs <- black_scholes(P)
bs <- rep(bs, 400)
lines(x,bs)

x <- 4:300
y <- c()
for (i in x){
  P <- parameterize(S=107, T_exp=0.5, n=i, r=0.035,K=100,sigma=0.12,choice=1,D_CF=c(0.25,10,0.5,10))
  y <-append(y, binomial_pricing(P,call_payoff)[[2]][1])
}
plot(x,y,xlab="Steps",ylab="Price (dollars)",main = "Binomial Convergence to BS Formula")
P <- bs_parameterize(S=107,T_exp=0.5,r=0.035,K=100,sigma=0.12,D_CF=c(0.25,10,0.5,10),put=F)
bs <- black_scholes_div(P)
lines(c(0,500),c(bs,bs),type='l',col="Blue")


# Binomial PDF Convergence
x <- c()
bs <- rep(4.0733,64)
y_1 <- c()
y_2 <- c()
y_3 <- c()
y_4 <- c()
y_5 <- c()
for (i in 1:64){
  k = 8*i
  x      <- append(x,k)
  P_CRR  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=k)
  y_1    <- append(y_1, binomial_pricing(P_CRR,payoff=put_payoff)[[2]][1])
  P_JR   <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=k)
  y_2    <- append(y_2, binomial_pricing(P_JR,payoff=put_payoff)[[2]][1])
  P_MCRR <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=k)
  y_3    <- append(y_3, binomial_pricing(P_MCRR,payoff=put_payoff)[[2]][1])
  P_TRG  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=5,n=k)
  y_5    <- append(y_5, binomial_pricing(P_TRG,payoff=put_payoff)[[2]][1])
}
plot(x,y_1,xlab="Steps",ylab="Price (dollars)", main="Put Convergence to Black-Scholes over Time")
lines(x,y_2,type="p",col="blue")
lines(x,y_3,type="p",col="green")
lines(x,y_5,type="p",col="red")
lines(x,bs)
x <- c()
y_1 <- c()
y_2 <- c()
y_3 <- c()
y_4 <- c()
y_5 <- c()
bs  <- rep(2.1334, 64)
for (i in 1:64){
  k = 8*i
  x      <- append(x,k)
  C_CRR  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=k)
  y_1    <- append(y_1, binomial_pricing(C_CRR,payoff=call_payoff)[[2]][1])
  C_JR   <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=k)
  y_2    <- append(y_2, binomial_pricing(C_JR,payoff=call_payoff)[[2]][1])
  C_MCRR <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=k)
  y_3    <- append(y_3, binomial_pricing(C_MCRR,payoff=call_payoff)[[2]][1])
  C_TRG  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=5,n=k)
  y_5    <- append(y_5, binomial_pricing(C_TRG,payoff=call_payoff)[[2]][1])
}
plot(x,y_1,xlab="Steps",ylab="Price (dollars)", main="Call Convergence to Black-Scholes over Time")
lines(x,y_2,type="p",col="blue")
lines(x,y_3,type="p",col="green")
lines(x,y_5,type="p",col="red")
lines(x,bs)

# Figure 10.4

P_fig10.4 <- parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=2,n=2,choice=1)
binomial_pricing(P_fig10.4,payoff=call_payoff)

# Figure 10.5

P_fig10.5 <- parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.5,payoff=call_payoff)

# Figure 10.6

P_fig10.6 <- parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.6,payoff=put_payoff)

# Figure 10.7

P_fig10.7 <- parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1,eur=F)
binomial_pricing(P_fig10.7,payoff=put_payoff)

# Figure 11.9
P_fig11.9 <- parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1,D = c(2/3, 5))
binomial_pricing(P_fig11.9,payoff=put_payoff)[[1]]

# Problem 7
square_payoff   <- function(S,K){
  return(S^2)
}
square_value    <- function(S,u,d,r,n){
  output <- (S^2)*(((exp(r))*(u+d) - u*d)^n)
  output <- output*exp(-r*n)
  return(output)
}
log_payoff      <- function(S,K){
  return(log(S))
}
log_value    <- function(S,u,d,r,n){
  q <- (exp(r) - d)/(u-d)
  output <- (log(S*(d^n))+n*q*log(u/d))
  output <- output*exp(-r*n)
  return(output)
}



#Create graphs to demonstrate functional relationship
T_exp = 1
n = 1
S = 50
K = 0
r = 0.05
choice = 0
u = 1.5
d = 0.75

#Square function: n vs price
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  n = i
  T_exp = i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,square_payoff)[[2]][1])
  y_2    <- append(y_2, square_value(S,u,d,r,n))
}
plot(x,y_1,xlab="Steps", ylab="Price (dollars)", main="Square Payoff: Price vs Number of Steps")
lines(x,y_2)

#Square function: r vs price
n = 1
T_exp = 1
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  r = 0.02*i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,square_payoff)[[2]][1])
  y_2    <- append(y_2, square_value(S,u,d,r,n))
}
plot(x,y_1,xlab="Interest (decimal)", ylab="Price (dollars)", main="Square Payoff: Price vs Interest Rate")
lines(x,y_2)

#Square function: S vs price
r = 0.05
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  S = 10*i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,square_payoff)[[2]][1])
  y_2    <- append(y_2, square_value(S,u,d,r,n))
}
plot(x,y_1, xlab="Asset Value (dollars)", ylab="Price (dollars)", main="Square Payoff: Price vs Asset Value")
lines(x,y_2)


T_exp = 1
n = 1
S = 50
K = 0
r = 0.05
choice = 0
u = 1.5
d = 0.75

#Log function: n vs price
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  n = i
  T_exp = i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,log_payoff)[[2]][1])
  y_2    <- append(y_2, log_value(S,u,d,r,n))
}
plot(x,y_1,xlab="Steps", ylab="Price (dollars)", main="Log Payoff: Price vs Number of Steps")
lines(x,y_2)

#Log function: r vs price
n = 1
T_exp = 1
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  r = 0.02*i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,log_payoff)[[2]][1])
  y_2    <- append(y_2, log_value(S,u,d,r,n))
}
plot(x,y_1,xlab="Interest (decimal)", ylab="Price (dollars)", main="Log Payoff: Price vs Interest Rate")
lines(x,y_2)

#Log function: S vs price
r = 0.05
x   <- c()
y_1 <- c()
y_2 <- c()
for (i in 1:50){
  S = 10*i
  x      <- append(x, i)
  P      <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  y_1    <- append(y_1, binomial_pricing(P,log_payoff)[[2]][1])
  y_2    <- append(y_2, log_value(S,u,d,r,n))
}
plot(x,y_1, xlab="Asset Value (dollars)", ylab="Price (dollars)", main="Log Payoff: Price vs Asset Value")
lines(x,y_2)



###############################################
###############################################
# Black Scholes
###############################################
###############################################


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

# McDonald Example 12.3
P <- bs_parameterize(S=41,K=40,sigma=0.3,r=0.08,T_exp=0.25,D_CF=c(3,1/12))
black_scholes_discrete_div(P)


# McDonald Example 12.4 - Currency

# For currency, the normal black_scholes() function works. with:
## In Parameters: S_0 = x_0 and delta = r_f (the foreign risk free rate)

P.Call <- list(S=1.25,sigma=.1,K=1.2,r=.01,T_exp=1,delta=0.03,put=F)
P.Put <-list(S=1.25,sigma=.1,K=1.2,r=.01,T_exp=1,delta=0.03,put=T)


###############################################
###############################################
# Volatility
###############################################
###############################################


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

# McDonald Example 12.10

impliedVolatility(S=50,K=45,r=0.08,T_exp=0.5,claimVal=8.07)

# Historical vs Implied Volatility Graphs

  #ADBE
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
legend(1.9,.55,cex=0.6,legend = c("Historical Volatility","Implied Volatility"),col=c("black","green"),lty=1:2)
  #TSCO
adobe.table <- read.csv(file.choose()) # Import the tesco DATA FILE
head(tesco.table)
tesco.prices <- tesco.table$Close

tesco.sigma <- ann_volatility(tesco.prices)

timeLengths <- seq(252/2,756,252/2)

tesco.sigmas <- c()

for(i in seq(30,756,1)){
  sigma <- ann_volatility(tail(tesco.prices,i))
  tesco.sigmas <- c(tesco.sigmas,sigma)
}

tesco.P <- list(S=452.13,K=450,r=0.0033,T_exp=54/252)
tesco.impVol     <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=54/252,claimVal=mean(c(36.15,37.65)))

plot(x=seq(30/252,756/252,1/252),y=tesco.sigmas,type = 'l',
     ylim = c(0.25,.55),
     main = 'tesco Historical and Implied Volatilities',
     xlab = 'Past Time Period in Years',
     ylab = 'Volatility',
)
lines(x=c(0,10),y=c(tesco.impVol,tesco.impVol),type='l',col='green',lty=2)
legend(1.9,.55,cex=0.6,legend = c("Historical Volatility","Implied Volatility"),col=c("black","green"),lty=1:2)
  #NYT
nyt.table <- read.csv(file.choose()) # Import the nyt DATA FILE
head(nyt.table)
nyt.prices <- nyt.table$Close

nyt.sigma <- ann_volatility(nyt.prices)

timeLengths <- seq(252/2,756,252/2)

nyt.sigmas <- c()

for(i in seq(30,756,1)){
  sigma <- ann_volatility(tail(nyt.prices,i))
  nyt.sigmas <- c(nyt.sigmas,sigma)
}

nyt.P <- list(S=452.13,K=450,r=0.0033,T_exp=54/252)
nyt.impVol     <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=54/252,claimVal=mean(c(36.15,37.65)))

plot(x=seq(30/252,756/252,1/252),y=nyt.sigmas,type = 'l',
     ylim = c(0.25,.55),
     main = 'nyt Historical and Implied Volatilities',
     xlab = 'Past Time Period in Years',
     ylab = 'Volatility',
)
lines(x=c(0,10),y=c(nyt.impVol,nyt.impVol),type='l',col='green',lty=2)
legend(1.9,.55,cex=0.6,legend = c("Historical Volatility","Implied Volatility"),col=c("black","green"),lty=1:2)
  #NOVN
novartis.table <- read.csv(file.choose()) # Import the novartis DATA FILE
head(novartis.table)
novartis.prices <- novartis.table$Close

novartis.sigma <- ann_volatility(novartis.prices)

timeLengths <- seq(252/2,756,252/2)

novartis.sigmas <- c()

for(i in seq(30,756,1)){
  sigma <- ann_volatility(tail(novartis.prices,i))
  novartis.sigmas <- c(novartis.sigmas,sigma)
}

novartis.P <- list(S=452.13,K=450,r=0.0033,T_exp=54/252)
novartis.impVol     <- impliedVolatility(S=452.13,K=450,r=0.0033,T_exp=54/252,claimVal=mean(c(36.15,37.65)))

plot(x=seq(30/252,756/252,1/252),y=novartis.sigmas,type = 'l',
     ylim = c(0.25,.55),
     main = 'novartis Historical and Implied Volatilities',
     xlab = 'Past Time Period in Years',
     ylab = 'Volatility',
)
lines(x=c(0,10),y=c(novartis.impVol,novartis.impVol),type='l',col='green',lty=2)
legend(1.9,.55,cex=0.6,legend = c("Historical Volatility","Implied Volatility"),col=c("black","green"),lty=1:2)

# TESCO Smile

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


###############################################
###############################################
# Greeks
###############################################
###############################################

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


#Novartis Graphs

Plots of the greeks with Novartis Parameters
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

###############################################
###############################################
# Exotics
###############################################
###############################################

# Standard Call Example

N = 10000
P.Call <- list(S=277.5,K=280,r=.004452,T_exp=52/252,put=F,
               sigma=0.24,delta=0)
payoffs <- c()
for(i in 1:N){
  priceRun <- blackscholesPriceRunSim(P.Call,c(0,52/252))
  S_T <- tail(priceRun,1)
  K   <- P.Call$K
  payoff = max(0,S_T-K)
  payoffs = c(payoffs,payoff)
}
price <- exp(-52/252*.004452)*mean(payoffs)
price
black_scholes(P.Call)

# McDonald Table 14.1

P_tbl14.1_bs = bs_parameterize(S = 40, K = 40, r = 0.08,
                               sigma = 0.3, delta = 0, T_exp = 1)
n_vec <- c(1,2,3,5,10)
names <- c("Avg Price Call", "Avg Price Put", 
           "Avg Strike Call", "Avg Strike Put")
bs_price <- matrix(, nrow = length(n_vec), ncol = 4)
colnames(bs_price) <- names
k <- 1
N <- 10000
for (i in n_vec){
  time_v <- 1:i
  time_v <- time_v/i
  time_v <- c(0, time_v)
  step_v <- 30*time_v
  step_v <- ceiling(step_v)
  
  C_strike_avg <- 0
  P_strike_avg <- 0
  C_under_avg  <- 0
  P_under_avg  <- 0
  for (j in 1:N){
    vals <- blackscholesPriceRunSim(P_tbl14.1_bs, time_v)
    mean <- exp(mean(log(vals[2:(i+1)])))
    C_strike <- pmax(0, vals[i+1] - mean)
    P_strike <- pmax(0, mean - vals[i+1])
    C_under  <- pmax(0, mean - 40)
    P_under  <- pmax(0, 40 - mean)
    C_strike_avg <- C_strike_avg + (C_strike/N)
    P_strike_avg <- P_strike_avg + (P_strike/N)
    C_under_avg <- C_under_avg + (C_under/N)
    P_under_avg <- P_under_avg + (P_under/N)
  }
  bs_price[k, 1] = C_under_avg
  bs_price[k, 2] = P_under_avg
  bs_price[k, 3] = C_strike_avg
  bs_price[k, 4] = P_strike_avg
  
  k <- k + 1
}
bs_price


#McDonald Table 14.3
P <- list(S=0.9,sigma=0.1,r=0.06,delta=0.03,T_exp=0.5, put = T)
N <- 1000
Strikes = c(0.8,0.9,1.0)

payoffs = list("aA"=c(),"bB"=c(),"cC"=c(),"dD"=c(),"eE"=c(),
               "fF"=c(),"gG"=c(),"hH"=c(),"iI"=c(),"jJ"=c(),"kK"=c(),
               "lL"=c(),"mM"=c(),"nN"=c(),"oO"=c())
prices  = rep(0,15)

for(i in 1:N){
  pricerun = blackscholesPriceRunSim(P,(0:126)/256)
  P$K = 0.8
  payoffaA = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payoffbB = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffcC = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffdD = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffeE = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  P$K = 0.9
  payofffF = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payoffgG = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffhH = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffiI = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffjJ = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  P$K = 1.0
  payoffkK = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.80)
  payofflL = barrierPayoff(P,pricerun,up=F,lockIn=T,barrier=0.85)
  payoffmM = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=0.95)
  payoffnN = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.00)
  payoffoO = barrierPayoff(P,pricerun,up=T,lockIn=F,barrier=1.05)
  
  payoffs$aA = c(payoffs$aA,payoffaA)
  payoffs$bB = c(payoffs$bB,payoffbB)
  payoffs$cC = c(payoffs$cC,payoffcC)
  payoffs$dD = c(payoffs$dD,payoffdD)
  payoffs$eE = c(payoffs$eE,payoffeE)
  payoffs$fF = c(payoffs$fF,payofffF)
  payoffs$gG = c(payoffs$gG,payoffgG)
  payoffs$hH = c(payoffs$hH,payoffhH)
  payoffs$iI = c(payoffs$iI,payoffiI)
  payoffs$jJ = c(payoffs$jJ,payoffjJ)
  payoffs$kK = c(payoffs$kK,payoffkK)
  payoffs$lL = c(payoffs$lL,payofflL)
  payoffs$mM = c(payoffs$mM,payoffmM)
  payoffs$nN = c(payoffs$nN,payoffnN)
  payoffs$oO = c(payoffs$oO,payoffoO)
}
for(i in 1:15){
  prices[i] = exp(-0.06*0.5)*mean(payoffs[[i]])
}
table = matrix(prices,3,5,byrow=T)
table


