###############################################
###############################################
# Homework 2
###############################################
###############################################


# Problem 5

P <- parameterize(S=4,K=2,r=1/3,T_exp=2,n=2)
payoff <- call_payoff
S_m <- c(4,0,0,6,2,0,9,3,1)
S_m <- matrix(S_m,3,3)
C_m <- matrix(0,P$n+1,P$n+1)
A_m <- matrix(0,P$n+1,P$n+1)
B_m <- matrix(0,P$n+1,P$n+1)
S_T <- S_m[,P$n+1]
C_m[,P$n+1] <- payoff(S_T,P$K)
solution <- solve_binomial_pricing_recombine(S_m,C_m,A_m,B_m,P$r_v,P$h,P$n,P$eur,payoff=payoff,P$K,P$delta,P$D_v)

# Standard Hedge: 
# Pu = Pd but Pu != P so this is not self financing

Pu = solution[[2]][4] - solution[[3]][1]*solution[[1]][4]
Pd = solution[[2]][5] - solution[[3]][1]*solution[[1]][5]
P  = solution[[2]][1] - solution[[3]][1]*solution[[1]][1]

# Hedge with cash:
# Pu = Pd = P = 0

Pu = solution[[2]][4] - solution[[3]][1]*solution[[1]][4] - solution[[4]][1]*exp(1/3)
Pd = solution[[2]][5] - solution[[3]][1]*solution[[1]][5] - solution[[4]][1]*exp(1/3)
P  = solution[[2]][1] - solution[[3]][1]*solution[[1]][1] - solution[[4]][1]



###############################################
###############################################
# Homework 3
###############################################
###############################################

# Problem 3
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
     xlab = "Strike Price\n",
     ylab = "Value of the Straddle")
xpoints <- c(K_min,K_min)
ypoints <- c(0,100)
lines(xpoints,ypoints,col="blue",lty = 2,type = 'l')
legend(100,25,cex=0.8,legend = c("Calculated Strike"),lty = 2,col='blue')



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
     xlab = "Volatility\n",
     ylab = "Claim Value of a Call Option")

puts <- BlackScholes(100,90,sigmas,.05,3/12,put=T)

plot(sigmas,puts,type='l',
     main = "Put Value v Volatility",
     sub  = "Parameters: S=100 K=90 r=.05 T=3/12",
     xlab = "Volatility\n",
     ylab = "Claim Value of a Put Option")
