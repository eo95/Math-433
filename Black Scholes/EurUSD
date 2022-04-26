eurusd.table <- read.csv(file.choose()) # Import the EURUSD DATA FILE
eurusd.prices <- eurusd.table$Close
eurusd.sigma <- ann_volatility(eurusd.prices[1:28])
P.Call <- list(S=1.0929,K=1.0946,r=.0033,sigma=eurusd.sigma,delta=-.00585,T_exp=41/252,put=F)
asianAveragePriceOptions_bs(P.Call)

N = 100000
payoffs = c()
for(i in 1:N){
  priceRun <- blackscholesPriceRunSim(P.Call,(0:41)/252)
  S_TRun <- tail(priceRun,1)
  payoff <- asianPayoff(P.Call,priceRun,S_Trun,avgAsStrike = F, algebraic = F)
  payoffs <- c(payoffs,payoff)
}
price <- exp(-41/252*.0033)*mean(payoffs)
price
