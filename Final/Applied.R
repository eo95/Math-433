###############################################
###############################################
# Historic Volatility
###############################################
###############################################

tesco.data <- read.csv(file.choose())
adobe.data <- read.csv(file.choose())
nyt.data <- read.csv(file.choose())
novartis.data <- read.csv(file.choose())

head(tesco.data)
summary(tesco.data)
tesco.sigma <- ann_volatility(as.numeric(tesco.data$Open[1:(252*0.5)]))

tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=270,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.P_bs <- bs_parameterize(S=277.5,r=.004452,T_exp=52/252,K=270,sigma=tesco.sigma,eur=F,put=F)
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]
black_scholes(tesco.P_bs)
tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=280,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.P_bs$K = 280
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]
black_scholes(tesco.P_bs)
tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=290,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.P_bs$K = 290
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]
black_scholes(tesco.P_bs)

head(adobe.data)
summary(adobe.data)
adobe.sigma <- ann_volatility(adobe.data$Open[1:(252*0.5)])


adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=440,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.P_bs <- bs_parameterize(S=452.13,r=.0033,T_exp=54/252,K=440,sigma=adobe.sigma,eur=F,put=F)
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]
black_scholes(adobe.P_bs)
adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=450,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.P_bs$K = 450
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]
black_scholes(adobe.P_bs)
adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=460,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.P_bs$K = 460
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]
black_scholes(adobe.P_bs)

head(nyt.data)
summary(nyt.data)
nyt.sigma <- ann_volatility(nyt.data$Open[1:(252*0.5)])

nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=45,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.P_bs <- bs_parameterize(S=46.01,r=.0021,T_exp=29/252,K=45,sigma=nyt.sigma,D_CF=c(22/252,.09),eur=F,put=F)
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]
black_scholes_div(nyt.P_bs)
nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=46,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.P_bs$K = 46
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]
black_scholes_div(nyt.P_bs)
nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=47,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.P_bs$K = 47
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]
black_scholes_div(nyt.P_bs)

head(novartis.data)
summary(novartis.data)
novartis.sigma <- ann_volatility(novartis.data$Open[1:(252*0.5)])

novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=76,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.P_bs <- bs_parameterize(S=77.64,r=-.0075,T_exp=53/252,K=76,sigma=novartis.sigma,D_CF=c(4/252,3.1),eur=F,put=F)
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]
black_scholes_div(novartis.P_bs)
novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=78,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.P_bs$K = 78
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]
black_scholes_div(novartis.P_bs)
novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=80,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.P_bs$K = 80
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]
black_scholes_div(novartis.P_bs)

adobe.sigma
adobe.model[[2]][1]

tesco.sigma
tesco.model[[2]][1]

novartis.sigma
novartis.model[[2]][1]

nyt.sigma
nyt.model[[2]][1]





###############################################
###############################################
# Implied Volatility
###############################################
###############################################

# TESCO ----

tesco.OTM.ImpVol <- impliedVolatility(S=277.5,K=290,r=.004452,T_exp=52/252,claimVal = mean(c(4.00,12.00)))
tesco.ITM.ImpVol <- impliedVolatility(S=277.5,K=270,r=.004452,T_exp=52/252,claimVal = mean(c(14.00,22.00)))

tesco.sigma = mean(c(tesco.OTM.ImpVol,tesco.ITM.ImpVol))

tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=280,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]


# ADOBE ----

adobe.OTM.ImpVol <- impliedVolatility(S=452.13,K=460,r=.0033,T_exp=54/252,claimVal = mean(c(31.35,32.65)))
adobe.ITM.ImpVol <- impliedVolatility(S=452.13,K=440,r=.0033,T_exp=54/252,claimVal = mean(c(41.55,43.00)))

adobe.sigma = mean(c(adobe.OTM.ImpVol,adobe.ITM.ImpVol))

adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=450,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]


# NYT ----

nyt.OTM.ImpVol <- impliedVolatility(S=46.01,K=47,r=.0021,T_exp=29/252,D_CF = c(22/252,.09),claimVal = mean(c(1.50,2.55)))
nyt.ITM.ImpVol <- impliedVolatility(S=46.01,K=46,r=.0021,T_exp=29/252,D_CF = c(22/252,.09),claimVal = mean(c(2.55,4.00)))

nyt.sigma = mean(c(nyt.OTM.ImpVol,nyt.ITM.ImpVol))

nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=25,K=46,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]

# NOVARTIS ----

novartis.OTM.ImpVol <- impliedVolatility(S=77.64,K=80,r=-.0075,T_exp=53/252,D_CF = c(4/252,3.1),claimVal = mean(c(1.28,1.87)))
novartis.ITM.ImpVol <- impliedVolatility(S=77.64,K=76,r=-.0075,T_exp=53/252,D_CF = c(4/252,3.1),claimVal = mean(c(2.63,3.54)))

novartis.sigma = mean(c(novartis.OTM.ImpVol,novartis.ITM.ImpVol))

novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=25,K=78,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]


###############################################
###############################################
# Currency Exchange Option
###############################################
###############################################

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
