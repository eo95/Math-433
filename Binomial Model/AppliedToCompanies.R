tesco.data <- read.csv(file.choose())
adobe.data <- read.csv(file.choose())
nyt.data <- read.csv(file.choose())
novartis.data <- read.csv(file.choose())

head(tesco.data)
summary(tesco.data)
tesco.sigma <- ann_volatility(as.numeric(tesco.data$Open[1:(252*0.5)]))

tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=270,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]
tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=280,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]
tesco.P <- parameterize(S=277.5,r=.004452,T_exp=52/252,n=25,K=290,sigma=tesco.sigma,choice=2,eur=F,recombine = F)
tesco.model <- binomial_pricing(tesco.P,payoff=call_payoff)
tesco.model[[2]][1]

head(adobe.data)
summary(adobe.data)
adobe.sigma <- ann_volatility(adobe.data$Open[1:(252*0.5)])


adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=440,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]
adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=450,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]
adobe.P <- parameterize(S=452.13,r=.0033,T_exp=54/252,n=25,K=460,sigma=adobe.sigma,choice=2,eur=F,recombine = F)
adobe.model <- binomial_pricing(adobe.P,payoff=call_payoff)
adobe.model[[2]][1]

head(nyt.data)
summary(nyt.data)
nyt.sigma <- ann_volatility(nyt.data$Open[1:(252*0.5)])

nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=45,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]
nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=46,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]
nyt.P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=20,K=47,sigma=nyt.sigma,choice=2,D_CF=c(22/252,.09),eur=F,recombine = F)
nyt.model <- binomial_pricing(nyt.P,payoff=call_payoff)
nyt.model[[2]][1]

head(novartis.data)
summary(novartis.data)
novartis.sigma <- ann_volatility(novartis.data$Open[1:(252*0.5)])

novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=76,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]
novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=78,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]
novartis.P <- parameterize(S=77.64,r=-.0075,T_exp=53/252,n=20,K=80,sigma=novartis.sigma,choice=2,D_CF=c(4/252,3.1),recombine = F)
novartis.model <- binomial_pricing(novartis.P,payoff=call_payoff)
novartis.model[[2]][1]

adobe.sigma
adobe.model[[2]][1]

tesco.sigma
tesco.model[[2]][1]

novartis.sigma
novartis.model[[2]][1]

nyt.sigma
nyt.model[[2]][1]

