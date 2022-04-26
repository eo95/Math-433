P <- parameterize(S=46.01,r=.0021,T_exp=29/252,n=25,K=46,sigma=0.29,choice=1,D_CF=c(22/252,.09),recombine=F)
binomial_pricing(P,call_payoff)[[2]][1]

P <- bs_parameterize(S=46.01,r=0.0021,T_exp=29/252,K=46,sigma=0.29,D_CF=c(.09,22/252))
black_scholes_div(P)
