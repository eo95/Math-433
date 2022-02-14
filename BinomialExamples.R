# Testing Functions with Exercises, Examples and Problems:

# Pole Exercises: ----

# Problem 2
# (a)
P_CRR.256  <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=256)
P_CRR.512  <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=512)
P_JR.256   <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=256)
P_JR.512   <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=512)
P_MCRR.256 <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=256)
P_MCRR.512 <- generate_P(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=512)

C_1 <- binomial_pricing(P_CRR.256,payoff=put_payoff)[[2]][1]
C_2 <- binomial_pricing(P_CRR.512,payoff=put_payoff)[[2]][1]
C_3 <- binomial_pricing(P_JR.256,payoff=put_payoff)[[2]][1]
C_4 <- binomial_pricing(P_JR.512,payoff=put_payoff)[[2]][1]
C_5 <- binomial_pricing(P_MCRR.256,payoff=put_payoff)[[2]][1]
C_6 <- binomial_pricing(P_MCRR.512,payoff=put_payoff)[[2]][1]

sol_2a <- matrix(c(C_1,C_2,C_3,C_4,C_5,C_6),3,2,byrow=T)
rownames(sol_2a) <- c("CRR","JR","MCRR")
colnames(sol_2a) <- c("N=256","N=512")

# (b)
P_CRR.256  <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=256)
P_CRR.512  <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=512)
P_JR.256   <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=256)
P_JR.512   <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=512)
P_MCRR.256 <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=256)
P_MCRR.512 <- generate_P(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=512)

C_1 <- binomial_pricing(P_CRR.256,payoff=call_payoff)[[2]][1]
C_2 <- binomial_pricing(P_CRR.512,payoff=call_payoff)[[2]][1]
C_3 <- binomial_pricing(P_JR.256,payoff=call_payoff)[[2]][1]
C_4 <- binomial_pricing(P_JR.512,payoff=call_payoff)[[2]][1]
C_5 <- binomial_pricing(P_MCRR.256,payoff=call_payoff)[[2]][1]
C_6 <- binomial_pricing(P_MCRR.512,payoff=call_payoff)[[2]][1]

sol_2b <- matrix(c(C_1,C_2,C_3,C_4,C_5,C_6),3,2,byrow=T)
rownames(sol_2b) <- c("CRR","JR","MCRR")
colnames(sol_2b) <- c("N=256","N=512")
sol_2a
sol_2b

# McDonald Examples: ----

# Figure 10.4

P_fig10.4 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=2,n=2,choice=1)
binomial_pricing(P_fig10.4,payoff=call_payoff)

# Figure 10.5

P_fig10.5 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.5,payoff=call_payoff)

# Figure 10.6

P_fig10.6 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1)
binomial_pricing(P_fig10.6,payoff=put_payoff)

# Figure 10.7

P_fig10.7 <- generate_P(S=41,K=40,sigma=0.3,r=0.08,T_exp=1,n=3,choice=1,eur=F)
binomial_pricing(P_fig10.7,payoff=put_payoff)

# Figure 10.8

P_fig10.8 <- generate_P(S=110,K=100,sigma=0.3,r=0.05,T_exp=1,n=3,choice=1,eur=F,delta=0.035)
binomial_pricing(P_fig10.8,payoff=call_payoff)

# McDonald Problems: ----

# Problem 10.1

P_10.1 <- generate_P(S=100,K=105,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)

#(a)
binomial_pricing(P_10.1,payoff=call_payoff)

#(b)
binomial_pricing(P_10.1,payoff=put_payoff)

# Problem 10.2
P_10.2 <- generate_P(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.2,payoff=call_payoff)[[2]][1]

# Problem 10.3

P_10.3 <- generate_P(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.3,payoff=put_payoff)[[2]][1]

# Problem 10.6

P_10.6 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.6,payoff=call_payoff)

# Problem 10.8

P_10.8 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.8,payoff=put_payoff)

# Problem 10.10

P_10.10 <- generate_P(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,eur=F,n=2)
binomial_pricing(P_10.10,payoff=call_payoff)
