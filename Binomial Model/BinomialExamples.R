# Testing Functions with Exercises, Examples and Problems:

# Pole Exercises: ----

# Convergence
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
  P_EQP  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=5,n=k)
  y_4    <- append(y_4, binomial_pricing(P_EQP,payoff=put_payoff)[[2]][1])
  P_TRG  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=6,n=k)
  y_5    <- append(y_5, binomial_pricing(P_TRG,payoff=put_payoff)[[2]][1])
}
plot(x,y_1, col="red")
lines(x,bs)
plot(x,y_2, col="green")
lines(x,bs)
plot(x,y_3, col="blue")
lines(x,bs)
plot(x,y_4, col="yellow")
lines(x,bs)
plot(x,y_5, col="black")
lines(x,bs)
x <- c()
y_1 <- c()
y_2 <- c()
y_3 <- c()
y_4 <- c()
y_5 <- c()
  for (i in 1:64){
    k = 8*i
    x      <- append(x,k)
    C_CRR  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=k)
    y_1    <- append(y_1, binomial_pricing(P_CRR,payoff=call_payoff)[[2]][1])
    C_JR   <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=k)
    y_2    <- append(y_2, binomial_pricing(P_JR,payoff=call_payoff)[[2]][1])
    C_MCRR <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=k)
    y_3    <- append(y_3, binomial_pricing(P_MCRR,payoff=call_payoff)[[2]][1])
    C_EQP  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=5,n=k)
    y_4    <- append(y_4, binomial_pricing(P_EQP,payoff=call_payoff)[[2]][1])
    C_TRG  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=6,n=k)
    y_5    <- append(y_4, binomial_pricing(P_TRG,payoff=call_payoff)[[2]][1])
  }
plot(x,y_1, col="red")
lines(x,bs)
plot(x,y_2, col="green")
lines(x,bs)
plot(x,y_3, col="blue")
lines(x,bs)
plot(x,y_4, col="yellow")
lines(x,bs)
plot(x,y_5, col="black")
lines(x,bs)

# Problem 2
# (a)
P_CRR.256  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=256)
P_CRR.512  <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=2,n=512)
P_JR.256   <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=256)
P_JR.512   <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=3,n=512)
P_MCRR.256 <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=256)
P_MCRR.512 <- parameterize(T_exp=1,S=5,K=10,r=0.12,sigma=0.5,choice=4,n=512)

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
P_CRR.256  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=256)
P_CRR.512  <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=2,n=512)
P_JR.256   <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=256)
P_JR.512   <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=3,n=512)
P_MCRR.256 <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=256)
P_MCRR.512 <- parameterize(T_exp=0.25,S=60,K=65,r=0.08,sigma=0.3,choice=4,n=512)

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

# Figure 10.8

P_fig10.8 <- parameterize(S=110,K=100,sigma=0.3,r=0.05,T_exp=1,n=3,choice=1,eur=F,delta=0.035)
binomial_pricing(P_fig10.8,payoff=call_payoff)

# McDonald Problems: ----

# Problem 10.1

P_10.1 <- parameterize(S=100,K=105,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)

#(a)
binomial_pricing(P_10.1,payoff=call_payoff)

#(b)
binomial_pricing(P_10.1,payoff=put_payoff)

# Problem 10.2
P_10.2 <- parameterize(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.2,payoff=call_payoff)[[2]][1]

# Problem 10.3

P_10.3 <- parameterize(S=100,K=95,r=0.08,T_exp=0.5,u=1.3,d=0.8,n=1)
binomial_pricing(P_10.3,payoff=put_payoff)[[2]][1]

# Problem 10.6

P_10.6 <- parameterize(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.6,payoff=call_payoff)

# Problem 10.8

P_10.8 <- parameterize(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,n=2)
binomial_pricing(P_10.8,payoff=put_payoff)

# Problem 10.10

P_10.10 <- parameterize(S=100,K=95,r=0.08,T_exp=1,u=1.3,d=0.8,eur=F,n=2)
binomial_pricing(P_10.10,payoff=call_payoff)







## Homework Problems


# Problem 5

P <- parameterize(S=4,K=2,r=1/3,T_exp=2,n=2)
payoff <- call_payoff
S_m <- c(4,0,0,6,2,0,9,3,1)
S_m <- matrix(S_m,3,3)
C_m <- matrix(0,P$n+1,P$n+1)
D_m <- matrix(0,P$n+1,P$n+1)
B_m <- matrix(0,P$n+1,P$n+1)
S_T <- S_m[,P$n+1]
C_m[,P$n+1] <- payoff(S_T,P$K)
solution <- solve_binomial_pricing_recombine(S_m,C_m,D_m,B_m,P$r_v,P$h,P$n,P$eur,payoff=payoff,P$K,P$delta,P$D_v)

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
  q <- (exp(r) - d)
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
plot(x,y_1)
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
plot(x,y_1)
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
plot(x,y_1)
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
plot(x,y_1)
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
plot(x,y_1)
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
plot(x,y_1)
lines(x,y_2)





# Run time --- 
T_exp = 1
n = 1
S = 50
K = 0
r = 0.05
choice = 0
u = 1.5
d = 0.75
x <- c()
y <- c()
for (i in 1:30){
  n <- i
  x <- append(x,n)
  P <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n)
  start_time <- Sys.time()
  append(y, binomial_pricing(P,call_payoff))
  end_time   <- Sys.time()
  y          <- append(y, as.numeric(end_time-start_time))
}
plot(x,y, col="blue")
x <- c()
y <- c()
for (i in 1:30){
  n <- i
  x <- append(x,n)
  P <- parameterize(T_exp=T_exp,S=S,K=K,r=r,choice=choice,u=u,d=d,n=n,recombine=F)
  start_time <- Sys.time()
  append(y, binomial_pricing(P,call_payoff))
  end_time   <- Sys.time()
  y         <- append(y, as.numeric(end_time-start_time))
}
lines(x,y, col="green")