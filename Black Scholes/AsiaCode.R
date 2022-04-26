asianAveragePriceOptions_bs <- function(P,N="Infinity"){
  # Uses Equations in Appendix 14.A to price an asian option.
  # Asian option calculated with Geometric Average as S_T, and at N evenly spaced time intervals in [0,T_exp]
  r     = P$r
  delta = P$delta
  sigma = P$sigma
  if(N == "Infinity"){
    deltaNew = (r+delta+sigma^2/6)/2
    sigmaNew = sigma * sqrt(1/3)
  } else {
    deltaNew = 0.5*(r*(N-1)/N+(delta+sigma^2/2)*(N+1)/N-sigma^2/N^2*(N+1)*(2*N+1)/(6))
    sigmaNew = sigma/N*sqrt((N+1)*(2*N+1)/(6))
  }
  P$delta = deltaNew
  P$sigma = sigmaNew
  return(black_scholes(P))
}

asianAverageStrikeOptions_bs <- function(P,N="Infinity"){
  # Uses Equations in Appendix 14.A to price an asian option.
  # Asian option calculated with Geometric Average as K, and at N evenly spaced time intervals in [0,T_exp]
  if(N=="Infinity"){N = 10^8}
  r        = P$r
  delta    = P$delta
  sigma    = P$sigma
  rNew     = 0.5*(r*(N-1)/N+(delta+sigma^2/2)*(N+1)/N-sigma^2/N^2*(N+1)*(2*N+1)/(6))
  corr     = 0.5*sqrt(6*(N+1)/(2*N+1))
  piece    = (N+1)*(2*N+1)/(6)/(N^2)
  sigmaNew = sigma*sqrt(1+piece-2*corr*sqrt(piece))
  strikeNew = P$S
  P$r     = rNew
  P$sigma = sigmaNew
  P$K     = strikeNew
  return(black_scholes(P))
}
