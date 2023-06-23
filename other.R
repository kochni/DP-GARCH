
## other stuff

packages_other = c(
  "dirmult"
)

lapply(packages_other, library, character.only = TRUE)

## random measures

gamma_dirichlet = function(theta0) {
  
  omega = rgamma(1, 1, theta0[1])
  alpha_beta = rdirichlet(1, c(1, theta0[2]))
  alpha = alpha_beta[1]
  beta = alpha_beta[2]
  
  coef = list("omega"=omega, "alpha"=alpha, "beta"=beta)
  return(coef)
}


## Chinese restaurant process

CRP = function(n, alpha) {
  
  K = 0     # no. of unique clusters
  z = c()   # cluster allocations
  N_k = c() # cluster sizes
  for (i in 1:n) {
    t = alpha / (alpha+i-1)   # Prob. of opening new cluster
    u = runif(1, 0, 1)
    if (u <= t) {
      z[i] = K + 1
    }
    else {
      p_k = N_k / (alpha+i-1)
      z[i] = sample(1:K, size=1, prob=p_k)
    }
    
    K = length(unique(z))
    for (k in 1:K) N_k[k] = sum(z==k)
    
  }
  
  K = length(unique(z))
  pi = N_k / n
  
  return(list("z"=z, "pi"=pi, "K"=K))
}


## Dirichlet process

DP = function(n, alpha, G0, G0_paras) {
  
  # get unique values & mixing coefficients from CRP
  crp = CRP(n=n, alpha=alpha)
  K = crp$K
  z = crp$z
  pi = crp$pi
  
  # draw the unique values iid from the base distribution
  theta = lapply(rep(list(paras), K), FUN=G0)
  
  return(list("K"=K, "pi"=pi, "z"=z, "theta"=theta))
  
}


## simulate from DP

A = DP(n=10000, alpha=1, G0=gamma_dirichlet, G0_paras=c(1,1))
X = 1
s = 1

vol = c()
vol_2 = c()
vol_4 = c()
for (i in 1:n) {
  r = A$z[i]
  omega = A$theta[[r]]$omega
  alpha = A$theta[[r]]$alpha
  beta = A$theta[[r]]$beta
  s2 = omega + alpha*X^2 + beta*s
  vol_2[i] = s2
  vol_4[i] = s2^2
  vol[i] = sqrt( s2 )
  
}
c("Kurt. incr.:" = round(mean(vol_4)/mean(vol_2)^2, 3))
c("vol"=mean(vol), "vol^2"=mean(vol_2), "vol^4"=mean(vol_4))





stickbreak = function(theta0) {
  
  omega = rgamma(1, theta0[1], 1)
  break1 = rbeta(1, theta0[2], 1)
  break2 = (1-alpha)*rbeta(1, theta0[2], 1)
  i = sample(1:2, size=1)
  alpha = c(break1, break2)[i]
  beta = c(break2, break1)[i]
  
  coef = list("omega"=omega, "alpha"=alpha, "beta"=beta)
  return(coef)
}
