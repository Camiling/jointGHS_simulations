library(BAS)
library(tailoredGlasso)

# Implement tau selection based on CCH distribution of the shrinkage coefficients 
# and solving Piironen & Vehtari effective degrees of freedom equal to prior expected number of edges equation

# Shows that it underselects edges, and does not give the shrinkage factor a horseshoe shape. 
# No plots are saved

H.gordy = function(p, q, r, s, nu, thet){
  # Compound confluent hypergeometric function
  nu^(-p)*exp(-s/nu)*BAS::phi1(q, r, p+q, s/nu, 1-thet)
}

E.CCH = function(p, q, r, s, nu, thet){
  # Find expectation (first moment) of CCH distribution
  res = p/(p+q) * H.gordy(p+1, q, r, s, nu, thet)/H.gordy(p, q, r, s, nu, thet)
  return(res)
}

E.shrinkagefactor = function(tau_sq, theta.hat, theta.diag.prior=1){
  # Find expectation of the shrinkage factor for a given neighbour j of node p 
  # tau is the global shrinkage parameter, theta.hat is the estimated prec matrix entry, theta.diag.prior is the p'th diagonal element
  # Must scale theta.hat as in formula
  E.CCH(1, 1/2, 1, (theta.hat/sqrt(theta.diag.prior))^2/2, 1, theta.diag.prior/tau_sq)
}

shinkagefactor.density = function(z, s, thet){
  # Use definition of CCH denstity by Gordy (1998)
  p = 1
  q = 1/2
  r = 1
  nu = 1
  res = z^(p-1)*(1-nu*z)^(q-1)*(thet + (1-thet)*nu*z)^(-r)*exp(-s*z)/(beta(p,q)*H.gordy(p,q,r,s,nu,thet))
  return(res)
}

E.meff.p = function(tau_sq, theta.vec, p, theta.diag.prior=1){
  # Find expected effective number of degrees of freedom through the shrinkage factors of all potential neighbours
  # Must supply theta_pj values for j=1,...,p-1 in theta.vec
  # theta.diag.prior is the p'th diagonal element of theta
  # Assumes scaled variables, so X'X has ones on the diagonal
  res = (p-1) - sum(sapply(theta.vec, FUN= function(s) E.shrinkagefactor(tau_sq, s, theta.diag.prior)))
  return(res)
}

E.meff.sum = function(tau_sq, theta, p, theta.diag.prior){
  sum(sapply(1:p, FUN = function(i) E.meff.p(tau_sq, theta[i,-i], p, theta.diag.prior[i])))
}

meff.solve = function(tau_sq, theta, p0, theta.diag.prior=NULL){
  # Note that tau cannot be zero, as we get zero division
  p = ncol(theta)
  if(is.null(theta.diag.prior)){
    theta.diag.prior = rep(1,p)
  }
  meff.vals = E.meff.sum(tau_sq, theta, p, theta.diag.prior)
  return(meff.vals-2*p0)
}



# Test on small data set 

p = 100
n=150

set.seed(123)
gg = huge.generator(n=n,d=p, graph='scale-free')
gg$omega[1:10,1:10]
theta.true = (as.matrix(gg$theta)+0)!=0
gg$sparsity
# 0.02

# Scale data first
dat = scale(gg$data)
# Use trick to invert
v0=10
sigma.0 = (t(dat)%*%dat + v0 * diag(p)) / (n-1)
theta.0 = as.matrix(Matrix::nearPD(solve(sigma.0))$mat)
# Prior expectation of the number of edges
p0=p

# Check how the average expected effective degrees of freedom change with tau
taus.sq = seq(0.0005,0.005, length.out = 40)
meff.vals = sapply(taus.sq, FUN = function(t) E.meff.sum(t, theta.0, p, rep(1, p)))
plot(taus.sq,meff.vals/2)
# Plot the prior guess of nonzero edges as a horizontal line. Our selected tau will be at the intersection
abline(h=p0)


tau.sq.meff = uniroot(meff.solve, c(0.0005,0.005), theta=theta.0, p0=p0, theta.diag.prior=rep(1,p))$root
tau.sq.meff
# 0.001036261

# See what value the OG Piironen & Vehtari approach would give
alpha.prior = gg$sparsity
tau.prior = alpha.prior/(1-alpha.prior)/sqrt(n)
tau.prior^2
# 4.164931e-06
# Very small, gives an empty graph

# Test selected value 
res.meff = fastGHS::fastGHS(dat, tau_sq = tau.sq.meff, fix_tau = T)
theta.meff = res.meff$theta
theta.meff[which(abs(theta.meff)<1e-5)] = 0 
sparsity(theta.meff) 
# 0
precision(theta.true,theta.meff!=0)

recall(theta.true,theta.meff!=0)


# Try with "true" values
taus.sq = seq(0.0005,0.005, length.out = 40)
meff.vals.true = sapply(taus.sq, FUN = function(t) E.meff.sum(t, gg$omega, p, diag(gg$omega)))
plot(taus.sq,meff.vals.true/2)
# Plot the prior guess of nonzero edges as a horizontal line. Our selected tau will be at the intersection
abline(h=p0)

# Still selects a too small value of tau => problem is not in estmiation of theta.


# Investigate larger values of tau
taus.sq.len = seq(0.0005,0.15, length.out = 40)
meff.vals.len = sapply(taus.sq.len, FUN = function(t) E.meff.sum(t, gg$omega, p, diag(gg$omega)))
plot(taus.sq.len,meff.vals.len/2)
# Plot the prior guess of nonzero edges as a horizontal line. Our selected tau will be at the intersection
abline(h=p0)

# E[meff|tau] gets very large for the right range of tau...


# Plot distribution of shrinkage factor, given an element of theta and tau_sq

z.vals = seq(0,1,length.out = 100)
tau_sq = 0.01
z.dens = shinkagefactor.density(z.vals, s=0.1^2/2, thet=1/tau_sq)
plot(z.vals, z.dens, type='l')

# Issue: does not have the horseshoe shape, for any tau or theta_ij => only large shrinkage. 




