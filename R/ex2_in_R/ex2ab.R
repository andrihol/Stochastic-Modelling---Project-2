rm(list = ls())

#rounding to avoid comparison errors in setdiff
theta = round(seq(from=0.25, to=0.50, by=0.005),3)
theta_obs = round(c(0.300, 0.350, 0.390, 0.410, 0.450),3)
theta_unkn = setdiff(theta, theta_obs)
mu = rep(0.5, length(theta_unkn))
y_obs = c(0.5, 0.32, 0.40, 0.35, 0.60)
mu_obs = rep(0.5, length(theta_obs))


covar = function(theta_i,theta_j)
{
  (1 +15*abs(theta_i-theta_j))*exp(-15*abs(theta_i-theta_j))*0.5^2
}

#constructing covariance matrices Sigma11, Sigma12, Sigma22, Sigma21
sigma11 = matrix(c(1:length(theta_unkn)^2)*0, nrow=length(theta_unkn), ncol = length(theta_unkn))
for(i in 1:length(theta_unkn))
{
  for (j in 1:length(theta_unkn))
  {
    sigma11[i,j] = covar(theta_unkn[i], theta_unkn[j])
  }
}

sigma12 = matrix(c(1:length(theta_unkn)*length(theta_obs))*0, nrow=length(theta_unkn), ncol=length(theta_obs))
for (i in 1:length(theta_unkn))
{
  for(j in 1:length(theta_obs))
  {
    sigma12[i,j] = covar(theta_unkn[i], theta_obs[j])
  }
}

sigma21 = matrix(c(1:length(theta_obs)*length(theta_unkn))*0, nrow=length(theta_obs), ncol=length(theta_unkn))
for(i in 1:length(theta_obs))
{
  for(j in 1:length(theta_unkn))
  {
    sigma21[i,j] = covar(theta_obs[i], theta_unkn[j])
  }
}

sigma22 = matrix(c(1:length(theta_obs)^2)*0, nrow=length(theta_obs), ncol=length(theta_obs))
for(i in 1:length(theta_obs))
{
  for(j in 1:length(theta_obs))
  {
    sigma22[i,j] = covar(theta_obs[i], theta_obs[j])
  }
}

sigma22_inv = solve(sigma22)

cond_exp = mu + sigma12%*%sigma22_inv%*%(y_obs - mu_obs)
cond_cov = sigma11 - sigma12%*%sigma22_inv%*%sigma21

plot(theta_unkn, cond_exp, main="Plot of prediction/conditional expectation on the observed data with 90% PIs", ylim = c(0,1))
points(theta_obs, y_obs, col="red")

#y_pred = rbind(cond_exp, y_obs)
#theta_pred = rbind(theta_unkn, theta_obs)

z = qnorm(0.95) # 5% quantile
for (i in 1:length(cond_exp))
{
  upper = cond_exp[i] + z*sqrt(cond_cov[i,i])
  lower = cond_exp[i] - z*sqrt(cond_cov[i,i])
  segments(theta_unkn[i], lower, theta_unkn[i], upper, col = "gray")
}
grid()

#exercise b

probabilities = c(1:length(cond_exp))*0
for (i in 1:length(cond_exp))
{
  probabilities[i] = pnorm((0.30-cond_exp[i])/sqrt(cond_cov[i,i]))
}
plot(theta_unkn, probabilities, main="P(Y(theta)<30 | Y_obs = y_obs) as function of theta", ylim=c(0,1), ylab="p", xlab="theta")
points(theta_obs, rep(1,length(theta_obs)), col="red")
grid()
