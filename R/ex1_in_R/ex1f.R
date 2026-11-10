rm(list = ls())

W_U = function(mu, lam, p)
{
  1/(mu-lam*p)
}

W_N = function(mu, lam, p)
{
  mu/((mu-lam)*(mu-p*lam))
}

p_vals = seq(from=0.001, to=0.999, length.out=100) #>0, <1?
mu = 1/10 #per min
lambda = 5/60 #per min
WU_vals = W_U(mu, lambda, p_vals)
WN_vals = W_N(mu, lambda, p_vals)
#plot(p_vals, WU_vals, type = 'l', col = 'red', xlab = 'p', ylab = 'W', main = 'Expected time spent in UCC by urgent and normal patients respectively, as function of p')
plot(NA, xlab = 'p', ylab='W', main = 'Expected time spent in UCC by urgent and normal patients respectively, as function of p',xlim = c(0,1), ylim = c(min(WN_vals, WU_vals), max(WN_vals, WU_vals)))
lines(p_vals, WN_vals, col = 'blue')
lines(p_vals, WU_vals, col = 'red')
legend('topleft', legend = c('Normal', 'Urgent'), col = c('blue', 'red'), lty = 1)
grid()

