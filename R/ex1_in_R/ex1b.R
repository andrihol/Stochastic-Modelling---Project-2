rm(list = ls())
#simulate stoch process for 50 days

#simulate with two competing processes. In each state, we jump to state i+1 if the first event is a birth, and to i-1 if the first event is a death

simulate_exp = function(rate) #inversion method
{
  u = runif(1)
  #find corresponding x
  x = - log(1-u)/rate
  x
}

simulate_N_days = function(N, X0, lambda, mu)
{
  X_arr = c(1: N*24*60)*0
  t = 0
  sojourn_times = c(1: N*24*60)*0 #element i = time spent in state X_arr[i]. R automaticaaly increases array size if an index larger than N is attempted modified
  X = X0
  idx = 1
  while(t < N*24*60) #minutes in N days
  {
    X_arr[idx] = X #savinvg current state
    #draw sojourn time until next event, s = min{T1, T2} where T1, T2 time until next birth and next death respectively, distributed like Exp(lambda + mu)
    if (X == 0)
    {
      s = simulate_exp(lambda)
      X = X+1
      
    }else{
      p = runif(1) #draw from uniform
      s = simulate_exp(lambda + mu)
      if (p <= lambda/(lambda + mu)) #prob of birth
      {
        #taking step
        X = X + 1
      } else {
        X = X - 1
      }
    }
    if(t+s > N*24*60) #checking if the soujorn time exeeds time limit of 50 days
    {
      s = N*24*60-t #if so only add the "allowed" time in state 
    }
    sojourn_times[idx] = s #save time spent in previous state (before the step was taken)
    idx = idx + 1
    t = t + s
  }
  list(X_arr, sojourn_times)
}

lam = 5/60 #per min
mu = 1/10 #per min
X0 = 0
N = 50
result = simulate_N_days(N, X0, lam, mu)
simulated_states = result[[1]]
times = result[[2]]

#Estimate of the long-run mean number of patients in the UCC: we use the last half of the simulated values for X in order to find the long-run mean
long_run_mean = function(N, X0, lambda, mu)
{
  simulations = simulate_N_days(N, X0, lam, mu)
  states = simulations[[1]]
  time_in_states = simulations[[2]]
  
  start = as.integer(length(states)/2) #index to start taking mean of
  #compute weighted mean (sum of xi*si / sum of si)
  weighted_states = states * time_in_states
  mean_num = sum(weighted_states[start:length(states)])/sum(time_in_states[start:length(states)])
  mean_num
}

L = long_run_mean(N, X0, lam, mu)
W = L/lam #Estimate of expected time a patient spends in the UCC (in minutes)


#repeat 30 times, check if normal distribution is a valid approximation (lrm = long run mean)
W_est_vals = c(1:30)*0
for(i in 1:30)
{
  X_lrm = long_run_mean(N, X0, lam, mu) # Estimate of L
  W_est = X_lrm/lam #little's law
  W_est_vals[i] = W_est
}
qqnorm(W_est_vals, main = "Normal Q-Q Plot for estimates of W (expected time a patient spends in the UCC)")
qqline(W_est_vals)
  
#compute average of the 30 lrm simulations 100 times to see if approximately normally distributed, if so: use CLT
W_est_averages = c(1:100)*0
for (n in 1:100)
{
  W_est_vals = c(1:30)*0
  for(i in 1:30)
  {
    X_lrm = long_run_mean(N, X0, lam, mu) # Estimate of L
    W_est = X_lrm/lam #little's law
    W_est_vals[i] = W_est
  }
  W_est_averages[n] = mean(W_est_vals)
}
qqnorm(W_est_averages, main = "Normal Q-Q Plot for average over 30 estimates of W")
qqline(W_est_averages) 
#so we can use CLT.


#CI:
#[W_bar-2.045*S_w/sqrt(30), W_bar + 2.045*S_w/sqrt(30)]
#compute W_bar
for(i in 1:30)
{
  X_lrm = long_run_mean(N, X0, lam, mu) # Estimate of L
  W_est = X_lrm/lam #little's law
  W_est_vals[i] = W_est
}
W_bar = mean(W_est_vals)

#compute S_w
S_w = 0
for (i in 1:30)
{
  S_w = S_w + (W_est_vals[i]-W_bar)^2
}
S_w = S_w / 29
S_w = sqrt(S_w)

#computing CI:
lower = W_bar - 2.045*S_w/sqrt(30)
upper = W_bar + 2.045*S_w/sqrt(30)
print(c("Approximate 95% CI for expected time a patient spends in the UCC (W): [", lower, ", ", upper, "]"))

#plotting one realization
realization = simulate_N_days(0.5, X0, lam, mu) #12 hours = 1/2 day
states = realization[[1]]
times = realization[[2]]
plot(NA, xlim=c(0,12*61), ylim=c(0,15), xlab="t (min)", ylab="X(t)", main = "One realization of {X(t):t≥0} for time 0-12 hours")
t = 0
for (i in 1:length(states)) {
  segments(t, states[i], t + times[i], states[i])
  t = t+times[i]
}
  
#verifying CI:

lam = 5/60 #per min
mu = 1/10 #per min
X0 = 0
N = 50
count = 0
W = 60
num = 1000

pb <- txtProgressBar(min = 1, max = num, style = 3)

for (n in 1:num) #running 1000 times to see if CI contains true W 95% of the times
{
  setTxtProgressBar(pb, n)
  for(i in 1:30)
  {
    X_lrm = long_run_mean(N, X0, lam, mu) # Estimate of L
    W_est = X_lrm/lam #little's law
    W_est_vals[i] = W_est
  }
  W_bar = mean(W_est_vals)
  
  #compute S_w
  S_w = 0
  for (i in 1:30)
  {
    S_w = S_w + (W_est_vals[i]-W_bar)^2
  }
  S_w = S_w / 29
  S_w = sqrt(S_w)
  
  #computing CI:
  lower = W_bar - 2.045*S_w/sqrt(30)
  upper = W_bar + 2.045*S_w/sqrt(30)
  
  if(W<= upper && W>= lower)
  {
    count = count + 1
  }
}
percentage = count/num * 100
print(c("True value W was contained in CI ", percentage, "% of the times."))

  
