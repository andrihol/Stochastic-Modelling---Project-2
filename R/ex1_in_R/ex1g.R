rm(list = ls())

simulate_exp = function(rate) #inversion method
{
  u = runif(1)
  #find corresponding x
  x = - log(1-u)/rate
  x
}

#X0 = 0 -> N0 = U0 = 0
simulate_UN_num_days = function(mu, lambda, p, num)
{
  realisations = matrix(c(1:2*100000)*0, nrow=2, ncol=100000) #2x100000 matrix of U(t) and N(t) realisations
  t = 0
  sojourn_times = c(1:100000)*0 #element j = time stayed in state realisations[i,j] before a new event occured
  N = 0
  U = 0
  idx = 1
  while(t < num*24*60)
  {
    realisations[1, idx] = U
    realisations[2, idx] = N
    if(N+U == 0) #no patients in system
    {
      s = simulate_exp(lambda) #only arrivals can occur
      u = runif(1) #decide whether U or N arrival
      if (u<=p) #urgent arrival
      {
        U = U+1
      }
      else #normal arrival
      {
        N = N+1 
      }
    }
    else if(U == 0) #N>0, normal patient in treatment
    {
      s = simulate_exp(lambda + mu) #time until first event, both births and deaths can occur
      u1 = runif(1) #to decide whether birth or death
      if (u1 <= lambda/(lambda + mu)) #arrival is first event
      {
        u2 = runif(1) #to decide whether U or N arrival
        if (u2<=p) #urgent arrival
        {
          U = U+1 #goes into treatment
        }
        else #normal arrival
        {
          N = N+1
        }
      }
      else #death first event
      {
        N = N-1 #only normal patients in UCC, so a normal patient must be the one to leave
      }
    }
    else if(N == 0) #U > 0, so an urgent patient is in treatment
    {
      s = simulate_exp(lambda + mu) #time until first event, both births and deaths can occur
      u1 = runif(1) #decide whether birth or death
      if (u1<= lambda/(lambda + mu)) #arrival occurs first
      {
        u2 = runif(1) #to decide whether U or N arrival
        if (u2<=p) #urgent arrival
        {
          U = U+1
        }
        else #normal arrival
        {
          N = N+1
        }
      }
      else #death occurs first
      {
        #since U>0 the person in treatment until now is urgent and therefore the person leaving must be urgent
        U = U-1
      }
    }
    else #both U, N >0, and since U>0 the patient in treatment is urgent
    {
      s = simulate_exp(lambda + mu) #time until first event, both births and deaths can occur
      u1 = runif(1) #decide whether birth or death is next event
      if (u1<=lambda/(lambda + mu)) #arrival occurs
      {
        u2 = runif(1) #decide whether U or N arrival
        if (u2<=p) #urgent arrival
        {
          U = U+1
        }
        else #normal arrival
        {
          N = N+1
        }
      }
      else #death
      {
        #U>0 so the person leaving must be urgent
        U = U-1
      }
    }
    
    #update stored values
    if(t+s > num*24*60) #checking if the soujorn time exeeds time limit of 50 days
    {
      s = num*24*60-t #if so only add the "allowed" time in state 
    }
    sojourn_times[idx] = s #save time spent in previous state (before the step was taken)
    idx = idx + 1
    t = t + s
  }
  #return results
  list(realisations[, 1:idx], sojourn_times[1:idx]) #remove unused spaces
}


#estimate Wn and Wu
lam = 5/60 #per min
mu = 1/10 #per min
num_days = 50
p = 0.80


#DEBUG
simulations = simulate_UN_num_days(mu, lam, p, num_days)
sims = simulations[[1]]
tims = simulations[[2]]

#Estimate of the long-run mean number of U and N patients in the UCC: we use the last half of the simulated values in order to find the long-run mean
long_run_mean = function(mu, lambda, p, num)
{
  simulations = simulate_UN_num_days(mu, lambda, p, num)
  states = simulations[[1]] #2D array for U (row 1) and N (row 2)
  time_in_states = simulations[[2]]
  
  start = as.integer(length(time_in_states)/2) #index to start taking mean of
  #compute weighted mean (sum of xi*si / sum of si)
  means = c(0,0)
  for (i in 1:2)
  {
    weighted_states = states[i,] * time_in_states #element-wise multiplication
    lr_mean = sum(weighted_states[start:length(time_in_states)])/sum(time_in_states[start:length(time_in_states)])
    means[i] = lr_mean
  }
  means
}

lrms = long_run_mean(mu, lam, p, num_days) #long run means from one run of 50 days
L_U = lrms[1]
L_N = lrms[2]
W_U = L_U/(lam*p)
W_N = L_N/((1-p)*lam)
cat("Estimated W_U = ", W_U , "\nEstimated W_N = ",W_N,  "\n(based on one simulation of 50 days)")

#repeat 30 times:
#store W_U and W_N values
WU_vals = c(1:30)*0
WN_vals = c(1:30)*0
for (i in 1:30)
{
  L = long_run_mean(mu, lam, p, num_days) #long run means from one run of 50 days
  WU_vals[i] = L[1]/(lam*p)
  WN_vals[i] = L[2]/((1-p)*lam)
}

WU_bar = mean(WU_vals)
WN_bar = mean(WN_vals)

#compute s_U and s_N:
s_U = 0
s_N = 0
for(i in 1:30)
{
  s_U = s_U + (WU_bar - WU_vals[i])^2
  s_N = s_N + (WN_bar - WN_vals[i])^2
}
s_U = sqrt(s_U/29)
s_N = sqrt(s_N/29)

lower_U = WU_bar - 2.045*s_U/sqrt(30)
upper_U = WU_bar + 2.045*s_U/sqrt(30)
lower_N = WN_bar - 2.045*s_N/sqrt(30)
upper_N = WN_bar + 2.045*s_N/sqrt(30)
#CI (approximation through normal approximation):
cat("95% CI for W_U: [", lower_U, ", ", upper_U, "]\n")
cat("95% CI for W_N: [", lower_N, ", ", upper_N, "]\n")


#verifying approximately 95% CI:
true_WU = 30 #min
true_WN = 180 #min
num_runs = 1000
count_U = 0
count_N = 0
pb <- txtProgressBar(min = 1, max = num_runs, style = 3)
for (m in 1:num_runs)
{
  setTxtProgressBar(pb, m)
  
  WU_vals = c(1:30)*0
  WN_vals = c(1:30)*0
  for (i in 1:30)
  {
    L = long_run_mean(mu, lam, p, num_days) #long run means from one run of 50 days
    WU_vals[i] = L[1]/(lam*p)
    WN_vals[i] = L[2]/((1-p)*lam)
  }
  
  WU_bar = mean(WU_vals)
  WN_bar = mean(WN_vals)
  
  #compute s_U and s_N:
  s_U = 0
  s_N = 0
  for(i in 1:30)
  {
    s_U = s_U + (WU_bar - WU_vals[i])^2
    s_N = s_N + (WN_bar - WN_vals[i])^2
  }
  s_U = sqrt(s_U/29)
  s_N = sqrt(s_N/29)
  
  lower_U = WU_bar - 2.045*s_U/sqrt(30)
  upper_U = WU_bar + 2.045*s_U/sqrt(30)
  lower_N = WN_bar - 2.045*s_N/sqrt(30)
  upper_N = WN_bar + 2.045*s_N/sqrt(30)
  if (true_WU >= lower_U && true_WU <= upper_U)
  {
    count_U = count_U + 1
  }
  if (true_WN >= lower_N && true_WN <= upper_N)
  {
    count_N = count_N + 1
  }
  
}
frac_U = count_U/num_runs
frac_N = count_N/num_runs
cat("Fraction of times where the true value W_U was contained in CI for W_U: ", frac_U, "\nFraction of times where the true value W_N was contained in CI for W_N: ", frac_N)
  
#plotting one joint realization of X(t) = U(t) + N(t) and U(t), for time 0-12hrs.
#plotting one realization
realization = simulate_UN_num_days(mu, lam, p, 0.5) #12 hours = 1/2 day
states = realization[[1]]
states_U = states[1,]
states_N = states[2,]
states_X = states_U + states_N
times = realization[[2]]
plot(NA, xlim=c(0,12*61), ylim=c(0,15), xlab="t (min)", ylab="number of patients", main = "One joint realization of {X(t):t≥0} and {U(t):t≥0} for time 0-12 hours")
t = 0
for (i in 1:length(times)) {
  segments(t, states_U[i], t + times[i], states_U[i], col = "red")
  segments(t, states_X[i], t + times[i], states_X[i], col = "black")
  t = t+times[i]
}
legend('topleft', legend = c('X(t)', 'U(t)'), col = c('black', 'red'), lty = 1)


