cost_function <- function(y){
  n = length(y)
  #log(sum((y - mean(y))**2)/n)
  #var(y)
  m = 0
  if (n==0) {
    C = 0
  }
  else{
    C=n * (log(2*pi) + log(mean(((y-m)**2)))+ 1)}
  return(C)
}

penalty_function <- function(n, params){
  beta = params  *  log(n)
  return(beta)
}

penalty_function2 <- function(data){
  n = length(data)
  return (2*var(data)*log(n))
}

minimisation <- function(data, F_comp,beta){
  t_star = length(data)
  #print(t_star)
  F_temp = cost_function(data)
  #print(F_temp)
  i_temp = 0
  for (i in (2:(t_star))) {
    #print(i)
    F_tau = F_comp[i-1] + cost_function(data[i:t_star])  + beta

    if (F_tau < F_temp) {
      F_temp = F_tau
      i_temp = i-1}
    #print(F_temp)
  }
  return(list(F_tau_star = F_temp, tau_prime = i_temp))
}

OP <- function(data,params=1,K=0){
  #"if cost function to be minus the log-likelihood then the constant K = 0
  #if we take -penalised log-likelihood then K =  penalisation factor."
  n = length(data)
  beta = penalty_function(n,params)
  #beta = penalty_function2(data)
  #print(beta)
  cp <- rep(0,n)
  F_comp <- rep(0,n)
  F_comp[1] = -beta
  #cp = c(Inf)
  for (tau in c(2:n)) {
    #print(tau)
    res_minimisation <- minimisation(data = data[1:tau],
                                     beta = beta,
                                     F_comp = F_comp)
    #print(res_minimisation)
    F_comp[tau] <- res_minimisation$F_tau_star
    cp[tau] <- res_minimisation$tau_prime
    #print(F_comp)
    #print(cp)
  }

  v <- cp[n]
  P <- cp[n]

  while (v > 1)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(unique(P))
}


minimisation_PELT <- function(data, t_star ,F_comp,beta, R ,K = 0){
  t_star = length(data)
  #print(t_star)
  F_temp = cost_function(data)
  #print(F_temp)
  F_tau_K = c()
  i_temp = 0
  for (i in c(R)) {
    #print(i)
    F_tau = F_comp[i-1] + cost_function(data[i:t_star])  + beta

    F_tau_K = c(F_tau_K,F_comp[i-1] + cost_function(data[i:t_star])  + K)

    if (F_tau < F_temp) {
      F_temp = F_tau
      i_temp = i-1}
    #print(F_temp)
  }
  R_new <- R[which(F_tau_K < F_temp)]

  return(list(F_tau_star = F_temp, tau_prime = i_temp, R_new = R_new))
}


PELT <- function(data,params=1,K=0){
  #"if cost function to be minus the log-likelihood then the constant K = 0
  #if we take -penalised log-likelihood then K =  penalisation factor."
  n = length(data)
  beta = penalty_function(n,params)
  #print(beta)
  cp <- rep(0,n)
  F_comp <- rep(0,n)
  F_comp[1] = -beta
  R = c(2)
  #cp = c(Inf)
  for (tau in c(2:n)) {
    #print(tau)
    res_minimisation <- minimisation_PELT(data = data[1:tau],
                                          beta = beta,
                                          F_comp = F_comp,
                                          R = R,
                                          K = 0)
    #print(res_minimisation)
    F_comp[tau] <- res_minimisation$F_tau_star
    cp[tau] <- res_minimisation$tau_prime
    R <- c(res_minimisation$R_new, tau)

    #print(F_comp)
    print(R)
  }
  v <- cp[n]
  P <- cp[n]

  while (v > 1)
  {
    P <- c(P, cp[v])
    v <- cp[v]
  }
  P <- rev(P)[-1]

  return(unique(P))
}
