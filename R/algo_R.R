
#générateur de données

data_generator <- function(n, positions, variance_values){
  library(gfpop)
  data = dataGenerator(n, positions, variance_values, type = 'variance')
  return(data)
}

#définition de la métrique css
css_stat <- function(y){
  stat = c()
  n = length(y)
  for(t in 1:n){
    a = sum(y[1:t]**2)/sum(y**2)
    b = t/n
    stat <- c(stat, sqrt(n/2)*abs(a-b))
  }
  return(stat)
}

#définition de la fonction BS
BS_R<- function(s,e,y,pen=1.358){
  if (e-s==1) {
    return(c())
  } else {
    css_stats <- css_stat(y[s:e])
    cp <- which.max(css_stats) + s
    message("cp: ",cp)
    message("pen: ",pen)
    message("max(css_stats): ",max(css_stats))
    if (max(css_stats)> pen) {
      return(c(binSeg(s,cp,y,pen),cp,binSeg(cp+1,e,y,pen)))
    } else {
      return(c())
    }
  }
  
}



#fonction coût
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

#fonction de minimisation
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

#fonction Optimal partitioning
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


#minimisation PELT

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

#fonction PELT
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


cpt_naive<- function(x,w,l) {
  
  library(gfpop)
  library(zoo)
  
  ftest_fstat<- function (x) {
    y=split(x,cut(seq_along(x),2, labels = FALSE))
    x1=y$'1'
    x2=y$'2'
    f=var.test(x=x2,y=x1)
    
    pval=f$p.value
    ftest=f$statistic
    re=data.frame(ftest,pval)
    
    return(re)
    
  } 
  
  window_ftest<-rollapply(x, width = w, FUN =ftest_fstat)
  ftests<- as.vector(window_ftest[,1])
  pvalues<- as.vector(window_ftest[,2])
  fstat <- ftests[ftests>l]
  biggest_values=sort(fstat,decreasing=TRUE)[1:w]
  
  
  new_vec=vector(mode='numeric')
  
  for (i in 1:w) {
    
    indexz=c(which(fstat==biggest_values[i]))
    new_vec=c(new_vec,indexz)
  }
  new_vec
  
  new_vecx=sort(new_vec,decreasing = F)
  
  deletx<- function(x,w) {
    
    n=length(x)-1
    
    new_vecu=vector(mode='numeric')
    
    for (o in 1:n) {
      
      if(abs((x[o])-(x[o+1])) > (w*2)) {
        
        
        a=c(x[o],x[o+1])
        
        
        
        
        new_vecu=c(new_vecu,a)
        
      }
      
      
    }
    return(new_vecu)
  }
  
  t=deletx(new_vecx,w)
  
  taux=c(t,max(new_vec))
  tau=sort(taux,decreasing = F)
  
  get_vec <- function(x, gap) {
    grp <- cumsum(c(TRUE, diff(x) > gap))
    i1 <- x[!!!ave(x, grp, FUN = function(i) length(i) > 1)]
    i2 <- unname(tapply(x, grp, function(i) mean(i[length(i) > 1])))
    
    return(c(i1, i2[!is.na(i2)]))
    
    
  }
  
  ptch=as.integer(get_vec(tau,(w*2)))
  
  return(ptch)
  
}
