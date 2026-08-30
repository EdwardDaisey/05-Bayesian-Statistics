HierarchicalData = function(nmbMcIterations)

{

  #### Getting the data
  
  Y.school.mathscore<-dget("http://www2.stat.duke.edu/~pdh10/FCBS/Inline/Y.school.mathscore")
  
  
  ### weakly informative priors
  nu_0<-1  ; 
  sigma2_0<-100
  eta_0<-1 ; 
  tau2_0<-100
  theta_0<-50 ; 
  gamma2_0<-25
  ###
  
  ### starting values
  # m<-length(Y.school.mathscore) 
  m <- 100;
  n<-sv<-ybar<-rep(NA,m) 
  for(j in 1:m) 
  { 
    schoolSubset = Y.school.mathscore[Y.school.mathscore[,"school"] == j,];
    ybar[j]<-mean(schoolSubset[,"mathscore"])
    sv[j]<-var(Y.school.mathscore[,"mathscore"])
    n[j]<-length(Y.school.mathscore[,"mathscore"]) 
  }
  temp_mu<-ybar
  sigma2<-mean(sv)
  
  temp_theta<-mean(temp_mu)
  tau2<-var(temp_mu)
  ###
  
  
  
  ##### setting up the Gibbs sampler
  temp_mu_storage<-matrix( nrow=nmbMcIterations,ncol=m)
  MST<-matrix( nrow=nmbMcIterations,ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:nmbMcIterations) 
  {
  
  	print(s);
    # sample new values of the mus
    for(j in 1:m) 
    {
      vmu<-1/(n[j]/sigma2+1/tau2)
      emu<-vmu*(ybar[j]*n[j]/sigma2+temp_theta/tau2)
      temp_mu[j]<-rnorm(1,emu,sqrt(vmu))
     }
  
    #sample new value of sigma2
     nun<-nu_0+sum(n)
     ss<-nu_0*sigma2_0;
     for(j in 1:m){
     		ss<-ss+sum((Y.school.mathscore[[j]]-temp_mu[j])^2)
     	}
     sigma2<-1/rgamma(1,nun/2,ss/2)
  
    #sample a new value of theta
     vtheta <- 1/(m/tau2+1/gamma2_0)
     etheta <- vtheta*(m*mean(temp_mu)/tau2 + theta_0/gamma2_0)
     temp_theta<-rnorm(1,etheta,sqrt(vtheta)) 
  
    # sample a new value of tau2
     etam<-eta_0+m
     ss<- eta_0*tau2_0 + sum( (temp_mu-temp_theta)^2 )
     tau2<-1/rgamma(1,etam/2,ss/2)
  
    #store results
     temp_mu_storage[s,]<-temp_mu
     MST[s,]<-c(temp_theta,sigma2,tau2)
  
  } 
  ###
  
  
  
  
  par(mfrow = c(2,2));
  hist(temp_mu_storage[,1], breaks = 50, xlab = "", main = "posterior of mu_1");
  hist(MST[,1], breaks = 50, xlab = "", main = "posterior of theta");
  hist(MST[,2], breaks = 50, xlab = "", main = "posterior of sigma2");
  hist(MST[,3], breaks = 50, xlab = "", main = "posterior of tau2");

}