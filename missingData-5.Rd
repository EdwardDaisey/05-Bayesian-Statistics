
missingData = function(nmbMcIterations)

{
  ##### get data
  library(MASS) ; data(Pima.tr)
  Y0<-Pima.tr[,2:5]
  Y<-Y0
  n<-dim(Y)[1]
  p<-dim(Y)[2]

  set.seed(1)
  O<-matrix(rbinom(n*p,1,.9),n,p)
  Y[O==0]<-NA


  #####
 


  print(Y);



  ############################


  ### prior parameters
  p<-dim(Y)[2]
  mu_0<-c(120,64,26,26)
  sd0<-(mu_0/2)
  Lambda_0<-matrix(.1,p,p) ; 
  diag(Lambda_0)<-1 ; 
  Lambda_0<-Lambda_0*outer(sd0,sd0)
  nu_0<-p+2 ; 
  Phi_0<-Lambda_0
  library(Rfast);
  ###

  ### starting values
  Sigma<-Phi_0
  Y.full<-Y
  O <- 1*(!is.na(Y));
  for(j in 1:p)
  {
    Y.full[is.na(Y.full[,j]),j]<-mean(Y.full[,j],na.rm=TRUE)
  }
  ###


  ### Gibbs sampler
  set.seed(1)
  results_of_mean = matrix(0, nmbMcIterations, 4);


  for(s in 1:nmbMcIterations)
  {

    if (s%%5 == 0) 
    {
	    print(paste("s = ", s));
	  }
  
    ###update mu
    ybar<-apply(Y.full,2,mean)
    Lambda_star <- solve( solve(Lambda_0) + n*solve(Sigma) )
    mu_star<-Lambda_star%*%( solve(Lambda_0)%*%mu_0 + n*solve(Sigma)%*%ybar )
    mu<-rmvnorm(1,mu_star,Lambda_star)
    ###
  
    ###update Sigma
    Sn<- Phi_0 + ( t(Y.full)-c(mu) )%*%t( t(Y.full)-c(mu) )
    Sigma<-solve( drop(rWishart(1, nu_0+n, solve(Sn)) ))
    ###
  
    ###update missing data
    for(i in 2:n)
    { 
      b <- ( O[i,]==0 )
      a <- ( O[i,]==1 )
      if (length(a[a == FALSE]) > 0)
      {
    	  iSa<- solve(Sigma[a,a])
    	  beta.j <- Sigma[b,a]%*%iSa
      	s2.j   <- Sigma[b,b] - Sigma[b,a]%*%iSa%*%Sigma[a,b]
      	mu.j<- mu[b] + beta.j%*%(t(Y.full[i,a])-mu[a])
      	Y.full[i,b] <- rmvnorm(1,mu.j,s2.j )
      }
    }
  
  
    ### save results
  
    results_of_mean[s,] = mu;
  }




  par(mfrow = c(2,2));
  hist(results_of_mean[,1], xlab = "", main = "Posterior of glu mean");
  hist(results_of_mean[,2], xlab = "", main = "Posterior of bp mean");
  hist(results_of_mean[,3], xlab = "", main = "Posterior of skin mean");
  hist(results_of_mean[,4], xlab = "", main = "Posterior of bmi mean");

}