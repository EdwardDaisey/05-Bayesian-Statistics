multivariateNormalBayes = function(nmbMcIterations)

{

	Y = structure(c(59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49, 72, 34, 
	+ 70, 34, 50, 41, 52, 60, 34, 28, 35, 77, 39, 46, 26, 38, 43, 68, 
	+ 86, 77, 60, 50, 59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33), .Dim = c(22L, 
	+ 2L), .Dimnames = list(NULL, c("pretest", "posttest")));


	library(mvtnorm);
	mu0<-c(50,50)
	L0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

	nu0<-4
	S0<-matrix( c(625,312.5,312.5,625),nrow=2,ncol=2)

	n<-dim(Y)[1] ; ybar<-apply(Y,2,mean)
	Sigma<-cov(Y) ; MEAN_POST<-SIGMA<-NULL
	set.seed(1)

	#### now looping through the MX Iterations

	for(s in 1:nmbMcIterations) 
	{
		print(paste("MC Iteration #", s));
  
  		###update mean
  		Ln<-solve( solve(L0) + n*solve(Sigma) )
  		mun<-Ln%*%( solve(L0)%*%mu0 + n*solve(Sigma)%*%ybar )
  		meanPost<-rmvnorm(1,mun,Ln)  
  		### 
   
	  	###update Sigma
  		Sn<- S0 + ( t(Y)-c(mun) )%*%t( t(Y)-c(mun) ) 
  		Sigma<-solve( drop(rWishart(1, nu0+n, solve(Sn))) )


  		### save/store results 
  		MEAN_POST<-rbind(MEAN_POST,meanPost) ; SIGMA<-rbind(SIGMA,c(Sigma))
  
	}


plot(MEAN_POST[,1], MEAN_POST[,2])

}
