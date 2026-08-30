treeRegression = function(n) {
	
	library(invgamma);
	library(MASS);
	beta0 = c(0,0,0);
	Sigma0 = 20*diag(3);
	nu0 = 5;
	sigma20 = 1;
	y = trees$Girth;
	X = cbind(rep(1,31), trees$Height, trees$Volume);
	sigma2 = sigma20;
	beta = beta0;
	betaSample = matrix(rep(0,3*n), 3,n);
	sigma2Sample = rep(0,n);
	
	for (j in 1:n) {
		
		Sigma_post = solve(solve(Sigma0) + t(X) %*% X/sigma2);
		mu_post = Sigma_post %*% (solve(Sigma0) %*% beta0 + t(X) %*% y/sigma2);
		sampledBeta = mvrnorm(1, mu_post, Sigma_post);
		betaSample[,j] = t(sampledBeta);
		
		
		SSR = t(y - X %*% sampledBeta) %*% (y - X%*% sampledBeta);
		postShape = (nu0 + 31)/2;
		postRate = (nu0*sigma20 + SSR)/2;
		sigma2 = rinvgamma(1,postShape, postRate);
		sigma2Sample[j] = sigma2;
		
		
		
	}
	
	
	par(mfrow= c(1,3));
	plot(betaSample[1,], betaSample[3,]);
	plot(betaSample[2,], betaSample[3,]);
	hist(sigma2Sample, breaks = 50, density = 50);
	
	
}
