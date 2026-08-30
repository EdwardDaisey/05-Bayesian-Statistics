#####    Metropolis-Hasting Algorithm in R

metHastings = function(n, sigma) {
	
	x = rep(0,n)
	j = 1
	acc = 0
	x[1] = 2;
	while (j < n) {
		
		######.  Sample candidate value
		x_cand = rnorm(1,x,sigma^2)
		
		######  Calculate acceptance probability
		acc_prob = acc_prob_calc(x_cand, x[j], sigma)
		
		
		######. Decide whether to accept or not
		acc_or_not = rbinom(1,1,acc_prob);
		
		
		#####.  if accept or not
		
		if (acc_or_not == 1) {
			acc = acc + 1;
			j = j + 1;
			x[j] = x_cand
		}
		else {
			j = j + 1
			x[j] = x[j-1]
		}	
	}
	acc_probability = acc/n;
	print(acc_probability);
	hist(x, breaks = 50, density = 50);
	
}




acc_prob_calc = function(candidate, former, sigma) {
	
	if (candidate < 0) {
		acc_probability = 0;
	}
	else {
		acc_probability = dens(candidate)*dnorm(former, candidate, sigma^2)/(dens(former)*dnorm(candidate, former, sigma^2));
	}
	acc_probability = min(1, acc_probability);
	return(acc_probability);
}






######## density is the noncentral chi-squared distribution

dens = function(x) {
	lambda = 3;
	p = 5;
	dens_evaluated_at_x = (1/2)*(x/lambda)^((p-2)/4)*(besselI(sqrt(3*x), (p-2)/2))*exp(-(p+x)/2)
	return(dens_evaluated_at_x);
}
