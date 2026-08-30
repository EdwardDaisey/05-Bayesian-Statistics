betaBinom_Bayesian = function(n, p, alpha, beta) {

	smpl = rbinom(1, n, p);
	print(paste("number of successes = ", smpl))
	xlimits = seq(0,1,.01);
	prior = dbeta(xlimits, alpha, beta);
	post_alpha = alpha + smpl;
	post_beta = beta + n - smpl;
	post = dbeta(xlimits, post_alpha, post_beta);
	plot(xlimits, prior, type = "l", ylim = c(0, max(post)))
	lines(xlimits, post, type = "l", col = "red")
}
