accReject = function(n) {
	
	
	simulatedValues = rep(0,n)
	simulationSuccesses = 0
	simulationTrials = 0
	M = 2;
	lambda = 1;
	
	while (simulationSuccesses < n) {
		
		simulationTrials = simulationTrials + 1
		
		##### generate candidate values of X
	
			x_cand = rexp(1,lambda)
	
		###### Generate Uniform
	
			unifRnd = runif(1,0, M*dexp(x_cand, lambda))
	
		##### Decide whether to accept or not
	
			if (unifRnd <  densToSampleFrom(x_cand))	{
				simulationSuccesses = simulationSuccesses + 1
				simulatedValues[simulationSuccesses] = x_cand
			}
	}
	
	accProb = simulationSuccesses/simulationTrials
	print(paste("acceptance probability = ", accProb))
	hist(simulatedValues, density = 50, breaks = 50)
}



densToSampleFrom = function(xcand) {
	
	densEval = xcand*2*exp(-xcand^2)
	return(densEval)
	
}
