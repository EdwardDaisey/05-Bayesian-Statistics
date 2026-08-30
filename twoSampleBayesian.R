twoSampleEqualVar = function(m, n, lambda, mu, phi, smplPost) {

x = rnorm(m, lambda, sqrt(phi))
y = rnorm(n, mu, sqrt(phi))
print(paste("sample mean = ", mean(x)))
print(paste("y mean = ", mean(y)))
smplVar = ((m-1)*var(x) + (n-1)*var(y))/(m + n - 2)
par(mfrow = c(3,2))
hist(x, density = 50, breaks = 50)
hist(y, density = 50, breaks = 50)
postPhi = smplVar*(n + m -2)*rinvchisq(smplPost, (n-1 + m - 1))
postLambda = rnorm(smplPost, mean(x), sqrt(postPhi/m))
postMu = rnorm(smplPost, mean(y), sqrt(postPhi/m))
hist(postPhi, density = 50, breaks = 50)
hist(postLambda, density = 50, breaks = 50)
hist(postMu, density = 50, breaks = 50)

}





twoSampleEqualVarInformPrior = function(m, n, lambda, mu, m0, n0, lambda0, mu0,  phi, v0, S0, smplPost) {

x = rnorm(m, lambda, sqrt(phi))
y = rnorm(n, mu, sqrt(phi))
print(paste("y mean = ", mean(y)))
print(paste("x mean = ", mean(x)))
n1 = n + n0
m1 = m + m0
smplVar = ((m-1)*var(x) + (n-1)*var(y))/(m + n - 2)
par(mfrow = c(3,2))
hist(x, density = 50, breaks = 50)
hist(y, density = 50, breaks = 50)
S1 = S0 + var(x)*(m-1) + var(y)*(n-1) + (1/m0 + 1/m)^(-1)*(mean(x) - lambda0)^2 + (1/n0 + 1/n)^(-1)*(mean(y) - mu0)^2
v1 = v0 + m + n
s1 = S1/v1
postPhi = S1*rinvchisq(smplPost, v1)
lambda1 = (m0*lambda0 + m*mean(x))/m1;
mu1 = (n0*mu0 + n*mean(y))/n1;
postLambda = rnorm(smplPost, lambda1, sqrt(postPhi/m1))
postMu = rnorm(smplPost, mu1, sqrt(postPhi/n1))
diff_T = rt(smplPost, v1)
postMeanDiff = diff_T*(s1*(1/m1 + 1/n1)^(1/2)) + (mean(x) - mean(y))
hist(postPhi, density = 50, breaks = 50)
hist(postLambda, density = 50, breaks = 50)
hist(postMu, density = 50, breaks = 50)
hist(postMeanDiff, density = 50, breaks = 50)

}








twoSampleUnequalVar = function(m, n, lambda, mu, m0, n0, lambda0, mu0,  phi, v0, S0, smplPost) {

x = rnorm(m, lambda, sqrt(phi))
y = rnorm(n, mu, sqrt(phi))
print(paste("y mean = ", mean(y)))
print(paste("x mean = ", mean(x)))
n1 = n + n0
m1 = m + m0
smplVar = ((m-1)*var(x) + (n-1)*var(y))/(m + n - 2)
par(mfrow = c(3,2))
hist(x, density = 50, breaks = 50)
hist(y, density = 50, breaks = 50)
sx = sd(x)
sy = sd(y)
theta = atan((sx/sqrt(m))/(sy/sqrt(n)))
vx = m-1
vy = n-1
diff_BF = rt(smplPost, vx)*sin(theta) - rt(smplPost, vy)*cos(theta)
postMeanDiff = diff_BF*(sqrt(sx^2/m + sy^2/n)) + (mean(x) - mean(y))

hist(postMeanDiff, density = 50, breaks = 50)


}