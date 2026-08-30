  gibbs.smpl.tnorm = function(n, theta1_0, theta2_0, theta3_0, k, x_1, x_2, x_3, phi) {
     smpl = matrix(0, 10000, 3)
     smpl[1,] = t(c(theta1_0, theta2_0, theta3_0))
     for (j in 2:n) { 
        # print(paste("j = ", j))
        smpl[j,1] = rtnorm(1,x_1,phi,0, smpl[(j-1),2])  
        smpl[j,2] = rtnorm(1,x_2,phi,smpl[(j),1], smpl[(j-1),3])
        smpl[j,3] = rtnorm(1,x_3,phi,smpl[(j),2],k)
     } 

     par(mfrow = c(2,2))
     plot(smpl[,1], smpl[,2], xlab = "theta_1", ylab = "theta_2")
     plot(smpl[,1], smpl[,3], xlab = "theta_1", ylab = "theta_3")
     plot(smpl[,2], smpl[,3], xlab = "theta_2", ylab = "theta_3")
}
