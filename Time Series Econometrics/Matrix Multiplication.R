D <- matrix(c(0,0.001,12.5,0,0,21.1,0,0,0), ncol=3)
A <- diag(3)-D
B <- matrix(c(0.35,0.004,10.6,-0.18,0.40,15.5,-0.0003,0.0003,0.98), ncol=3)
E <- matrix(c(0,0,1))

irf0 <- solve(A) %*% E
irf0

irf1 <- solve(A) %*% B %*% irf0
irf1

irf2 <- solve(A) %*% B %*% irf1
irf2


