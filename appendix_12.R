nrow(x);ncol(x);dim(x)
x = matrix(rnorm(24), 4, 6); x
x[c(2, 1), ]
x[, c(1, 3)]
x[2, 1]
x[x[, 1]>0, 1]
x[, 1]>0
sum(x[,1]>0)
sum(x[,1]<=0)
x[,-c(1,3)]
diag(x)
diag(1:5)
diag(5)
diag(rep(5,5))
x[-2,-c(1,3)]
x[x[,1]>0&x[,3]<=1,1]
x[x[,2]>0|x[,1]<.51, 1]
apply(x, 1, mean)
apply(x, 2, sum)
x = matrix(rnorm(24), 4, 6); x
x[lower.tri(x)] = 0; x
x[upper.tri(x)] = 0; x
