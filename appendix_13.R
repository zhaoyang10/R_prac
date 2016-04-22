x = array(runif(24), c(4, 3, 2));x
is.matrix(x)
dim(x)
x = array(1:24, c(4, 3, 2));x
x[c(1,3), ,]
x = array(1 : 24, c(4, 3, 2));x 
apply(x, 1, mean)
apply(x, 1:2, sum)
apply(x, c(1, 3), prod)
