x = sample(1:100, 12);x
all(x > 0); all(x!=0);any(x>0);(1:10)[x>0]
diff(x)
x
diff(x, lag = 2)
x = matrix(1 : 20, 4, 5);x
x = matrix(1 : 20, 4, 5, byrow = T);x
t(x)
x = matrix(sample(1:100, 20), 4, 5);x
(2*x)
(x + 5)
y = matrix(sample(1:100, 20), 5, 4);y
(x + t(y))
(z = x%*%y)
z1 = solve(z);z1
z1%*%z
round(z1%*%z, 13)
b = solve(z, 1:4);b
