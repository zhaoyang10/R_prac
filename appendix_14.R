x = matrix(1:20, 5, 4);x
sweep(x, 1, 1:5, "*")
sweep(x, 1, 1:4, "+")
x*1:5
(x = matrix(sample(1:100, 24), 6, 4))
(x1 = scale(x))
x2 = scale(x, scale = F); x2
x3 = scale(x, scale = T, center = F); x3
round(apply(x1, 2, mean), 14)
apply(x1, 2, sd)
round(apply(x2, 2, mean), 14)
apply(x2, 2, sd)
round(apply(x3, 2, mean), 14)
apply(x3, 2, sd)
round(apply(x, 2, mean), 14)
apply(x, 2, sd)
