x = 1 : 100
(x = 1:100)
set.seed(0)
sample(1:10, 3)
z=sample(1:200000, 10000)
z[1:10]
y=c(1, 3, 7, 3, 4, 2)
z[y]
(z = sample(x, 100, rep = T))
(z1 = unique(z))
length(z1)
xz=setdiff(x, z)
sort(union(xz, z))
setequal(union(xz, z), x)
intersect(1:10, 7:50)
sample(1:100, 20, prob = 1 : 100)
