x = scan()
#1.5 2.6 3.7 2.1 8.9 12 -1.2 -4
x = c(1.5, 2.6, 3.7, 2.1, 8.9, 12, -1.2, -4)
w = read.table(file.choose(), header = T)

setwd("E:\\Code\\github\\R_prac")
(x = rnorm(20))
write(x, "test.txt")
y = scan("test.txt")
y
y = iris
y[1:5, ]
str(y)
write.table(y, "test.txt", row.names = F)
w = read.table("test.txt", header = T)
str(w)
write.csv(y, "test.csv")
v = read.csv("test.csv")
str(v)
data = read.table("clipboard")
