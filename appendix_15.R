airquality
complete.cases(airquality)
which(complete.cases(airquality) == FALSE)
sum(complete.cases(airquality))
na.omit(airquality)
x = 1 : 10; x[12] = 3
x1 = append(x, 77, after = 5); x1
cbind(1:5, rnorm(5))
rbind(1:5, rnorm(5))
x = rbind(1 : 5, runif(5), runif(5), 1 : 5, 7 : 11); x
x[!duplicated(x),]
unique(x)
