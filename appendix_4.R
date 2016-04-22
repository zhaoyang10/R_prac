ncol(cars)
nrow(cars)
cars
lm(dist~speed, data=cars)
cars$qspeed = cut(cars$speed, breaks = quantile(cars$speed), include.lowest = TRUE)
names(cars)
