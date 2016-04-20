loadPath = "E:\\datasets\\time-series\\R_prac\\TimeSeriesDatasets_130207\\"
dataset = "NZRainfall.csv"
data = paste(loadPath, dataset, sep="")
w = read.csv(data)
x = ts(w[,-1], start = c(2000,1), freq = 12)
plot(x, ylab = "Monthly rainfall (mm)", type = "o", pch = 16, nc = 1, main = "")
title("New Zealand Rainfall")