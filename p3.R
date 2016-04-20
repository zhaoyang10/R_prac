loadPath = "E:\\datasets\\time-series\\R_prac\\"
dataset = "divorces-1855-2011.xls"
data = paste(loadPath, dataset, sep = "")
library(xlsx)
w <- read.xlsx(data, 1)
x = ts(w[,-1], start = c(1855), freq = 1)
plot(x, ylab = "Number of Marriages and Divorces", main = "", type = "o", pch = 16, plot.type = "single", lty=1:2)
title("Marriage and Divorce in Scotland")
legend("topleft", c("Marriage", "Divorce"), lty = 1:2)