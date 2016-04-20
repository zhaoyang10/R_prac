library(MASS)
library(cluster)
.search.vanilla <- search()
detach.everything <- function () {
  for (i in setdiff( search(), .search.vanilla )) {
    detach(pos = match(i, search()))
  }
}
detach.everything()

png(file="g1.png", width=600, height=600)
  library(Zelig)
  d <- data.frame(
    y = crabs$sp,
    x1 = crabs$FL,
    x2 = crabs$RW
  )
  r <- zelig( y ~ x1 + x2, model="probit", data=d )
  summary(r)
  op <- par(mfrow=c(2,2), mar=c(4,4,3,2))
  plot(r)
  par(op)
dev.off()

png(file="g2.png", width=600, height=600)
  x <- rnorm(100)
  hist(x, col = "light blue")
dev.off()

png(file="g3.png", width=600, height=600)
  N <- 100
  x <- rnorm(N)
  y <- x + rnorm(N)
  plot(y ~ x)
dev.off()

png(file="g4.png", width=600, height=600)
  N <- 100
  x <- rnorm(N)
  y <- x + rnorm(N)
  plot(y ~ x)
  abline( lm(y ~ x), col = "red" )
dev.off()
detach.everything()

png(file="g5.png", width=600, height=600)
  data(EuStockMarkets)
  x <- EuStockMarkets
  # We aren't interested in the spot prices, but in the returns
  # return[i] = ( price[i] - price[i-1] ) / price[i-1]
  y <- apply(x, 2, function (x) { diff(x)/x[-length(x)] })
  # We normalize the data
  z <- apply(y, 2, function (x) { (x-mean(x))/sd(x) })
  # A single time series
  r <- z[,1]
  # The runs
  f <- factor(cumsum(abs(diff(sign(r))))/2)
  r <- r[-1]
  accumulated.return <- tapply(r, f, sum)
  trend.number <- table(f)
  boxplot(abs(accumulated.return) ~ trend.number, col='pink',
          main="Accumulated return")
dev.off()

png(file="g6.png", width=600, height=600)
  boxplot(abs(accumulated.return)/trend.number ~ trend.number,
          col='pink', main="Average return")
dev.off()

png(file="g7.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    r <- z[,i]
    f <- factor(cumsum(abs(diff(sign(r))))/2)
    r <- r[-1]
    accumulated.return <- tapply(r, f, sum)
    trend.number <- table(f)
    boxplot(abs(accumulated.return) ~ trend.number, col='pink')
  }
  par(op)    
dev.off()

png(file="g8.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    r <- z[,i]
    f <- factor(cumsum(abs(diff(sign(r))))/2)
    r <- r[-1]
    accumulated.return <- tapply(r, f, sum)
    trend.number <- table(f)
    boxplot(abs(accumulated.return)/trend.number ~ trend.number, col='pink')
  }
  par(op)    
dev.off()

png(file="g9.png", width=600, height=600)
  data(volcano)
  M <- volcano
  n <- dim(M)[1]
  m <- dim(M)[2]
  M1 <- M [1:(n-1),] [,1:(m-1)]
  M2 <- M [2:n,] [,1:(m-1)]
  M3 <- M [1:(n-1),] [,2:m]
  M4 <- M [2:n,] [,2:m]
  # Overlapping quadripixels
  M0 <- (M1+M2+M3+M4)/4

  # Non-overlapping quadripixels
  nn <- floor((n-1)/2)
  mm <- floor((m-1)/2)
  M00 <- M0 [2*(1:nn),] [,2*(1:mm)]

  op <- par(mfrow=c(2,2))
  image(M, main="Initial image")
  image(M0, main="Overlapping Quadripixels")
  image(M00, main="Non Overlapping Quadripixels")
  par(op)
dev.off()

png(file="g10.png", width=600, height=600)
  n <- 100
  m <- matrix(runif(2*n),nc=2)
  library(ape)
  r <- mst(dist(m)) # The incidence matrix (of the minimum spanning
                    # tree of the points)
  plot(m)
  n <- dim(r)[1]
  w <- which(r!=0)
  i <- as.vector(row(r))[w]
  j <- as.vector(col(r))[w]
  segments( m[i,1], m[i,2], m[j,1], m[j,2], col='red' )
dev.off()

png(file="g11.png", width=600, height=600)
  several.times <- function (n, f, ...) {
    for (i in 1:n) {
      f(...)
    }
  }
  matrix.multiplication <- function (s) {
    A <- matrix(1:(s*s), nr=s, nc=s) 
    B <- matrix(1:(s*s), nr=s, nc=s) 
    C <- A %*% B
  }
  v <- NULL
  for (i in 2:10) {
    v <- append(
      v, 
      system.time( 
        several.times( 
          10000, 
          matrix.multiplication, 
          i 
        )
      ) [1]
    )
  }
  plot(v, type = 'b', pch = 15, 
       main = "Matrix product computation time")
dev.off()

png(file="g12.png", width=600, height=600)
  x <- seq(0,4, length=100)
  y <- sqrt(x)
  plot(y~x, type='l', lwd=3, main=expression(y == sqrt(x)) )
dev.off()
detach.everything()

png(file="g13.png", width=600, height=600)
  n <- 1000
  k <- 20
  p <- 3
  i <- sample(1:p, n, replace=TRUE)
  x <- 10 * matrix(rnorm(p*k), nr=p, nc=k)
  x <- x[i,] + matrix(rnorm(n*k), nr=n, nc=k)
  L1L2 <- function (x) {
    cbind(L1 = apply(x, 1, mean),
          L2 = apply(x, 1, sd))
  }
  plot(L1L2(x), col=i)
dev.off()

png(file="g14.png", width=600, height=600)
  x <- crabs$FL
  y <- crabs$CL  # The two vectors need not 
                 # have the same length
  op <- par(mfrow=c(2,1))
  hist(x, col="light blue", xlim=c(0,50))
  hist(y, col="light blue", xlim=c(0,50))
  par(op)
dev.off()

png(file="g15.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  hist( (x - mean(x)) / sd(x), 
        col = "light blue",
        xlim = c(-3, 3) )
  hist( (y - mean(y)) / sd(y), 
        col = "light blue",
        xlim = c(-3, 3) )
  par(op)
dev.off()

png(file="g16.png", width=600, height=600)
  N <- 50           # Sample size
  set.seed(2)
  x1 <- runif(N)    # Uniform distribution
  x2 <- rt(N,2)     # Fat-tailed distribution
  x3 <- rexp(N)     # Skewed distribution
  x4 <- c(x2,20)    # Outlier (not that uncommon,
                    # with fat-tailed distributions)
  f <- function (x, ...) {
    x <- (x - mean(x)) / sd(x)
    N <- length(x)
    hist( x,
          col = "light blue",
          xlim = c(-3, 3),
          ylim = c(0, .8),
          probability = TRUE,
          ...
        )
    lines(density(x), 
          col = "red", lwd = 3)
    rug(x)
  }
  op <- par(mfrow=c(2,2))
  f(x1, main = "Uniform distribution")
  f(x2, main = "Fat-tailed distribution")
  f(x3, main = "Skewed distribution")
  f(x4, main = "As above, with one outlier")  
  par(op)
dev.off()

png(file="g17.png", width=600, height=600)
  x <- read.csv("2006-08-27_pe.csv")
  op <- par(mfrow=c(1,2))
  plot(p ~ eps, data=x, main="Before")
  plot(log(p) ~ log(eps), data=x, main="After")
  par(op)
dev.off()

png(file="g18.png", width=600, height=600)
  f <- function (x, main, FUN) { 
    hist(x, 
         col = "light blue", 
         probability = TRUE, 
         main = paste(main, "(before)"),
         xlab = "")
    lines(density(x), col = "red", lwd = 3)
    rug(x)
    x <- FUN(x)
    hist(x, 
         col = "light blue", 
         probability = TRUE, 
         main = paste(main, "(after)"),
         xlab = "")
    lines(density(x), col = "red", lwd = 3)
    rug(x)
  }
  op <- par(mfrow=c(2,2))
  f(x3, 
    main="Skewed distribution", 
    FUN = log)
  f(x2, 
    main="Fat tailed distribution", 
    FUN = function (x) {  # If you have an idea of the 
                          # distribution followed by 
                          # your variable, you can use 
                          # that distribution to get a 
                          # p-value (i.e., a number between
                          # 0 and 1: just apply the inverse
                          # of the cumulative distribution 
                          # function -- in R, it is called
                          # the p-function of the 
                          # distribution) then apply the 
                          # gaussian cumulative distribution 
                          # function (in R, it is called the 
                          # quantile function or the 
                          # q-function).
            qnorm(pcauchy(x)) 
          } 
    )
  par(op)
dev.off()

png(file="g19.png", width=600, height=800)
  uniformize <- function (x) { # This could be called 
                               # "forceful uniformization".
                               # More about it when we introduce
                               # the notion of copula.
    x <- rank(x, 
              na.last = "keep", 
              ties.method = "average")
    n <- sum(!is.na(x))
    x / (n + 1)
  }
  normalize <- function (x) {
    qnorm(uniformize(x))
  }
  op <- par(mfrow=c(4,2))
  f(x1, FUN = normalize, main = "Uniform distribution")
  f(x3, FUN = normalize, main = "Skewed distribution")
  f(x2, FUN = normalize, main = "Fat-tailed distribution")
  f(x4, FUN = normalize, main = "Idem with one outlier")
  par(op)
dev.off()

png(file="g20.png", width=600, height=600)
  library(e1071) # For the "skewness" and "kurtosis" functions
  n <- 1000
  x <- rnorm(n)
  op <- par(mar=c(3,3,4,2)+.1)
  hist(x, col="light blue", probability=TRUE, 
       main=paste("skewness =", round(skewness(x), digits=2)),
       xlab="", ylab="")
  lines(density(x), col="red", lwd=3)
  par(op)
dev.off()

png(file="g21.png", width=600, height=600)
  x <- rexp(n)
  op <- par(mar=c(3,3,4,2)+.1)
  hist(x, col="light blue", probability=TRUE, 
       main=paste("skewness =", round(skewness(x), digits=2)),
       xlab="", ylab="")
  lines(density(x), col="red", lwd=3)
  par(op)
dev.off()

png(file="g22.png", width=600, height=600)
  x <- -rexp(n)
  op <- par(mar=c(3,3,4,2)+.1)
  hist(x, col="light blue", probability=TRUE, 
       main=paste("skewness =", round(skewness(x), digits=2)),
       xlab="", ylab="")
  lines(density(x), col="red", lwd=3)
  par(op)
dev.off()

png(file="g23.png", width=600, height=600)
  library(e1071) # For the "skewness" and "kurtosis" functions
  n <- 1000
  x <- rnorm(n)
  qqnorm(x, main=paste("kurtosis =", round(kurtosis(x), digits=2),
                       "(gaussian)"))
  qqline(x, col="red")
  op <- par(fig=c(.02,.5,.5,.98), new=TRUE)
  hist(x, probability=T,
       col="light blue", xlab="", ylab="", main="", axes=F)
  lines(density(x), col="red", lwd=2)
  box()
  par(op)
dev.off()

png(file="g24.png", width=600, height=600)
  set.seed(1)
  x <- rt(n, df=4)
  qqnorm(x, main=paste("kurtosis =", round(kurtosis(x), digits=2),
                       "(T, df=4)"))
  qqline(x, col="red")
  op <- par(fig=c(.02,.5,.5,.98), new=TRUE)
  hist(x, probability=T,
       col="light blue", xlab="", ylab="", main="", axes=F)
  lines(density(x), col="red", lwd=2)
  box()
  par(op)
dev.off()

png(file="g25.png", width=600, height=600)
  x <- runif(n)
  qqnorm(x, main=paste("kurtosis =", round(kurtosis(x), digits=2),
                       "(uniform)"))
  qqline(x, col="red")
  op <- par(fig=c(.02,.5,.5,.98), new=TRUE)
  hist(x, probability=T,
       col="light blue", xlab="", ylab="", main="", axes=F)
  lines(density(x), col="red", lwd=2)
  box()
  par(op)
dev.off()

png(file="g26.png", width=600, height=600)
  op <- par(mfrow=c(2,2), mar=c(3,2,2,2)+.1)
  data(EuStockMarkets)
  x <- EuStockMarkets
  # We aren't interested in the spot prices, but in the returns
  # return[i] = ( price[i] - price[i-1] ) / price[i-1]
  y <- apply(x, 2, function (x) { diff(x)/x[-length(x)] })
  # We normalize the data
  z <- apply(y, 2, function (x) { (x-mean(x))/sd(x) })
  for (i in 1:4) {
    d <- density(z[,i])
    plot(d$x,log(d$y),ylim=c(-5,1),xlim=c(-5,5))
    curve(log(dnorm(x)),col='red',add=T)
    mtext(colnames(x)[i], line=-1.5, font=2)
  }
  par(op)
  mtext("Are stock returns gaussian?", line=3, font=2)
dev.off()

png(file="g27.png", width=600, height=600)
  n <- dim(z)[1]
  N <- 2000    # Two thousand samples of the same size
  m <- matrix(rnorm(n*N), nc=N, nr=n)
  a <- apply(m^3,2,mean)
  b <- apply(z^3,2,mean)
  op <- par(mar=c(3,3,4,1)+.1)
  hist(a, col='light blue', xlim=range(c(a,b)),
       main="Third moment (skewness)",
       xlab="", ylab="")
  h <- rep(.2*par("usr")[3] + .8*par("usr")[4], length(b))
  points(b, h, type='h', col='red',lwd=3)
  points(b, h, col='red', lwd=3)
  text(b, h, names(b), pos=3)
  par(op)
dev.off()

png(file="g28.png", width=600, height=600)
  n <- dim(z)[1]
  N <- 2000
  m <- matrix(rnorm(n*N), nc=N, nr=n)
  a <- apply(m^4,2,mean) - 3
  b <- apply(z^4,2,mean) - 3
  op <- par(mar=c(3,3,4,1)+.1)
  hist(a, col='light blue', xlim=range(c(a,b)),
       main="Expected kurtosis distribution and observed values",
       xlab="", ylab="")
  h <- rep(.2*par("usr")[3] + .8*par("usr")[4], length(b))
  points(b, h, type='h', col='red',lwd=3)
  points(b, h, col='red', lwd=3)
  text(b, h, names(b), pos=3)
  par(op)
dev.off()

png(file="g29.png", width=600, height=600)
  data(EuStockMarkets)
  x <- EuStockMarkets
  y <- apply(x, 2, function (x) { diff(x)/x[-length(x)] })

  library(lmomco)
  n <- dim(z)[1]
  N <- 200
  m <- matrix(rnorm(n*N), nc=N, nr=n)

  # We normalize the data in the same way
  f <- function (x) { 
    r <- lmom.ub(x)
    (x - r$L1) / r$L2
  }
  z <- apply(y, 2, f)
  m <- apply(m, 2, f)

  a <- apply(m, 2, function (x) lmom.ub(x)$TAU3)
  b <- apply(z, 2, function (x) lmom.ub(x)$TAU3)
  op <- par(mar=c(3,3,4,1)+.1)
  hist(a, col='light blue', xlim=range(c(a,b)),
       main="Expected L-skewness distribution and observed values",
       xlab="", ylab="")
  h <- rep(.2*par("usr")[3] + .8*par("usr")[4], length(b))
  points(b, h, type='h', col='red',lwd=3)
  points(b, h, col='red', lwd=3)
  text(b, h, names(b), pos=3)
  par(op)
dev.off()

png(file="g30.png", width=600, height=600)
  a <- apply(m, 2, function (x) lmom.ub(x)$TAU4)
  b <- apply(z, 2, function (x) lmom.ub(x)$TAU4)
  op <- par(mar=c(3,3,4,1)+.1)
  hist(a, col='light blue', xlim=range(c(a,b)),
       main="Expected L-kurtosis distribution and observed values",
       xlab="", ylab="")
  h <- rep(.2*par("usr")[3] + .8*par("usr")[4], length(b))
  points(b, h, type='h', col='red',lwd=3)
  points(b, h, col='red', lwd=3)
  text(b, h, names(b), pos=3)
  par(op)
dev.off()

png(file="g31.png", width=600, height=200)
  data(faithful)
  stripchart(faithful$eruptions, main="The \"stripchart\" function")
dev.off()

png(file="g32.png", width=600, height=200)
  # Only horizontal noise
  stripchart(faithful$eruptions, jitter=TRUE,
             main="jittered scatterplot")
dev.off()

png(file="g33.png", width=600, height=200)
  stripchart(faithful$eruptions, method='jitter',
             main="jittered scatterplot")
dev.off()

png(file="g34.png", width=600, height=200)
  my.jittered.stripchart <- function (x) {
    x.name <- deparse(substitute(x))
    n <- length(x)
    plot( runif(n) ~ x, xlab=x.name, ylab='noise',
          main="jittered scatterplot" )
  }
  my.jittered.stripchart(faithful$eruptions)
dev.off()

png(file="g35.png", width=600, height=200)
  my.jittered.stripchart <- function (x) {
    x.name <- deparse(substitute(x))
    n <- length(x)
    x <- x + diff(range(x))*.05* (-.5+runif(n))
    plot( runif(n) ~ x, 
          xlab=paste("jittered", x.name), ylab='noise',
          main="jittered scatterplot" )
  }
  my.jittered.stripchart(faithful$eruptions)
dev.off()

png(file="g36.png", width=600, height=600)
  op <- par(mar=c(3,4,2,2)+.1)
  plot( sort( faithful$eruptions ),
        xlab = "" 
      )
  par(op)
dev.off()

png(file="g37.png", width=600, height=600)
  op <- par(mar=c(3,4,2,2)+.1)
  plot(sort(faithful$eruptions), xlab="")
  rug(faithful$eruptions, side=2)
  par(op)
dev.off()

png(file="g38.png", width=600, height=600)
  op <- par(mar=c(3,4,2,2)+.1)
  x <- round( rnorm(100), digits=1 )
  plot(sort(x))
  rug(jitter(x), side=2)
  par(op)
dev.off()

png(file="g39.png", width=600, height=600)
  cumulated.frequencies <- function (x, main="") {
    x.name <- deparse(substitute(x))
    n <- length(x)
    plot( 1:n ~ sort(x), 
          xlab = x.name, 
          ylab = 'Cumulated frequencies',
          main = main
        )
  }    
  cumulated.frequencies(faithful$eruptions,
                        main = "Eruption lengths")
dev.off()

png(file="g40.png", width=600, height=800)
  data(islands)
  dotchart(islands, main="Island area")
dev.off()

png(file="g41.png", width=600, height=800)
  dotchart(sort(log(islands)), 
           main="Island area (logarithmic scale)")
dev.off()

png(file="g42.png", width=800, height=600)
  op <- par(mfcol=c(2,4), mar=c(2,2,1,1)+.1)
  do.it <- function (x) {
    hist(x, probability=T, col='light blue',
         xlab="", ylab="", main="", axes=F)
    axis(1)
    lines(density(x), col='red', lwd=3)
    x <- sort(x)
    q <- ppoints(length(x))
    plot(q~x, type='l',
         xlab="", ylab="", main="")
    abline(h=c(.25,.5,.75), lty=3, lwd=3, col='blue')
  }
  n <- 200
  do.it(rnorm(n))
  do.it(rlnorm(n))
  do.it(-rlnorm(n))
  do.it(rnorm(n, c(-5,5)))
  par(op)
dev.off()

png(file="g43.png", width=600, height=600)
  N <- 2000
  x <- rnorm(N)
  op <- par(mar=c(0,0,0,0), oma=c(0,0,0,0)+.1)
  layout(matrix(c(1,1,1,2), nc=1))
  y <- ppoints( length(x) )
  plot(sort(x), y, type="l", lwd=3,
       xlab="", ylab="", main="")
  abline(h=c(0,.25,.5,.75,1), lty=3)
  abline(v = quantile(x), col = "blue", lwd = 3, lty=2)
  points(quantile(x), c(0,.25,.5,.75,1), lwd=10, col="blue")
  boxplot(x, horizontal = TRUE, col = "pink", lwd=5)  
  abline(v = quantile(x), col = "blue", lwd = 3, lty=2)
  par(new=T)
  boxplot(x, horizontal = TRUE, col = "pink", lwd=5)  
  par(op)
dev.off()

png(file="g44.png", width=600, height=600)
  boxplot(faithful$eruptions, range=0)
dev.off()

png(file="g45.png", width=600, height=200)
  boxplot(faithful$eruptions, range=0, horizontal=TRUE)
dev.off()

png(file="g46.png", width=600, height=600)
  op <- par(mfrow=c(1,2), mar=c(3,2,4,2)+.1)
  do.it <- function (x, xlab="", ylab="", main="") {
    d <- density(x)
    plot(d, type='l', xlab=xlab, ylab=ylab, main=main)
    q <- quantile(x)
    do.it <- function (i, col) {
      x <- d$x[i]
      y <- d$y[i]
      polygon( c(x,rev(x)), c(rep(0,length(x)),rev(y)), border=NA, col=col )
    }
    do.it(d$x <= q[2], 'red')
    do.it(q[2] <= d$x & d$x <= q[3], 'green')
    do.it(q[3] <= d$x & d$x <= q[4], 'blue')
    do.it(d$x >= q[4], 'yellow')
    lines(d, lwd=3)
  }
  do.it( rnorm(2000), main="Gaussian" )
  do.it( rexp(200), main="Exponential" )
  par(op)
  mtext("Quartiles", side=3, line=3, font=2, cex=1.2)
dev.off()

png(file="g47.png", width=600, height=200)
  boxplot(faithful$eruptions, horizontal = TRUE,
          main = "No outliers")
dev.off()

png(file="g48.png", width=600, height=200)
  # There are outliers, they might bring trouble, 
  # but it is normal, it is not pathological
  boxplot(rnorm(500), horizontal = TRUE,
          main = "Normal outliers")
dev.off()

png(file="g49.png", width=600, height=200)
  x <- c(rnorm(30),20)
  x <- sample(x, length(x))
  boxplot( x, horizontal = TRUE,
           main = "An outlier" )
dev.off()

png(file="g50.png", width=600, height=200)
  library(boot)
  data(aml)
  boxplot( aml$time, horizontal = TRUE,
           main = "An outlier" )
dev.off()

png(file="g51.png", width=600, height=200)
  data(attenu)
  boxplot(attenu$dist, horizontal = TRUE,
          main = "Non gaussian (asymetric) data")
dev.off()

png(file="g52.png", width=600, height=200)
  data(attenu)
  boxplot(log(attenu$dist), horizontal = TRUE,
          main = "Transformed variable")
dev.off()

png(file="g53.png", width=600, height=300)
  y <- c(rnorm(10+100+1000+10000+100000))
  x <- c(rep(1,10), rep(2,100), rep(3,1000), rep(4,10000), rep(5,100000))
  x <- factor(x)
  plot(y~x, 
       horizontal = TRUE, 
       col = "pink",
       las = 1,
       xlab = "", ylab = "",
       main = "The larger the sample, the more outliers")
dev.off()

png(file="g54.png", width=600, height=200)
  boxplot(faithful$eruptions, 
          notch = TRUE,
          horizontal = TRUE,
          main = "Confidence interval on the median...")
dev.off()

png(file="g55.png", width=600, height=300)
  library(boot)
  data(breslow)
  boxplot(breslow$n, 
          notch = TRUE,
          horizontal = TRUE, 
          col = "pink",
          main = "...that goes beyond the quartiles")
dev.off()

png(file="g56.png", width=600, height=200)
  boxplot(faithful$eruptions, 
          horizontal = TRUE,
          col = "pink")
  rug(faithful$eruption, 
      ticksize = .2)
dev.off()

png(file="g57.png", width=600, height=600)
  hist(faithful$eruptions)
dev.off()

png(file="g58.png", width=600, height=600)
  hist(faithful$eruptions, breaks=20, col="light blue")
dev.off()

png(file="g59.png", width=600, height=600)
  op <- par(mfrow=c(2,1), mar=c(2,2,2,1)+.1)
  hist(faithful$eruptions, breaks=seq(1,6,.5), 
       col='light blue',
       xlab="", ylab="", main="")
  hist(faithful$eruptions, breaks=.25+seq(1,6,.5), 
       col='light blue',
       xlab="", ylab="", main="")
  par(op)
  mtext("Is the first peak symetric or not?", 
        side=3, line=2.5, font=2.5, size=1.5)
dev.off()

png(file="g60.png", width=600, height=600)
  hist(faithful$eruptions, 
       probability=TRUE, breaks=20, col="light blue",
       xlab="", ylab="",
       main="Histogram and density estimation")
  points(density(faithful$eruptions, bw=.1), type='l', 
         col='red', lwd=3)
dev.off()

png(file="g61.png", width=600, height=600)
  hist(faithful$eruptions, 
       probability=TRUE, breaks=20, col="light blue",
       xlab="", ylab="",
       main="Histogram and density estimation")
  points(density(faithful$eruptions, bw=1),  type='l', 
         lwd=3, col='black')
  points(density(faithful$eruptions, bw=.5), type='l', 
         lwd=3, col='blue')
  points(density(faithful$eruptions, bw=.3), type='l', 
         lwd=3, col='green')
  points(density(faithful$eruptions, bw=.1), type='l', 
         lwd=3, col='red')
dev.off()

png(file="g62.png", width=600, height=600)
  hist(faithful$eruptions, 
       probability=TRUE, breaks=20, col="light blue",
       main="")
  rug(faithful$eruptions)
  points(density(faithful$eruptions, bw=.1), type='l', lwd=3, col='red')
  f <- function(x) {
    dnorm(x, 
          mean=mean(faithful$eruptions), 
          sd=sd(faithful$eruptions),
    ) 
  }
  curve(f, add=T, col="red", lwd=3, lty=2)
dev.off()

png(file="g63.png", width=600, height=600)
  symetry.plot <- function (x0, 
                            main="Symetry plot",
                            breaks="Sturges", ...) {
    x <- x0[ !is.na(x0) ]
    x <- sort(x)
    x <- abs(x - median(x))
    n <- length(x)
    nn <- ceiling(n/2)
    plot( x[n:(n-nn+1)] ~ x[1:nn] ,
          xlab='Distance below median',
          ylab='Distance above median', 
          main=main,
          ...)
    abline(0,1, col="blue", lwd=3)  
    op <- par(fig=c(.02,.5,.5,.98), new=TRUE)
    hist(x0, probability=T, breaks=breaks,
         col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x0), col="red", lwd=2)
    box()
    par(op)
  }

  symetry.plot(rnorm(500), 
               main="Symetry plot (gaussian distribution)")
dev.off()

png(file="g64.png", width=600, height=600)
  symetry.plot(rexp(500),
               main="Symetry plot (exponential distribution)")

dev.off()

png(file="g65.png", width=600, height=600)
  symetry.plot(-rexp(500),
               main="Symetry plot (negative skewness)")

dev.off()

png(file="g66.png", width=600, height=600)
  symetry.plot(rexp(500),
               main="Symetry plot, logarithmic scales)")

dev.off()

png(file="g67.png", width=600, height=600)
  symetry.plot(faithful$eruptions, breaks=20)
dev.off()

png(file="g68.png", width=600, height=600)
  symetry.plot.2 <- function (x, N=1000, 
                              pch=".", cex=1, ...) {
    x <- x[ !is.na(x) ]
    x <- sort(x)
    x <- abs(x - median(x))
    n <- length(x)
    nn <- ceiling(n/2)
    plot( x[n:(n-nn+1)] ~ x[1:nn] ,
          xlab='Distance below median',
          ylab='Distance above median', 
          ...)
    for (i in 1:N) {
      y <- sort( rnorm(n) )
      y <- abs(y - median(y))
      m <- ceiling(n/2)
      points( y[n:(n-m+1)] ~ y[1:m], 
              pch=pch, cex=cex, col='red' )
    }
    points(x[n:(n-nn+1)] ~ x[1:nn] , ...)
    abline(0,1, col="blue", lwd=3)  
  }

  n <- 100
  symetry.plot.2( rnorm(n), pch='.', lwd=3, 
                  main=paste("Symetry plot: gaussian,", n, "observations"))
dev.off()

png(file="g69.png", width=600, height=600)
  n <- 10
  symetry.plot.2( rnorm(n), pch=15, lwd=3, type="b", cex=.5,
                  main=paste("Symetry plot: gaussian,", n, "observations"))
dev.off()

png(file="g70.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g71.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g72.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g73.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g74.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g75.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g76.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g77.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2, lwd = 3)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"CAC"])
dev.off()

png(file="g78.png", width=600, height=600)
  robust.symetry.plot <- function (x, 
                                   N = max(ceiling(1000/length(x)),2),
                                   alpha = .05, 
                                   xlab = "Distance below the median",
                                   ylab = "Distance above the median", 
                                   main = "Symetry plot",
                                   ...) {
    cat(N, "\n")
    # The symetry plot
    x <- x[!is.na(x)]
    n <- length(x)
    nn <- ceiling(n/2)
    x <- sort(x)
    d <- abs(x - median(x))       # Distance to the median
    plot( d[1:nn], d[n:(n-nn+1)], 
          xlab = xlab, ylab = ylab, 
          main = main,
          ... )
    # The symetry plot of resampled, symetric data
    y <- c(x, 2 * median(x) - x)  # We symetrize the data
    X <- Y <- rep(NA, N * nn)
    for (i in 1:N) {
      a <- sort(sample(y, n))
      a <- abs(a - median(a))
      j <- ((i-1) * nn + 1) : (i * nn) 
      X[j] <- a[1:nn]
      Y[j] <- a[n:(n-nn+1)]
    }
    points(X, Y, col="red")
    points( d[1:nn], d[n:(n-nn+1)], ...) 
    # The 5% confidence interval stemming from the resampled data
    require(quantreg)
    for (tau in c(alpha, 1-alpha)) {
      r <- lprq(X, Y,
                h = bw.nrd0(x),  # See ?density
                tau = tau)
      lines(r$xx, r$fv, col = "blue", lwd = 3)
    }
    abline(0, 1, col = "blue", lty = 2)
    # The histogram, in a corner
    op <- par(fig = if (skewness(x)>0) 
                      c(.02,.5,.5,.98)       # Top left corner
                    else c(.5,.98,.02,.5),   # Bottom right
              new = TRUE)
    hist(x, probability=T,
        col="light blue", xlab="", ylab="", main="", axes=F)
    lines(density(x), col="red", lwd=2)
    box()
    par(op)
  }
  robust.symetry.plot(EuStockMarkets[,"FTSE"])
dev.off()

png(file="g79.png", width=600, height=600)
  robust.symetry.plot(rnorm(100), N=100, pch=16)
dev.off()

png(file="g80.png", width=600, height=600)
  data(airquality)
  x <- airquality[,4]
  hist(x, probability=TRUE, breaks=20, col="light blue")
  rug(jitter(x, 5))
  points(density(x), type='l', lwd=3, col='red')
  f <- function(t) {
    dnorm(t, mean=mean(x), sd=sd(x) ) 
  }
  curve(f, add=T, col="red", lwd=3, lty=2)
dev.off()

png(file="g81.png", width=600, height=600)
  x <- airquality[,4]
  qqnorm(x)
  qqline(x,
         col="red", lwd=3)
dev.off()

png(file="g82.png", width=600, height=600)
  y <- rnorm(100)
  qqnorm(y, main="Gaussian random variable")
  qqline(y, 
         col="red", lwd=3)
dev.off()

png(file="g83.png", width=600, height=600)
  y <- rnorm(100)^2
  qqnorm(y, main="Non gaussian variable")
  qqline(y,
         col="red", lwd=3)
dev.off()

png(file="g84.png", width=600, height=600)
  my.qqnorm <- function (x, N=1000, ...) {
    op <- par()
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    s <- sd(x)
    print("a")
    qqnorm(x, axes=F, ...)
    for (i in 1:N) {
      par(new=T)
      qqnorm(rnorm(n, mean=m, sd=s), col='red', pch='.', 
             axes=F, xlab='', ylab='', main='')
    }
    par(new=T)
    qqnorm(x, ...)
    qqline(x, col='blue', lwd=3)
    par(op)
  }
  my.qqnorm(rnorm(100), 
            main = "QQplot: Gaussian distribution")
dev.off()

png(file="g85.png", width=600, height=600)
  my.qqnorm(runif(100), 
            main = "uniform distribution")
dev.off()

png(file="g86.png", width=600, height=600)
  my.qqnorm(exp(rnorm(100)), 
            main = 'log-normal distribution')
dev.off()

png(file="g87.png", width=600, height=600)
  my.qqnorm(c(rnorm(50), 5+rnorm(50)), 
            main = 'bimodal distribution')
dev.off()

png(file="g88.png", width=600, height=600)
  my.qqnorm(c(rnorm(50), 20+rnorm(50)), 
            main = 'two remote peaks')
dev.off()

png(file="g89.png", width=600, height=600)
  x <- rnorm(100)
  x <- x + x^3
  my.qqnorm(x, main = 'fat tails')
dev.off()

png(file="g90.png", width=600, height=600)
  y <- exp(rnorm(100))
  qqnorm(y, 
         main = '(1) Log-normal distribution')
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g91.png", width=600, height=600)
  y <- rnorm(100)^2
  qqnorm(y, ylim = c(-2,2), 
         main = "(2) Square of a gaussian variable")
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g92.png", width=600, height=600)
  y <- -exp(rnorm(100))
  qqnorm(y, ylim = c(-2,2), 
         main = "(3) Opposite of a log-normal variable")
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g93.png", width=600, height=600)
  y <- runif(100, min=-1, max=1)
  qqnorm(y, ylim = c(-2,2), 
         main = '(4) Uniform distribution')
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g94.png", width=600, height=600)
  y <- rnorm(10000)^3
  qqnorm(y, ylim = c(-2,2), 
         main = "(5) Cube of a gaussian r.v.")
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g95.png", width=600, height=600)
  y <- c(rnorm(50), 5+rnorm(50))
  qqnorm(y, 
         main = '(6) Two peaks')
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g96.png", width=600, height=600)
  y <- c(rnorm(50), 20+rnorm(50))
  qqnorm(y, 
         main = '(7) Two peaks, farther away')
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g97.png", width=600, height=600)
  y <- sample(seq(0,1,.1), 100, replace=T)
  qqnorm(y, 
         main = '(7) Discrete distribution')
  qqline(y, 
         col = 'red', lwd = 3)
dev.off()

png(file="g98.png", width=600, height=600)
  x <- seq(from=0, to=2, length=100)
  y <- exp(x)-1
  plot( y ~ x, type = 'l', col = 'red',
        xlim = c(-2,2), ylim = c(-2,2), 
        xlab = "Theoretical (gaussian) quantiles", 
        ylab = "Sample quantiles")        
  lines( x~y, type='l', col='green')
  x <- -x
  y <- -y
  lines( y~x, type='l', col='blue', )
  lines( x~y, type='l', col='cyan')
  abline(0,1)
  legend( -2, 2,
          c( "less concentrated on the right",  
             "more concentrates on the right",
             "less concentrated on the left",
             "more concentrated on the left"
           ),
          lwd=3,
          col=c("red", "green", "blue", "cyan")
        )
  title(main="Reading a qqplot")
dev.off()

png(file="g99.png", width=600, height=600)
  op <- par()
  layout( matrix( c(2,2,1,1), 2, 2, byrow=T ),
          c(1,1), c(1,6),
        )
  # The plot
  n <- 100
  y <- rnorm(n)
  x <- qnorm(ppoints(n))[order(order(y))]
  par(mar=c(5.1,4.1,0,2.1))
  plot( y ~ x, col = "blue", 
        xlab = "Theoretical (gaussian) quantiles", 
        ylab = "Sample quantiles" )
  y1 <- scale( rnorm(n)^2 )
  x <- qnorm(ppoints(n))[order(order(y1))]
  lines(y1~x, type="p", col="red")
  y2 <- scale( -rnorm(n)^2 )
  x <- qnorm(ppoints(n))[order(order(y2))]
  lines(y2~x, type="p", col="green")
  abline(0,1)

  # The legend
  par(bty='n', ann=F)
  g <- seq(0,1, length=10)
  e <- g^2
  f <- sqrt(g)
  h <- c( rep(1,length(e)), rep(2,length(f)), rep(3,length(g)) )
  par(mar=c(0,4.1,1,0))
  boxplot( c(e,f,g) ~ h, horizontal=T, 
           border=c("red", "green", "blue"),
           col="white", # Something prettier?
           xaxt='n',
           yaxt='n', 
           )
  title(main="Reading a qqplot")
  par(op)
dev.off()

png(file="g100.png", width=600, height=600)
  y <- rnorm(100)^2
  y <- scale(x)
  y <- sort(x)
  x <- qnorm( seq(0,1,length=length(y)) )
  plot(y~x)
  abline(0,1)
dev.off()

png(file="g101.png", width=600, height=600)
  qq <- function (y, ylim, quantiles=qnorm,
      main = "Q-Q Plot", xlab = "Theoretical Quantiles",
      ylab = "Sample Quantiles", plot.it = TRUE, ...)
  {
    y <- y[!is.na(y)]
    if (0 == (n <- length(y)))
      stop("y is empty")
    if (missing(ylim))
      ylim <- range(y)
    x <- quantiles(ppoints(n))[order(order(y))]
    if (plot.it)
      plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim,
            ...)
    # From qqline
    y <- quantile(y, c(0.25, 0.75))
    x <- quantiles(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1] - slope * x[1]
    abline(int, slope, ...)
    invisible(list(x = x, y = y))
  }

  y <- runif(100)
  qq(y, quantiles=qunif)
dev.off()

png(file="g102.png", width=600, height=600)
  two.point.line <- function (x1,y1,x2,y2, ...) {
    a1 <- (y2-y1)/(x2-x1)
    a0 <- y1 - a1 * x1
    abline(a0,a1,...)
  }
  trended.probability.plot <- function (x, q=qnorm) {
    n <- length(x)
    plot( sort(x) ~ q(ppoints(n)), 
          xlab='theoretical quantiles', 
          ylab='sample quantiles')
    two.point.line(q(.25), quantile(x,.25), 
                   q(.75), quantile(x,.75), col='red')
  }
  detrended.probability.plot <- function (x, q=qnorm,
                                          xlab="", ylab="") {
    n <- length(x)
    x <- sort(x)
    x1 <- q(.25)
    y1 <- quantile(x,.25)
    x2 <- q(.75)
    y2 <- quantile(x,.75)
    a1 <- (y2-y1)/(x2-x1)
    a0 <- y1 - a1 * x1
    u <- q(ppoints(n))
    x <- x - (a0 + a1 * u)
    plot(x ~ u, 
         xlab=xlab, ylab=ylab)
    abline(h=0, col='red')
  }

  op <- par(mfrow = c(3,2), mar = c(2,2,2,2) + .1)
  x <- runif(20)
  trended.probability.plot(x)
  detrended.probability.plot(x)
  x <- runif(500)
  trended.probability.plot(x)
  detrended.probability.plot(x)
  trended.probability.plot(x, qunif)
  detrended.probability.plot(x,qunif)
  par(op)
  mtext("Detrended quantile-quantile plots", 
        side=3, line=3, font=2, size=1.5)
dev.off()

png(file="g103.png", width=600, height=600)
  xy <- matrix(c( 0, 0,
                 .2, .9, 
                 .3, .95, 
                 .5, .99,
                  1, 1), byrow = T, nc = 2)
  plot(xy, type = 'b', pch = 15,
       main = "Conventration curve",
       xlab = "patients",
       ylab = "expenses")
  polygon(xy, border=F, col='pink')
  lines(xy, type='b', pch=15)
  abline(0,1,lty=2)
dev.off()

png(file="g104.png", width=600, height=600)
  x <- c(0,20,2720,3000)/3000
  y <- c(0,100000,100100,100110)/100110
  plot(x,y, type='b', pch=15,
       xlab = "Genes", ylab = "ARNm",
       main = "Conventration curve")
  polygon(x,y, border=F, col='pink')
  lines(x,y, type='b', pch=15)
  abline(0,1,lty=2)
dev.off()

png(file="g105.png", width=600, height=600)
  library(ineq)
  op <- par(mfrow=c(3,3), mar=c(2,3,3,2)+.1, oma=c(0,0,2,0))
  n <- 500
  set.seed(1)
  plot(Lc(runif(n,0,1)),   
       main="uniform on [0,1]", col='red',
       xlab="", ylab="")
  do.it <- function (x, main="", xlab="", ylab="") { 
    plot(Lc(x), col = "red",
         main=main, xlab=xlab, ylab=ylab)
  }
  do.it(runif(n,0,10),       main="uniform on [0,10]")
  do.it(runif(n,10,11),      main="uniform on [10,11]")
  do.it(rlnorm(n),          main="log-normal")
  do.it(rlnorm(n,0,4),      main="log-normal, wider")
  do.it(abs(rcauchy(n,1)),  main="half-Cauchy")
  do.it(abs(rnorm(n,1)),    main="half-Gaussian")
  do.it(rpois(n,1),         main="Poisson with mean 1")
  do.it(rpois(n,10),        main="Poisson with mean 10")
  par(op)
  mtext("Gini concentration curves", side=3, line=3, 
        font=2, cex=1.5)
dev.off()

png(file="g106.png", width=600, height=400)
  data(esoph)
  dotchart(table(esoph$agegp))
  mtext("Misleading plot", side=3, line=2.5, font=2, cex=1.2)
  mtext("The origin is not on the plot", side=3, line=1)
dev.off()

png(file="g107.png", width=600, height=600)
  barplot(table(esoph$agegp))
dev.off()

png(file="g108.png", width=600, height=600)
  hist(as.numeric(esoph$agegp), 
       breaks=seq(.5,.5+length(levels(esoph$agegp)),step=1),
       col='light blue')
dev.off()

png(file="g109.png", width=600, height=200)
  boxplot(as.numeric(esoph$agegp), 
          horizontal = T, col = "pink")
dev.off()

png(file="g110.png", width=600, height=200)
  stripchart(jitter(as.numeric(esoph$agegp),2), method='jitter')
dev.off()

png(file="g111.png", width=600, height=600)
  plot(table(esoph$agegp), type='b', pch=7)
dev.off()

png(file="g112.png", width=600, height=600)
  data(HairEyeColor)
  x <- apply(HairEyeColor, 2, sum)
  barplot(x)
  title(main="Column plot")
dev.off()

png(file="g113.png", width=600, height=600)
  barplot(x, col = 1, density = c(3,7,11,20), 
          angle = c(45,-45,45,-45))
  title(main = "Column plot")
dev.off()

png(file="g114.png", width=200, height=600)
  x <- apply(HairEyeColor, 2, sum)
  barplot(as.matrix(x), legend.text = TRUE)
  title("Bar plot")
dev.off()

png(file="g115.png", width=600, height=200)
  barplot(as.matrix(x), 
          horiz = TRUE, 
          col = rainbow(length(x)), 
          legend.text = TRUE)
  title(main = "Bar plot")
dev.off()

png(file="g116.png", width=400, height=600)
  op <- par(no.readonly=TRUE)
  par(mar=c(5,4,4,7)+.1)
  barplot(as.matrix(x))
  title("Bar plot, with legend")
  par(xpd=TRUE)  # Do not clip to the drawing area
  lambda <- .025
  legend(par("usr")[2], 
         par("usr")[4],
         names(x),
         fill = grey.colors(length(x)),         
         xjust = 0, yjust = 1
        )
  par(op)
dev.off()

png(file="g117.png", width=600, height=200)
  op <- par(no.readonly=TRUE)
  par(mar=c(3,1,4,7)+.1)
  barplot(as.matrix(x), 
          horiz = TRUE, 
          col = rainbow(length(x)))
  title(main = "Bar plot, with legend")
  par(xpd=TRUE)  # Do not clip to the drawing area
  lambda <- .05
  legend((1+lambda)*par("usr")[2] - lambda*par("usr")[1], 
         par("usr")[4],
         names(x),
         fill = rainbow(length(x)),         
         xjust = 0, yjust = 1
        )
  par(op)
dev.off()

png(file="g118.png", width=600, height=600)
  data(attenu)
  op <- par(las=2) # Write the labels perpendicularly to the axes
  barplot(table(attenu$event))
  title(main="Column plot")
  par(op)
dev.off()

png(file="g119.png", width=600, height=600)
  op <- par(las=2)
  barplot(rev(sort(table(attenu$event))))
  title(main="Pareto Plot")
  par(op)
dev.off()

png(file="g120.png", width=600, height=600)
  # I cannot seem to manage to do it with 
  # the "barplot" function...
  pareto <- function (x, main = "", ylab = "Value") {
    op <- par(mar = c(5, 4, 4, 5) + 0.1,
              las = 2)
    if( ! inherits(x, "table") ) {
      x <- table(x)
    }
    x <- rev(sort(x))
    plot( x, type = 'h', axes = F, lwd = 16,
          xlab = "", ylab = ylab, main = main )
    axis(2)
    points( x, type = 'h', lwd = 12, 
            col = heat.colors(length(x)) )
    y <- cumsum(x)/sum(x)
    par(new = T)
    plot(y, type = "b", lwd = 3, pch = 7, 
         axes = FALSE, 
         xlab='', ylab='', main='')
    points(y, type = 'h')
    axis(4)
    par(las=0)
    mtext("Cumulated frequency", side=4, line=3)
    print(names(x))
    axis(1, at=1:length(x), labels=names(x))
    par(op)
  }
  pareto(attenu$event)    
  title(main="Pareto plot with cumulated frequencies")
dev.off()

png(file="g121.png", width=600, height=600)
  x <- apply(HairEyeColor, 2, sum)
  pie(x)
  title(main="Pie chart")
dev.off()

png(file="g122.png", width=600, height=800)
  op <- par(mfrow=c(4,2), mar=c(2,4,2,2))

  # Barchart (1 bar)
  set.seed(1)
  x <- rlnorm(6)
  barplot(as.matrix(x), 
          xlim = c(-2,3),
          main = "Barchart")

  # Barchart with an added dimension (stacked area chart) (p.173)
  y <- matrix(rnorm(60), nc=6)
  y <- apply(y, 2, cumsum)
  y <- exp(y/5)
  stacked_area_chart <- function (y, axes = TRUE, ...) {
    stopifnot(all(y>=0))
    y <- t(apply(y, 1, cumsum))
    plot.new()
    plot.window(xlim = c(1,nrow(y)), 
                ylim = range(y) + .1*c(-1,1)*diff(range(y)))
    for (i in ncol(y):1) {
      polygon(c(1,1:nrow(y),nrow(y)), 
              c(0,y[,i],0),
              col=i, border=NA)
      lines(1:nrow(y), y[,i], lwd=3)
    }
    if (axes) {
      axis(1)
      axis(2)
    }
    box()
  }
  stacked_area_chart(y, axes = FALSE)
  title(main = "Barchart with an added dimension",
        sub = "Stacked area chart")

  # Pie chart
  pie(x, 
      col = 1:length(x), 
      labels = LETTERS[1:length(x)],
      main = "Pie chart")

  # Annular chart
  annular_chart <- function (x, r1=1, r2=2) {
    stopifnot(x>=0, r1 >= 0, r2 > 0, r1 < r2)
    x <- cumsum(x) / sum(x)
    x <- c(0,x)
    plot.new()
    plot.window(xlim = c(-1.1,1.1)*r2,
                ylim = c(-1.1,1.1)*r2)
    for (i in 2:length(x)) {
      theta <- 2*pi*seq(x[i-1], x[i], length=100)
      polygon( c(r1 * cos(theta), r2 * cos(rev(theta))),
               c(r1 * sin(theta), r2 * sin(rev(theta))),
               col = i )
    }
  }
  annular_chart(x)
  title("Annular chart")

  # Pie chart
  pie(x, 
      col = 1:length(x), 
      labels = LETTERS[1:length(x)],
      main = "From bad...")

  # Concentrical chart
  # Grid graphics would be better for this: they would
  # help you enforce orthonormal coordinates, and thus
  # circular circles...
  circular_pie <- function (x, ...) {
    stopifnot(is.vector(x), 
              all(x >= 0),
              length(x) >= 1)
    plot.new()
    radii <- sqrt(cumsum(x)) # The areas are 
                             # proportional to the 
                             # inital x
    plot.window(xlim = max(radii)*c(-1.1,1.1),
                ylim = max(radii)*c(-1.1,1.1) )
    theta <- seq(0, 2*pi, length=100)[-1]
    x <- cos(theta)
    y <- sin(theta)
    for (i in length(x):1) {
      polygon(radii[i] * x, radii[i] * y,
              col = i, border = NA)
      lines(radii[i] * x, radii[i] * y)
    }                
  }
  circular_pie(x)
  title("...to worse")

  # barchart (several bars)
  xx <- sample(x)
  barplot(cbind("1" = x, "2" = xx), 
          space = 1, 
          xlim = c(0,5),
          col = 1:length(x),
          main = "Barchart with several bars")

  # Several annular charts p.212
  annular_chart_ <- function (x, r1=1, r2=2) {
    stopifnot(x>=0, r1 >= 0, r2 > 0, r1 < r2)
    x <- cumsum(x) / sum(x)
    x <- c(0,x)
    for (i in 2:length(x)) {
      theta <- 2*pi*seq(x[i-1], x[i], length=100)
      polygon( c(r1 * cos(theta), r2 * cos(rev(theta))),
               c(r1 * sin(theta), r2 * sin(rev(theta))),
               col = i )
    }
  }
  two_annular_charts <- function (x, y, 
                                 r1=1, r2=1.9, 
                                 r3=2, r4=2.9) {
    plot.new()
    plot.window(xlim = c(-1.1,1.1)*r4,
                ylim = c(-1.1,1.1)*r4)
    annular_chart_(x, r1, r2)
    annular_chart_(y, r3, r4)
  }
  two_annular_charts(x, xx)
  title("Two annular charts")

  par(op)
dev.off()

png(file="g123.png", width=600, height=600)
  library(ape)
  example(plot.ancestral)
dev.off()

png(file="g124.png", width=600, height=600)
  example(plot.phylo)
dev.off()

png(file="g125.png", width=600, height=600)
  ##
  ## barplot2D(area, colour)
  ##
  ##
  ## The algorithm is not that obvious.
  ##   - Start with a rectangle, representing 100%, to be filled
  ##     by other rectangles.
  ##   - Try to put the first rectangle on the left
  ##   - If it too elongated, try to put two rectangles, on
  ##     top of each other, on the left
  ##   - Go on, until you are satisfied
  ##   - When you have put those rectangles, proceed with the
  ##     remaining of the large rectangle container.
  ## More precisely, we choose the number of rectables to
  ## stack so as to minimize the following penalty:
  ##   penalty for the first rectangle in the stack + penalty for the last
  ## where the penalty of a rectangle is
  ##   ratio - 1.1
  ## where "ratio" is the ratio of the longer side by the smaller.
  ##
  ## Arguments:
  ##  area:      vector, containing positive number (NAs are discarded),
  ##             to be used as the area of the rectangles
  ##  colour:    vector (same length) of strings containing the colours
  ##             You can create it with "rgb", or "cm.colors".
  ##  threshold: The maximum acceptable aspect ratio of the rectangles
  ##  width, height: Dimensions of the initial rectangle.
  ##                 I suggest to plot the picture in a rectangular
  ##                 device, e.g.,
  ##                   pdf(width=6, height=4)
  ##                 but to tell the function that this rectangle is
  ##                 actually a square, i.e.,
  ##                   barplot2D(area, colour, width=1, height=1)
  ##                 so that the cells be horizontal
  ##                 rectangles: you get more space to add
  ##                 labels
  ##
  ## Returns:
  ##   A matrix, one row per cell, containing the x- and
  ##   y-coordinates of the corners of all the cells (first
  ##   eight columns), and the coordinates of the center of
  ##   those cells (last two columns).
  ##   The rows are in one-to-one correspondance with the
  ##   elements of the "area" vector: if there were missing
  ##   values, we have rows of missing values.
  ##   The row names are the same as the names of the "area"
  ##   vector, in the same order.
  ##

  barplot2D <- function (area, colour, 
                         threshold=1.1, 
                         width=1, height=1) {
    stopifnot(is.vector(area), is.vector(colour),
              length(area) == length(colour),
              !all(is.na(area)))
    if (is.null(names(area))) {
      names(area) <- as.character(1:length(area))
    }
    area0 <- area
    if (any(is.na(area))) {
      warning("Discarding NAs")
      i <- which(!is.na(area))
      area <- area[i]
      colour <- colour[i]
    }
    stopifnot(all(area>=0), sum(area)>0)
    i <- order(-area)
    area <- area[i]
    colour <- colour[i]
    n <- length(area)
    res <- matrix(NA, nr=n, nc=8)
    colnames(res) <- as.vector(t(outer(LETTERS[1:4], 1:2, paste, sep="")))
    rownames(res) <- names(area)
    A <- c(0,height)
    B <- c(0,0)
    C <- c(width,0)
    D <- c(width,height)
    plot.new()
    plot.window(xlim=c(0,1), ylim=c(0,1))
    i <- 1
    while (i <= n) {
      lambda <- cumsum(area[i:n]) / sum(area[i:n])
      mu <- area[i]   / cumsum(area[i:n])
      nu <- area[i:n] / cumsum(area[i:n])
      penalty1 <- mu * sum(abs(A-B)) / ( lambda * sum(abs(A-D)) )
      penalty1 <- ifelse(penalty1 <= threshold, 0, penalty1 - threshold)
      penalty2 <- lambda * sum(abs(A-D)) / ( nu * sum(abs(A-B)) )
      penalty2 <- ifelse(penalty2 <= threshold, 0, penalty2 - threshold)
      j <- which.min(penalty1 + penalty2)[1] + i - 1
      cat(i, " => ", j, "\n")
      lambda <- sum(area[i:j]) / sum(area[i:n])
      A1 <- A
      B1 <- B
      C1 <- (1-lambda) * B + lambda * C
      D1 <- (1-lambda) * A + lambda * D
      AA <- C1
      BB <- C
      CC <- D
      DD <- D1
      while (i <= j) {
        lambda <- area[i] / sum(area[i:j])
        B2 <- (1-lambda) * A1 + lambda * B1
        C2 <- (1-lambda) * D1 + lambda * C1
        polygon(rbind(A1, B2, C2, D1), col=colour[i])
        res[i,] <- c(A1, B2, C2, D1)
        A1 <- B2 
        D1 <- C2
        i <- i + 1
      }
      A <- AA
      B <- BB
      C <- CC
      D <- DD
    } # Main loop
    res0 <- matrix(NA, nr=length(area0), nc=10)
    colnames(res0) <- c(colnames(res), "x", "y")
    rownames(res0) <- names(area0)
    res0[ names(area), 1:8] <- res  
    res0[, "x"] <- apply(res0[,c("A1","B1","C1","D1")],1,mean)
    res0[, "y"] <- apply(res0[,c("A2","B2","C2","D2")],1,mean)
    invisible(res0)
  }

  N <- 20
  area <- rlnorm(N)
  names(area) <- LETTERS[1:N]
  value <- rt(N, df=4)
  # Difficult part: compute the colours...
  colour <- cm.colors(255)[
    1 + round(
      254 * (value - min(value, na.rm = TRUE)) /
      diff(range(value, na.rm = TRUE))
    )
  ]
  r <- barplot2D(area, colour)
  title("2-dimensional barplot")
  # Add the labels
  text(r[,"x"], r[,"y"], names(area), cex=.8)
dev.off()

png(file="g126.png", width=600, height=600)
  library(portfolio)
  example(map.market)
dev.off()

png(file="g127.png", width=600, height=600)
  olap <- function (x, i) {
    # Project (drill-up?) a data cube
    y <- x <- apply(x, i, sum)
    if (length(i) > 1) {
      y <- as.vector(x)
      n <- dimnames(x)
      m <- n[[1]]
      for (i in (1:length(dim(x)))[-1]) {
        m <- outer(m, n[[i]], paste)
      }
      names(y) <- m
    }
    y
  }
  col1 <- c("red", "green", "blue", "brown")
  col2 <- c("red", "light coral",
            "green", "light green",
            "blue", "light blue",
            "brown", "rosy brown")
  col3 <- col2[c(1,2,1,2,3,4,3,4,5,6,5,6,7,8,7,8)]
  op <- par(mfrow=c(3,1), mar=c(8,4,0,2), oma=c(0,0,2,0), las=2)
  barplot(olap(Titanic,1),   space=0, col=col1)
  barplot(olap(Titanic,2:1), space=0, col=col2)
  barplot(olap(Titanic,3:1), space=0, col=col3)
  par(op)
  mtext("Region tree", font = 2, line = 3)
dev.off()

png(file="g128.png", width=600, height=600)
  x1 <- olap(Titanic,3:1)
  x2 <- rep(olap(Titanic,2:1), each=dim(Titanic)[3])
  x3 <- rep(olap(Titanic,1),   each=prod(dim(Titanic)[2:3]))
  x4 <- rep(sum(Titanic),      each=prod(dim(Titanic)[1:3]))
  op <- par(mar=c(8,4,4,2))
  barplot(x4, names.arg="", axes = FALSE, col = "light coral")
  barplot(x3, names.arg="", axes = FALSE, col = "light green",  add = TRUE)
  barplot(x2, names.arg="", axes = FALSE, col = "light blue",   add = TRUE)
  barplot(x1, las=2,        axes = FALSE, col = "yellow", add = TRUE)
  mtext("TempleMVV Plot", line=2, font=2, cex=1.2)
  par(op)
dev.off()

png(file="g129.png", width=600, height=200)
  x <- apply(HairEyeColor, 2, sum)
  dotchart(x, main="dotchart")
dev.off()

png(file="g130.png", width=600, height=600)
  library(MASS)  # For the Cars93 data set
  dotchart(table(Cars93$Manufacturer))
dev.off()

png(file="g131.png", width=600, height=1100)
  library(nlme)
  data(Milk)
  dotchart(table(Milk$Cow))
dev.off()

png(file="g132.png", width=600, height=600)
  data(cars)
  plot(cars$dist ~ cars$speed, 
       xlab = "Speed (mph)", 
       ylab = "Stopping distance (ft)", 
       las = 1)
  title(main = "Point cloud")
dev.off()

png(file="g133.png", width=600, height=600)
  plot(cars$dist ~ cars$speed, 
       xlab = "Speed (mph)", 
       ylab = "Stopping distance (ft)", 
       las = 1)
  title(main = "cars data")
  rug(side=1, jitter(cars$speed, 5))
  rug(side=2, jitter(cars$dist, 20))
dev.off()

png(file="g134.png", width=600, height=600)
  op <- par()
  layout( matrix( c(2,1,0,3), 2, 2, byrow=T ),
          c(1,6), c(4,1),
        )

  par(mar=c(1,1,5,2))
  plot(cars$dist ~ cars$speed, 
       xlab='', ylab='',
       las = 1)
  rug(side=1, jitter(cars$speed, 5) )
  rug(side=2, jitter(cars$dist, 20) )
  title(main = "cars data")

  par(mar=c(1,2,5,1))
  boxplot(cars$dist, axes=F)
  title(ylab='Stopping distance (ft)', line=0)

  par(mar=c(5,1,1,2))
  boxplot(cars$speed, horizontal=T, axes=F)
  title(xlab='Speed (mph)', line=1)

  par(op)
dev.off()

png(file="g135.png", width=600, height=600)
  plot(dist ~ speed, data = cars,
       main = "\"cars\" data and regression line")
  abline(lm( dist ~ speed, data = cars), 
         col = 'red')
dev.off()

png(file="g136.png", width=600, height=600)
  plot(cars, 
       xlab = "Speed (mph)", 
       ylab = "Stopping distance (ft)",
       las = 1)
  # lines(loess(dist ~ speed, data=cars), 
  #       col = "red") # Didn't that use to work?
  r <- loess(dist ~ speed, data=cars)
  lines(r$x, r$fitted, col="red")
  title(main = "\"cars\" data and loess curve")
dev.off()

png(file="g137.png", width=600, height=600)
  plot(cars, 
       xlab = "Speed (mph)", 
       ylab = "Stopping distance (ft)",
       las = 1)
  lines(lowess(cars$speed, cars$dist, 
               f = 2/3, iter = 3), 
        col = "red")
  title(main = "\"cars\" data and lowess curve")
dev.off()

png(file="g138.png", width=600, height=600)
  x <- c(15, 9, 75, 90, 1, 1, 11, 5, 9, 8, 33, 11, 11, 
         20, 14, 13, 10, 28, 33, 21, 24, 25, 11, 33)
  # I tried to produce the same with the "stars" 
  # function, with no success.
  clock.plot <- function (x, col = rainbow(n), ...) {
    if( min(x)<0 ) x <- x - min(x)
    if( max(x)>1 ) x <- x/max(x)
    n <- length(x)
    if(is.null(names(x))) names(x) <- 0:(n-1)
    m <- 1.05
    plot(0, 
         type = 'n', # do not plot anything
         xlim = c(-m,m), ylim = c(-m,m), 
         axes = F, xlab = '', ylab = '', ...)
    a <- pi/2 - 2*pi/200*0:200
    polygon( cos(a), sin(a) )
    v <- .02
    a <- pi/2 - 2*pi/n*0:n
    segments( (1+v)*cos(a), (1+v)*sin(a), 
              (1-v)*cos(a), (1-v)*sin(a) )
    segments( cos(a), sin(a), 
              0, 0, 
              col = 'light grey', lty = 3) 
    ca <- -2*pi/n*(0:50)/50
    for (i in 1:n) {
      a <- pi/2 - 2*pi/n*(i-1)
      b <- pi/2 - 2*pi/n*i
      polygon( c(0, x[i]*cos(a+ca), 0),
               c(0, x[i]*sin(a+ca), 0),
               col=col[i] )
      v <- .1
      text((1+v)*cos(a), (1+v)*sin(a), names(x)[i])
    }
  }
  clock.plot(x, 
    main = "Number of visitors to a web site for each hour of the day")
dev.off()

png(file="g139.png", width=600, height=600)
  library(plotrix)
  clock24.plot(x,
               line.col = "blue", 
               lwd = 10)
  # See also polar.plot, radial.plot
dev.off()

png(file="g140.png", width=600, height=600)
  library(circular)
  rose.diag(x)
  # x <- as.circular(rep( 2*pi / 24 * (0:23), x ))
  detach("package:circular")   # redefines "var"...
dev.off()

png(file="g141.png", width=600, height=600)
  # Polar plot to spot seasonal patterns
  x <- as.vector(UKgas)
  n <- length(x)
  theta <- seq(0, by=2*pi/4, length=n)
  plot(x * cos(theta), x * sin(theta),
       type = "l",
       xlab = "", ylab = "",
       main = "UK gas consumption")
  abline(h=0, v=0, col="grey")
  abline(0, 1, col="grey")
  abline(0, -1, col="grey")
  circle <- function (x, y, r, N=100, ...) {
    theta <- seq(0, 2*pi, length=N+1)
    lines(x + r * cos(theta), y + r * sin(theta), ...)
  }
  circle(0,0, 250, col="grey")
  circle(0,0, 500, col="grey")
  circle(0,0, 750, col="grey")
  circle(0,0, 1000, col="grey")
  circle(0,0, 1250, col="grey")
  segments( x[-n] * cos(theta[-n]),
            x[-n] * sin(theta[-n]),
            x[-1] * cos(theta[-1]),
            x[-1] * sin(theta[-1]),
            col = terrain.colors(length(x)),
            lwd = 3)
  text(par("usr")[2], 0, "Winter", adj=c(1,0))
  text(0, par("usr")[4], "Spring", adj=c(0,1))
  text(par("usr")[1], 0, "Summer", adj=c(0,0))
  text(0, par("usr")[3], "Autumn", adj=c(0,0))
  legend("topright", legend = c(1960, 1973, 1986),
         fill = terrain.colors(3))
dev.off()

png(file="g142.png", width=600, height=600)
  conformal_plot <- function (x, y, ...) {
    # To be used when y is thought to be a periodic function of x, 
    # with period 2pi.
    z <- y + 1i * x
    z <- exp(z)
    x <- Re(z)
    y <- Im(z)
    plot(x, y, ...)
  }
  conformal_abline <- function (h=NULL, v=NULL, a=NULL, b=NULL, ...) {
    if (!is.null(a) | ! is.null(b)) {
      stop("Do not set a or b but only h or v")
    }
    if (!is.null(h)) {
      theta <- seq(0, 2*pi, length=200)
      for (i in 1:length(h)) {
        rho <- exp(h[i])
        lines(rho * cos(theta), rho * sin(theta), type = "l", ...)
      }
    }
    if (!is.null(v)) {
      rho <- sqrt(2) * max(abs(par("usr")))
      segments(0, 0, rho * cos(v), rho * sin(v), ...)
    }
  }

  op <- par(mar=c(1,1,3,1))
  x <- as.vector(sunspots)
  conformal_plot(2 * pi * seq(from=0, by=1/(11*12), length=length(x)),
                 x / 100, 
                 type = "l",
                 lwd = 2,
                 col = "blue",
                 xlab = "", ylab = "",
                 main = "Sunspots after conformal transformation")
  conformal_abline(h=seq(0,3, by=.25), col="grey")
  conformal_abline(v = seq(0, 2*pi, length=12), ## 11 years...
                   col = "grey")
  par(op)
dev.off()

png(file="g143.png", width=600, height=600)
  library(lattice)
  y <- cars$dist
  x <- cars$speed
  # vitesse <- shingle(x, co.intervals(x, number=6))
  vitesse <- equal.count(x)
  histogram(~ y | vitesse)
dev.off()

png(file="g144.png", width=600, height=600)
  bwplot(~ y | vitesse, layout=c(1,6))
dev.off()

png(file="g145.png", width=600, height=600)
  densityplot(~ y | vitesse, aspect='xy')
dev.off()

png(file="g146.png", width=600, height=600)
  densityplot(~ y | vitesse, layout=c(1,6))
dev.off()

png(file="g147.png", width=600, height=300)
  y <- cars$dist
  x <- cars$speed
  q <- quantile(x)
  o1 <- x<q[2]
  o2 <- q[2]<x & x<q[3]
  o3 <- q[3]<x & x<q[4]
  o4 <- q[4]<x 
  dotchart(c(median(y[o1]), median(y[o2]), 
             median(y[o3]), median(y[o4])),
           labels = as.character(1:4),
           xlab = "speed", ylab = "distance",
           main = "Before I knew lattice plots")
dev.off()

png(file="g148.png", width=600, height=600)
  my.dotchart <- function (y,x,...) {
    x <- as.matrix(x)
    z <- NULL
    for (i in 1:dim(x)[2]) {
      q <- quantile(x[,i])
      for (j in 1:4) {
        z <- append(z, median(y[ 
          q[j] <= x[,i] & x[,i] <= q[j+1] 
        ]))
        names(z)[length(z)] <- 
          paste(colnames(x)[i], as.character(j))
      }
    }
    dotchart(z, ...)
  }
  my.dotchart(y, x, xlab = "speed", ylab = "distance",
           main = "Before I knew lattice plots")
dev.off()

png(file="g149.png", width=600, height=600)
  my.dotchart <- function (y,x,...) {
    x <- as.matrix(x)
    z <- NULL
    for (i in 1:dim(x)[2]) {
      q <- quantile(x[,i])
      for (j in 1:4) {
        ya <- y[ q[j] <= x[,i] & x[,i] <= q[j+1] ]
        z <- rbind(z, quantile(ya))
        rownames(z)[dim(z)[1]] <- 
          paste(colnames(x)[i], as.character(j))
      }
    }
    dotchart(t(z), ...)
  }
  my.dotchart(y, x, xlab = "speed", ylab = "distance",
           main = "Before I knew lattice plots")
dev.off()

png(file="g150.png", width=600, height=300)
  my.dotchart <- function (y,x,...) {
    x <- as.matrix(x)
    z <- NULL
    for (i in 1:dim(x)[2]) {
      q <- quantile(x[,i])
      for (j in 1:4) {
        ya <- y[ q[j] <= x[,i] & x[,i] <= q[j+1] ]
        z <- rbind(z, quantile(ya))
        rownames(z)[dim(z)[1]] <- 
          paste(colnames(x)[i], as.character(j))
      }
    }
    xmax <- max(z)
    xmin <- min(z)
    n <- dim(z)[1]
    plot( z[,3], 1:n, xlim = c(xmin,xmax), ylim = c(1,n), 
          axes=F, frame.plot = TRUE, pch = '.', 
          ... )
    axis(1)
    axis(2, at=1:n, las=1)
    abline( h=1:n, lty=3 )
    # median
    points( z[,3], 1:n, pch=16, cex=3 )
    # quartiles
    segments( z[,2], 1:n, z[,4], 1:n, lwd=7 )
    # min and max
    segments( z[,1], 1:n, z[,5], 1:n )
  }
  my.dotchart(y,x, xlab="speed", ylab="distance",
              main = "Before I knew lattice plots")
dev.off()

png(file="g151.png", width=600, height=600)
  plot(cars)
  polygon( cars[chull(cars),], col="pink", lwd=3)
  points(cars)
dev.off()

png(file="g152.png", width=600, height=600)
  draw.ellipse <- function (
      x, y = NULL, 
      N = 100,
      method = lines, 
      ...
    ) {
    if (is.null(y)) {
      y <- x[,2]
      x <- x[,1]
    }
    centre <- c(mean(x), mean(y))
    m <- matrix(c(var(x),cov(x,y),
                  cov(x,y),var(y)),
                nr=2,nc=2)
    e <- eigen(m)
    r <- sqrt(e$values)
    v <- e$vectors
    theta <- seq(0,2*pi, length=N)
    x <- centre[1] + r[1]*v[1,1]*cos(theta) +
         r[2]*v[1,2]*sin(theta)
    y <- centre[2] + r[1]*v[2,1]*cos(theta) +
         r[2]*v[2,2]*sin(theta)
    method(x,y,...)
  }
  plot(cars)
  draw.ellipse(cars, col="blue", lwd=3)
dev.off()

png(file="g153.png", width=600, height=600)
  library(chplot)
  data(hdr)
  x <- hdr$age
  y <- log(hdr$income)
  library(MASS)
  z <- kde2d(x,y, n=50)
  image(z, main = "Density estimation")
dev.off()

png(file="g154.png", width=600, height=600)
  contour(z, 
          col = "red", drawlabels = FALSE,
          main = "Density estimation: contour plot")
dev.off()

png(file="g155.png", width=600, height=600)
  i <- sample(1:length(x), 1000)
  plot(jitter(x[i]), y[i])
  contour(z, 
          col = "red", lwd = 3, drawlabels = FALSE, 
          add = TRUE,
          main = "Density estimation: contour plot")
dev.off()

png(file="g156.png", width=600, height=600)
  persp(z, main = "Density estimation: perspective plot")
dev.off()

png(file="g157.png", width=600, height=600)
  persp(z, 
        phi = 45, theta = 30, 
        xlab = "age", ylab = "income", zlab = "density",
        main = "Density estimation: perspective plot")
dev.off()

png(file="g158.png", width=600, height=600)
  op <- par(mar=c(0,0,2,0)+.1)
  persp(z, phi = 45, theta = 30, 
        xlab = "age", ylab = "income", zlab = "density", 
        col = "yellow", shade = .5, border = NA,
        main = "Density estimation: perspective plot")
  par(op)
dev.off()

png(file="g159.png", width=600, height=600)
  library(chplot)
  data(hdr)
  x <- hdr$age
  y <- log(hdr$income)
  FUNPerFractile <- function (x, y, N, FUN=mean, ...) {
    y <- cut(y, 
             breaks = quantile(y, seq(0, 1, length = N+1), 
                               na.rm=T), 
             labels = FALSE)
    a <- 1:N
    b <- tapply(x, y, FUN, ...)
    data.frame(a,b)
  }
  MeanPerFractilePlot <- function (x, y, N=20, ...) {
    plot(FUNPerFractile(x,y,N=N,FUN=mean,na.rm=T),...)
  }
  MeanPerFractilePlot(x, y, type = "b", lwd = 3, 
                      xlab = "age fractiles", 
                      ylab = "mean income",
                      main = "Mean-per-fractile plot")
dev.off()

png(file="g160.png", width=600, height=600)
  MedianPerFractilePlot <- function (x, y, N=20, ...) {
    plot(FUNPerFractile(x,y,N=N,FUN=median,na.rm=T),...)
  }
  MedianPerFractilePlot(x,y, type="b", lwd=3,
                      xlab = "age fractiles", 
                      ylab = "median income",
                      main = "Median-per-fractile plot")
dev.off()

png(file="g161.png", width=600, height=600)
  r <- loess(y~x)
  o <- order(x)
  plot( r$x[o], r$fitted[o], type = "l",
        xlab = "age", ylab = "income",
        main = "loess curve, without the points" )
dev.off()

png(file="g162.png", width=600, height=600)
  r <- loess(y~x)
  o <- order(x)
  plot( jitter(x, amount = .5), y, pch = ".",
        xlab = "age", ylab = "income",
       main = "Loess curve")
  lines(r$x[o], r$fitted[o], col="blue", lwd=3)
  r <- kde2d(x,y)
  contour(r, drawlabels=F, col="red", lwd=3, add=T)
dev.off()

png(file="g163.png", width=600, height=600)
  data(InsectSprays)
  boxplot(count ~ spray, 
          data = InsectSprays,
          xlab = "Type of spray", 
          ylab = "Insect count",
          main = "InsectSprays data", 
          varwidth = TRUE, 
          col = "lightgray")
dev.off()

png(file="g164.png", width=600, height=400)
  my.dotchart <- function (y,x,...) {
    x <- data.frame(x)
    z <- NULL
    cn <- NULL
    for (i in 1:dim(x)[2]) {
      if( is.numeric(x[,i]) ) {
        q <- quantile(x[,i])
        for (j in 1:4) {
          ya <- y[ q[j] <= x[,i] & x[,i] <= q[j+1] ]
          z <- rbind(z, quantile(ya))
          cn <- append(cn, paste(colnames(x)[i], as.character(j)))
        }
      } else {
        for (j in levels(x[,i])) {
          ya <- y[ x[,i] == j ]
          z <- rbind(z, quantile(ya))
          cn <- append(cn, paste(colnames(x)[i], as.character(j)))
        }
      }
    }
    xmax <- max(z)
    xmin <- min(z)
    n <- dim(z)[1]
    plot( z[,3], 1:n, 
          xlim=c(xmin,xmax), ylim=c(1,n), 
          axes=F, frame.plot=T, pch='.', ... )
    axis(1)
    axis(2, at=1:n, labels=cn, las=1)
    abline( h=1:n, lty=3 )
    # median
    points( z[,3], 1:n, pch=16, cex=3 )
    # quartiles
    segments( z[,2], 1:n, z[,4], 1:n, lwd=7 )
    # min and max
    segments( z[,1], 1:n, z[,5], 1:n )
  }
  spray <- InsectSprays$spray
  y <- InsectSprays$count
  my.dotchart(y,spray, xlab="count", ylab="spray")
dev.off()

png(file="g165.png", width=600, height=600)
  # (This package used to be called UsingR)
  library(UsingR)
  n <- 1000
  k <- 10
  x <- factor(sample(1:5, n, replace=T))
  m <- rnorm(k,sd=2)
  s <- sample( c(rep(1,k-1),2) )
  y <- rnorm(n, m[x], s[x])
  simple.violinplot(y~x, col='pink')
  detach("package:UsingR")
dev.off()

png(file="g166.png", width=600, height=600)
  library(vioplot)
  vioplot(y[x=="1"], y[x=="2"], y[x=="3"], 
          y[x=="4"], y[x=="5"])
  title( main = "vioplot" )
  # The following does not work because the function 
  # wants its first argument to be called "x": it was 
  # defined as function(x,...) instead of function(...).
  # do.call("vioplot", tapply(y, x, function (x) x))
dev.off()

png(file="g167.png", width=600, height=600)
  library(lattice)
  bwplot( y ~ x, 
          panel = panel.violin,
          main = "panel.violin" )
dev.off()

png(file="g168.png", width=600, height=300)
  f <- function (x, N=20) {
    plot.new()
    plot.window( xlim = range(x), 
                 ylim = c(0,1) )
    q <- quantile(x, seq(0, 1, by=1/N))
    segments(q, 0, q, 1)
    lines(c(min(x), min(x), max(x), max(x), min(x)),
            c(0, 1, 1, 0, 0))
  }
  x <- rnorm(1000)
  f(x)
dev.off()

png(file="g169.png", width=600, height=300)
  f <- function (x, N=20) {
    plot.new()
    plot.window( xlim=range(x), ylim=c(-.6,.6) )
    q <- quantile(x, seq(0,1, by=1/N))
    for (i in 1:N) {
      y <- if (i <= N/2) (i-1)/N else (N-i)/N
      lines( c(q[i], q[i], q[i+1], q[i+1], q[i]),
             c(y, -y, -y, y, y) )
    }
  }
  f(x, N=100)
dev.off()

png(file="g170.png", width=600, height=600)
  op <- par(mfrow=c(3,1), mar=c(2,2,3,2))

  n <- length(x)
  plot(sort(x), (1:n)/n, 
       type = "l",
       main = "Cumulative distribution function")

  a <- sort(x)
  b <- (1:length(x))/length(x)
  plot(a, b, 
       type = "l",
       main = "We reverse its second half")
  k <- ceiling(n/2)
  lines(a, c( b[1:(k-1)], (1-b)[k:n] ), 
        col = "blue", lwd = 3)

  plot.new()
  plot.window( xlim=range(x), ylim=c(-.6, .6) )
  lines(a, c( b[1:(k-1)], (1-b)[k:n] ), 
        col="blue", lwd=3,)
  lines(a, -c( b[1:(k-1)], (1-b)[k:n] ), 
        col="blue", lwd=3)
  axis(1)
  title("We symetrize it to get the box-percentile plot")
  abline(h=0, lty=3)
  par(op)
dev.off()

png(file="g171.png", width=600, height=600)
  library(Hmisc)
  bpplot(tapply(y, x, function (x) x))
dev.off()

png(file="g172.png", width=600, height=600)
  library(Hmisc)
  bpplt()    # This is the documentation
  title(main = "bpplt()")
dev.off()

png(file="g173.png", width=600, height=600)
  bwplot( ~ y | x, 
          panel = panel.bpplot, 
          main = "panel.bpplot",
          layout = c(1,5) )
dev.off()

png(file="g174.png", width=600, height=600)
  bpplot( faithful$waiting,
          main = "Box-precentile plot of bimodal data" )
dev.off()

png(file="g175.png", width=400, height=600)
  library(hdrcde)
  hdr.boxplot(rnorm(1000), col = "pink",
              main = "Highest Density Region Plot")
dev.off()

png(file="g176.png", width=400, height=600)
  hdr.boxplot(faithful$waiting, 
              col = "pink",
              main = "Highest Density Region Plot")
dev.off()

png(file="g177.png", width=600, height=600)
  stripchart(InsectSprays$count ~ InsectSprays$spray, 
             method = 'jitter')
dev.off()

png(file="g178.png", width=600, height=600)
  data(iris)
  plot(iris[1:4], 
       pch = 21, 
       bg = c("red", "green", "blue")[
         as.numeric(iris$Species)
       ])
dev.off()

png(file="g179.png", width=600, height=600)
  a <- InsectSprays$count
  b <- rnorm(length(a))
  plot(b ~ a, 
       pch = 21, 
       bg = c("red", "green", "blue", 
            "cyan", "yellow", "black")
         [as.numeric(InsectSprays$spray)],
       main = "1-dimensional scatter plot",
       xlab = "Number of insects",
       ylab = "")
dev.off()

png(file="g180.png", width=600, height=600)
  a <- as.vector(t(iris[1]))
  b <- rnorm(length(a))
  plot(b ~ a, 
       pch = 21, 
       bg = c("red", "green", "blue")[
         as.numeric(iris$Species)
       ],
       main = "1-dimensional scatter plot",
       xlab = "Number of insects",
       ylab = "")
dev.off()

png(file="g181.png", width=600, height=600)
  do.it <- function (v, ...) {
    n <- 100
    y <- sample( 1:3, n, replace=T )
    a <- runif(1)
    b <- runif(1)
    c <- runif(1)
    x <- ifelse( y==1, a+v*rnorm(n), 
         ifelse( y==2, b+v*rnorm(n), c+v*rnorm(n) ))
    r <- rnorm(n)
    plot( r ~ x, 
          pch = 21, 
          bg = c('red', 'green', 'blue')[y],
          ... )
  }
  do.it(.1, main = "1-dimensional scatterplot")
dev.off()

png(file="g182.png", width=600, height=600)
  do.it(.05, main = "1-dimensional scatterplot")
dev.off()

png(file="g183.png", width=400, height=800)
  hists <- function (x, y, ...) {
    y <- factor(y)
    n <- length(levels(y))
    op <- par( mfcol=c(n,1), mar=c(2,4,1,1) )    
    b <- hist(x, ..., plot=F)$breaks
    for (l in levels(y)){
      hist(x[y==l], breaks=b, probability=T, ylim=c(0,.3), 
           main="", ylab=l, col='lightblue', xlab="", ...)
      points(density(x[y==l]), type='l', lwd=3, col='red')
    }
    par(op)
  }
  hists(InsectSprays$count, InsectSprays$spray)
dev.off()

png(file="g184.png", width=600, height=600)
  library(lattice)
  histogram( ~ count | spray, data=InsectSprays)
dev.off()

png(file="g185.png", width=600, height=600)
  densityplot( ~ count | spray, data = InsectSprays )
dev.off()

png(file="g186.png", width=600, height=600)
  bwplot( ~ count | spray, data = InsectSprays )
dev.off()

png(file="g187.png", width=600, height=600)
  bwplot( ~ count | spray, data = InsectSprays, layout=c(1,6) )
dev.off()

png(file="g188.png", width=600, height=600)
  data(HairEyeColor)
  a <- as.table( apply(HairEyeColor, c(1,2), sum) )
  barplot(a, legend.text = attr(a, "dimnames")$Hair)
dev.off()

png(file="g189.png", width=600, height=600)
  barplot(a, 
          beside = TRUE, 
          legend.text = attr(a, "dimnames")$Hair)
dev.off()

png(file="g190.png", width=600, height=600)
  barplot(t(a), 
          legend.text = attr(a, "dimnames")$Eye)
dev.off()

png(file="g191.png", width=600, height=600)
  barplot(t(a), 
          beside = TRUE, 
          legend.text = attr(a, "dimnames")$Eye)
dev.off()

png(file="g192.png", width=600, height=600)
  b <- a / apply(a, 1, sum)
  barplot(t(b))
dev.off()

png(file="g193.png", width=600, height=600)
  c <- t( t(a) / apply(a, 2, sum) )
  barplot(c)
dev.off()

png(file="g194.png", width=600, height=600)
  plot(a, main = "Mosaic plot")
dev.off()

png(file="g195.png", width=600, height=600)
  plot(t(a), main = "Mosaic plot")
dev.off()

png(file="g196.png", width=600, height=600)
  plot(a, 
       col = heat.colors(dim(a)[2]),
       main = "Mosaic plot")
dev.off()

png(file="g197.png", width=600, height=600)
  plot(a, 
       color = TRUE,
       main = "Mosaic plot")
dev.off()

png(file="g198.png", width=600, height=600)
  plot(a, 
       shade = TRUE,
       main = "Mosaic plot")
dev.off()

png(file="g199.png", width=600, height=600)
  plot(t(a), 
       shade = TRUE,
       main = "Mosaic plot")
dev.off()

png(file="g200.png", width=600, height=600)
  data(HairEyeColor)
  a <- apply(HairEyeColor, c(1,2) , sum)
  qualplot <- function (a) {
    matplot( row(a), a, 
             type = 'l', axes = FALSE,
             col = 1:dim(a)[2]+1,
             lty = 1:dim(a)[2],
             lwd=3,
             xlab = names(dimnames(a))[1], 
             ylab = names(dimnames(a))[2] )
    axis(1, 1:dim(a)[1], row.names(a))
    axis(2)
    legend(1, max(a), row.names(t(a)),
           lwd = 3, cex = 1.5,
           col = 1:dim(a)[2]+1, 
           lty = 1:dim(a)[2])
  }
  # For interactive use
  qualplots <- function (a) {
    op <- par(ask=TRUE)
    qualplot(a)
    qualplot(t(a))
    par(op)
  }
  qualplot(a)
dev.off()

png(file="g201.png", width=600, height=600)
  qualplot(t(a))
dev.off()

png(file="g202.png", width=600, height=600)
  qualplotfreq <- function (a) {
    a <- t( t(a) / apply(a,2,sum) )
    qualplot(a)
  }
  qualplotsfreq <- function (a) {
    op <- par(ask=TRUE)
    qualplotfreq(a)
    qualplotfreq(t(a))
    par(op)
  }
  qualplotfreq(a)
dev.off()

png(file="g203.png", width=600, height=600)
  qualplotfreq(t(a))
dev.off()

png(file="g204.png", width=600, height=600)
  data(bacteria, package="MASS")
  fourfoldplot( table(bacteria$y, bacteria$ap) )
dev.off()

png(file="g205.png", width=600, height=600)
  fourfoldplot( table(bacteria$y, 
                      bacteria$ap, 
                      bacteria$week) )
dev.off()

png(file="g206.png", width=600, height=600)
  n <- 50
  x <- rnorm(n)
  y <- rnorm(n)
  z <- rnorm(n)
  my.renorm <- function (z) {
    z <- abs(z)
    z <- 10*z/max(z)
    z
  }
  z <- my.renorm(z)
  op <- par(mar = c(3,2,4,2)+.1)
  plot(x, y, cex = z,
       xlab = "", ylab = "", 
       main = "Bubble plot")
dev.off()

png(file="g207.png", width=600, height=600)
  plot(x, y, cex = z, 
       pch = 16, col = 'red',
       xlab = "", ylab = "", 
       main = "Bubble plot")
  points(x, y, cex = z)
dev.off()

png(file="g208.png", width=600, height=600)
  u <- sample(c('red','green','blue'),n,replace=T)
  plot(x, y, cex = z, col = u,
       pch = 16,
       xlab = "", ylab = "", 
       main = "Bubble plot")
  points(x, y, cex = z)
dev.off()

png(file="g209.png", width=600, height=600)
  z2 <- rnorm(n)
  z2 <- my.renorm(z2)
  plot(x, y, cex = z,
       xlab = "", ylab = "", 
       main = "Bubble plot")
  points(x, y, cex = z2, col = 'red')
dev.off()

png(file="g210.png", width=600, height=600)
  # Other renormalization (if there is no zero)
  my.renorm <- function (z) { 
    z <- (z-min(z)) / (max(z)-min(z))
    z <- 1+9*z
    z
  }
  z <- my.renorm(z)
  z2 <- my.renorm(z2)
  plot(x, y, cex = z,
       xlab = "", ylab = "", 
       main = "Bubble plot")
  points(x, y, cex = z2, col = 'red')
dev.off()

png(file="g211.png", width=600, height=600)
  n <- 50
  x <- runif(n)
  y <- runif(n)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  z3 <- rnorm(n)
  z4 <- rnorm(n)
  z5 <- rnorm(n)
  stars( data.frame(z1,z2,z3,z4,z5), location=cbind(x,y), 
         labels=NULL, len=1/sqrt(n)/2,
         main = "Star plot" )
dev.off()

png(file="g212.png", width=600, height=600)
  v <- .2
  n <- 50
  x <- runif(n)
  y <- runif(n)
  z1 <- x+y+v*rnorm(n)
  z2 <- x*y+v*rnorm(n)
  z3 <- x^2 + y^2 + v*rlnorm(n)
  stars( data.frame(z1,z2,z3), 
         location = cbind(x,y), 
         labels = NULL, 
         len = 1/sqrt(n)/2, 
         axes = TRUE,
         draw.segments = TRUE, 
         col.segments = 1:5,
         main = "Star plot" )
dev.off()

png(file="g213.png", width=600, height=600)
  n <- 10
  d <- data.frame(y1 = abs(rnorm(n)),
                  y2 = abs(rnorm(n)),
                  y3 = abs(rnorm(n)),
                  y4 = abs(rnorm(n)),
                  y5 = abs(rnorm(n))
                 )
  matplot(d, 
          type = 'l',
          ylab = "", 
          main = "Matplot")
dev.off()

png(file="g214.png", width=600, height=600)
  barplot(t(as.matrix(d)))
dev.off()

png(file="g215.png", width=600, height=600)
  line.chart <- function (d, 
                          xlab = "", ylab = "", 
                          main = "") {
    m <- d
    m <- t(apply(m,1,cumsum))
    #print(m)
    n1 <- dim(m)[1]
    n2 <- dim(m)[2]
    col <- rainbow(n)
    plot.new()
    plot.window(xlim = c(1, n1), 
                ylim = c(min(m), max(m)))
    axis(1)
    axis(2)
    title(xlab = xlab, ylab = ylab,
          main = main)
    for (i in n2:1) {
      polygon(c(1:n1,n1,1), c(m[,i],0,0), 
              col = col[i], 
              border = 0)
    }
    for (i in n2:1) {
      lines(m[,i], lwd = 2)
    }
  }
  line.chart(d, main = "Linechart")
dev.off()

png(file="g216.png", width=600, height=600)
  data(LifeCycleSavings)
  plot(LifeCycleSavings)
dev.off()

png(file="g217.png", width=600, height=600)
  panel.hist <- function(x, ...) {
    usr <- par("usr"); 
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, ...)
  }
  # Correlation coefficient
  my.panel.smooth <- 
  function (x, y, 
            col = par("col"), 
            bg = NA, 
            pch = par("pch"),
            cex = 1, 
            col.smooth = "red", 
            span = 2/3, 
            iter = 3, ...) {
    points(x, y, 
           pch = pch, col = col, 
           bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok))
      lines(lowess(x[ok], y[ok], 
                   f = span, 
                   iter = iter), 
            col = col.smooth,
            ...)
    usr <- par('usr')
    text( (usr[1]+usr[2])/2, (usr[3]+9*usr[4])/10, 
          floor(100*cor(x,y))/100,
          col='blue', cex=3, adj=c(.5,1) )
  }
  pairs(LifeCycleSavings,
        diag.panel  = panel.hist,
        upper.panel = panel.smooth,
        lower.panel = my.panel.smooth,
        gap = 0)
dev.off()

png(file="g218.png", width=600, height=600)
  cor.plot <- function (x, 
                        xlab = "", ylab = "", 
                        main = "") {
    n <- dim(x)[1]
    m <- dim(x)[2]
    N <- 1000
    col = topo.colors(N)
    plot(NA, xlim = c(0,1.2), ylim = c(-1,0),
         xlab = xlab, ylab = ylab, main = main)
    for (i in 1:n) {
      for (j in 1:m) {
        polygon( c((j-1)/m, (j-1)/m, j/m, j/m),
                 -c((i-1)/m, i/m, i/m, (i-1)/m),
                 col = col[ N*(x[i,j]+1)/2 ] )
      }
    }
    for (i in 1:N) {
      polygon( c(1.1, 1.1, 1.2, 1.2),
               -c((i-1)/N, i/N, i/N, (i-1)/N),
               col = col[N-i+1],
               border = NA )
    }
    # Exercice: add a legend
  }

  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x <- rnorm(n)
  x4 <- x + .1*rnorm(n)
  x5 <- x + .1*rnorm(n)
  y <- 1 + x1 + x4 + rnorm(n)
  d <- data.frame(y,x1,x2,x3,x4,x5)
  op <- par(mar=c(3,3,4,2)+.1)
  cor.plot(cor(d), main = "Correlation plot")
  par(op)
dev.off()

png(file="g219.png", width=600, height=600)
  library(sma)
  plot.cor(cor(d),
           labels = colnames(d),
           main = "plot.cor (in the \"sma\" package)")
dev.off()

png(file="g220.png", width=600, height=600)
  library(ellipse)
  plotcorr(cor(d), main = "plotcorr (in the \"ellipse\" package)")
dev.off()

png(file="g221.png", width=600, height=600)
  uniformize <- function (x) {
    x <- rank(x, na.last="keep")
    x <- (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
    x
  }
  scagnostic_contour <- function (x, y, ..., FUN = median) {
    x <- uniformize(x)
    y <- uniformize(y)
    require(MASS) # For kde2d()
    r <- kde2d(x, y, ...)
    r$z > FUN(r$z)
  }
  translate <- function (x,i,j,zero=0) {
    n <- dim(x)[1]
    m <- dim(x)[2]
    while (i>0) {
      x <- rbind( rep(zero,m), x[-n,] )
      i <- i - 1
    }
    while (i<0) {
      x <- rbind( x[-1,], rep(zero,m) )
      i <- i + 1
    }
    while (j>0) {
      x <- cbind( rep(zero,n), x[,-m] )
      j <- j - 1
    }
    while (j<0) {
      x <- cbind( x[,-1], rep(zero,n) )
      j <- j + 1
    }
    x
  }
  scagnostic_perimeter <- function (x, y, ...) {
    z <- scagnostic_contour(x, y, ...)
    zz <- z | 
          translate(z,1,0)  | translate(z,0,1)  | 
          translate(z,-1,0) | translate(z,0,-1)
    sum(zz & ! z)
  }
  scagnostic_area <- function (x, y, ...) {
    z <- scagnostic_contour(x, y, ..., FUN = mean)
    sum(z) / length(z)
  }
  connected_components <- function (x) {
    stopifnot(is.matrix(x), is.logical(x))
    m <- dim(x)[1]
    n <- dim(x)[2]
    x <- rbind( rep(FALSE, n+2),
                cbind( rep(FALSE, m), x, rep(FALSE, m) ),
                rep(FALSE, n+2))
    x[ is.na(x) ] <- FALSE
    # Assign a label to each pixel, so that pixels with the same 
    # label be in the same connected component -- but pixels in the 
    # same connected component may have different labels.
    current_label <- 0
    result <- ifelse(x, 0, 0)
    equivalences <- list()
    for (i in 1 + 1:m) {
      for (j in 1 + 1:n) {
        if (x[i,j]) {
          number_of_neighbours <- x[i-1,j-1] + x[i-1,j] + x[i-1,j+1] + x[i,j-1]
          labels <- c( result[i-1,j-1], result[i-1,j],
                       result[i-1,j+1], result[i,j-1] )
          labels <- unique(labels[ labels > 0 ])
          neighbour_label <- max(0,labels)
          if (number_of_neighbours == 0) {
            current_label <- current_label + 1
            result[i,j] <- current_label
          } else if (length(labels) == 1) {
            result[i,j] <- neighbour_label
          } else {
            result[i,j] <- neighbour_label
            equivalences <- append(equivalences, list(labels))
          }
        }
      }
    }
    # Build the matrix containing the equivalences between those labels
    # We just have the matrix of a (non-equivalence) relation: we compute
    # the equivalence relation it generates.
    E <- matrix(FALSE, nr=current_label, nc=current_label)
    for (e in equivalences) {
      stopifnot( length(e) > 1 )
      for (i in e) {
        for (j in e) {
          if (i != j) {
            E[i,j] <- TRUE
          }
        }
      }
    }
    E <- E | t(E)
    diag(E) <- TRUE
    for (k in 1:current_label) {
      E <- E | (E %*% E > 0)
    }
    stopifnot( E == E | (E %*% E > 0) )
    # Find the equivalence classes, i.e., the unique rows of this matrix
    E <- apply(E, 2, function (x) min(which(x)))
    # Finally, label the equivalence classes
    for (i in 1:current_label) {
      result[ result == i ] <- E[i]
    }
    result
  }
  connected_components_TEST <- function () {
    n <- 100
    x <- matrix(NA, nr=n, nc=n)
    x <- abs(col(x) - (n/3)) < n/8 & abs(row(x) - n/3) < n/8
    x <- x | ( (col(x) - 2*n/3)^2 + (row(x) - 2*n/3)^2 < (n/8)^2 )
    image(!x)
    image(-connected_components(x))    
  }
  scagnostic_modality <- function (x, y, ...) {
    z <- scagnostic_contour(x, y, ...)
    z <- connected_components(z)
    max(z)
  }
  scagnostic_slope <- function (x,y) {
    x <- uniformize(x)
    y <- uniformize(y)
    pc1 <- prcomp(cbind(x,y))$rotation[,1]
    pc1[2] / pc1[1]
  }
  scagnostic_sphericity <- function (x,y) {
    x <- uniformize(x)
    y <- uniformize(y)
    # Ratio of the eigenvalues of the PCA
    # For a spherical cloud of points, the slope
    # is not well defined, but this ratio is close to 1.
    ev <- prcomp(cbind(x,y))$sdev
    ev[1] / ev[2]
  }
  scagnostic_curvature <- function (x,y) {
    x <- uniformize(x)
    y <- uniformize(y)
    require(pcurve)
    # BUG: pcurve() starts a new plot by fiddling with par() --
    # it also fails to set it back to what it was...
    par <- function (...) { }
    r <- NULL
    try(
    r <- pcurve(cbind(x,y),
                start = "pca", # Defaults to CA,
                               # which only works with count data...
                plot.pca = FALSE,
                plot.true = FALSE,
                plot.init = FALSE, 
                plot.segs = FALSE, 
                plot.resp = FALSE,
                plot.cov = FALSE,
                use.loc = FALSE)
    )
    if (is.null(r)) return(0)
    X <- r$s[,1:2]               # The principal curve
    n <- dim(X)[1]
    V <- X[2:n,] - X[1:(n-1),]
    V <- V / sqrt(V[,1]^2 + V[,2]^2) # The direction of the principal 
                                     # curve, at each point on it
    C <- apply( V[1:(n-2),] * V[2:(n-1),], 1, sum )
    C <- acos(C)                 # The angles
    sum(abs(C)) / pi
  }
  scagnostic_distance <- function (x,y) {
    i <- is.finite(x) & is.finite(y)
    if (length(i) < 2) {
      return(NA)
    }
    x <- uniformize(x)[i]
    y <- uniformize(y)[i]
    d <- as.matrix(dist(cbind(x,y)))
    diag(d) <- Inf
    d <- apply(d, 2, min)   # Nearest neighbour distance
    mean(d)
  }
  scagnostics <- function (
    x, 
    functions = list(
      Perimeter  = scagnostic_perimeter,
      Area       = scagnostic_area,
      Modality   = scagnostic_modality,
      Slope      = scagnostic_slope,
      Sphericity = scagnostic_sphericity,
      Curvature  = scagnostic_curvature,
      "Nearest neighbour distance" = scagnostic_distance
    )
  ) {
    stopifnot( is.matrix(x) || is.data.frame(x) )
    number_of_variables   <- dim(x)[2]
    number_of_scagnostics <- length(functions)
    res <- array(NA, dim=c(number_of_variables,
                           number_of_variables,
                           number_of_scagnostics))
    dimnames(res) <- list(
      Variable1 = colnames(x),
      Variable2 = colnames(x),
      Scagnostic = names(functions)
    )
    for (i in 1:number_of_variables) {
      for (j in 1:number_of_variables) {
        if (i != j) {
          for (k in 1:number_of_scagnostics) {
            res[i,j,k] <- functions[[k]] (x[,i], x[,j])
          }
        }
      }
    }
    class(res) <- "scagnostics"
    res
  }
  plot.scagnostics <- function (x, FUN=pairs, ...) {
    stopifnot(inherits(x, "scagnostics"))
    y <- apply(x, 3, as.vector)
    colnames(y) <- dimnames(x)[[3]]
    rownames(y) <- outer(dimnames(x)[[1]], dimnames(x)[[2]], paste, sep="-")
    FUN(y, ...)
  }

  pairs(USJudgeRatings, gap=0)
dev.off()

png(file="g222.png", width=600, height=600)
  plot(scagnostics(USJudgeRatings), gap=0)
dev.off()

png(file="g223.png", width=600, height=600)
  x <- Harman74.cor[[1]]
  pairs(x, gap=0)
dev.off()

png(file="g224.png", width=600, height=600)
  plot(scagnostics(x), gap=0)
dev.off()

png(file="g225.png", width=600, height=600)
  image(t(as.matrix(USJudgeRatings)))
dev.off()

png(file="g226.png", width=600, height=600)
  # This uses cluster analysis
  heatmap(as.matrix(USJudgeRatings))
dev.off()

png(file="g227.png", width=600, height=600)
  my.dotchart(LifeCycleSavings[,1], LifeCycleSavings[,-1],
              xlab='savings', ylab='')
dev.off()

png(file="g228.png", width=600, height=600)
  to.factor.vector <- function (x, number = 4) {
    resultat <- NULL
    intervalles <- co.intervals(x, number, 
                                overlap = 0)
    for (i in 1:number) {
      if ( i == 1 ) {
        intervalles[i,1] = min(x)
      } else {
        intervalles[i,1] <- intervalles[i-1,2]
      }
      if( i == number ) {
        intervalles[i,2] <- max(x)
      }
    }
    for (valeur in x) {
      r <- NA
      for (i in 1:number) {
        if( valeur >= intervalles[i,1] & 
            valeur <= intervalles[i,2] )
          r <- i
      }
      resultat <- append(resultat, r)
    }
    factor(resultat, levels = 1:number)
  }
  to.factor <- function (x, number = 4) {
    if(is.vector(x)) 
      r <- to.factor.vector(x, number)
    else {
      r <- NULL
      for (v in x) {
        a <- to.factor.vector(v)
        if( is.null(r) ) 
          r <- data.frame(a)
        else 
          r <- data.frame(r,a)
      }
      names(r) <- names(x)
    }
    r
  }
  x <- to.factor(LifeCycleSavings[,-1])
  y <- LifeCycleSavings[,1]
  y <- as.vector(matrix(y, 
                        nr = length(y), 
                        ncol = dim(x)[2]))
  for (i in names(x)) {
    levels(x[[i]]) <- paste(i, levels(x[[i]]))
  }
  col <- gl( dim(x)[2], length(levels(x[,1])), 
             labels = rainbow( dim(x)[2] ))
  col <- as.vector(col)
  x <- factor(as.vector(as.matrix(x)))
  boxplot(y ~ x, 
          horizontal = TRUE, 
          las = 1, 
          col = col,
          main = "Boxplot for each quartile")
dev.off()

png(file="g229.png", width=600, height=300)
  bwplot( ~ LifeCycleSavings[,1] | 
            equal.count(LifeCycleSavings[,2], number=4),
          layout=c(1,4) )
dev.off()

png(file="g230.png", width=600, height=300)
  bwplot( ~ LifeCycleSavings[,1] | 
            equal.count(LifeCycleSavings[,3], number=4),
          layout=c(1,4) )
dev.off()

png(file="g231.png", width=600, height=300)
  bwplot( ~ LifeCycleSavings[,1] | 
            equal.count(LifeCycleSavings[,4], number=4),
          layout=c(1,4) )
dev.off()

png(file="g232.png", width=600, height=300)
  bwplot( ~ LifeCycleSavings[,1] | 
            equal.count(LifeCycleSavings[,5], number=4),
          layout=c(1,4) )
dev.off()

png(file="g233.png", width=600, height=600)
  data(mtcars)
  stars(mtcars[, 1:7], 
        key.loc = c(14, 2),
        main = "Motor Trend Cars : stars(*, full = FALSE)", 
        full = FALSE)
dev.off()

png(file="g234.png", width=600, height=600)
  stars(mtcars[, 1:7], 
        key.loc = c(14, 1.5),
        main = "Motor Trend Cars : full stars()",
        flip.labels = FALSE)
dev.off()

png(file="g235.png", width=600, height=600)
  palette(rainbow(12, s = 0.6, v = 0.75))
  stars(mtcars[, 1:7], 
        len = 0.8, 
        key.loc = c(12, 1.5),
        main = "Motor Trend Cars", 
        draw.segments = TRUE)
dev.off()

png(file="g236.png", width=600, height=600)
  stars(mtcars[, 1:7], 
        locations = c(0,0), 
        radius = FALSE,
        key.loc=c(0,0), 
        main="Motor Trend Cars", 
        lty = 2)
dev.off()

png(file="g237.png", width=600, height=600)
  library(circular)
  rose.diag(mtcars[,5])
dev.off()

png(file="g238.png", width=600, height=600)
  rose.diag(mtcars)
dev.off()

png(file="g239.png", width=600, height=600)
  # From the manual
  x <- seq(-10, 10, length=50)
  y <- x
  f <- function(x,y) {
    r <- sqrt(x^2+y^2)
    10 * sin(r)/r
  }
  z <- outer(x, y, f)
  z[is.na(z)] <- 1
  persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
        shade=.5,
        xlab = "X", ylab = "Y", zlab = "Z")
dev.off()

png(file="g240.png", width=600, height=600)
  # From the manual
  data(volcano)
  z <- 2 * volcano        # Exaggerate the relief
  x <- 10 * (1:nrow(z))   # 10-meter spacing (S to N)
  y <- 10 * (1:ncol(z))   # 10-meter spacing (E to W)
  persp(x, y, z, 
        theta = 120, phi = 15, 
        scale = FALSE, axes = FALSE)
  # See also the other examples in 
  #   demo(persp)
dev.off()

png(file="g241.png", width=600, height=600)
  # From the manual
  data("volcano")
  rx <- range(x <- 10*1:nrow(volcano))
  ry <- range(y <- 10*1:ncol(volcano))
  ry <- ry + c(-1,1) * (diff(rx) - diff(ry))/2
  tcol <- terrain.colors(12)
  op <- par(pty = "s", bg = "lightcyan")
  plot(x = 0, y = 0,
       type = "n", 
       xlim = rx, ylim = ry, 
       xlab = "", ylab = "")
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], 
       col = tcol[8], 
       border = "red")
  contour(x, y, volcano, 
          col = tcol[2], 
          lty = "solid", 
          add = TRUE,
          vfont = c("sans serif", "plain"))
  title("A Topographic Map of Maunga Whau", font = 4)
  abline(h = 200*0:4, v = 200*0:4, 
         col = "lightgray", 
         lty = 2, 
         lwd = 0.1)
  par(op)
dev.off()

png(file="g242.png", width=600, height=600)
  # From the manual
  data(volcano)
  x <- 10*(1:nrow(volcano))
  y <- 10*(1:ncol(volcano))
  image(x, y, volcano, 
        col = terrain.colors(100),
        axes = FALSE,
        xlab = "", ylab = "")
  contour(x, y, volcano, 
          levels = seq(90, 200, by=5), 
          add = TRUE, 
          col = "peru")
  axis(1, at = seq(100, 800, by = 100))
  axis(2, at = seq(100, 600, by = 100))
  box()
  title(main = "Maunga Whau Volcano", font.main = 4)
dev.off()

png(file="g243.png", width=600, height=600)
  data(volcano)
  x <- 10*(1:nrow(volcano))
  x <- rep(x, ncol(volcano))
  y <- 10*(1:ncol(volcano))
  y <- rep(y, each=nrow(volcano))
  z <- as.vector(volcano)
  wireframe( z ~ x * y )
dev.off()

png(file="g244.png", width=600, height=600)
  cloud( z ~ x * y )
dev.off()

png(file="g245.png", width=600, height=600)
  data(iris)
  print(cloud(Sepal.Length ~ Petal.Length * Petal.Width,
              data = iris, cex = .8,
              groups = Species,
              subpanel = panel.superpose,
              main = "Stereo",
              screen = list(z = 20, x = -70, y = 3)),
        split = c(1,1,2,1), more = TRUE)
  print(cloud(Sepal.Length ~ Petal.Length * Petal.Width,
              data = iris, cex = .8,
              groups = Species,
              subpanel = panel.superpose,
              main = "Stereo",
              screen = list(z = 20, x = -70, y = 0)),
        split = c(2,1,2,1))
dev.off()

png(file="g246.png", width=600, height=800)
  z <- matrix(rnorm(24),nr=4)
  library(akima) # non-free
  r <- interp( as.vector(row(z)),
               as.vector(col(z)),
               as.vector(z),
               seq(1, dim(z)[1], length=500),
               seq(1, dim(z)[2], length=500) )
  op <- par(mfrow=c(2,1))
  image(t(z), main="Data to be smoothed or interpolated")
  box()
  image(t(r$z), main="Linear interpolation")
  box()
  par(op)
dev.off()

png(file="g247.png", width=600, height=800)
  library(fields)
  loc <- make.surface.grid(list( seq(1,dim(z)[1],length=500),
                                 seq(1,dim(z)[2],length=500) ))
  r <- interp.surface(
    list(x=1:dim(z)[1], y=1:dim(z)[2], z=z),
    loc
  )
  op <- par(mfrow=c(2,1))
  image.plot(z, main="Raw data")
  image.plot(as.surface(loc,r), main="Linear interpolation")
  par(op)
dev.off()

png(file="g248.png", width=600, height=600)
  # You may not want to interpolate, but rather to smooth
  # (the initial data set need not be on a grid)
  # Also check the Tps() function in the fields package
  library(tgp)
  r <- interp.loess( as.vector(row(z)),
                     as.vector(col(z)),
                     as.vector(z),
                     gridlen = 500 )
  image(t(r$z), main="Loess 2-dimensional smoothing")
dev.off()

png(file="g249.png", width=600, height=600)
  library(fields)
  data(lennon)
  x <- lennon[201:240,201:240]
  op <- par(mfrow=c(2,1))
  image(x, main="image()")
  image.plot(x, main="image.plot()")
  par(op)
dev.off()

png(file="g250.png", width=600, height=600)
  library(RColorBrewer)
  display.brewer.all(type="div")
  title(main="RColorBrewer: diverging palettes (i.e., with a zero)")
dev.off()

png(file="g251.png", width=600, height=600)
  display.brewer.all(type="seq")
  title(main="RColorBrewer: sequential palettes")
dev.off()

png(file="g252.png", width=600, height=600)
  display.brewer.all(type="qual")
  title(main="RColorBrewer: qualitative palettes")
dev.off()

png(file="g253.png", width=600, height=600)
  breaks <- function (x, N) {
    x <- as.vector(x)
    x <- x[ !is.na(x) ]
    if (length(x) == 0) {
      return( rep(NA, N) )
    }
    if (N %% 2 == 0) {
      if (all(x >= 0)) {
        res <- c(rep(0, N/2), seq(0, max(x), length=N/2+1))
      } else if (all(x <= 0)) {
        res <- c(seq(min(x), 0, length=N/2+1), rep(0,N/2))
      } else {
        res <- c(seq(min(x), 0, length=N/2+1),
                 seq(0, max(x), length=N/2+1)[-1])
      }
    } else {
      if (all(x >= 0)) {
        res <- c(rep(0,length=(N+1)/2), seq(0, max(x), length=(N+1)/2))
      } else if (all(x <= 0)) {
        res <- c(seq(min(x), 0, length=(N+1)/2), rep(0, length=(N+1)/2))
      } else {
        res <- c(seq(min(x), 0, length=N+1) [seq(1, N, by=2)],
                 seq(0, max(x), length=N+1) [seq(2, N+1, by=2)])
      }
    }
    stopifnot( length(res) == N+1 )
    stopifnot( res == sort(res) )
    stopifnot( all(x <= max(res)), all(x >= min(res)) )
    res
  }

  breaks(  0:10,  5) == c(0,0,0,    0,5,10)
  breaks(-(0:10), 5) == c(-10,-5,0, 0,0,0)
  breaks(-20:10,  5)  == c(-20, -12, -4, 2, 6, 10)
  breaks(   0:9,  6)  == c(0,0,0, 0, 3,6,9)
  breaks(-(0:9),  6)  == c(-9,-6,-3, 0, 0,0,0)
  breaks( -30:9,  6) == c(-30,-20,-10,0,3,6,9)

  # Example from the "fields" manual
  data(ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
  # Remove the missing values
  i <- !is.na(y)
  y <- y[i]
  x <- x[i,]
  # The residuals of a regression
  r <- Tps(x,y)
  z <- residuals(r)
  # Put those residuals on a regular grid
  # We cannot use interp.surface(): it assumes that the data is regular
  library(tgp)
  r <- interp.loess(x[,1], x[,2], z, gridlen=500)

  # I wanted an example with skewed data: residuals tend to be symetric...
  op <- par(mfrow=c(2,2))
  image(r)
  image.plot(r)
  image.plot(r,
             breaks=breaks(r$z, 9), # Fine for the plot, but no for the legend...
             col=rev(brewer.pal(9, "RdBu")))
  par(op)
dev.off()

png(file="g254.png", width=600, height=600)
  n <- 100
  m <- matrix( rnorm(5*n)+c(1,-1,3,0,2), 
               nr = n, nc = 5, byrow = TRUE )
  matplot(1:5, t(m), type = 'l',
          xlab = "", ylab = "")
  title(main = "Parallel plot: Homogeneous cloud")
dev.off()

png(file="g255.png", width=600, height=600)
  n <- 50
  k <- 2
  m <- matrix( rnorm(5*k*n) + 
                 runif(5*k, min = -10, max = 10), 
               nr = n, nc = 5, byrow = TRUE )
  matplot(1:5, t(m), type = 'l',
          xlab = "", ylab = "")
  title(main = "Parallel plot: two clusters")
dev.off()

png(file="g256.png", width=600, height=600)
  n <- 50
  k <- 5
  m <- matrix( rnorm(5*k*n) + 
                 runif(5*k, min = -10, max = 10), 
               nr = n, nc = 5, byrow = TRUE )
  matplot(1:5, t(m), type = 'l',
          xlab = "", ylab = "")
  title(main = "Parallel plot, 5 clusters")
dev.off()

png(file="g257.png", width=600, height=600)
  matplot(1:5, t(princomp(m)$scores), type = 'l')
  title(main = "Idem, after PCA")
dev.off()

png(file="g258.png", width=600, height=600)
  matplot(1:5, t(m), type = 'l')
  title(main = "Point cloud with less visible clusters")
dev.off()

png(file="g259.png", width=600, height=600)
  library(lattice)
  parallel(as.data.frame(m))
dev.off()

png(file="g260.png", width=600, height=600)
  polar_parallel_plot <- function (d, col = par("fg"),
                                   type = "l", lty = 1, ...) {
    d <- as.matrix(d)
    d <- apply(d, 2, function (x) .5 + (x - min(x)) / (max(x) - min(x)))
    theta <- (col(d) - 1) / ncol(d) * 2 * pi
    d <- cbind(d, d[,1])
    theta <- cbind(theta, theta[,1])
    matplot( t(d * cos(theta)), 
             t(d * sin(theta)),
             col = col, 
             type = type, lty = lty, ...,
             axes = FALSE, xlab = "", ylab = "" )
    segments(rep(0,ncol(theta)), rep(0, ncol(theta)),
             1.5 * cos(theta[1,]), 1.5 * sin(theta[1,]))
    if (! is.null(colnames(d))) {
      text(1.5 * cos(theta[1,-ncol(theta)]),
           1.5 * sin(theta[1,-ncol(theta)]),
           colnames(d)[-ncol(d)])
    }
  }
  op <- par(mar=c(0,0,0,0))
  polar_parallel_plot(iris[1:4], col = as.numeric(iris$Species))
  par(op)
dev.off()

png(file="g261.png", width=600, height=600)
  parallel(~iris[1:4], groups = Species, iris)
dev.off()

png(file="g262.png", width=600, height=600)
  parallel(~iris[c(2,4,1,3)], groups= Species, iris)
dev.off()

png(file="g263.png", width=600, height=600)
  x <- seq(-pi, pi, length=100)
  y <- apply(as.matrix(iris[,1:4]),
             1,
             function (u) u[1] + u[2] * cos(x) +
                                 u[3] * sin(x) + u[4] * cos(2*x))
  matplot(x, y,
          type = "l",
          lty = 1,
          col = as.numeric(iris[,5]),
          xlab = "", ylab = "",
          main = "Fourier (Andrew) curves")
dev.off()

png(file="g264.png", width=600, height=600)
  matplot(y * cos(x), y * sin(x),
          type = "l",
          lty = 1,
          col = as.numeric(iris[,5]),
          xlab = "", ylab = "",
          main = "Fourier blob")
dev.off()

png(file="g265.png", width=600, height=600)
  library(TeachingDemos)
  faces(longley[1:9,], main="Macro-economic data")
dev.off()

png(file="g266.png", width=600, height=600)
  library(MASS)
  data(Skye)

  ternary <- function(X, pch = par("pch"), lcex = 1,
                      add = FALSE, ord = 1:3, ...)
  {
    X <- as.matrix(X)
    if(any(X) < 0) stop("X must be non-negative")
    s <- drop(X %*% rep(1, ncol(X)))
    if(any(s<=0)) stop("each row of X must have a positive sum")
    if(max(abs(s-1)) > 1e-6) {
      warning("row(s) of X will be rescaled")
      X <- X / s
    }
    X <- X[, ord]
    s3 <- sqrt(1/3)
    if(!add)
    {
      oldpty <- par("pty")
      on.exit(par(pty=oldpty))
      par(pty="s")
      plot(c(-s3, s3), c(0.5-s3, 0.5+s3), type="n", axes=FALSE,
           xlab="", ylab="")
      polygon(c(0, -s3, s3), c(1, 0, 0), density=0)
      lab <- NULL
      if(!is.null(dn <- dimnames(X))) lab <- dn[[2]]
      if(length(lab) < 3) lab <- as.character(1:3)
      eps <- 0.05 * lcex
      text(c(0, s3+eps*0.7, -s3-eps*0.7),
           c(1+eps, -0.1*eps, -0.1*eps), lab, cex=lcex)
    }
    points((X[,2] - X[,3])*s3, X[,1], ...)
  }

  ternary(Skye/100, ord=c(1,3,2))
dev.off()

png(file="g267.png", width=600, height=600)
  tri <-
  function(a, f, m, symb = 2, grid = F, ...)
  {
    ta <- paste(substitute(a))
    tf <- paste(substitute(f))
    tm <- paste(substitute(m))

    tot <- 100/(a + f +m)
    b <- f * tot
    y <- b * .878
    x <- m * tot + b/2
    par(pty = "s")
    oldcol <- par("col")
    plot(x, y, axes = F, xlab = "", ylab = "", 
         xlim = c(-10, 110), ylim= c(-10, 110), type = "n", ...)
    points(x,y,pch=symb)
    par(col = oldcol)
    trigrid(grid)
    text(-5, -5, ta)
    text(105, -5, tm)
    text(50, 93, tf)
    par(pty = "m")
    invisible()
  }

  trigrid  <-
  function(grid = F)
  {
    lines(c(0, 50, 100, 0), c(0, 87.8, 0, 0)) #draw frame
    if(!grid) {
      for(i in 1:4 * 20) {
        lines(c(i, i - 1), c(0, 2 * .878)) #side a-c (base)
        lines(c(i, i + 1), c(0, 2 * .878))
        T.j <- i/2 #side a-b (left)
        lines(c(T.j, T.j + 2), c(i * .878, i * .878))
        lines(c(T.j, T.j + 1), c(i * .878, (i - 2) * .878))
        T.j <- 100 - i/2 #side b-c (right)
        lines(c(T.j, T.j - 2), c(i * .878, i * .878))
        lines(c(T.j, T.j - 1), c(i * .878, (i - 2) * .878))
      }
    } else {
      for(i in 1:4 * 20) {
        # draw dotted grid
        lines(c(i, i/2), c(0, i * .878), lty = 4, col = 3)
        lines(c(i, (50 + i/2)), c(0, .878 * (100 - i)), lty = 4, col = 3)
        lines(c(i/2, (100 - i/2)), c(i * .878, i * .878), lty = 4, col = 3)
      }
      par(lty = 1, col = 1)
    }
  }

  # some random data in three variables
  c1<- runif(5, 10, 20)
  c2<- runif(5, 1, 5)
  c3 <- runif(5, 15, 25)
  # basic plot
  tri(c1,c2,c3)
dev.off()

png(file="g268.png", width=600, height=600)
  # plot with different symbols and a grid
  tri(c1,c2,c3, symb=7, grid=T) 
dev.off()

png(file="g269.png", width=600, height=600)
  histogram( ~ Sepal.Length | Species, data = iris, 
             layout = c(1,3) )
dev.off()

png(file="g270.png", width=600, height=600)
  xyplot( Sepal.Length ~ Sepal.Width | Species, data = iris,
          layout = c(1,3) )
dev.off()

png(file="g271.png", width=600, height=600)
  xyplot( Sepal.Length ~ Sepal.Width, group = Species, data = iris,
          panel = function (x, y, groups, ...) {
            panel.superpose(x, y, groups = groups, ...)
            groups <- as.factor(groups)
            for (i in seq(along=levels(groups))) {
              g <- levels(groups)[i]
              panel.lmline( x[groups == g], y[groups == g], 
                            col = trellis.par.get("superpose.line")$col[i] )
            }
          }
        )
dev.off()

png(file="g272.png", width=600, height=600)
  xyplot( Sepal.Length ~ Sepal.Width, group = Species, data = iris,
          panel = function (x, y, groups, ...) {
            panel.superpose(x, y, groups = groups, ...)
            groups <- as.factor(groups)
            for (i in seq(along=levels(groups))) {
              g <- levels(groups)[i]
              panel.loess( x[groups == g], y[groups == g], 
                            col = trellis.par.get("superpose.line")$col[i] )
            }
          }
        )
dev.off()

png(file="g273.png", width=600, height=600)
  data(iris)
  plot(iris[1:4], pch=21, 
       bg=c("red", "green", "blue")[as.numeric(iris$Species)])
dev.off()

png(file="g274.png", width=600, height=600)
  a <- rnorm(10)
  b <- 1+ rnorm(10)
  c <- 1+ rnorm(10)
  d <- rnorm(10)
  x <- c(a,b,c,d)
  y <- factor(c( rep("A",20), rep("B",20)))
  z <- factor(c( rep("U",10), rep("V",20), rep("U",10) ))
  op <- par(mfrow=c(2,2))
  plot(x~y)
  plot(x~z)
  plot(x[z=="U"] ~ y[z=="U"], border="red", ylim=c(min(x),max(x)))
  plot(x[z=="V"] ~ y[z=="V"], border="blue", add=T)
  plot(x[y=="A"] ~ z[y=="A"], border="red", ylim=c(min(x),max(x)))
  plot(x[y=="B"] ~ z[y=="B"], border="blue", add=T)
  par(op)
dev.off()

png(file="g275.png", width=600, height=600)
  l <- rep("",length(x))
  for (i in 1:length(x)){
    l[i] <- paste(y[i],z[i])
  }
  l <- factor(l)
  boxplot(x~l)
dev.off()

png(file="g276.png", width=600, height=600)
  # l is a 2-element list
  myplot1 <- function (x, l, ...) {
    t <- tapply(x,l,mean)
    l1 <- levels(l[[1]])
    l2 <- levels(l[[2]])
    matplot(t,
            type='l', lty=1, col=1:length(l2),
            axes=F, ...)
    axis(1, 1:2, l1)
    axis(2)
    lim <- par("usr")
    legend(lim[1] + .05*(lim[2]-lim[1]), lim[4],
           l2, lwd=1, lty=1, col=1:length(l2) )
  }
  op <- par(mfrow=c(1,2))
  myplot1( x, list(y,z), ylim=c(0,2), ylab = "" )
  myplot1( x, list(z,y), ylim=c(0,2), ylab = "" )
  par(op)
dev.off()

png(file="g277.png", width=600, height=600)
  myplot3 <- function (x, l, ...) {
    l1 <- levels(l[[1]])
    l2 <- levels(l[[2]])
    t0 <- tapply(x,l,min)
    t1 <- tapply(x,l,function(x)quantile(x,.25))
    t2 <- tapply(x,l,median)
    t3 <- tapply(x,l,function(x)quantile(x,.75))
    t4 <- tapply(x,l,max)
    matplot(cbind(t0,t1,t2,t3,t4),
            type='l', 
            lty=c(rep(3,length(l2)), rep(2,length(l2)), 
                  rep(1,length(l2)), rep(2,length(l2)),
                  rep(3,length(l2)) ),
            col=1:length(l2),
            axes=F, ...)
    axis(1, 1:2, l1)
    axis(2)
    lim <- par("usr")
    legend(lim[1] + .05*(lim[2]-lim[1]), lim[4],
           l2, lwd=1, lty=1, col=1:length(l2) )
  }
  op <- par(mfrow=c(1,2))
  myplot3( x, list(y,z), ylab = "" )
  myplot3( x, list(z,y), ylab = "" )
  par(op)
dev.off()

png(file="g278.png", width=600, height=600)
  shaded.pie <- function (...) {
    pie(...)
    op <- par(new=T)
    a <- seq(0,2*pi,length=100)
    for (i in (256:64)/256) {
      r <- .8-.1*(1-i)
      polygon( .1+r*cos(a), -.2+r*sin(a), border=NA, col=rgb(i,i,i))
    }
    par(new=T)
    pie(...)
    par(op)
  }
  x <- rpois(10,5)
  x <- x[x>0]
  shaded.pie(x)
dev.off()
detach.everything()

png(file="g279.png", width=600, height=600)
  library(MASS)
  data(beav1)
  plot(beav1$temp ~ beav1$time)
dev.off()

png(file="g280.png", width=600, height=600)
  x <- beav1$time
  y <- beav1$temp
  o <- order(x)
  x <- x[o]
  y <- y[o]
  plot(y ~ x, 
       type = "l",
       xlab = "Time",
       ylab = "Temperature", 
       main = "The \"plot\" function, with type=\"l\"")
dev.off()

png(file="g281.png", width=600, height=600)
  plot(y ~ x, 
       type = "b",
       lwd = 3,
       xlim = c(0, 400),
       xlab = "Time",
       ylab = "Temperature", 
       main = "The \"plot\" function, with type=\"b\"")
dev.off()

png(file="g282.png", width=600, height=600)
  x <- as.matrix( EuStockMarkets[1:50,] )
  matplot(x,                 # By default: not lines, 
          main = "matplot",  # but unconnected coloured numbers
          xlab = "",
          ylab = "")
dev.off()

png(file="g283.png", width=600, height=600)
  matplot(x, 
          type = "l",    # Lines -- but I am not happy 
          lty = 1,       # with the axes
          xlab = "",
          ylab = "",
          main = "matplot")
dev.off()

png(file="g284.png", width=600, height=600)
  x <- as.matrix( EuStockMarkets )
  matplot(time(EuStockMarkets), 
          x, 
          log = "y",
          type = 'l', 
          lty = 1, 
          ylab = "Closing price", 
          xlab = "Date", 
          main = "matplot",
          axes = FALSE)
  axis(1)
  axis(2)
  box()
dev.off()

png(file="g285.png", width=600, height=600)
  pairs(longley)
dev.off()

png(file="g286.png", width=600, height=600)
  pairs(longley, 
        gap=0,
        diag.panel = function (x, ...) {
          par(new = TRUE)
          hist(x, 
               col = "light blue", 
               probability = TRUE, 
               axes = FALSE, 
               main = "")
          lines(density(x), 
                col = "red", 
                lwd = 3)
          rug(x)
        })
dev.off()

png(file="g287.png", width=600, height=200)
  stripchart(longley$Unemployed)  
dev.off()

png(file="g288.png", width=600, height=600)
  hist(longley$Unemployed)
dev.off()

png(file="g289.png", width=600, height=600)
  hist(longley$Unemployed,
       probability = TRUE,    # Change the vertical units, 
                              # to overlay a density estimation
       col = "light blue")
  lines(density(longley$Unemployed),
        col = "red",
        lwd = 3)
dev.off()

png(file="g290.png", width=200, height=600)
  boxplot(longley$Unemployed)
dev.off()

png(file="g291.png", width=600, height=200)
  boxplot(longley$Unemployed,
          horizontal = TRUE,
          col = "pink",
          main = "Box-and-whiskers plot (boxplot)")
dev.off()

png(file="g292.png", width=600, height=600)
  data(InsectSprays)
  boxplot(count ~ spray, 
          data = InsectSprays, 
          col = "pink",
          xlab = "Spray", 
          ylab = "Count",
          main = "Insect sprays")
dev.off()

png(file="g293.png", width=600, height=600)
  boxplot(count ~ spray, 
          data = InsectSprays, 
          col = "pink", 
          horizontal = TRUE, 
          las = 1,            # Horizontal labels
          xlab = "Count", 
          ylab = "Spray",
          main = "Insect sprays")
dev.off()

png(file="g294.png", width=600, height=600)
  N <- 50
  x <- seq(-1, 1, length=N)
  y <- seq(-1, 1, length=N)
  xx <- matrix(x, nr=N, nc=N)
  yy <- matrix(y, nr=N, nc=N, byrow=TRUE)
  z <- 1 / (1 + xx^2 + (yy + .2 * sin(10*yy))^2)
  contour(x, y, z,
          main = "Contour plot")
dev.off()

png(file="g295.png", width=600, height=600)
  image(z)
dev.off()

png(file="g296.png", width=600, height=600)
  image(x, y, z,
        xlab = "",
        ylab = "")
  contour(x, y, z, lwd=3, add=TRUE)
dev.off()

png(file="g297.png", width=600, height=600)
  persp(z)
dev.off()

png(file="g298.png", width=600, height=600)
  op <- par(mar=c(0,0,3,0)+.1)
  persp(x, y, z, 
        theta = 45, phi = 30, 
        shade = .5, 
        col = rainbow(N), 
        border = "green",
        main = "perspective plot, theta=45, phi=30")
  par(op)
dev.off()

png(file="g299.png", width=600, height=600)
  # From the manual: the sinc function
  x <- seq(-10, 10, length= 30)
  y <- x
  f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
  z <- outer(x, y, f)
  z[is.na(z)] <- 1
  op <- par(bg = "white", mar=c(0,2,3,0)+.1)
  persp(x, y, z, 
        theta = 30, phi = 30, 
        expand = 0.5, 
        col = "lightblue",
        ltheta = 120, 
        shade = 0.75, 
        ticktype = "detailed",
        xlab = "X", ylab = "Y", zlab = "Sinc(r)",
        main = "The sinc function"
  )
  par(op)
dev.off()

png(file="g300.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- 1 - x^2 + .3*rnorm(n)
  plot(y ~ x, 
       xlab = 'X axis', 
       ylab = "Y axis", 
       main = "Title")
dev.off()

png(file="g301.png", width=600, height=600)
  plot(y ~ x, 
       xlab = "", 
       ylab = "", 
       main = "")
  title(main = "Title", 
        xlab = "X axis", 
        ylab = "Y axis")
dev.off()

png(file="g302.png", width=600, height=600)
  set.seed(1)
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  box()
  N <- 50
  text(  
    runif(N), runif(N), 
    sample(    # Random words...
      scan("/usr/share/dict/cracklib-small", character(0)), 
      N
    ) 
  )
dev.off()

png(file="g303.png", width=600, height=600)
  N <- 200
  x <- runif(N, -4, 4)
  y <- sin(x) + .5 * rnorm(N)
  plot(x, y, 
       xlab = "", ylab = "", 
       main = paste("The \"mtext\" function",
                    paste(rep(" ", 60), collapse="")))
  mtext("Line 0", 3, line=0)
  mtext("Line 1", 3, line=1)
  mtext("Line 2", 3, line=2)
  mtext("Line 3", 3, line=3)
dev.off()

png(file="g304.png", width=600, height=600)
  N <- 200
  x <- runif(N, -4, 4)
  y <- sin(x) + .5 * rnorm(N)
  plot(x, y, xlab="", ylab="", main="")
  mtext("Subtitle", 3, line=.8)
  mtext("Title",    3, line=2, cex=1.5)
  mtext("X axis", 1,          line=2.5, cex=1.5)
  mtext("X axis subtitle", 1, line=3.7)
dev.off()

png(file="g305.png", width=600, height=600)
  N  <- 200
  x  <- seq(-4,4, length=N)
  y1 <- sin(x)
  y2 <- cos(x)
  op <- par(mar=c(5,4,4,4)) # Add some space in the right margin
                            # The default is c(5,4,4,2) + .1
  xlim <- range(x)
  ylim <- c(-1.1, 1.1)
  plot(x, y1, col="blue", type="l", 
       xlim=xlim, ylim=ylim,
       axes=F, xlab="", ylab="", main="Title")
  axis(1)
  axis(2, col="blue")
  par(new=TRUE)
  plot(x, y2, col="red", type="l", 
       xlim=xlim, ylim=ylim,
       axes=F, xlab="", ylab="", main="")
  axis(4, col="red")
  mtext("First Y axis",  2, line=2, col="blue", cex=1.2)
  mtext("Second Y axis", 4, line=2, col="red",  cex=1.2)
dev.off()

png(file="g306.png", width=600, height=600)
  x <- seq(-5,5,length=200)
  y <- sqrt(1+x^2)
  plot(y~x, type='l',
       ylab=expression( sqrt(1+x^2) ))
  title(main=expression(
    "graph of the function f"(x) == sqrt(1+x^2)
  ))
dev.off()

png(file="g307.png", width=600, height=600)
  x <- seq(-5,5,length=200)
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    y <- sqrt(i+x^2)
    plot(y ~ x, 
         type = 'l', 
         ylim = c(0,6),
         ylab = substitute( 
           expression( sqrt(i+x^2) ), 
           list(i=i) 
         ))
    title(main = substitute(
      "graph of the function f"(x) == sqrt(i+x^2),
      list(i=i)))
  }
  par(op)
dev.off()

png(file="g308.png", width=600, height=600)
  # From the manual
  plot(1:10, 1:10, main = "text(...) examples\n~~~~~~~~~~~~~~",
       sub = "R is GNU ©, but not ® ...")
  mtext("«ISO-accents»: ± éè øØ å<Å æ<Æ", side=3)
  points(c(6,2), c(2,1), pch = 3, cex = 4, col = "red")
  text(6, 2, "the text is CENTERED around (x,y) = (6,2) by default",
       cex = .8)
  text(2, 1, "or Left/Bottom - JUSTIFIED at (2,1) by `adj = c(0,0)'",
       adj = c(0,0))
  text(4, 9, expression(hat(beta) == (X^t * X)^{-1} * X^t * y))
  text(4, 8.4, "expression(hat(beta) == (X^t * X)^{-1} * X^t * y)", cex = .75)
  text(4, 7, expression(bar(x) == sum(frac(x[i], n), i==1, n)))
dev.off()

png(file="g309.png", width=600, height=600)
  # From the manual
  plot(1:9, type="n", axes=FALSE, frame=TRUE, ylab="",
       main= "example(Japanese)", xlab= "using Hershey fonts")
  par(cex=3)
  Vf <- c("serif", "plain")
  text(4, 2, "\\#J2438\\#J2421\\#J2451\\#J2473", vfont = Vf)
  text(4, 4, "\\#J2538\\#J2521\\#J2551\\#J2573", vfont = Vf)
  text(4, 6, "\\#J467c\\#J4b5c", vfont = Vf)
  text(4, 8, "Japan", vfont = Vf)
  par(cex=1)
  text(8, 2, "Hiragana")
  text(8, 4, "Katakana")
  text(8, 6, "Kanji")
  text(8, 8, "English")
dev.off()

png(file="g310.png", width=600, height=600)
  # Other example from the manual 
  # (it also contains katakana and kanji tables)
  make.table <- function(nr, nc) {
    opar <- par(mar=rep(0, 4), pty="s")
    plot(c(0, nc*(10%/%nc) + 1), c(0, -(nr + 1)),
         type="n", xlab="", ylab="", axes=FALSE)
    invisible(opar)
  }

  get.r <- function(i, nr)   i %% nr + 1
  get.c <- function(i, nr)   i %/% nr + 1
  Esc2 <- function(str)      paste("\\", str, sep="")

  draw.title <- function(title, nc)
    text((nc*(10%/%nc) + 1)/2, 0, title, font=2)

  draw.vf.cell <- function(typeface, fontindex, string, i, nr, raw.string=NULL) {
    r <- get.r(i, nr)
    c <- get.c(i, nr)
    x0 <- 2*(c - 1)
    if (is.null(raw.string)) raw.string <- Esc2(string)
    text(x0 + 1.1, -r, raw.string, col="grey")
    text(x0 + 2,   -r, string, vfont=c(typeface, fontindex))
    rect(x0 +  .5, -(r - .5), x0 + 2.5, -(r + .5), border="grey")
  }

  draw.vf.cell2 <- function(string, alt, i, nr) {
    r <- get.r(i, nr)
    c <- get.c(i, nr)
    x0 <- 3*(c - 1)
    text(x0 + 1.1, -r, Esc2(string <- Esc2(string)), col="grey")
    text(x0 + 2.2, -r, Esc2(Esc2(alt)), col="grey", cex=.6)
    text(x0 + 3,   -r, string, vfont=c("serif", "plain"))
    rect(x0 +  .5, -(r - .5), x0 + 3.5, -(r + .5), border="grey")
  }

  tf <- "serif"
  fi <- "plain"
  nr <- 25
  nc <- 4
  oldpar <- make.table(nr, nc)
  index <- 0
  digits <- c(0:9,"a","b","c","d","e","f")
  draw.title("Hiragana : \\\\#J24nn", nc)
  for (i in 2:7) {
    for (j in 1:16) {
      if (!((i == 2 && j == 1) || (i == 7 && j > 4))) {
        draw.vf.cell(tf, fi, paste("\\#J24", i, digits[j], sep=""),
                     index, nr)
        index <- index + 1
      }
    }
  }
dev.off()

png(file="g311.png", width=600, height=600)
  plot(runif(5), ylim=c(0,1), type='l')
  for (i in c('red', 'blue', 'green')) {
    lines( runif(5), col=i )
  }
  title(main="Lines in various colours")
dev.off()

png(file="g312.png", width=600, height=600)
  plot(runif(5), ylim=c(0,1), type='n')
  for (i in 5:1) {
    lines( runif(5), col=i, lwd=i )
  }
  title(main = "Varying the line thickness")
dev.off()

png(file="g313.png", width=600, height=800)
  op <- par(mfrow=c(3,2))
  plot(runif(5), type = 'p', 
       main = "plot type 'p' (points)")
  plot(runif(5), type = 'l', 
       main = "plot type 'l' (lines)")
  plot(runif(5), type = 'b', 
       main = "plot type 'b' (both points and lines)")
  plot(runif(5), type = 's', 
       main = "plot type 's' (stair steps)")
  plot(runif(5), type = 'h', 
       main = "plot type 'h' (histogram)")
  plot(runif(5), type = 'n', 
       main = "plot type 'n' (no plot)")
  par(op)
dev.off()

png(file="g314.png", width=600, height=800)
  op <- par(mfrow=c(3,2), mar=c(3,1,5,1))
  plot(runif(5), lty = 1,
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 1 (default, solid)")
  plot(runif(5), lty = 2,
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 2 (dashed)")
  plot(runif(5), lty = 3,
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 3 (dotted)")
  plot(runif(5), lty = "dotdash",
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 4 (dot, dash)")
  plot(runif(5), lty = "longdash",
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 5 (longdash)")
  plot(runif(5), lty = "twodash",
       axes = FALSE, type = "l", lwd = 3, 
       main = "lty = 6 (twodash)")
  par(op)
dev.off()

png(file="g315.png", width=600, height=600)
  # You can also cook up your own line type
  # by providing the length of each segment and 
  # each space
  op <- par(mfrow=c(2,2), mar=c(3,1,5,1))
  for (lty in c("42", "14", "8222", "82624222")) { 
    plot(runif(5), lty = lty,
       axes = FALSE, type = "l", lwd = 3, 
       main = paste("lty =", lty))
  }
  par(op)
dev.off()

png(file="g316.png", width=600, height=600)
  op <- par(mar=c(1,1,4,1)+.1)
  plot(0,0, 
       xlim = c(1,5), ylim = c(-.5,4), 
       axes = F, 
       xlab = '', ylab = '',
       main = "Available symbols")
  for (i in 0:4) {
    for (j in 1:5) {
      n <- 5*i+j
      points(j, i, 
             pch = n, 
             cex = 3)
      text(j,i-.25, as.character(n))
    }
  }
  par(op)
dev.off()

png(file="g317.png", width=600, height=600)
  hist(longley$Unemployed, density=3, angle=45)
dev.off()

png(file="g318.png", width=600, height=600)
  op <- par(mfrow = c(2, 2))
  for (i in 1:4) 
    plot(runif(20), runif(20), 
         main=paste("random plot (",i,")",sep=''))
  par(op)
  mtext("Four plots, without enough room for this title", 
         side=3, font=2, cex=2, col='red')
dev.off()

png(file="g319.png", width=600, height=600)
  op <- par(mfrow = c(2, 2), 
            oma = c(0,0,3,0)   # Outer margins
           )
  for (i in 1:4) 
    plot(runif(20), runif(20), 
         main=paste("random plot (",i,")",sep=''))
  par(op)
  mtext("Four plots, with some room for this title", 
        side=3, line=1.5, font=2, cex=2, col='red')
dev.off()

png(file="g320.png", width=600, height=600)
  op <- par(mfrow = c(2, 2), 
            oma = c(0,0,3,0),
            mar = c(3,3,4,1) + .1     # Margins
           )
  for (i in 1:4) 
    plot(runif(20), runif(20), 
         xlab = "", ylab = "",
         main=paste("random plot (",i,")",sep=''))
  par(op)
  mtext("Title", 
        side=3, line=1.5, font=2, cex=2, col='red')
  par(op)
dev.off()

png(file="g321.png", width=600, height=600)
  n <- 20
  x <- rnorm(n)
  y <- x^2 - 1 + .3*rnorm(n)
  plot(y ~ x,
       main = "The \"fig\" graphic parameter")  
  op <- par()
  for (i in 2:10) {
    done <- FALSE
    while(!done) {
      a <- c( sort(runif(2,0,1)), 
              sort(runif(2,0,1)) )
      par(fig=a, new=T)
      r <- try(plot(runif(5), type='l', col=i))
      done <- !inherits(r, "try-error")
    }
  }
  par(op)
dev.off()

png(file="g322.png", width=600, height=600)
  n <- 1000
  x <- rt(n, df=10)
  hist( x, 
        col = "light blue",
        probability = "TRUE",
        ylim = c(0, 1.2*max(density(x)$y)))
  lines(density(x),
        col = "red",
        lwd = 3)
  op <- par(fig = c(.02,.4,.5,.98), 
            new = TRUE)
  qqnorm(x, 
         xlab = "", ylab = "", main = "",
         axes = FALSE)
  qqline(x, col = "red", lwd = 2)
  box(lwd=2)
  par(op)
dev.off()

png(file="g323.png", width=600, height=600)
  op <- par(oma = c(0,0,3,0))
  layout(matrix(c(1, 1, 1,
                  2, 3, 4,
                  2, 3, 4), nr = 3, byrow = TRUE))
  hist( rnorm(n), col = "light blue" )
  hist( rnorm(n), col = "light blue" )
  hist( rnorm(n), col = "light blue" )
  hist( rnorm(n), col = "light blue" )
  mtext("The \"layout\" function",
        side = 3, outer = TRUE, 
        font = 2, size = 1.2)
  par(op)
dev.off()

png(file="g324.png", width=600, height=600)
  random.plot <- function () {
    N <- 200
    f <- sample(list(rnorm, 
                     function (x) { rt(x, df=2) }, 
                     rlnorm, 
                     runif), 
                1) [[1]]
    x <- f(N)
    hist(x, col="light blue", main="", xlab="", ylab="", axes=F)
    axis(1)
  }
  op <- par(bg="white", mar=c(2.5,2,1,2))
  split.screen(c(2,1))
  split.screen(c(1,3), screen = 2)
  screen(1); random.plot()
  #screen(2); random.plot() # Screen 2 was split into three screens: 3, 4, 5
  screen(3); random.plot()
  screen(4); random.plot()
  screen(5); random.plot()
  close.screen(all=TRUE)
  par(op)
dev.off()

png(file="g325.png", width=600, height=600)
  plot(runif(5), runif(5), 
       xlim = c(0,1), ylim = c(0,1))
  points(runif(5), runif(5), 
         col = 'orange', pch = 16, cex = 3)
  lines(runif(5), runif(5), 
        col = 'red')
  segments(runif(5), runif(5), runif(5), runif(5), 
           col = 'blue')
  title(main = "Overlaying points, segments, lines...")
dev.off()

png(file="g326.png", width=600, height=600)
  my.col <- function (f, g, xmin, xmax, col, N=200,
                      xlab="", ylab="", main="") {
    x <- seq(xmin, xmax, length = N)
    fx <- f(x)
    gx <- g(x)
    plot(0, 0, type = 'n', 
         xlim = c(xmin,xmax),
         ylim = c( min(fx,gx), max(fx,gx) ),
         xlab = xlab, ylab = ylab, main = main)
    polygon( c(x,rev(x)), c(fx,rev(gx)), 
             col = 'red', border = 0 )
    lines(x, fx, lwd = 3)
    lines(x, gx, lwd = 3)
  }
  op <- par(mar=c(3,3,4,1)+.1)
  my.col( function(x) x^2, function(x) x^2+10*sin(x), 
          -6, 6,
          main = "The \"polygon\" function")
  par(op)
dev.off()

png(file="g327.png", width=600, height=600)
  # From the manual
  ch.col <- c("rainbow(n, start=.7, end=.1)", 
              "heat.colors(n)",
              "terrain.colors(n)", 
              "topo.colors(n)", 
              "cm.colors(n)")
  n <- 16
  nt <- length(ch.col)
  i <- 1:n
  j <- n/nt
  d <- j/6
  dy <- 2*d
  plot(i, i+d, 
       type="n", 
       yaxt="n", 
       ylab="", 
       main=paste("color palettes;  n=",n))
  for (k in 1:nt) {
    rect(i-.5, (k-1)*j+ dy, i+.4,  k*j, 
         col = eval(parse(text=ch.col[k])))
    text(2*j,  k * j +dy/4, ch.col[k])
  }
dev.off()

png(file="g328.png", width=600, height=600)
  x <- seq(-6,6,length=200)
  y <- sin(x)
  z <- cos(x)
  plot(y~x, type='l', lwd=3, 
       ylab='', xlab='angle', main="Trigonometric functions")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  lines(z~x, type='l', lwd=3, col='red')
  legend(-6,-1, yjust=0,
         c("Sine", "Cosine"),
         lwd=3, lty=1, col=c(par('fg'), 'red'),
        )
dev.off()

png(file="g329.png", width=600, height=600)
  plot(y~x, type='l', lwd=3, 
       ylab='', xlab='angle', main="Trigonometric functions")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  lines(z~x, type='l', lwd=3, col='red')
  legend("bottomleft", 
         c("Sine", "Cosine"),
         lwd=3, lty=1, col=c(par('fg'), 'red'),
        )
dev.off()

png(file="g330.png", width=600, height=600)
  plot(y~x, type='l', lwd=3, 
       ylab='', xlab='angle', main="Trigonometric functions")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  lines(z~x, type='l', lwd=3, col='red')
  legend("bottomleft", 
         c("Sine", "Cosine"),
         inset = c(.03, .03),
         lwd=3, lty=1, col=c(par('fg'), 'red'),
        )
dev.off()

png(file="g331.png", width=600, height=600)
  op <- par(no.readonly=TRUE)
  plot(y~x, type='l', lwd=3, 
       ylab='', xlab='angle', main="Trigonometric functions")
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  lines(z~x, type='l', lwd=3, col='red')
  par(xpd=TRUE)  # Do not clip to the drawing area
  lambda <- .025
  legend(par("usr")[1], 
         (1 + lambda) * par("usr")[4] - lambda * par("usr")[3],
         c("Sine", "Cosine"),
         xjust = 0, yjust = 0,
         lwd=3, lty=1, col=c(par('fg'), 'red'),
        )
  par(op)
dev.off()

png(file="g332.png", width=600, height=600)
  library(grid)
  grid.show.viewport(viewport(x=0.6, y=0.6,
                     w=unit(1, "inches"), h=unit(1, "inches")))
dev.off()

png(file="g333.png", width=600, height=600)
  grid.show.layout(grid.layout(4,2,
                   heights=unit(rep(1, 4),
                                c("lines", "lines", "lines", "null")),
                   widths=unit(c(1, 1), "inches")))
dev.off()

png(file="g334.png", width=600, height=600)
  dessine <- function () {
    push.viewport(viewport(w = 0.9, h = 0.9, 
                           xscale=c(-.1,1.1), yscale=c(-.1,1.1)))
    grid.rect(gp=gpar(fill=rgb(.5,.5,0)))
    grid.points( runif(50), runif(50) )
    pop.viewport()
  }
  grid.newpage()
  grid.rect(gp=gpar(fill=rgb(.3,.3,.3)))
  push.viewport(viewport(layout=grid.layout(2, 2)))
  for (i in 1:2) {
    for (j in 1:2) {
      push.viewport(viewport(layout.pos.col=i,
                             layout.pos.row=j))
      dessine()
      pop.viewport()
    }
  }
  pop.viewport()
dev.off()

png(file="g335.png", width=600, height=600)
  grid.multipanel(vp=viewport(0.5, 0.5, 0.8, 0.8))
dev.off()

png(file="g336.png", width=600, height=600)
  do.it <- function (x=runif(100), y=runif(100), 
                     a=.9, b=.1, 
                     col1=rgb(0,.3,0), col2=rgb(1,1,0)) {
    xscale <- range(x) + c(-1,1)*.05
    yscale <- range(y) + c(-1,1)*.05
    grid.newpage()
    grid.rect(gp=gpar(fill=col1, col=col1))
    w1 <- a - b/2
    w2 <- 1 - a - b/2
    c1 <- b/3 + w1/2
    c2 <- a + b/6 + w2/2
    vp1 <- viewport(x=c1, y=c1, width=w1, height=w1,
                    xscale=xscale, yscale=yscale)
    push.viewport(vp1)
    grid.rect(gp=gpar(fill=col2, col=col2))
    grid.points(x,y)
    pop.viewport()
    vp2 <- viewport(x=c1, y=c2, width=w1, height=w2,
                    xscale=xscale, yscale=c(0,1))
    push.viewport(vp2)
    grid.rect(gp=gpar(fill=col2, col=col2))
    grid.points(x,rep(.5,length(x)))
    pop.viewport()
    vp3 <- viewport(x=c2, y=c1, width=w2, height=w1,
                    xscale=c(0,1), yscale=yscale)
    push.viewport(vp3)
    grid.rect(gp=gpar(fill=col2, col=col2))
    grid.points(rep(.5,length(y)),y)
    pop.viewport()
  }
  do.it()
dev.off()

png(file="g337.png", width=600, height=600)
  data(quakes)
  library(lattice)
  Depth <- equal.count(quakes$depth, number=8, overlap=.1)
  xyplot(lat ~ long | Depth, data = quakes)
dev.off()

png(file="g338.png", width=600, height=600)
  plot(lat ~ long, data=quakes)
dev.off()

png(file="g339.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  plot(lat ~ long, data=quakes)
  plot(lat ~ -depth, data=quakes)
  plot(-depth ~ long, data=quakes)
  par(op)
dev.off()

png(file="g340.png", width=600, height=600)
  library(mva)
  biplot(princomp(quakes[1:3]))
dev.off()

png(file="g341.png", width=600, height=600)
  pairs( princomp(quakes[1:3])$scores )
dev.off()

png(file="g342.png", width=600, height=600)
  library(scatterplot3d)
  scatterplot3d(quakes[,1:3],
                highlight.3d = TRUE,
                pch = 20)
dev.off()

png(file="g343.png", width=600, height=600)
  data(barley)
  barchart(yield ~ variety | year * site, data=barley)
dev.off()

png(file="g344.png", width=600, height=600)
  barchart(yield ~ variety | year * site, data = barley,
           ylab = "Barley Yield (bushels/acre)",
           scales = list(x = list(0, abbreviate = TRUE,
                         minlength = 5)))
dev.off()

png(file="g345.png", width=600, height=600)
  dotplot(yield ~ variety | year * site, data = barley)
dev.off()

png(file="g346.png", width=400, height=1000)
  dotplot(variety ~ yield | year * site, data = barley)
dev.off()

png(file="g347.png", width=600, height=600)
  library(nlme)
  data(bdf)
  d <- data.frame( iq=bdf$IQ.perf, sex=bdf$sex, den=bdf$denomina )
  d <- d[1:100,] 
  bwplot( ~ d$iq | d$sex + d$den )
dev.off()

png(file="g348.png", width=600, height=600)
  histogram( ~ d$iq | d$sex + d$den )
dev.off()

png(file="g349.png", width=600, height=600)
  densityplot( ~ d$iq | d$sex + d$den )
dev.off()

png(file="g350.png", width=600, height=600)
  d <- data.frame( x=bdf$aritPOST, y=bdf$sex, z=equal.count(bdf$langPRET) )
  bwplot( ~ x | y + z, data=d )
dev.off()

png(file="g351.png", width=600, height=600)
  histogram( ~ x | y + z, data=d )
dev.off()

png(file="g352.png", width=600, height=600)
  densityplot( ~ x | y + z, data=d )
dev.off()

png(file="g353.png", width=600, height=600)
  d <- data.frame( x= (bdf$IQ.perf>11), y=bdf$sex, z=bdf$denomina )
  d <- as.data.frame(table(d))
  barchart( Freq ~ x | y * z, data=d )
dev.off()

png(file="g354.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  y <- x^3+rnorm(n)
  plot1 <- xyplot(y~x)
  plot2 <- bwplot(x)
  # Beware, the order is xmin, ymin, xmax, ymax
  print(plot1, position=c(0, .2, 1, 1),  more=T)
  print(plot2, position=c(0, 0,  1, .2), more=F)
dev.off()

png(file="g355.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  y <- x^4+rnorm(n)
  k <- .7
  op <- par(mar=c(0,0,0,0))
  # Attention : l'ordre est xmin, xmax, ymin, ymax
  par(fig=c(0,k,0,k))
  plot(y~x)
  par(fig=c(0,k,k,1), new=T)
  boxplot(x, horizontal=T)  
  par(fig=c(k,1,0,k), new=T)
  boxplot(y, horizontal=F)  
  par(op)
dev.off()

png(file="g356.png", width=800, height=800)
  show.settings()
dev.off()

png(file="g357.png", width=600, height=600)
  y <- cars$dist
  x <- cars$speed
  vitesse <- shingle(x, co.intervals(x, number=6))
  histogram(~ x | vitesse, type = "density",
            panel = function(x, ...) {
              ps <- trellis.par.get('plot.symbol')
              nps <- ps
              nps$cex <- 1
              trellis.par.set('plot.symbol', nps)
              panel.histogram(x, ...)
              panel.densityplot(x, col = 'brown', lwd=3) 
              panel.xyplot(x = jitter(x), y = rep(0, length(x)), col='brown' )
              panel.mathdensity(dmath = dnorm,
                                args = list(mean=mean(x),sd=sd(x)),
                                lwd=3, lty=2, col='white')
              trellis.par.set('plot.symbol', ps)
            })
dev.off()

png(file="g358.png", width=600, height=600)
  data(sunspot.year)
  sunspot <- sunspot.year[20 + 1:37]
  xyplot(sunspot ~ 1:37 ,type = "l",
         scales = list(y = list(log = TRUE)),
         sub = "log scales")
dev.off()

png(file="g359.png", width=600, height=600)
  xyplot(sunspot ~ 1:37 ,type = "l", aspect="xy",
         scales = list(y = list(log = TRUE)),
         sub = "log scales")
dev.off()
detach.everything()

png(file="g360.png", width=600, height=600)
  data(USArrests)
  p <- prcomp(USArrests)
  biplot(p)
dev.off()

png(file="g361.png", width=600, height=600)
  plot(p)
dev.off()

png(file="g362.png", width=600, height=600)
  a <- seq(0,2*pi,length=100)
  plot( cos(a), sin(a), 
        type = 'l', lty = 3,
        xlab = 'comp 1', ylab = 'comp 2', 
        main = "Correlation circle")
  v <- t(p$rotation)[1:2,]
  arrows(0,0, v[1,], v[2,], col='red')
  text(v[1,], v[2,],colnames(v))
dev.off()

png(file="g363.png", width=600, height=600)
  # Copy-pasted with the help of the "deparse" command:
  #   cat( deparse(x), file='foobar')
  notes <- matrix( c(15, NA, 7, 15, 11, 7, 7, 8, 11, 11, 13,
  6, 14, 19, 9, 8, 6, NA, 7, 14, 11, 13, 16, 10, 18, 7, 7,
  NA, 11, NA, NA, 6, 15, 5, 11, 7, 3, NA, 3, 1, 10, 1, 1,
  18, 13, 2, 2, 0, 7, 9, 13, NA, 19, 0, 17, 8, 2, 9, 2, 5,
  12, 0, 8, 12, 8, 4, 8, 0, 5, 5.5, 1, 12, 4, 13, 5, 11, 6,
  0, 7, 8, 11, 9, 9, 9, 14, 8, 5, 8, 5, 5, 12, 6, 16.5,
  13.5, 15, 3, 10.5, 1.5, 10.5, 9, 15, 7.5, 12, 13.5, 4.5,
  13.5, 13.5, 6, 12, 7.5, 9, 6, 13.5, 13.5, 15, 13.5, 6, NA,
  13.5, 4.5, 14, NA, 14, 14, 14, 8, 16, NA, 6, 6, 12, NA, 7,
  15, 13, 17, 18, 5, 14, 17, 17, 13, NA, NA, 16, 14, 18, 13,
  17, 17, 8, 4, 16, 16, 16, 10, 15, 8, 10, 13, 12, 14, 8,
  19, 7, 7, 9, 8, 15, 16, 8, 7, 12, 5, 11, 17, 13, 13, 7,
  12, 15, 8, 17, 16, 16, 6, 7, 11, 15, 15, 19, 12, 15, 16,
  13, 19, 14, 4, 13, 13, 19, 11, 15, 7, 20, 16, 10, 12, 16,
  14, 0, 0, 11, 9, 4, 10, 0, 0, 5, 11, 12, 7, 12, 17, NA, 6,
  6, 9, 7, 0, 7, NA, 15, 3, 20, 11, 10, 13, 0, 0, 6, 1, 5,
  6, 5, 4, 2, 0, 8, 9, NA, 0, 11, 11, 0, 7, 0, NA, NA, 7, 0,
  NA, NA, 6, 9, 6, 4, 5, 4, 3 ), nrow=30)
  notes <- data.frame(notes)
  # These are not the real names
  row.names(notes) <- 
    c("Anouilh", "Balzac", "Camus", "Dostoievski",
    "Eschyle", "Fenelon", "Giraudoux", "Homer",
    "Ionesco", "Jarry", "Kant", "La Fontaine", "Marivaux",
    "Nerval", "Ossian", "Platon", "Quevedo", "Racine",
    "Shakespeare", "Terence", "Updike", "Voltaire",
    "Whitman", "X", "Yourcenar", "Zola", "27", "28", "29",
    "30")
  attr(notes, "names") <- c("C1", "DM1", "C2", "DS1", "DM2",
                            "C3", "DM3", "DM4", "DS2")
  notes <- as.matrix(notes)
  notes <- t(t(notes) - apply(notes, 2, mean, na.rm=T))
  # Get rid of NAs
  notes[ is.na(notes) ] <- 0
  # plots
  plot(princomp(notes))
dev.off()

png(file="g364.png", width=600, height=600)
  biplot(princomp(notes))
dev.off()

png(file="g365.png", width=600, height=600)
  pairs(princomp(notes)$scores, gap=0)
dev.off()

png(file="g366.png", width=600, height=600)
  pairs(princomp(notes)$scores[,1:3])
dev.off()

png(file="g367.png", width=600, height=600)
  p <- princomp(notes)
  pairs( rbind(p$scores, p$loadings)[,1:3], 
         col=c(rep(1,p$n.obs),rep(2, length(p$center))),
         pch=c(rep(1,p$n.obs),rep(3, length(p$center))),
       )
dev.off()

png(file="g368.png", width=600, height=600)
  library(lattice)
  splom(as.data.frame(
    princomp(notes)$scores[,1:3]
  ))
dev.off()

png(file="g369.png", width=600, height=600)
  my.acp <- function (x) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    # Translation, to use linear algebra
    centre <- apply(x, 2, mean)
    x <- x - matrix(centre, nr=n, nc=p, byrow=T)
    # diagonalizations, base changes
    e1 <- eigen( t(x) %*% x, symmetric=T )
    e2 <- eigen( x %*% t(x), symmetric=T )
    variables <- t(e2$vectors) %*% x
    subjects <- t(e1$vectors) %*% t(x)
    # The vectors we want are the columns of the 
    # above matrices. To draw them, with the "pairs"
    # function, we have to transpose them.
    variables <- t(variables)
    subjects <- t(subjects)
    eigen.values <- e1$values
    # Plot
    plot( subjects[,1:2], 
          xlim=c( min(c(subjects[,1],-subjects[,1])), 
                  max(c(subjects[,1],-subjects[,1])) ),
          ylim=c( min(c(subjects[,2],-subjects[,2])), 
                  max(c(subjects[,2],-subjects[,2])) ),
          xlab='', ylab='', frame.plot=F )
    par(new=T)
    plot( variables[,1:2], col='red',
          xlim=c( min(c(variables[,1],-variables[,1])), 
                  max(c(variables[,1],-variables[,1])) ),
          ylim=c( min(c(variables[,2],-variables[,2])), 
                  max(c(variables[,2],-variables[,2])) ),
          axes=F, xlab='', ylab='', pch='.')
    axis(3, col='red')
    axis(4, col='red')
    arrows(0,0,variables[,1],variables[,2],col='red')
    # Return the data
    invisible(list(data=x, centre=centre, subjects=subjects, 
                   variables=variables, eigen.values=eigen.values))
  }

  n <- 20
  p <- 5
  x <- matrix( rnorm(p*n), nr=n, nc=p )
  my.acp(x)
  title(main="ACP by hand")
dev.off()

png(file="g370.png", width=600, height=600)
  biplot(princomp(x))
dev.off()

png(file="g371.png", width=600, height=600)
  b <- read.table('ling.txt')
  names(b) <- c(letters[1:26], 'language')
  a <- b[,1:26]
  a <- a/apply(a,1,sum)
  biplot(princomp(a))
dev.off()

png(file="g372.png", width=600, height=600)
  plot(hclust(dist(a)))
dev.off()

png(file="g373.png", width=600, height=600)
  kmeans.plot <- function (data, n=2, iter.max=10) {
    k <- kmeans(data,n,iter.max)
    p <- princomp(data)
    u <- p$loadings
    # The observations
    x <- (t(u) %*% t(data))[1:2,]
    x <- t(x)
    # The centers
    plot(x, col=k$cluster, pch=3, lwd=3)
    c <- (t(u) %*% t(k$center))[1:2,]
    c <- t(c)
    points(c, col = 1:n, pch=7, lwd=3)
    # A segment joining each observation to its group center
    for (i in 1:n) {
      for (j in (1:length(data[,1]))[k$cluster==i]) {
        segments( x[j,1], x[j,2], c[i,1], c[i,2], col=i )
      }
    }
    text( x[,1], x[,2], attr(x, "dimnames")[[1]] )
  }
  kmeans.plot(a,2)
dev.off()

png(file="g374.png", width=600, height=600)
  n <- 100
  v <- .1
  a <- rcauchy(n)
  b <- rcauchy(n)
  c <- rcauchy(n)
  d <- data.frame( x1 =  a+b+c+v*rcauchy(n),
                   x2 =  a+b-c+v*rcauchy(n),
                   x3 =  a-b+c+v*rcauchy(n),
                   x4 = -a+b+c+v*rcauchy(n),
                   x5 =  a-b-c+v*rcauchy(n),
                   x6 = -a+b-c+v*rcauchy(n) )
  biplot(princomp(d))
dev.off()

png(file="g375.png", width=600, height=600)
  rank.and.normalize.vector <- function (x) {
    x <- (rank(x)-.5)/length(x)
    x <- qnorm(x)
  }
  rank.and.normalize <- function (x) {
    if( is.vector(x) )
      return( rank.and.normalize.vector(x) )
    if( is.data.frame(x) ) {
      d <- NULL
      for (v in x) {
        if( is.null(d) )
          d <- data.frame( rank.and.normalize(v) )
        else 
          d <- data.frame(d, rank.and.normalize(v))
      }
      names(d) <- names(x)
      return(d)
    }
    stop("Data type not handled")
  }
  biplot(princomp(apply(d,2,rank.and.normalize)))
dev.off()

png(file="g376.png", width=600, height=600)
  pairs( princomp(d)$scores )
dev.off()

png(file="g377.png", width=600, height=600)
  pairs( princomp(apply(d,2,rank.and.normalize))$scores )
dev.off()

png(file="g378.png", width=600, height=600)
  N <- 1000
  X <- matrix(runif(2*N, -1, 1), nc=2)
  plot(X)
dev.off()

png(file="g379.png", width=600, height=600)
  M <- matrix(rnorm(4), nc=2)
  Y <- X %*% M
  plot(Y)
dev.off()

png(file="g380.png", width=600, height=600)
  plot(Y)
  p <- prcomp(Y)$rotation
  abline(0, p[2,1] / p[1,1], col="red", lwd=3)
  abline(0, -p[1,1] / p[2,1], col="red", lwd=3)
  abline(0, M[1,2]/M[1,1], col="blue", lwd=3)  
  abline(0, M[2,2]/M[2,1], col="blue", lwd=3)  
dev.off()

png(file="g381.png", width=600, height=600)
  op <- par(mfrow=c(2,2), mar=c(1,1,1,1))
  for (i in 1:4) {
    N <- 1000
    X <- matrix(runif(2*N, -1, 1), nc=2)
    M <- matrix(rnorm(4), nc=2)
    Y <- X %*% M
    plot(Y, xlab="", ylab="")
    p <- prcomp(Y)$rotation
    abline(0, p[2,1] / p[1,1], col="red", lwd=3)
    abline(0, -p[1,1] / p[2,1], col="red", lwd=3)
    abline(0, M[1,2]/M[1,1], col="blue", lwd=3)  
    abline(0, M[2,2]/M[2,1], col="blue", lwd=3)  
  }
  par(op)
dev.off()

png(file="g382.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  hist(X[,1], col="light blue")
  hist(X[,2], col="light blue")
  hist(Y[,1], col="light blue")
  hist(Y[,2], col="light blue")  
  par(op)
dev.off()

png(file="g383.png", width=600, height=600)
  library(e1071)
  r <- ica(Y,.1)
  plot(r$projection)
dev.off()

png(file="g384.png", width=600, height=600)
  library(mlica)
  ICA <- function (x,...) {
     prPCA <- PriorNormPCA(x);
     prNCP <- proposeNCP(prPCA,0.1);
     mlica(prNCP,...)
  }

  set.seed(1) # It sometimes crashes...
  N <- 1000
  X <- matrix(runif(2*N, -1, 1), nc=2)
  M <- matrix(rnorm(4), nc=2)
  Y <- X %*% M
  r <- ICA(Y)
  plot(r$S)
dev.off()

png(file="g385.png", width=600, height=600)
  n <- 100
  v <- .1
  # (Almost) planar data
  x <- rnorm(n)
  y <- x*x + v*rnorm(n)
  z <- v*rnorm(n)
  d <- cbind(x,y,z)
  # A rotation
  random.rotation.matrix <- function (n=3) {
    m <- NULL
    for (i in 1:n) {
      x <- rnorm(n)
      x <- x / sqrt(sum(x*x))
      y <- rep(0,n)
      if (i>1) for (j in 1:(i-1)) {
        y <- y + sum( x * m[,j] ) * m[,j]
      }
      x <- x - y
      x <- x / sqrt(sum(x*x))
      m <- cbind(m, x)
    }
    m
  }
  m <- random.rotation.matrix(3)
  d <- t( m %*% t(d) )
  pairs(d)
  title(main="The points lie in a plane")
dev.off()

png(file="g386.png", width=600, height=600)
  pairs(cmdscale(dist(d),3))
  title(main="MDS")
dev.off()

png(file="g387.png", width=600, height=600)
  pairs(princomp(d)$scores)
  title(main="Principal Component Analysis")
dev.off()

png(file="g388.png", width=600, height=600)
  # Data
  n <- 200  # Number of patients, number of columns
  k <- 10   # Dimension of the ambient space
  nb.points <- 5
  p <- matrix( 5*rnorm(nb.points*k), nr=k )
  barycentre <- function (x, n) {
    # Add number between the values of x in order to get a length n vector
    i <- seq(1,length(x)-.001,length=n)
    j <- floor(i)
    l <- i-j
    (1-l)*x[j] + l*x[j+1]
  }
  m <- apply(p, 1, barycentre, n)
  data.broken.line <- t(m)
  data.noisy.broken.line <- data.broken.line + rnorm(length(data.broken.line))
  library(splines)
  barycentre2 <- function (y,n) {
    m <- length(y)
    x <- 1:m
    r <- interpSpline(x,y)    
    #r <- lm( y ~ bs(x, knots=m) )
    predict(r, seq(1,m,length=n))$y
  }
  data.curve <- apply(p, 1, barycentre2, n)
  data.curve <- t(data.curve)
  data.noisy.curve <- data.curve + rnorm(length(data.curve))
  data.real <- read.table("Tla_z.txt", sep=",")
  r <- prcomp(t(data.real))
  data.real.3d <- r$x[,1:3]

  library(cluster)
  library(ape)
  mst.of.classification <- function (x, k=6, ...) {
    x <- t(x) 
    x <- t( t(x) - apply(x,2,mean) )
    r <- prcomp(x)
    y <- r$x
    u <- r$rotation
    r <- kmeans(y,k)
    z <- r$centers
    m <- mst(dist(z))
    plot(y[,1:2], ...)
    points(z[,1:2], col='red', pch=15)
    w <- which(m!=0)
    i <- as.vector(row(m))[w]
    j <- as.vector(col(m))[w]
    segments( z[i,1], z[i,2], z[j,1], z[j,2], col='red' )
  }
  set.seed(1)
  mst.of.classification(data.curve, 6)
dev.off()

png(file="g389.png", width=600, height=600)
  mst.of.classification(data.curve, 6)
dev.off()

png(file="g390.png", width=600, height=600)
  op <- par(mfrow=c(2,2),mar=.1+c(0,0,0,0))
  for (k in c(4,6,10,15)) {
    mst.of.classification(data.curve, k, axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g391.png", width=600, height=600)
  op <- par(mfrow=c(2,2),mar=.1+c(0,0,0,0))
  for (k in c(4,6,10,15)) {
    mst.of.classification(data.noisy.curve, k, axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g392.png", width=600, height=600)
  op <- par(mfrow=c(2,2),mar=.1+c(0,0,0,0))
  for (k in c(4,6,10,15)) {
    mst.of.classification(data.broken.line, k, axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g393.png", width=600, height=600)
  op <- par(mfrow=c(2,2),mar=.1+c(0,0,0,0))
  for (k in c(4,6,10,15)) {
    mst.of.classification(data.noisy.broken.line, k, axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g394.png", width=600, height=600)
  op <- par(mfrow=c(2,2),mar=.1+c(0,0,0,0))
  for (k in c(4,6,10,15)) {
    mst.of.classification(data.real, k, axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g395.png", width=600, height=600)
  op <- par(mfrow=c(3,3),mar=.1+c(0,0,0,0))
  for (k in c(4:6)) {
    for (i in 1:3) {
      mst.of.classification(data.real, k, axes=F)
      box()
    }
  }
  par(op)
dev.off()

png(file="g396.png", width=600, height=600)
  op <- par(mfrow=c(3,3),mar=.1+c(0,0,0,0))
  for (k in c(7:9)) {
    for (i in 1:3) {
      mst.of.classification(data.real, k, axes=F)
      box()
    }
  }
  par(op)
dev.off()

png(file="g397.png", width=600, height=600)
  op <- par(mfrow=c(3,3),mar=.1+c(0,0,0,0))
  for (k in c(10:12)) {
    for (i in 1:3) {
      mst.of.classification(data.real, k, axes=F)
      box()
    }
  }
  par(op)
dev.off()

png(file="g398.png", width=600, height=600)
  op <- par(mfrow=c(3,5),mar=.1+c(0,0,0,0))
  for (k in c(13:15)) {
    for (i in 1:3) {
      mst.of.classification(data.real, k, axes=F)
      box()
    }
  }
  par(op)
dev.off()

png(file="g399.png", width=600, height=600)
  library(ape)
  my.plot.mst <- function (d) {
    r <- mst(dist(t(d)))
    d <- prcomp(t(d))$x[,1:2]
    plot(d)
    n <- dim(r)[1]
    w <- which(r!=0)
    i <- as.vector(row(r))[w]
    j <- as.vector(col(r))[w]
    segments( d[i,1], d[i,2], d[j,1], d[j,2], col='red' )
  }
  my.plot.mst(data.broken.line)
dev.off()

png(file="g400.png", width=600, height=600)
  my.plot.mst(data.noisy.broken.line)
dev.off()

png(file="g401.png", width=600, height=600)
  my.plot.mst(data.curve)
dev.off()

png(file="g402.png", width=600, height=600)
  my.plot.mst(data.noisy.curve)
dev.off()

png(file="g403.png", width=600, height=600)
  my.plot.mst(data.real)
dev.off()

png(file="g404.png", width=600, height=600)
  my.plot.mst(t(data.real.3d))
dev.off()

png(file="g405.png", width=600, height=600)
  # Gives the list of the oaths from the branching nodes to the leaves
  chemins.vers.les.feuilles <- function (G) {
    nodes <- which(apply(G,2,sum)>2)
    leaves <- which(apply(G,2,sum)==1)
    res <- list()
    for (a in nodes) {
      for (b in which(G[a,]>0)) {
        if (! b %in% nodes) {
          res <- append(res,list(c(a,b)))
        }
      }
    }
    chemins.vers.les.feuilles.suite(G, nodes, leaves, res)
  }
  # Last coordinate of a vector
  end1 <- function (x) {
    n <- length(x)
    x[n]
  }
  # Last two coordinates of a vector
  end2 <- function (x) {
    n <- length(x)
    x[c(n-1,n)]
  }
  chemins.vers.les.feuilles.suite <- function (G, nodes, leaves, res) {
    new <- list()
    done <- T
    for (ch in res) {
      if( end1(ch) %in% nodes ) {
        # Pass
      } else if ( end1(ch) %in% leaves ) {
        new <- append(new, list(ch))
      } else {
        done <- F
        a <- end2(ch)[1]
        b <- end2(ch)[2]
        for (x in which(G[b,]>0)) {
          if( x != a ){
            new <- append(new, list(c( ch, x )))
          }
        }
      }
    }
    if(done) {
      return(new)
    } else {
      return(chemins.vers.les.feuilles.suite(G,nodes,leaves,new))
    }
  } 

  G <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,0, 0,1,0,0), nr=4)
  chemins.vers.les.feuilles(G)

  # Compute the length of a path
  longueur.chemin <- function (chemin, d) {
    d <- as.matrix(d)
    n <- length(chemin)
    i <- chemin[ 1:(n-1) ]
    j <- chemin[ 2:n ]
    if (n==2) {
      d[i,j]
    } else {
      sum(diag(d[i,][,j]))
    }
  }

  G <- matrix(c(0,1,0,0, 1,0,1,1, 0,1,0,0, 0,1,0,0), nr=4)
  d <- matrix(c(0,2,4,3, 2,0,2,1, 4,2,0,3, 3,1,3,0), nr=4)
  chemins <- chemins.vers.les.feuilles(G)
  chemins
  l <- sapply(chemins, longueur.chemin, d)
  l
  stopifnot( l == c(2,2,1) )

  elague <- function (G0, d0) {
    d0 <- as.matrix(d0)
    G <- G0
    d <- d0
    res <- 1:dim(d)[1]
    chemins <- chemins.vers.les.feuilles(G)
    while (length(chemins)>0) {
      longueurs <- sapply(chemins, longueur.chemin, d)
      # Number of the shortest path
      i <- which( longueurs == min(longueurs) )[1]
      cat(paste("Removing", paste(res[chemins[[i]]],collapse=' '), "length =", longueurs[i],"\n"))
      # Nodes to remove
      j <- chemins[[i]] [-1]
      res <- res[-j]
      G <- G[-j,][,-j]
      d <- d[-j,][,-j]
      cat(paste("Removing", paste(j), "\n" ))
      cat(paste("", paste(res, collapse=' '), "\n"))
      chemins <- chemins.vers.les.feuilles(G)      
    }
    res
  }

  library(ape)
  my.plot.mst <- function (x0) {
    cat("Plotting the points\n")
    x <- prcomp(t(x0))$x[,1:2]
    plot(x)
    cat("Computing the distance matrix\n")
    d <- dist(t(x0))
    cat("Computing the MST (Minimum Spanning Tree)\n")
    G <- mst(d)
    cat("Plotting the MST\n")
    n <- dim(G)[1]
    w <- which(G!=0)
    i <- as.vector(row(G))[w]
    j <- as.vector(col(G))[w]
    segments( x[i,1], x[i,2], x[j,1], x[j,2], col='red' )
    cat("Pruning the tree\n")
    k <- elague(G,d)
    cat("Plotting the pruned tree\n")
    x <- x[k,]
    G <- G[k,][,k]
    n <- dim(G)[1]
    w <- which(G!=0)
    i <- as.vector(row(G))[w]
    j <- as.vector(col(G))[w]
    segments( x[i,1], x[i,2], x[j,1], x[j,2], col='red', lwd=3 )        
  }

  my.plot.mst(data.noisy.broken.line)
dev.off()

png(file="g406.png", width=600, height=600)
  my.plot.mst(data.noisy.curve)
dev.off()

png(file="g407.png", width=600, height=600)
  my.plot.mst(data.real)
dev.off()

png(file="g408.png", width=600, height=600)
  my.plot.mst(t(data.real.3d))
dev.off()

png(file="g409.png", width=600, height=600)
  # k: each point is linked to its k nearest neighbors
  # eps: each point is linked to all its neighbors within a radius eps

  isomap.incidence.matrix <- function (d, eps=NA, k=NA) {
    stopifnot(xor( is.na(eps), is.na(k) ))
    d <- as.matrix(d)
    if(!is.na(eps)) {
      im <- d <= eps
    } else {
      im <- apply(d,1,rank) <= k+1
      diag(im) <- F
    }
    im | t(im)
  }
  plot.graph <- function (im,x,y=NULL, ...) {
    if(is.null(y)) {
      y <- x[,2]
      x <- x[,1]
    }
    plot(x,y, ...)
    k <- which(  as.vector(im)  )
    i <- as.vector(col(im))[ k ]
    j <- as.vector(row(im))[ k ]
    segments( x[i], y[i], x[j], y[j], col='red' )
  }

  d <- dist(t(data.noisy.curve))
  r <- princomp(t(data.noisy.curve))
  x <- r$scores[,1]
  y <- r$scores[,2]

  plot.graph(isomap.incidence.matrix(d, k=5), x, y)
dev.off()

png(file="g410.png", width=600, height=600)
  plot.graph(isomap.incidence.matrix(d, eps=quantile(as.vector(d), .05)), 
             x, y)
dev.off()

png(file="g411.png", width=600, height=600)
  isomap.incidence.matrix <- function (d, eps=NA, k=NA) {
    stopifnot(xor( is.na(eps), is.na(k) ))
    d <- as.matrix(d)
    if(!is.na(eps)) {
      im <- d <= eps
    } else {
      im <- apply(d,1,rank) <= k+1
      diag(im) <- F
    }
    im | t(im) | mst(d)
  }
  plot.graph(isomap.incidence.matrix(d, eps=quantile(as.vector(d), .05)), 
             x, y)
dev.off()

png(file="g412.png", width=600, height=600)
  inf <- function (x,y) { ifelse(x<y,x,y) }
  isomap.distance <- function (im, d) {
    d <- as.matrix(d)
    n <- dim(d)[1]
    dd <- ifelse(im, d, Inf)
    for (k in 1:n) {
      dd <- inf(dd, matrix(dd[,k],nr=n,nc=n) + matrix(dd[k,],nr=n,nc=n,byrow=T))
    }
    dd
  }

  isomap <- function (x, d=dist(x), eps=NA, k=NA) {
    im <- isomap.incidence.matrix(d, eps, k)
    dd <- isomap.distance(im,d)
    r <- list(x,d,incidence.matrix=im,distance=dd)
    class(r) <- "isomap"
    r    
  }

  r <- isomap(t(data.noisy.curve), k=5)
  xy <- cmdscale(r$distance,2)   # long: around 30 seconds
  plot.graph(r$incidence.matrix, xy)
dev.off()

png(file="g413.png", width=600, height=600)
  plot.graph(r$incidence.matrix, xy, ylim=range(xy))
dev.off()

png(file="g414.png", width=600, height=600)
  library(pixmap)
  x <- read.pnm("photo1.ppm")
  d <- cbind( as.vector(x@red),
              as.vector(x@green),
              as.vector(x@blue) )
  m <- apply(d,2,mean)
  d <- t(  t(d) - m )
  s <- apply(d,2,sd)
  d <- t(  t(d) / s )
  library(som)
  r1 <- som(d,5,5)
  plot(r1)
dev.off()

png(file="g415.png", width=600, height=600)
  x <- r1$code.sum$x
  y <- r1$code.sum$y
  n <- r1$code.sum$nobs
  co <- r1$code   # Is it in the same order that x, y and n?
  co <- t( t(co) * s + m )
  plot(x, y, 
       pch=15,
       cex=5,      
       col=rgb(co[,1], co[,2], co[,3])
      )
dev.off()

png(file="g416.png", width=600, height=600)
  x <- r1$code.sum$x
  y <- r1$code.sum$y
  n <- r1$code.sum$nobs
  co <- r1$code   # Is it in the same order that x, y and n?
  co <- t( t(co) * s + m )
  plot(x, y, 
       pch=15,
       cex=5*n/max(n),      
       col=rgb(co[,1], co[,2], co[,3])
      )
dev.off()

png(file="g417.png", width=600, height=600)
  library(class)
  r2 <- SOM(d)
  plot(r2)
dev.off()

png(file="g418.png", width=600, height=600)
  x <- r2$grid$pts[,1]
  y <- r2$grid$pts[,2]
  n <- 1   # Where???
  co <- r2$codes
  co <- t( t(co) * s + m )
  plot(x, y, 
       pch=15,
       cex=5*n/max(n),
       col=rgb(co[,1], co[,2], co[,3])
      )
dev.off()

png(file="g419.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    r2 <- SOM(d) 
    plot(r2)
  }
  par(op)
dev.off()

png(file="g420.png", width=600, height=600)
  x <- r1$code.sum$x
  y <- r1$code.sum$y
  v <- r1$code

  op <- par(mfrow=c(2,2), mar=c(1,1,1,1))
  for (k in 1:dim(v)[2]) {
    m <- matrix(NA, nr=max(x), nc=max(y))
    for (i in 1:length(x)) {
      m[ x[i], y[i] ] <- v[i,k]
    }
    image(m, col=rainbow(255), axes=F)
  }
  par(op)
dev.off()

png(file="g421.png", width=600, height=600)
  x <- r2$grid$pts[,1]
  y <- r2$grid$pts[,2]
  v <- r2$codes

  op <- par(mfrow=c(2,2), mar=c(1,1,1,1))
  for (k in 1:dim(v)[2]) {
    m <- matrix(NA, nr=max(x), nc=max(y))
    for (i in 1:length(x)) {
      m[ x[i], y[i] ] <- v[i,k]
    }
    image(m, col=rainbow(255), axes=F)
  }
  par(op)
dev.off()

png(file="g422.png", width=600, height=600)
  library(MASS)
  data(HairEyeColor)
  x <- HairEyeColor[,,1]+HairEyeColor[,,2]
  biplot(corresp(x, nf = 2))
dev.off()

png(file="g423.png", width=600, height=600)
  biplot(corresp(t(x), nf = 2))
dev.off()

png(file="g424.png", width=600, height=600)
  # ???
  plot(corresp(x, nf=1))
dev.off()

png(file="g425.png", width=600, height=600)
  n <- 100
  m <- matrix(sample(c(T,F),n^2,replace=T), nr=n, nc=n)
  biplot(corresp(m, nf=2), main="Correspondance Analysis of Random Data")
dev.off()

png(file="g426.png", width=600, height=600)
  vp <- corresp(m, nf=100)$cor
  plot(vp, ylim=c(0,max(vp)), type='l', 
       main="Correspondance Analysis of Random Data")
dev.off()

png(file="g427.png", width=600, height=600)
  n <- 100
  x <- matrix(1:n, nr=n, nc=n, byrow=F)
  y <- matrix(1:n, nr=n, nc=n, byrow=T)
  m <- abs(x-y) <= n/10
  biplot(corresp(m, nf=2), 
         main='Correspondance Analysis of "Band" Data')
dev.off()

png(file="g428.png", width=600, height=600)
  vp <- corresp(m, nf=100)$cor
  plot(vp, ylim=c(0,max(vp)), type='l', 
         main='Correspondance Analysis of "Band" Data')
dev.off()

png(file="g429.png", width=600, height=600)
  n <- 100
  x <- matrix(1:n, nr=n, nc=n, byrow=F)
  y <- matrix(1:n, nr=n, nc=n, byrow=T)
  m <- abs(x-y) <= n/10
  plot.boolean.matrix <- function (m) { # Voir aussi levelplot
    nx <- dim(m)[1]
    ny <- dim(m)[2]
    x <- matrix(1:nx, nr=nx, nc=ny, byrow=F)
    y <- matrix(1:ny, nr=nx, nc=ny, byrow=T)
    plot( as.vector(x)[ as.vector(m) ], as.vector(y)[ as.vector(m) ], pch=16 )
  }
  plot.boolean.matrix(m)
dev.off()

png(file="g430.png", width=600, height=600)
  ox <- sample(1:n, n, replace=F)
  oy <- sample(1:n, n, replace=F)
  reorder.matrix <- function (m,ox,oy) {
    m <- m[ox,]
    m <- m[,oy]
    m
  }
  m2 <- reorder.matrix(m,ox,oy)
  plot.boolean.matrix(m2)
dev.off()

png(file="g431.png", width=600, height=600)
  a <- corresp(m2)
  o1 <- order(a$rscore)
  o2 <- order(a$cscore)
  m3 <- reorder.matrix(m2,o1,o2)
  plot.boolean.matrix(m3)
dev.off()

png(file="g432.png", width=600, height=600)
  n <- 100
  p <- .05
  done <- F 
  while( !done ){
    # We often get singular matrices
    m2 <- matrix( sample(c(F,T), n*n, replace=T, prob=c(1-p, p)), nr=n, nc=n )
    done <- det(m2) != 0
  }
  plot.boolean.matrix(m2)
dev.off()

png(file="g433.png", width=600, height=600)
  a <- corresp(m2)
  o1 <- order(a$rscore)
  o2 <- order(a$cscore)
  m3 <- reorder.matrix(m2,o1,o2)
  plot.boolean.matrix(m3)
dev.off()

png(file="g434.png", width=600, height=600)
  # Correspondance Analysis
  my.ac <- function (x) {
    if(any(x<0))
      stop("Expecting a contingency matrix -- with no negative entries")
    x <- x/sum(x)
    nr <- dim(x)[1]
    nc <- dim(x)[2]
    marges.lignes <- apply(x,2,sum)
    marges.colonnes <- apply(x,1,sum)
    profils.lignes <- x / matrix(marges.lignes, nr=nr, nc=nc, byrow=F)
    profils.colonnes <- x / matrix(marges.colonnes, nr=nr, nc=nc, byrow=T)
    # Do not forget to center the matrix: we compute the frequency matrix
    # we would have if the variables were independant and we take the difference.
    x <- x - outer(marges.colonnes, marges.lignes)
    e1 <- eigen( t(x) %*% diag(1/marges.colonnes) %*% x %*% diag(1/marges.lignes) )
    e2 <- eigen( x %*% diag(1/marges.lignes) %*% t(x) %*% diag(1/marges.colonnes) )
    v.col <- solve( e2$vectors, x )
    v.row <- solve( e1$vectors, t(x) )
    v.col <- t(v.col)
    v.row <- t(v.row)
    if(nr<nc)
      valeurs.propres <- e1$values
    else
      valeurs.propres <- e2$values
    # Dessin
    plot( v.row[,1:2], 
          xlab='', ylab='', frame.plot=F )
    par(new=T)
    plot( v.col[,1:2], col='red',
          axes=F, xlab='', ylab='', pch='+')
    axis(3, col='red')
    axis(4, col='red')
    # Return the data
    invisible(list(donnees=x, colonnes=v.col, lignes=v.row, 
                   valeurs.propres=valeurs.propres))
  }    

  nr <- 3
  nc <- 5
  x <- matrix(rpois(nr*nc,10), nr=nr, nc=nc)
  my.ac(x)
dev.off()

png(file="g435.png", width=600, height=600)
  plot(corresp(x,nf=2))
dev.off()

png(file="g436.png", width=600, height=600)
  library(MASS)
  data(farms)
  farms.mca <- mca(farms, abbrev=TRUE)
  farms.mca
  plot(farms.mca)
dev.off()

png(file="g437.png", width=600, height=600)
  # Not pretty
  my.table.to.data.frame <- function (a) {
    r <- NULL
    d <- as.data.frame.table(a)
    n1 <- dim(d)[1]
    n2 <- dim(d)[2]-1
    for (i in 1:n1) {
      for (j in 1:(d[i,n2+1])) {
        r <- rbind(r, d[i,1:n2])
      }
      row.names(r) <- 1:dim(r)[1]
    }
    r
  }
  r <- my.table.to.data.frame(HairEyeColor)
  plot(mca(r))
dev.off()

png(file="g438.png", width=1000, height=500)
  x <- HairEyeColor[,,1]+HairEyeColor[,,2]
  op <- par(mfcol=c(1,2))
  biplot(corresp(x, nf = 2), 
         main="Simple Correspondance Analysis")  
  plot(mca(my.table.to.data.frame(x)), rows=F, 
       main="Multiple Correspondance Analysis")
  par(op)
dev.off()

png(file="g439.png", width=600, height=600)
  # Multiple Correspondance Analysis
  tableau.disjonctif.vecteur <- function (x) {
    y <- matrix(0, nr=length(x), nc=length(levels(x)))
    for (i in 1:length(x)) {
      y[i, as.numeric(x[i])] <- 1
    }
    y
  }
  tableau.disjonctif <- function (x) {
    if( is.vector(x) )
      y <- tableau.disjonctif.vecteur(x)
    else {
      y <- NULL
      y.names <- NULL
      for (i in 1:length(x)) {
        y <- cbind(y, tableau.disjonctif.vecteur(x[,i]))
        y.names <- c( y.names, paste(names(x)[i], levels(x[,i]), sep='') )
      }
    }
    colnames(y) <- y.names
    y
  }
  my.acm <- function (x, garder.un=F) {
    # x is a data.frame that contains only factors
    # y is a boolean matrix (it only contains 0s and 1s)
    y <- tableau.disjonctif(x)
    # Number of observations
    n <- dim(y)[1]
    # Number of variables
    s <- length(x)
    # Number of columns in the disjunctive table
    p <- dim(y)[2]
    # The matrix and the weights
    F <- y/(n*s)
    Dp <- diag(t(y)%*%y) / (n*s)
    Dn <- rep(1/n,n)
    # Let us perform the analysis
    # Do NOT forget to remove 1 as an eigenvalue!!!
    # (it comes from the fact that we did not center the matrix)
    e1 <- eigen( t(F) %*% diag(1/Dn) %*% F %*% diag(1/Dp) )
    e2 <- eigen( F %*% diag(1/Dp) %*% t(F) %*% diag(1/Dn) )
    variables <- t(e2$vectors) %*% F
    individus <- t(e1$vectors) %*% t(F)
    variables <- t(variables)
    individus <- t(individus)
    valeurs.propres <- e1$values
    if( !garder.un ) {
      variables <- variables[,-1]
      individus <- individus[,-1]
      valeurs.propres <- valeurs.propres[-1]
    }
    plot( jitter(individus[,2], factor=5) ~ jitter(individus[,1], factor=5), 
          xlab='', ylab='', frame.plot=F )
    par(new=T)
    plot( variables[,1:2], col='red',
          axes=FALSE, xlab='', ylab='', type='n' )
    text( variables[,1:2], colnames(y), col='red')
    axis(3, col='red')
    axis(4, col='red')
    j <- 1
    col = rainbow(s)
    for (i in 1:s) {
      jj <- j + length(levels(x[,i])) - 1
      print( paste(j,jj) )
      lines( variables[j:jj,1], variables[j:jj,2], col=col[i] )
      j <- jj+1
    }
    invisible(list( donnees=x, variables=variables, individus=individus, 
                    valeurs.propres=valeurs.propres ))
  }
  
  random.data.1 <- function () {
    n <- 100
    m <- 3
    l <- c(3,2,5)
    x <- NULL    
    for (i in 1:m) {
      v <- factor( sample(1:l[i], n, replace=T), levels=1:l[i] )
      if( is.null(x) )
        x <- data.frame(v)
      else 
        x <- data.frame(x, v)
    }
    names(x) <- LETTERS[1:m]
    x
  }
  r <- NULL
  while(is.null(r)) {
    x1 <- random.data.1()
    try(r <- my.acm(x1))
  }
dev.off()

png(file="g440.png", width=600, height=600)
  my.random.data.2 <- function () {
    n <- 100
    m <- 3
    l <- c(3,2,3)
    x <- NULL    
    for (i in 1:m) {
      v <- factor( sample(1:l[i], n, replace=T), levels=1:l[i] )
      if( is.null(x) )
        x <- data.frame(v)
      else 
        x <- data.frame(x, v)
    }
    x[,3] <- factor( ifelse( runif(n)>.8, x[,1], x[,3] ), levels=1:l[1])
    names(x) <- LETTERS[1:m]
    x
  }
  r <- NULL
  while(is.null(r)) {
    x2 <- my.random.data.2()
    try(r <- my.acm(x2))
  }
dev.off()

png(file="g441.png", width=600, height=600)
  library(MASS)
  plot(mca(x1))
dev.off()

png(file="g442.png", width=600, height=600)
  plot(mca(x2))
dev.off()

png(file="g443.png", width=600, height=600)
  my.acm(x1, garder.un=T)
dev.off()

png(file="g444.png", width=600, height=600)
  my.acm(x2, garder.un=T)
dev.off()

png(file="g445.png", width=600, height=600)
  tableau.disjonctif.vecteur <- function (x) {
    if( is.factor(x) ){
      y <- matrix(0, nr=length(x), nc=length(levels(x)))
      for (i in 1:length(x)) {
      y[i, as.numeric(x[i])] <- 1
      }
      return(y)
    } else 
    return(x)
  }
  tableau.disjonctif <- function (x) {
    if( is.vector(x) )
      y <- tableau.disjonctif.vecteur(x)
    else {
      y <- NULL
      y.names <- NULL
      for (i in 1:length(x)) {
        y <- cbind(y, tableau.disjonctif.vecteur(x[,i]))
        if( is.factor(x[,i]) )
          y.names <- c( y.names, paste(names(x)[i], levels(x[,i]), sep='') )
        else 
          y.names <- c(y.names, names(x)[i])
      }
    }
    colnames(y) <- y.names
    y
  }
  my.am <- function (x) {
    y <- tableau.disjonctif(x)
    # Number of observations
    n <- dim(y)[1]
    # Number of variables
    s <- length(x)
    # Number of columns of the disjunctive table
    p <- dim(y)[2]

    # The matrix and the weights
    # Each variable weight is 1, but for the qualitative variables
    # that are "split" into several "subvariables", it is the sum of
    # those "subvariable" weights that is 1.
    Dp <- diag(t(y)%*%y) / n
    j <- 1
    for(i in 1:s){
      if( is.factor(x[,i]) )
        j <- j + length(levels(x[,i]))
      else {
        Dp[j] <- 1
        j <- j+1
      }      
    }
    Dn <- rep(1,n)
    # It if not a good idea to center the variable...
    f <- ne.pas.centrer(y)
    # We perform the analysis...
    e1 <- eigen( t(f) %*% diag(1/Dn) %*% f %*% diag(1/Dp) )
    e2 <- eigen( f %*% diag(1/Dp) %*% t(f) %*% diag(1/Dn) )
    variables <- t(e2$vectors) %*% f
    individus <- t(e1$vectors) %*% t(f)
    variables <- t(variables)
    individus <- t(individus)
    valeurs.propres <- e1$values
    # Sometimes, because of rounding errors, the matrix becomes
    # non-diagonalizable in R. We then end up with complex eigenvalues
    # and (worse) complex eigenvectors. That is why it is a better
    # idea to use the SVD, that always yields real results. 
    if( any(Im(variables)!=0) | any(Im(individus)!=0) | 
        any(Im(valeurs.propres)!=0) ){
      warning("Matrix not diagonalizable on R!!!")
      variables <- Re(variables)
      individus <- Re(individus)
      valeurs.propres <- Re(valeurs.propres)
    }
    plot( jitter(individus[,2], factor=5) ~ jitter(individus[,1], factor=5), 
          xlab='', ylab='', frame.plot=F )
    par(new=T)
    plot( variables[,1:2], col='red',
          axes=FALSE, xlab='', ylab='', type='n' )
    text( variables[,1:2], colnames(y), col='red')
    axis(3, col='red')
    axis(4, col='red')
    col = rainbow(s)
    j <- 1
    for (i in 1:s) {
      if( is.factor(x[,i]) ){
        jj <- j + length(levels(x[,i])) - 1
        print( paste(j,jj) )
        lines( variables[j:jj,1], variables[j:jj,2], col=col[i] )
        j <- jj+1
      } else {
        arrows(0,0,variables[j,1],variables[j,2],col=col[i])
        j <- j+1
      }
    }
  }
  ne.pas.centrer <- function (y) { y }
  
  n <- 500
  m <- 3
  p <- 2
  l <- c(3,2,5)
  x <- NULL    
  for (i in 1:m) {
    v <- factor( sample(1:l[i], n, replace=T), levels=1:l[i] )
    if( is.null(x) )
      x <- data.frame(v)
    else 
      x <- data.frame(x, v)
  }
  x <- cbind( x, matrix( rnorm(n*p), nr=n, nc=p ) )
  names(x) <- LETTERS[1:(m+p)]
  x1 <- x
  my.am(x1)
dev.off()

png(file="g446.png", width=600, height=600)
  n <- 500
  m <- 3
  p <- 2
  l <- c(3,2,5)
  x <- NULL    
  for (i in 1:m) {
    v <- factor( sample(1:l[i], n, replace=T), levels=1:l[i] )
    if( is.null(x) )
      x <- data.frame(v)
    else 
      x <- data.frame(x, v)
  }
  x <- cbind( x, matrix( rnorm(n*p), nr=n, nc=p ) )
  names(x) <- LETTERS[1:(m+p)]
  x[,3] <- factor( ifelse( runif(n)>.8, x[,1], x[,3] ), levels=1:l[1])
  x[,5] <- scale( ifelse( runif(n)>.9, as.numeric(x[,1]), as.numeric(x[,2]) )) + .1*rnorm(n)
  x2 <- x
  my.am(x2)
dev.off()

png(file="g447.png", width=600, height=600)
  ne.pas.centrer <- function (y) {
    centre <- apply(y, 2, mean)
    y - matrix(centre, nr=n, nc=dim(y)[2], byrow=T)
  }
  my.am(x1)
dev.off()

png(file="g448.png", width=600, height=600)
  my.am(x2)
dev.off()

png(file="g449.png", width=600, height=600)
  to.factor.numeric.vector <- function (x, number) {
    resultat <- NULL
    intervalles <- co.intervals(x,number,overlap=0)
    for (i in 1:number) {
      if( i==1 ) intervalles[i,1] = min(x)
      else
        intervalles[i,1] <- intervalles[i-1,2]
      if( i==number )
        intervalles[i,2] <- max(x)
    }
    for (valeur in x) {
      r <- NA
      for (i in 1:number) {
        if( valeur >= intervalles[i,1] & valeur <= intervalles[i,2] )
          r <- i
      }
      resultat <- append(resultat, r)
    }
    factor(resultat, levels=1:number)
  }
  to.factor.vector <- function (x, number) {
    if( is.factor(x) ) 
      return(x)
    else
      return(to.factor.numeric.vector(x,number))
  }
  to.factor <- function (x, number=4 ) {
    y <- NULL
    for (a in x) {
      aa <- to.factor.vector(a, number)
      if( is.null(y) )
        y <- data.frame(aa)
      else 
        y <- data.frame(y,aa)
    }
    names(y) <- names(x)
    y
  }
  my.am(to.factor(x1))
dev.off()

png(file="g450.png", width=600, height=600)
  my.am(to.factor(x1, number=3))
dev.off()

png(file="g451.png", width=600, height=600)
  my.am(to.factor(x2))
dev.off()

png(file="g452.png", width=600, height=600)
  my.am(to.factor(x2, number=3))
dev.off()

png(file="g453.png", width=600, height=600)
  library(MASS)
  n <- 100
  k <- 5
  x1 <- runif(k,-5,5) + rnorm(k*n)
  x2 <- runif(k,-5,5) + rnorm(k*n)
  x3 <- runif(k,-5,5) + rnorm(k*n)
  x4 <- runif(k,-5,5) + rnorm(k*n)
  x5 <- runif(k,-5,5) + rnorm(k*n)
  y <- factor(rep(1:5,n))
  plot(lda(y~x1+x2+x3+x4+x5))
dev.off()

png(file="g454.png", width=600, height=600)
  plot(lda(y~x1+x2+x3+x4+x5), dimen=2) 
dev.off()

png(file="g455.png", width=600, height=600)
  op <- par(mar=c(4,4,0,2)+.1)
  plot(lda(y~x1))
  par(op)
dev.off()

png(file="g456.png", width=600, height=600)
  plot(lda(y~x1+x2+x3+x4+x5)$svd, type='h')
dev.off()

png(file="g457.png", width=600, height=600)
  # Data
  n <- 200  # Number of patients, number of columns
  k <- 10   # Dimension of the ambient space
  nb.points <- 5
  p <- matrix( 5*rnorm(nb.points*k), nr=k )
  library(splines)
  barycentre2 <- function (y,n) {
    m <- length(y)
    x <- 1:m
    r <- interpSpline(x,y)    
    #r <- lm( y ~ bs(x, knots=m) )
    predict(r, seq(1,m,length=n))$y
  }
  data.curve <- apply(p, 1, barycentre2, n)
  data.curve <- t(data.curve)
  pairs(t(data.curve))
dev.off()

png(file="g458.png", width=600, height=600)
  m <- data.curve
  mm <- apply(m, 2, function (x) { x %o% x } )
  r <- princomp(t(mm))
  plot(r)
dev.off()

png(file="g459.png", width=600, height=600)
  pairs(r$scores[,1:5])
dev.off()

png(file="g460.png", width=600, height=600)
  # Degree 1 kernel ("noyau" is the French word for "kernel")
  noyau1 <- function (x,y) { sum(x*y) }
  m <- data.curve
  m <- t(t(m) - apply(m,2,mean))
  k <- dim(m)[1]
  wrapper <- function(x, y, my.fun, ...) {
    sapply(seq(along=x), FUN = function(i) my.fun(x[i], y[i], ...))
  }
  mm <- outer(1:k, 1:k, wrapper, function (i,j) { noyau1(m[,i],m[,j]) })

  # Degree 2 kernel
  noyau2 <- function (x,y) {
    a <- x*y
    n <- length(a)
    i <- gl(n,1,n^2)
    j <- gl(n,n,n^2)
    i <- as.numeric(i)
    j <- as.numeric(j)
    w <- which( i <= j & j <= k )
    i <- i[w]
    j <- j[w]
    sum(a[i]*a[j])
  }
  stopifnot( noyau2(1:2,2:1) == 12 )

  # Degree 3 kernel
  noyau3 <- function (x,y) {
    a <- x*y
    n <- length(a)
    i <- gl(n,1,n^3)
    j <- gl(n,n,n^3)
    k <- gl(n,n^2,n^3)
    i <- as.numeric(i)
    j <- as.numeric(j)
    k <- as.numeric(k)
    w <- which( i <= j & j <= k )
    i <- i[w]
    j <- j[w]
    k <- k[w]
    sum(a[i]*a[j]*a[k])
  }
  stopifnot( noyau3(1:2,2:1) == 32 )
  
  wrapper <- function(x, y, my.fun, ...) {
    sapply(seq(along=x), FUN = function(i) my.fun(x[i], y[i], ...))
  }
  k <- dim(m)[1]
  mm <- outer(1:k, 1:k, FUN=wrapper, my.fun = function (i,j) { noyau3(m[,i],m[,j]) })

  r <- princomp(covmat=mm)
  plot(r)
dev.off()

png(file="g461.png", width=600, height=600)
  #pairs(r$scores[,1:5])
  plot((t(m) %*% r$loadings) [,1:2])
dev.off()

png(file="g462.png", width=600, height=600)
  data.noisy.curve <- data.curve + rnorm(length(data.curve))
  k <- dim(data.noisy.curve)[1]
  mm <- outer(1:k, 1:k, FUN=wrapper, my.fun = function (i,j) { 
    noyau3(data.noisy.curve[,i],data.noisy.curve[,j]) 
  })
  r <- princomp(t(data.noisy.curve), covmat=mm)
  plot(r)
dev.off()

png(file="g463.png", width=600, height=600)
  pairs((t(m) %*% r$loadings) [,1:5])
dev.off()

png(file="g464.png", width=600, height=600)
  plot((t(m) %*% r$loadings) [,1:2])
  lines((t(m) %*% r$loadings) [,1:2], col='red')
dev.off()

png(file="g465.png", width=600, height=600)
  x <- 1:dim(data.noisy.curve)[2]
  y1 <- (t(m) %*% r$loadings) [,1]
  y2 <- (t(m) %*% r$loadings) [,2]
  r1 <- loess(y1~x)
  r2 <- loess(y2~x)
  plot(y1,y2)
  lines(y1,y2, col='red')
  lines(predict(r1), predict(r2), col='blue', lwd=3)
dev.off()

png(file="g466.png", width=600, height=600)
  noyau <- function (x,y) {
    noyau1(x,y) + noyau2(x,y) + noyau3(x,y)
  }
  mm <- outer(1:k, 1:k, FUN=wrapper, my.fun = function (i,j) { 
    noyau(data.curve[,i],data.curve[,j]) 
  })
  r <- princomp(covmat=mm)
  plot(r)
dev.off()

png(file="g467.png", width=600, height=600)
  pairs((t(m) %*% r$loadings) [,1:5])
dev.off()

png(file="g468.png", width=600, height=600)
  perceptron_learn <- function (input, output,
                                max_iterations = 100) {
    stopifnot( is.matrix(input),
               is.matrix(output),
               is.numeric(input),
               is.logical(output) )
    stopifnot( dim(input)[2] == dim(output)[2] )
    input <- rbind(input, rep(1, N)) # Biases
    N <- dim(input)[2]
    dim_input  <- dim(input)[1]
    dim_output <- dim(output)[1]
    W <- matrix(rnorm(dim_input) * rnorm(dim_output),
                nc = dim_input, nr = dim_output)
    hardlim <- function (x) ifelse( x > 0, 1, 0 )
    finished <- FALSE
    remaining_iterations <- max_iterations
    while (! finished) {
      W_previous <- W
      for (i in 1:N) {
        forecast <- hardlim( W %*% input[,i] )
        error <- output[,i] - forecast
        W <- W + error %*% input[,i]
      }
      remaining_iterations <- remaining_iterations - 1
      finished <- remaining_iterations <= 0 | 
                  sum(abs(W - W_previous)) < 1e-6
    }
    attr(W, "converged") <- all(W == W_previous)
    dimnames(W) <- list(output=NULL, input=NULL)
    class(W) <- "perceptron"
    W
  }
  perceptron_predict <- function (W, input) {
    input <- as.matrix(input)
    N <- dim(input)[2]
    input <- rbind(input, rep(1, N)) # Biases
    stopifnot( dim(W)[2] == dim(input)[1] )
    W %*% input
  }

  # Test
  k <- 2   # BUG: It does not work with other values...
  N <- 50
  set.seed(1)
  centers <- matrix(rnorm(2*k), nc=2)
  output <- sample(1:k, N, replace=T)
  input <- t(matrix(
    rnorm(k*N, mean=centers[output,], sd=.1), 
             nc=2))
  output <- t(output == 1)
  w <- perceptron_learn(input, output)

  plot(t(input), col=1+output,
       xlab="", ylab="", main="Perceptron learning")
  abline(-w[3]/w[2], -w[1]/w[2], lwd=3, lty=2)
dev.off()

png(file="g469.png", width=600, height=600)
  library(MASS)
  z <- output[1,]
  x <- input[1,]
  y <- input[2,]
  r <- lda( z ~ x + y )
  n <- 200
  x <- rep(seq(-2,2,length=n), each=n)
  y <- rep(seq(-2,2,length=n), n)
  z <- predict(r, data.frame(x, y))
  z <- c(rgb(.7,.7,.7), rgb(1,.7,.7))[ as.numeric(z$class) ]
  plot(x, y, col=z, pch=15,
       xlab="", ylab="", 
       main="Non-converging perceptron and LDA")
  points(t(input), col = 1 + output)
  abline(-w[3]/w[2], -w[1]/w[2], lwd=3, lty=2)
dev.off()

png(file="g470.png", width=600, height=600)
  library(e1071)
  y <- as.factor(output[1,])
  x <- t(input)
  r <- svm(x, y)
  n <- 200
  x <- cbind( rep(seq(-2,2,length=n), each=n),
              rep(seq(-2,2,length=n), n)       )
  z <- predict(r, x)
  z <- c(rgb(.7,.7,.7), rgb(1,.7,.7))[ as.numeric(z) ]
  plot(x, col=z, pch=15,
       xlab="", ylab="", 
       main="Non-converging perceptron and SVM")
  points(t(input), col = 1 + output)
  abline(-w[3]/w[2], -w[1]/w[2], lwd=3, lty=2)
dev.off()

png(file="g471.png", width=600, height=600)
  n <- 200  # Number of patients, number of columns
  k <- 10   # Dimension of the ambient space
  nb.points <- 5
  p <- matrix( 5*rnorm(nb.points*k), nr=k )
  barycentre <- function (x, n) {
    # Add number between the values of x in order to get a length n vector
    i <- seq(1,length(x)-.001,length=n)
    j <- floor(i)
    l <- i-j
    (1-l)*x[j] + l*x[j+1]
  }
  m <- apply(p, 1, barycentre, n)
  data.broken.line <- t(m)
  pairs(t(data.broken.line))
dev.off()

png(file="g472.png", width=600, height=600)
  plot(princomp(t(data.broken.line)))
dev.off()

png(file="g473.png", width=600, height=600)
  pairs(princomp(t(data.broken.line))$scores[,1:5])
dev.off()

png(file="g474.png", width=600, height=600)
  data.noisy.broken.line <- data.broken.line + rnorm(length(data.broken.line))
  pairs(princomp(t(data.noisy.broken.line))$scores[,1:5])
dev.off()

png(file="g475.png", width=600, height=600)
  library(splines)
  barycentre2 <- function (y,n) {
    m <- length(y)
    x <- 1:m
    r <- interpSpline(x,y)    
    #r <- lm( y ~ bs(x, knots=m) )
    predict(r, seq(1,m,length=n))$y
  }

  k <- 5
  y <- sample(1:100,k)
  x <- seq(1,k,length=100)
  plot(barycentre2(y,100) ~ x)
  lines(y, col='red', lwd=3, lty=2)
dev.off()

png(file="g476.png", width=600, height=600)
  data.curve <- apply(p, 1, barycentre2, n)
  data.curve <- t(data.curve)
  pairs(t(data.curve))
dev.off()

png(file="g477.png", width=600, height=600)
  plot(princomp(t(data.curve)))
dev.off()

png(file="g478.png", width=600, height=600)
  pairs(princomp(t(data.curve))$scores[,1:5])
dev.off()

png(file="g479.png", width=600, height=600)
  data.noisy.curve <- data.curve + rnorm(length(data.curve))
  pairs(t(data.noisy.curve))
dev.off()

png(file="g480.png", width=600, height=600)
  plot(princomp(t(data.noisy.curve)))
dev.off()

png(file="g481.png", width=600, height=600)
  pairs(princomp(t(data.noisy.curve))$scores[,1:5])
dev.off()

png(file="g482.png", width=600, height=600)
  random.rotation.matrix <- function (n=3) {
    m <- NULL
    for (i in 1:n) {
      x <- rnorm(n)
      x <- x / sqrt(sum(x*x))
      y <- rep(0,n)
      if (i>1) for (j in 1:(i-1)) {
        y <- y + sum( x * m[,j] ) * m[,j]
      }
      x <- x - y
      x <- x / sqrt(sum(x*x))
      m <- cbind(m, x)
    }
    m
  }

  n <- 200
  k <- 10
  x <- seq(0,2*pi,length=n)
  data.circle <- matrix(0, nr=n, nc=k)
  data.circle[,1] <- cos(x)    
  data.circle[,2] <- sin(x)
  data.circle <- data.circle %*% random.rotation.matrix(k)
  data.circle <- t( t(data.circle) + rnorm(k) )
  pairs(data.circle[,1:3])
dev.off()

png(file="g483.png", width=600, height=600)
  data.circle <- data.circle + .1*rnorm(n*k)
  pairs(data.circle[,1:3])
dev.off()

png(file="g484.png", width=600, height=600)
  pairs(princomp(data.circle)$scores[,1:3])
dev.off()

png(file="g485.png", width=600, height=600)
  data.real <- read.table("Tla_z.txt", sep=",")
  data.real.group <- factor(substr(names(data.real),0,1))
  r <- prcomp(t(data.real))
  plot(r$sdev, type='h')
dev.off()

png(file="g486.png", width=600, height=600)
  data.real.3d <- r$x[,1:3]
  pairs(data.real.3d, pch=16, col=as.numeric(data.real.group))
dev.off()

png(file="g487.png", width=600, height=600)
  draw.ellipse <- function (
      x, y=NULL, N=100,
      method=lines, ...
    ) {
    if (is.null(y)) {
      y <- x[,2]
      x <- x[,1]
    }
    centre <- c(mean(x), mean(y))
    m <- matrix(c(var(x),cov(x,y),
                  cov(x,y),var(y)),
                nr=2,nc=2)
    e <- eigen(m)
    r <- sqrt(e$values)
    v <- e$vectors
    theta <- seq(0,2*pi, length=N)
    x <- centre[1] + r[1]*v[1,1]*cos(theta) +
         r[2]*v[1,2]*sin(theta)
    y <- centre[2] + r[1]*v[2,1]*cos(theta) +
         r[2]*v[2,2]*sin(theta)
    method(x,y,...)
  }
  draw.star <- function (x, y=NULL, ...) {
    if (is.null(y)) {
      y <- x[,2]
      x <- x[,1]
    }
    d <- cbind(x,y)
    m <- apply(d, 2, mean)
    segments(m[1],m[2],x,y,...)
  }
  my.plot <- function (
      d, f=rep(1,dim(d)[1]),
      col=rainbow(length(levels(f))),
      variables=NULL, legend=T, legend.position=1,
      draw=draw.ellipse, ...) {
    xlim <- range(d[,1])
    ylim <- range(d[,2])
    if(!is.null(variables)){
      xlim <- range(xlim, variables[,1])
      ylim <- range(ylim, variables[,2])
    }
    plot(d, col=col[as.numeric(f)], pch=16,
         xlim=xlim, ylim=ylim, ...)
    for (i in 1:length(levels(f))) {
      try(
        draw(d[ as.numeric(f)==i, ], col=col[i])
      )
    }
    if(!is.null(variables)){
      arrows(0,0,variables[,1],variables[,2])
      text(1.05*variables,rownames(variables))
    }
    abline(h=0,lty=3)
    abline(v=0,lty=3)
    if(legend) {
      if(legend.position==1) {
        l=c( par('usr')[1],par('usr')[4], 0, 1 )
      } else if (legend.position==2) {
        l=c( par('usr')[2],par('usr')[4], 1, 1 )
      } else if (legend.position==3) {
        l=c( par('usr')[1],par('usr')[3], 0, 0 )
      } else if (legend.position==4) {
        l=c( par('usr')[2],par('usr')[3], 1, 0 )
      } else {
        l=c( mean(par('usr')[1:2]),
             mean(par('usr')[1:2]), .5, .5 )
      }
      legend(l[1], l[2], xjust=l[3], yjust=l[4],
             levels(f),
             col=col,
             lty=1,lwd=3)
    }
  }
  my.plot(data.real.3d[,1:2], data.real.group)
dev.off()

png(file="g488.png", width=600, height=600)
  op <- par(mfrow=c(3,3), mar=.1+c(0,0,0,0))
  plot.new();
  my.plot(data.real.3d[,c(2,1)], data.real.group, 
          xlab='',ylab='',axes=F,legend=F)
  box()
  my.plot(data.real.3d[,c(3,1)], data.real.group, 
          xlab='',ylab='',axes=F,legend=F)
  box()
  
  my.plot(data.real.3d[,c(1,2)], data.real.group, 
          draw=draw.star, 
          xlab='',ylab='',axes=F,legend=F)
  box()
  plot.new()
  my.plot(data.real.3d[,c(3,2)], data.real.group, 
          xlab='',ylab='',axes=F,legend=F)
  box()

  my.plot(data.real.3d[,c(1,2)], data.real.group, 
          draw=draw.star, 
          xlab='',ylab='',axes=F,legend=F)
  box()
  my.plot(data.real.3d[,c(2,3)], data.real.group, 
          draw=draw.star, 
          xlab='',ylab='',axes=F,legend=F)
  box()
  plot.new()
  par(op)
dev.off()

png(file="g489.png", width=600, height=600)
  set.seed(66327)
  random.rotation.matrix <- function (n=3) {
    m <- NULL
    for (i in 1:n) {
      x <- rnorm(n)
      x <- x / sqrt(sum(x*x))
      y <- rep(0,n)
      if (i>1) for (j in 1:(i-1)) {
        y <- y + sum( x * m[,j] ) * m[,j]
      }
      x <- x - y
      x <- x / sqrt(sum(x*x))
      m <- cbind(m, x)
    }
    m
  }
  op <- par(mfrow=c(3,3), mar=.1+c(0,0,0,0))
  for (i in 1:9) {
    #plot( (data.real.3d %*% random.rotation.matrix(3))[,1:2],
    #      pch=16, col=as.numeric(data.real.group),
    #      xlab='', ylab='', axes=F )
    my.plot((data.real.3d %*% random.rotation.matrix(3))[,1:2],
            data.real.group, 
            draw=draw.ellipse, 
            xlab='',ylab='',axes=F,legend=F)
    box()
  }
  par(op)
dev.off()

png(file="g490.png", width=600, height=600)
  longueur <- function (d,o) {
    n <- length(o)
    sum(diag( d [o[1:(n-1)],] [,o[2:n]] ))
  }
  tsp.descent <- function (d, N=1000) {
    # d: distance matrix
    d <- as.matrix(d)
    n <- dim(d)[1]
    o <- sample(1:n)
    v <- longueur(d,o)
    k <- 0
    res <- list()
    k.res <- c()
    l.res <- c(v)
    while (k<N) {
      i <- sample(1:n, 2)
      oo <- o
      oo[ i[1] ] <- o[ i[2] ]
      oo[ i[2] ] <- o[ i[1] ]
      w <- longueur(d,oo)
      if (w<v) {
        v <- w
        o <- oo
        res <- append(res, list(o))
        k.res <- append(k.res, k)
        l.res <- append(l.res, v)
        k <- 0
      } else {
        k <- k+1
      }
    }
    list(o=o, details=res, k=k.res, longueur=v, longueurs=l.res)
  }

  n <- 100
  x <- matrix(runif(2*n), nr=n)
  r <- tsp.descent(dist(x))
  o <- r$o
  plot(x)
  lines(x[o,], col='red')
dev.off()

png(file="g491.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  n <- length(r$details)  
  for (i in floor(c(1, n/4, n/2, n))) {
    plot(x, main=paste("n =",i))
    lines(x[r$details[[i]],], col='red')
  }
  par(op)
dev.off()

png(file="g492.png", width=600, height=600)
  plot(r$k)
dev.off()

png(file="g493.png", width=600, height=600)
  n <- 100
  x <- seq(0,6, length=n)
  x <- sample(x)
  x <- cbind(cos(x), sin(x))
  r <- tsp.descent(dist(x), N=10000)
  op <- par(mfrow=c(2,2))
  n <- length(r$details)  
  for (i in floor(c(1, n/4, n/2, n))) {
    plot(x, main=paste("n =",i))
    lines(x[r$details[[i]],], col='red')
  }
  par(op)
dev.off()

png(file="g494.png", width=600, height=600)
  plot(r$k)
dev.off()

png(file="g495.png", width=600, height=600)
  # Very long...
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    r <- tsp.descent(dist(x), N=1000)
    plot(x, main=signif(r$longueur))
    lines( x[r$o,], col='red', type='l' )
  }
  par(op)
dev.off()

png(file="g496.png", width=600, height=600)
  tsp.recuit <- function (d, N=1000, taux=.99) {
    # d: distance matrix
    d <- as.matrix(d)
    n <- dim(d)[1]
    o <- sample(1:n)
    v <- longueur(d,o)
    # Computing the initial temperature
    dE <- NULL
    j <- 0
    for (i in 1:100) {
      i <- sample(1:n, 2)
      oo <- o
      oo[ i[1] ] <- o[ i[2] ]
      oo[ i[2] ] <- o[ i[1] ]
      w <- longueur(d,oo)
      if (w<v) {
        dE <- append(dE,abs(w-v))
        j <- j+1        
      }
    }
    print(dE)
    T <- - max(dE) / log(.8)
    print(T)
    # Initialisations
    k <- 0
    res <- list()
    k.res <- c()
    l.res <- c(v)
    t.res <- c(T)
    p.res <- c()
    while (k<N) {
      i <- sample(1:n, 2)
      oo <- o
      oo[ i[1] ] <- o[ i[2] ]
      oo[ i[2] ] <- o[ i[1] ]
      w <- longueur(d,oo)
      p.res <- append(p.res, exp((v-w)/T))
      if ( runif(1) < exp((v-w)/T) ) {
        v <- w
        o <- oo
        res <- append(res, list(o))
        k.res <- append(k.res, k)
        l.res <- append(l.res, v)
        k <- 0
        T <- T*taux
        t.res <- append(t.res, T)
      } else {
        k <- k+1
      }
    }
    list(o=o, details=res, k=k.res, longueur=v, longueurs=l.res, T=t.res, p=p.res)
  }

  n <- 100
  x <- seq(0,5, length=n)
  x <- sample(x)
  x <- cbind(cos(x), sin(x))
  r <- tsp.recuit(dist(x), N=1000, taux=.995)
  op <- par(mfrow=c(2,2))
  n <- length(r$details)  
  for (i in floor(c(1, n/4, n/2, n))) {
    plot(x, main=paste("n =",i))
    lines(x[r$details[[i]],], col='red')
  }
  par(op)
dev.off()

png(file="g497.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  plot(r$k)
  plot(r$longueurs, main="Lengths")
  plot(r$T, main="Temperature")
  plot(r$p, ylim=0:1, main="Probabilities")
  par(op)
dev.off()

png(file="g498.png", width=600, height=600)
  print.tsp <- function (x, d=dist(x), f="") {
    ## BEWARE: The documentation states "All explicit [weight] data is
    ##         integral"
    ## If the best tour has length 0, the problem may come from there.
    d <- round(as.dist(d))
    n <- dim(x)[1]
    cat("TYPE : TSP\n", file=f) # Creating the file
    cat("DIMENSION : ", file=f, append=T)
    cat(n, file=f, append=T)
    cat("\n", file=f, append=T)
    cat("EDGE_WEIGHT_TYPE : EXPLICIT\n", file=f, append=T)
    cat("EDGE_WEIGHT_FORMAT : UPPER_ROW\n", file=f, append=T)
    cat("DISPLAY_DATA_TYPE : NO_DISPLAY\n", file=f, append=T)
    cat("EDGE_WEIGHT_SECTION\n", file=f, append=T)
    cat(paste(as.vector(d), collapse=" "), file=f, append=T)    
    cat("\n", file=f, append=T)
    cat("EOF\n", file=f, append=T)
  }
  tsp.plot.linkern <- function (x, d=dist(x), ...) {
    print.tsp(x, d, f="tmp.tsp")
    system(paste("./linkern", "-o", "tmp.tsp_result", "tmp.tsp"))
    ij <- 1+read.table("tmp.tsp_result", skip=1)[,1:2]
    xy <- prcomp(x)$x[,1:2]
    plot(xy, ...)
    segments( xy[,1][ ij[,1] ], 
              xy[,2][ ij[,1] ], 
              xy[,1][ ij[,2] ], 
              xy[,2][ ij[,2] ],
              col='red' )
    ij
  }

  op <- par(mfrow=c(3,3), mar=.1+c(0,0,0,0))
  for (i in 1:9) {  
    n <- 50
    k <- 3
    x <- matrix(rnorm(n*k), nr=n, nc=k)
    tsp.plot.linkern(x, d=round(100*dist(x)), axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g499.png", width=600, height=600)
  tsp.plot.concorde <- function (x, d=dist(x), plot=T, ...) {
    print.tsp(x, d, f="tmp.tsp")
    system(paste("./concorde", "tmp.tsp"))
    i <- scan("tmp.sol")[-1]
    ij <- 1+matrix(c(i, i[-1], i[1]), nc=2)
    xy <- prcomp(x)$x[,1:2]
    if (plot) {
      plot(xy, ...)
      n <- dim(ij)[1]
      segments( xy[,1][ ij[,1] ], 
                xy[,2][ ij[,1] ], 
                xy[,1][ ij[,2] ], 
                xy[,2][ ij[,2] ],
                col=rainbow(n) )
    }
    ij
  }
  op <- par(mfrow=c(3,3), mar=.1+c(0,0,0,0))
  for (i in 1:9) {  
    n <- 50
    k <- 3
    x <- matrix(rnorm(n*k), nr=n, nc=k)
    tsp.plot.concorde(x, d=round(100*dist(x)), axes=F)
    box()
  }
  par(op)
dev.off()

png(file="g500.png", width=600, height=600)
  tsp.plot.concorde(t(data.broken.line))
dev.off()

png(file="g501.png", width=600, height=600)
  tsp.plot.concorde(t(data.noisy.broken.line))
dev.off()

png(file="g502.png", width=600, height=600)
  x <- t(data.noisy.broken.line)
  ij <- tsp.plot.concorde(x, plot=F)
  xy <- prcomp(x)$x[,1:2]
  y <- xy[ij[,1],]
  x <- 1:dim(ij)[1]
  plot(xy)
  x1 <- predict(loess(y[,1]~x, span=.2))
  x2 <- predict(loess(y[,2]~x, span=.2))
  n <- length(x1)
  segments(x1[-n], x2[-n], x1[-1], x2[-1], col=rainbow(n-1), lwd=3)
dev.off()

png(file="g503.png", width=600, height=600)
  print.open.tsp <- function (x, d=dist(x), f="") {
    ## BEWARE: The documentation states "All explicit [weight] data is
    ##         integral"
    ## If the best tour has length 0, the problem may come from there.
    d <- as.matrix(d)
    d <- cbind(0,rbind(0,d))
    d <- round(as.dist(d))
    n <- dim(x)[1] + 1
    cat("TYPE : TSP\n", file=f) # Creating the file
    cat("DIMENSION : ", file=f, append=T)
    cat(n, file=f, append=T)
    cat("\n", file=f, append=T)
    cat("EDGE_WEIGHT_TYPE : EXPLICIT\n", file=f, append=T)
    cat("EDGE_WEIGHT_FORMAT : UPPER_ROW\n", file=f, append=T)
    cat("DISPLAY_DATA_TYPE : NO_DISPLAY\n", file=f, append=T)
    cat("EDGE_WEIGHT_SECTION\n", file=f, append=T)
    cat(paste(as.vector(d), collapse=" "), file=f, append=T)    
    cat("\n", file=f, append=T)
    cat("EOF\n", file=f, append=T)
  }
  tsp.plot.concorde.open <- function (x, d=dist(x), plot=T, smooth=F, span=.2,...) {
    print.open.tsp(x, d, f="tmp.tsp")
    system(paste("./concorde", "tmp.tsp"))
    i <- scan("tmp.sol")[-1]  # Remove the number of
    i <- i[-1]                # Remove the first node: 0
    ij <- matrix(c(i[-length(i)], i[-1] ), nc=2)
    xy <- prcomp(x)$x[,1:2]
    if (plot) {
      if (smooth) {
        y <- xy[ij[,1],]
        x <- 1:dim(ij)[1]
        plot(xy, ...)
        x1 <- predict(loess(y[,1]~x, span=span))
        x2 <- predict(loess(y[,2]~x, span=span))
        n <- length(x1)
        segments(x1[-n], x2[-n], x1[-1], x2[-1], col=rainbow(n-1), lwd=3)
      } else {
        plot(xy, ...)
        n <- dim(ij)[1]
        segments( xy[,1][ ij[,1] ], 
                  xy[,2][ ij[,1] ], 
                  xy[,1][ ij[,2] ], 
                  xy[,2][ ij[,2] ],
                  col=rainbow(n) )
      }
    }
    ij
  }
  x <- t(data.noisy.broken.line)
  tsp.plot.concorde.open(x)
dev.off()

png(file="g504.png", width=600, height=600)
  tsp.plot.concorde.open(x, smooth=T)
dev.off()

png(file="g505.png", width=600, height=600)
  tsp.plot.concorde.open(t(data.curve), smooth=T)
dev.off()

png(file="g506.png", width=600, height=600)
  tsp.plot.concorde.open(t(data.noisy.curve), smooth=T)
dev.off()

png(file="g507.png", width=600, height=600)
  tsp.plot.concorde.open(t(data.real), smooth=F)
dev.off()

png(file="g508.png", width=600, height=600)
  tsp.plot.concorde.open(t(data.real), smooth=T)
dev.off()

png(file="g509.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (s in c(.2, .3, .4, .5)) {
    tsp.plot.concorde.open(t(data.real), smooth=T, span=s)
  }
  par(op)
dev.off()

png(file="g510.png", width=600, height=600)
  tsp.plot.concorde.open(data.real.3d, smooth=F)
dev.off()

png(file="g511.png", width=600, height=600)
  tsp.plot.concorde.open(data.real.3d, smooth=T)
dev.off()

png(file="g512.png", width=600, height=600)
  tsp.plot.concorde.open(data.real.3d, smooth=T, span=.5)
dev.off()
detach.everything()

png(file="g513.png", width=600, height=600)
  x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
             matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
  cl <- kmeans(x, 2, 20)
  plot(x, col = cl$cluster, pch=3, lwd=1)
  points(cl$centers, col = 1:2, pch = 7, lwd=3)
  segments( x[cl$cluster==1,][,1], x[cl$cluster==1,][,2], 
            cl$centers[1,1], cl$centers[1,2])
  segments( x[cl$cluster==2,][,1], x[cl$cluster==2,][,2], 
            cl$centers[2,1], cl$centers[2,2], 
            col=2)
dev.off()

png(file="g514.png", width=600, height=600)
  # This is not what we want...
  load(file="data_notes.Rda") # reads in the "notes" matrix, defined
                              # in a preceding chapter.
  cl <- kmeans(notes, 6, 20)
  plot(notes, col = cl$cluster, pch=3, lwd=3)
  points(cl$centers, col = 1:6, pch=7, lwd=3)
dev.off()

png(file="g515.png", width=600, height=600)
  n <- 6
  cl <- kmeans(notes, n, 20)
  p <- princomp(notes)
  u <- p$loadings
  x <- (t(u) %*% t(notes))[1:2,]
  x <- t(x)
  plot(x, col=cl$cluster, pch=3, lwd=3)
  c <- (t(u) %*% t(cl$center))[1:2,]
  c <- t(c)
  points(c, col = 1:n, pch=7, lwd=3)
  for (i in 1:n) {
    print(paste("Cluster", i))
    for (j in (1:length(notes[,1]))[cl$cluster==i]) {
      print(paste("Point", j))
      segments( x[j,1], x[j,2], c[i,1], c[i,2], col=i )
    }
  }
  text( x[,1], x[,2], attr(x, "dimnames")[[1]] )
dev.off()

png(file="g516.png", width=600, height=600)
  library(cluster)
  clusplot( notes, pam(notes, 6)$clustering )
dev.off()

png(file="g517.png", width=600, height=600)
  clusplot( notes, pam(notes, 2)$clustering )
dev.off()

png(file="g518.png", width=600, height=600)
  clusplot( daisy(notes), pam(notes,2)$clustering, diss=T )
dev.off()

png(file="g519.png", width=600, height=600)
  clusplot( daisy(notes), pam(notes,6)$clustering, diss=T )
dev.off()

png(file="g520.png", width=600, height=600)
  clusplot(fanny(notes,2))
dev.off()

png(file="g521.png", width=600, height=600)
  clusplot(fanny(notes,10))
dev.off()

png(file="g522.png", width=600, height=600)
  data(USArrests)
  hc <- hclust(dist(USArrests), "ave")
  plot(hc)
dev.off()

png(file="g523.png", width=600, height=600)
  plot(hc, hang = -1)
dev.off()

png(file="g524.png", width=600, height=600)
  plot(hclust(dist(notes)))
dev.off()

png(file="g525.png", width=600, height=600)
  p <- princomp(notes)
  u <- p$loadings
  x <- (t(u) %*% t(notes))[1:2,]
  x <- t(x)
  plot(x, col="red")
  c <- hclust(dist(notes))$merge
  y <- NULL
  n <- NULL
  for (i in (1:(length(notes[,1])-1))) {
    print(paste("Step", i))
    if( c[i,1]>0 ){
      a <- y[ c[i,1], ]
      a.n <- n[ c[i,1] ]
    } else {
      a <- x[ -c[i,1], ]
      a.n <- 1
    }
    if( c[i,2]>0 ){
      b <- y[ c[i,2], ]
      b.n <- n[ c[i,2] ]
    } else {
      b <- x[ -c[i,2], ]
      b.n <- 1
    }
    n <- append(n, a.n+b.n)
    m <- ( a.n * a + b.n * b )/( a.n + b.n )
    y <- rbind(y, m)
    segments( m[1], m[2], a[1], a[2], col="red" )
    segments( m[1], m[2], b[1], b[2], col="red" )
    if( i> length(notes[,1])-1-6 ){
      op=par(ps=30)
      text( m[1], m[2], paste(length(notes[,1])-i), col="red" )
      par(op)
    }
  }
  text( x[,1], x[,2], attr(x, "dimnames")[[1]] )
dev.off()

png(file="g526.png", width=600, height=600)
  do.it <- function (x, k=3, main="") {
    if (inherits(x, "dist")) {
      d <- x
      x <- sammon(d)$points
    } else {
      d <- dist(x)      
      if (!is.vector(x)) {
        x <- sammon(d)$points
      }
    }
    op <- par(mfrow = c(2,2), 
              mar = c(3,2,4,2)+.1,
              oma = c(0,0,0,0))
    if (is.vector(x)) {
      hist(x, 
           col="light blue",
           xlab="", ylab="", main="Data")
    } else {
      plot(x, main="Data")
    }
    hist(as.vector(d), 
         col="light blue",
         xlab="", ylab="", main="Distances")
    plot(hclust(d), labels=FALSE, 
         main = "Hierarchical clustering")
    plot(silhouette(cutree(hclust(d),k), d),
         main = "Silhouette")
    par(op)
    mtext(main, line=3, font=2, cex=1.2)
  }
  do.it(runif(100), k=2,
        main="Random gaussian data")
dev.off()

png(file="g527.png", width=600, height=600)
  set.seed(1)
  x <- rt(100,2)
  do.it(x, main="Random data with fat tails")
dev.off()

png(file="g528.png", width=600, height=600)
  x <- matrix(rnorm(10000), nr=100, nc=100)
  do.it(x, main="High-dimensional random gaussian data")
dev.off()

png(file="g529.png", width=600, height=600)
  set.seed(1)
  x <- matrix(rt(10000, df=4), nr=100, nc=100)
  do.it(x, k=5, main="High-dimensional, fat tails")
dev.off()

png(file="g530.png", width=600, height=600)
  N <- 100
  x <- rnorm(N, sample(c(-2,2), N, replace=T))
  do.it(x, k=2, main="Two real clusters")
dev.off()

png(file="g531.png", width=600, height=600)
  N <- 100  # Number of observations
  n <- 10   # Dimension
  k <- 3    # Number of clusters
  i <- sample( 1:k, N, replace=T )
  mu <- matrix(runif(n*k, -2, 2), nr=k, nc=n)
  x <- matrix(rnorm(N*n), nc=n) + mu[i,]   # One subject per row
  do.it(x, k=3, main="Three cluster, higher dimension")
dev.off()

png(file="g532.png", width=600, height=600)
  N <- 1000
  d <- matrix(rnorm(2*N), nc=2)
  par(mar=c(2,2,4,2))
  plot(d, xlim=c(-2,2), ylim=c(-2,2),
       axes = FALSE,
       xlab="", ylab="",
       main="Contour plot of the Euclidian distance to the origin")
  box()
  abline(h=0, v=0, lty=3)
  x <- seq(min(d[,1]), max(d[,1]), length=100)
  y <- seq(min(d[,2]), max(d[,2]), length=100)
  z <- outer(x, y, function (x,y) sqrt(x^2 + y^2))
  contour(x, y, z, 
          add = TRUE,
          col = "blue", lwd = 3)
dev.off()

png(file="g533.png", width=600, height=600)
  N <- 5000
  library(MASS)
  d <- mvrnorm(N, mu=c(0,0), Sigma=matrix(c(1,.8,.8,1),2))
  par(mar=c(2,2,4,2))
  plot(d, xlim=c(-2,2), ylim=c(-2,2),
       axes = FALSE,
       xlab="", ylab="",
       main="Contour plot of the Euclidian distance to the origin")
  box()
  abline(h=0, v=0, lty=3)
  x <- seq(min(d[,1]), max(d[,1]), length=100)
  y <- seq(min(d[,2]), max(d[,2]), length=100)
  z <- outer(x, y, function (x,y) sqrt(x^2 + y^2))
  contour(x, y, z, 
          add = TRUE,
          col = "blue", lwd = 1)
  points(sqrt(2)*c(1,-1), sqrt(2)*c(1,1), 
         lwd = 10, cex = 5, pch = 4, col = "blue")
  text(sqrt(2)*c(1,-1), sqrt(2)*c(1,1) + .1, 
       c("2", "1"),
       pos = 3, adj = .5,
       font = 2, cex = 3, col = "blue")
dev.off()

png(file="g534.png", width=600, height=600)
  N <- 1000
  V <- matrix(c(1,.8,.8,1),2)
  d <- mvrnorm(N, mu=c(0,0), Sigma=V)
  par(mar=c(2,2,4,2))
  plot(d, xlim=c(-2,2), ylim=c(-2,2),
       axes = FALSE,
       xlab="", ylab="",
       main="Contour plot of the Mahalanobis distance to the origin")
  box()
  abline(h=0, v=0, lty=3)
  x <- seq(min(d[,1]), max(d[,1]), length=100)
  y <- seq(min(d[,2]), max(d[,2]), length=100)
  z <- outer(x, y, function (x, y) {
    sqrt(apply(rbind(x,y) * solve(V, rbind(x,y)), 2, sum))
  } )   # BUG: MEMORY PROBLEM...
  contour(x, y, z, 
          add = TRUE,
          col = "blue", lwd = 3)
dev.off()

png(file="g535.png", width=600, height=600)
  N <- 1000
  V <- matrix(c(1,.8,.8,1),2)
  d <- mvrnorm(N, mu=c(0,0), Sigma=V)
  par(mar=c(2,2,4,2))
  plot(d, xlim=c(-2,2), ylim=c(-2,2),
       axes = FALSE,
       xlab="", ylab="",
       main="Contour plot of the Tracking Error distance to the origin")
  box()
  abline(h=0, v=0, lty=3)
  x <- seq(min(d[,1]), max(d[,1]), length=100)
  y <- seq(min(d[,2]), max(d[,2]), length=100)
  z <- outer(x, y, function (x, y) {
    sqrt(apply(rbind(x,y) * (V %*% rbind(x,y)), 2, sum))
  } )
  contour(x, y, z, 
          add = TRUE,
          col = "blue", lwd = 3)
dev.off()

png(file="g536.png", width=600, height=600)
  ## Do the same with centroid clustering and squared Euclidean distance,
  ## cut the tree into ten clusters and reconstruct the upper part of the
  ## tree from the cluster centers.
  hc <- hclust(dist(USArrests)^2, "cen")
  memb <- cutree(hc, k = 10)
  cent <- NULL
  for(k in 1:10){
    cent <- rbind(cent, colMeans(USArrests[memb == k, , drop = FALSE]))
  }
  hc1 <- hclust(dist(cent)^2, method = "cen", members = table(memb))
  opar <- par(mfrow = c(1, 2))
  plot(hc,  labels = FALSE, hang = -1, main = "Original Tree")
dev.off()

png(file="g537.png", width=600, height=600)
  plot(hc1, labels = FALSE, hang = -1, main = "Re-start from 10 clusters")
dev.off()

png(file="g538.png", width=600, height=600)
  n <- 1000
  x <- runif(n,-1,1)
  y <- ifelse(runif(n)>.5,-.1,.1) + .02*rnorm(n)
  d <- data.frame(x=x, y=y)
  plot(d,type='p', xlim=c(-1,1), ylim=c(-1,1))
dev.off()

png(file="g539.png", width=600, height=600)
  test.kmeans <- function (d, ...) {
    cl <- kmeans(d,2)
    plot(d, col=cl$cluster, main="kmeans", ...)
    points(cl$centers, col=1:2, pch=7, lwd=3)
  }
  test.kmeans(d, xlim=c(-1,1), ylim=c(-1,1))
dev.off()

png(file="g540.png", width=600, height=600)
  test.hclust <- function (d, ...) {
    hc <- hclust(dist(d))
    remplir <- function (m, i, res=NULL) {
      if(i<0) {
        return( c(res, -i) )
      } else {
        return( c(res, remplir(m, m[i,1], NULL), remplir(m, m[i,2], NULL) ) )
      }
    }
    a <- remplir(hc$merge, hc$merge[n-1,1])
    b <- remplir(hc$merge, hc$merge[n-1,2])
    co <- rep(1,n)
    co[b] <- 2
    plot(d, col=co, main="hclust", ...)
  }
  test.hclust(d, xlim=c(-1,1), ylim=c(-1,1))
dev.off()

png(file="g541.png", width=600, height=600)
  get.sample <- function (n=1000, p=.7) {
    x1 <- rnorm(n)
    y1 <- rnorm(n)
    r2 <- 7+rnorm(n)
    t2 <- runif(n,0,2*pi)
    x2 <- r2*cos(t2)
    y2 <- r2*sin(t2)
    r <- runif(n)>p
    x <- ifelse(r,x1,x2)
    y <- ifelse(r,y1,y2)
    d <- data.frame(x=x, y=y)
    d
  }
  d <- get.sample()
  plot(d,type='p', xlim=c(-10,10), ylim=c(-10,10))
dev.off()

png(file="g542.png", width=600, height=600)
  test.kmeans(d, xlim=c(-10,10), ylim=c(-10,10))
dev.off()

png(file="g543.png", width=600, height=600)
  test.hclust(d, xlim=c(-10,10), ylim=c(-10,10))
dev.off()

png(file="g544.png", width=600, height=600)
  y <- abs(y)
  d <- data.frame(x=x, y=y)
  plot(d,type='p', xlim=c(-10,10), ylim=c(0,10))
dev.off()

png(file="g545.png", width=600, height=600)
  test.kmeans(d)
dev.off()

png(file="g546.png", width=600, height=600)
  test.hclust(d)
dev.off()

png(file="g547.png", width=600, height=600)
  d <- get.sample()
  x <- d$x; y <- d$y
  d <- data.frame(x=x, y=y, xx=x*x, yy=y*y, xy=x*y)
  test.kmeans(d)
dev.off()

png(file="g548.png", width=600, height=600)
  test.hclust(d)
dev.off()

png(file="g549.png", width=600, height=600)
  d <- data.frame(x=x, y=y, xx=x*x, yy=y*y, xy=x*y, xpy=x*x+y*y)
  test.kmeans(d)
dev.off()

png(file="g550.png", width=600, height=600)
  test.hclust(d)
dev.off()

png(file="g551.png", width=600, height=600)
  library(KernSmooth)
  r <- bkde2D(d, bandwidth=c(.5,.5))
  persp(r$fhat)
dev.off()

png(file="g552.png", width=600, height=600)
  n <- length(r$x1)
  plot( matrix(r$x1, nr=n, nc=n, byrow=F), 
        matrix(r$x2, nr=n, nc=n, byrow=T), 
        col=r$fhat>.001 )
dev.off()

png(file="g553.png", width=600, height=600)
  density.classification.plot <- function (x,y,d.lim=.5,n.lim=5) {
    n <- length(x)
    # Distance computation
    a <- matrix(x, nr=n, nc=n, byrow=F) - matrix(x, nr=n, nc=n, byrow=T) 
    b <- matrix(y, nr=n, nc=n, byrow=F) - matrix(y, nr=n, nc=n, byrow=T) 
    a <- a*a + b*b
    # Which observations are dense (i.e., in a cluster)?
    b <- apply(a<d.lim, 1, sum)>=n.lim
    plot(x, y, col=b)
    points(x,y,pch='.')
  }
  density.classification.plot(d$x,d$y)
dev.off()

png(file="g554.png", width=600, height=600)
  # Beware: the following code is very recursive -- but, by default, 
  # R limits the size of the function calls stack to 500 elements.
  # We first increment this value.
  options(expressions=10000)

  density.classification.plot <- function (x,y,d.lim=.5,n.lim=5) {
    n <- length(x)
    # Distance computations
    a <- matrix(x, nr=n, nc=n, byrow=F) - matrix(x, nr=n, nc=n, byrow=T) 
    b <- matrix(y, nr=n, nc=n, byrow=F) - matrix(y, nr=n, nc=n, byrow=T) 
    a <- a*a + b*b
    # Which observations are dense (in a cluster)?
    b <- apply(a<d.lim, 1, sum)>=n.lim
    # We sort the observations
    cl <- rep(0,n)
    m <- 1
    numerote <- function (i,co,cl) {
      print(paste(co, i))
      for (j in (1:n)[ a[i,]<d.lim & b & cl==0 ]) {
        #print(paste("  ",j))
        cl[j] <- co
        try( cl <- numerote(j,co,cl) )   # Too recursive...
      }
      cl
    }
    for (i in 1:n) {
      if (b[i]) { # Are we in a cluster?
        # Cluster number
        if (cl[i] == 0) {
          co <- m
          cl[i] <- co
          m <- m+1
        } else {
          co <- cl[i]          
        }
        # We number the nearby points
        #print(co)
        cl <- numerote(i,co,cl)
      }
    }
    plot(x, y, col=cl)
    points(x,y,pch='.')
  }
  density.classification.plot(d$x,d$y)
dev.off()

png(file="g555.png", width=600, height=600)
  density.classification.plot <- function (x,y,d.lim=.5,n.lim=5, ...) {
    n <- length(x)
    # Distance computation
    a <- matrix(x, nr=n, nc=n, byrow=F) - matrix(x, nr=n, nc=n, byrow=T) 
    b <- matrix(y, nr=n, nc=n, byrow=F) - matrix(y, nr=n, nc=n, byrow=T) 
    a <- a*a + b*b
    # Which observations are dense (in a cluster)?
    b <- apply(a<d.lim, 1, sum)>=n.lim
    # We sort the observations
    cl <- rep(0,n)
    m <- 0
    for (i in 1:n) {
      if (b[i]) { # Are we in a cluster?
        # Cluster number
        if (cl[i] == 0) {
          m <- m+1
          co <- m
          cl[i] <- co
          print(paste("Processing cluster", co))
          done <- F 
          while (!done) {
            done <- T
            for (j in (1:n)[cl==co]) {
              l <- (1:n)[ a[j,]<d.lim & b & cl==0 ]
              if( length(l)>0 ) {
                done <- F
                for (k in l) { 
                  cl[k] <- co 
                  #print(k)
                }
              }
            }
            
          }
        } else {
          # Already processed cluster: pass
        }
      }
    }
    plot(x, y, col=cl, ...)
    points(x,y,pch='.')
  }
  density.classification.plot(d$x,d$y)
dev.off()

png(file="g556.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (d.lim in c(.2,1,2)) {
    for (n.lim in c(3,5,10)) {
      density.classification.plot(d$x,d$y,d.lim,n.lim, 
        main=paste("d.lim = ",d.lim, ",  n.lim = ",n.lim, sep=''))
    }
  }
  par(op)
dev.off()
detach.everything()

png(file="g557.png", width=600, height=200)
  n <- 100
  x <- sample(c(-1,1), n, replace=T)
  plot(x, type='h', main="Bernoulli variables")
dev.off()

png(file="g558.png", width=600, height=600)
  n <- 1000
  x <- sample(c(-1,1), n, replace=T)
  plot(cumsum(x), type='l',
       main="Cumulated sums of Bernoulli variables")
dev.off()

png(file="g559.png", width=600, height=200)
  n <- 100
  x <- sample(c(-1,1), n, replace=T, prob=c(.2,.8))
  plot(x, type='h',
       main="Bernoulli variables, different probabilities")
dev.off()

png(file="g560.png", width=600, height=600)
  n <- 200
  x <- sample(c(-1,1), n, replace=T, prob=c(.2,.8))
  plot(cumsum(x), type='l',
       main="Cummulative sums of Bernoulli variables")
dev.off()

png(file="g561.png", width=600, height=200)
  n <- 200
  x <- runif(n)
  x <- x>.3
  plot(x, type='h', main="Bernoulli variables")
dev.off()

png(file="g562.png", width=600, height=600)
  N <- 10000
  n <- 20
  p <- .5
  x <- rep(0,N)
  for (i in 1:N) {
    x[i] <- sum(runif(n)<p)
  }
  hist(x, 
       col='light blue',
       main="Simulating a binomial law")
dev.off()

png(file="g563.png", width=600, height=600)
  N <- 1000
  n <- 10
  p <- .5
  x <- rbinom(N,n,p)
  hist(x, 
       xlim = c(min(x), max(x)), 
       probability = TRUE, 
       nclass = max(x) - min(x) + 1, 
       col = 'lightblue',
       main = 'Binomial distribution, n=10, p=.5')
  lines(density(x, bw=1), col = 'red', lwd = 3)
dev.off()

png(file="g564.png", width=600, height=600)
  N <- 100000
  n <- 100
  p <- .5
  x <- rbinom(N,n,p)
  hist(x, 
       xlim = c(min(x), max(x)), 
       probability = TRUE, 
       nclass = max(x) - min(x) + 1, 
       col = 'lightblue',
       main = 'Binomial distribution, n=100, p=.5')
  lines(density(x,bw=1), col = 'red', lwd = 3)
dev.off()

png(file="g565.png", width=600, height=600)
  p <- .9
  x <- rbinom(N,n,p)
  hist(x, 
       xlim = c(min(x), max(x)), 
       probability = TRUE, 
       nclass = max(x) - min(x) + 1, 
       col = 'lightblue',
       main = 'Binomial distribution, n=100, p=.9')
  lines(density(x,bw=1), col = 'red', lwd = 3)
dev.off()

png(file="g566.png", width=600, height=600)
  N <- 10000
  n <- 100
  p <- .5
  x <- NULL
  for (i in 1:N) {
    x <- append(x, sum(sample( c(1,0), 
                               n, 
                               replace = TRUE, 
                               prob = c(p, 1-p) )))
  }
  hist(x, 
       xlim = c(min(x), max(x)), 
       probability = TRUE, 
       nclass = max(x) - min(x) + 1, 
       col = 'lightblue',
       main = 'Binomial distribution, n=100, p=.5')
  lines(density(x,bw=1), col = 'red', lwd = 3)
dev.off()

png(file="g567.png", width=600, height=600)
  N <- 10000
  n <- 5
  urn <- c(rep(1,15),rep(0,5))
  x <- NULL
  for (i in 1:N) {
    x <- append(x, sum(sample( urn, n, replace=F )))
  }
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Hypergeometric distribution, n=20, p=.75; k=5')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g568.png", width=600, height=600)
  N <- 10000
  n <- 5
  x <- rhyper(N, 15, 5, 5)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Hypergeometric distribution, n=20, p=.75, k=5')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g569.png", width=600, height=600)
  N <- 10000
  n <- 5
  x <- rhyper(N, 300, 100, 100)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Hypergeometric distribution, n=400, p=.75, k=100')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g570.png", width=600, height=600)
  N <- 10000
  x <- rpois(N, 1)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Poisson distribution, lambda=1')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g571.png", width=600, height=600)
  N <- 10000
  x <- rpois(N, 3)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Poisson distribution, lambda=3')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g572.png", width=600, height=600)
  N <- 10000
  x <- rpois(N, 5)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Poisson distribution, lambda=5')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g573.png", width=600, height=600)
  N <- 10000
  x <- rpois(N, 20)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Poisson distribution, lambda=20')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g574.png", width=600, height=600)
  my.rgeom <- function (N, p) {
    bernoulli <- sample( c(0,1), N, replace=T, prob=c(1-p, p) )
    diff(c(0, which(bernoulli == 1))) - 1
  }
  hist( my.rgeom(10000, .5), col="light blue",
        main="Geometric distribution" )
dev.off()

png(file="g575.png", width=600, height=600)
  N <- 10000
  x <- rgeom(N, .5)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Geometric distribution, p=.5')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g576.png", width=600, height=600)
  N <- 10000
  x <- rgeom(N, .1)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Geometric distribution, p=.1')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g577.png", width=600, height=600)
  N <- 10000
  x <- rgeom(N, .01)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=20,
       col='lightblue',
       main='Geometric distribution, p=.01')
  lines(density(x), col='red', lwd=3)
dev.off()

png(file="g578.png", width=600, height=600)
  N <- 100000
  x <- rnbinom(N, 10, .25)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='Negative binomial distribution, n=10, p=.25')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g579.png", width=600, height=600)
  N <- 10000
  x <- rnbinom(N, 10, .5)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='negative binomial distribution, n=10, p=.5')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g580.png", width=600, height=600)
  N <- 10000
  x <- rnbinom(N, 10, .75)
  hist(x, 
       xlim=c(min(x),max(x)), probability=T, nclass=max(x)-min(x)+1, 
       col='lightblue',
       main='negative binomial distribution, n=10, p=.75')
  lines(density(x,bw=1), col='red', lwd=3)
dev.off()

png(file="g581.png", width=600, height=600)
  curve(dexp(x), xlim=c(0,10), col='red', lwd=3,
        main='Exponential Probability Distribution Function')
dev.off()

png(file="g582.png", width=600, height=600)
  n <- 1000
  x <- rexp(n)
  hist(x, probability=T,
       col='light blue', main='Exponential Distribution')
  lines(density(x), col='red', lwd=3)
  curve(dexp(x), xlim=c(0,10), col='red', lwd=3, lty=2,
        add=T)
dev.off()

png(file="g583.png", width=600, height=600)
  limite.centrale <- function (r=runif, m=.5, s=1/sqrt(12), n=c(1,3,10,30), N=1000) {
    for (i in n) {
      x <- matrix(r(i*N),nc=i)
      x <- ( apply(x, 1, sum) - i*m )/(sqrt(i)*s)
      hist(x, col='light blue', probability=T, main=paste("n =",i), 
           ylim=c(0,max(.4, density(x)$y)))
      lines(density(x), col='red', lwd=3)
      curve(dnorm(x), col='blue', lwd=3, lty=3, add=T)
      if( N>100 ) {
        rug(sample(x,100))
      } else {
        rug(x)
      }
    }
  }
  op <- par(mfrow=c(2,2))
  limite.centrale()
  par(op)
dev.off()

png(file="g584.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  limite.centrale(rexp, m=1, s=1)
  par(op)
dev.off()

png(file="g585.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  limite.centrale(rexp, m=1, s=1)
  par(op)
dev.off()

png(file="g586.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  limite.centrale(function (n) { rnorm(n, sample(c(-3,3),n,replace=T)) }, 
                  m=0, s=sqrt(10), n=c(1,2,3,10))
  par(op)
dev.off()

png(file="g587.png", width=600, height=600)
  curve(dnorm(x), xlim=c(-3,3), col='red', lwd=3)
  title(main='Gaussian Probability Distribution Function')
dev.off()

png(file="g588.png", width=600, height=600)
  curve(pnorm(x), xlim=c(-3,3), col='red', lwd=3)
  title(main='Cumulative gaussian distribution function')
dev.off()

png(file="g589.png", width=600, height=600)
  curve(qnorm(x), xlim=c(0,1), col='red', lwd=3)
  title(main='Gaussian quantiles function')
dev.off()

png(file="g590.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  hist(x, probability=T, col='light blue', main='Gaussian Distribution')
  lines(density(x), col='red', lwd=3)
  curve(dnorm(x), add=T, col='red', lty=2, lwd=3)
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('sample density', 'theoretical density'),
         lwd=2, lty=c(1,2),
         col='red')
dev.off()

png(file="g591.png", width=600, height=600)
  curve(dchisq(x,1), xlim=c(0,5), col='red', lwd=3)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  title(main="Chi2, one degree of freedom")
dev.off()

png(file="g592.png", width=600, height=600)
  curve(dchisq(x,1), xlim=c(0,10), ylim=c(0,.6), col='red', lwd=3)
  curve(dchisq(x,2), add=T, col='green', lwd=3)
  curve(dchisq(x,3), add=T, col='blue', lwd=3)
  curve(dchisq(x,5), add=T, col='orange', lwd=3)
  abline(h=0,lty=3)
  abline(v=0,lty=3)
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('df=1', 'df=2', 'df=3', 'df=5'),
         lwd=3,
         lty=1,
         col=c('red', 'green', 'blue', 'orange')
        )
  title(main='Chi^2 Distributions')
dev.off()

png(file="g593.png", width=600, height=600)
  curve( dt(x,1), xlim=c(-3,3), ylim=c(0,.4), col='red', lwd=2 )
  curve( dt(x,2), add=T, col='blue', lwd=2 )
  curve( dt(x,5), add=T, col='green', lwd=2 )
  curve( dt(x,10), add=T, col='orange', lwd=2 )
  curve( dnorm(x), add=T, lwd=3, lty=3 )
  title(main="Student T distributions")
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('df=1', 'df=2', 'df=5', 'df=10', 'Gaussian distribution'),
         lwd=c(2,2,2,2,2), 
         lty=c(1,1,1,1,3),
         col=c('red', 'blue', 'green', 'orange', par("fg")))
dev.off()

png(file="g594.png", width=600, height=600)
  curve(df(x,1,1), xlim=c(0,2), ylim=c(0,.8), lty=2)
  curve(df(x,3,1), add=T)
  curve(df(x,6,1), add=T, lwd=3)
  curve(df(x,3,3), add=T, col='red')
  curve(df(x,6,3), add=T, lwd=3, col='red')
  curve(df(x,3,6), add=T, col='blue')
  curve(df(x,6,6), add=T, lwd=3, col='blue')
  title(main="Fisher's F")
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('df=(1,1)', 'df=(3,1)', 'df=(6,1)', 
           'df=(3,3)', 'df=(6,3)', 
           'df=(3,6)', 'df=(6,6)'),
         lwd=c(1,1,3,1,3,1,3),
         lty=c(2,1,1,1,1,1,1),
         col=c(par("fg"), par("fg"), par("fg"), 'red', 'red', 'blue', 'blue'))
dev.off()

png(file="g595.png", width=600, height=600)
  curve(dlnorm(x), xlim=c(-.2,5), lwd=3,
        main="Log-normal distribution")
dev.off()

png(file="g596.png", width=600, height=600)
  N <- 100                         # Number of arrows
  alpha <- runif(N, -pi/2, pi/2)   # Direction of the arrow
  x <- tan(alpha)                  # Arrow impact
  plot.new()
  plot.window(xlim=c(-5, 5), ylim=c(-1.1, 2))
  segments( 0, -1,     # Position of the Bowman
            x, 0   )   # Impact
  d <- density(x)
  lines(d$x, 5*d$y, col="red", lwd=3 )
  box()
  abline(h=0)
  title(main="The bowman's distribution (Cauchy)")
  # Exercise: turn this into an animation...
dev.off()

png(file="g597.png", width=600, height=600)
  N <- 10000
  x <- tan(runif(N, -pi/2, pi/2))
  xlim <- qcauchy(2/N)
  xlim <- c(xlim, -xlim)
  plot(qcauchy(ppoints(N)),  sort(x),
       xlim=xlim, ylim=xlim,
       main="The bowman's distribution and Cauchy's")
dev.off()

png(file="g598.png", width=600, height=600)
  curve(dcauchy(x),xlim=c(-5,5), ylim=c(0,.5), lwd=3)
  curve(dnorm(x), add=T, col='red', lty=2)
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('Cauchy distribution', 'Gaussian distribution'),
         lwd=c(3,1),
         lty=c(1,2),
         col=c(par("fg"), 'red'))
dev.off()

png(file="g599.png", width=600, height=600)
  N <- 100
  x <- cumsum(rnorm(N))
  y <- cumsum(rnorm(N))
  plot(x, y, 
       type = "o", pch = 16, lwd = 2,
       xlab = "", ylab = "", 
       axes = FALSE,
       main = "Brownian Motion")
  box()
dev.off()

png(file="g600.png", width=600, height=600)
  set.seed(1)
  x <- cumsum(rt(N, df=2))
  y <- cumsum(rt(N, df=2))
  plot(x, y, 
       type = "o", pch = 16, lwd = 2,
       xlab = "", ylab = "", 
       axes = FALSE,
       main = "Levy flight")
  box()
dev.off()

png(file="g601.png", width=600, height=600)
  N <- 10000
  m <- c(-2,0,2)    # Means
  p <- c(.3,.4,.3)  # Probabilities
  s <- c(1, 1, 1)   # Standard deviations
  x <- cbind( rnorm(N, m[1], s[1]), 
              rnorm(N, m[2], s[2]), 
              rnorm(N, m[3], s[3]) )
  a <- sample(1:3, N, prob=p, replace=TRUE)
  y <- x[ 1:N + N*(a-1) ]

  qqnorm(y, 
         main="Gaussian QQ-plot of a mixture of gaussians")
  qqline(y, col="red", lwd=3)
dev.off()

png(file="g602.png", width=600, height=600)
  hist(y, col="light blue", probability=TRUE,
       ylim=c(0,.25),
       main="Mixture of gaussians")
  curve(dnorm(x, mean=mean(y), sd=sd(y)),
        add=TRUE, col="red", lwd=3, lty=2)
  lines(density(x), col="red", lwd=3)
dev.off()

png(file="g603.png", width=600, height=600)
  curve( p[2] * dnorm(x, mean=m[2], sd=s[2]),
         col = "green", lwd = 3, 
         xlim = c(-5,5),
         main = "The three gaussian distributions in our mixture",
         xlab = "", ylab = "")
  curve( p[1] * dnorm(x, mean=m[1], sd=s[1]),
         col="red", lwd=3, add=TRUE)
  curve( p[3] * dnorm(x, mean=m[3], sd=s[3]),
         col="blue", lwd=3, add=TRUE)
dev.off()

png(file="g604.png", width=600, height=600)
  n <- 200
  x <- seq(-5, 5, length=n)
  y3 <- p[1] * dnorm(x, mean=m[1], sd=s[1]) +
        p[2] * dnorm(x, mean=m[2], sd=s[2]) +
        p[3] * dnorm(x, mean=m[3], sd=s[3])
  y2 <- p[1] * dnorm(x, mean=m[1], sd=s[1]) +
        p[2] * dnorm(x, mean=m[2], sd=s[2])
  y1 <- p[1] * dnorm(x, mean=m[1], sd=s[1])
  plot.new()
  plot.window(xlim=range(x), ylim=range(0,y1,y2,y3))
  polygon(c(x[1],x,x[n]), c(0,y3,0), col="blue",  border=NA)
  polygon(c(x[1],x,x[n]), c(0,y2,0), col="green", border=NA)
  polygon(c(x[1],x,x[n]), c(0,y1,0), col="red",   border=NA)
  lines(x, y1, lwd=3)
  lines(x, y2, lwd=3)
  lines(x, y3, lwd=3)
  box()
  axis(1)
  title("Mixture of gaussians")
dev.off()

png(file="g605.png", width=600, height=600)
  q <- function (p, a=1, b=0, c=0, d=0) {
    a * qnorm(p) + b + c * p + d * p^2
  }
  N <- 10000
  x <- runif(N)
  op <- par(mfrow=c(2,3))
  for (a in c(-2,2)) {
    y <- q(x, 1, -a/2, a)
    hist(y, 
         xlab = "",
         main = "Gaussian-polynomial quantile mixture",
         col = "light blue", 
         probability = TRUE)
    curve(q(x, 1, -a/2, a), 
          lwd = 3, 
          xlim = c(0,1), ylim = c(-3,3),
          xlab = "",
          main = "Quantile function")
    abline(h = 0, 
           v = c(0, .5, 1), 
           lty = 3)
    curve(qnorm(x), lty = 2, add = TRUE)
    qqnorm(y)
    qqline(y, col="red")
  }
  par(op)
dev.off()

png(file="g606.png", width=600, height=600)
  library(fBasics)
  # alpha=2 is the gaussian distribution
  # alpha<2 distributions have fat tails
  x <- rsymstb(10000, alpha=1.5)
  y <- x[ abs(x) < 10 ]
  hist(y, 
       ylim = c(0, .4),
       col = "light blue",
       probability = TRUE,
       xlab = "", ylab = "",
       main = "A stable distribution (alpha=1.5)")
  lines(density(y), col = "red", lwd = 3)
  curve(dnorm(x),   col = "red", lwd = 3, lty=2, 
                    add = TRUE)  
dev.off()

png(file="g607.png", width=600, height=600)
  set.seed(1)
  N <- 1e7
  x <- sample(c(-1,+1), N, replace = TRUE)
  x <- cumsum(x)          # Random walk
  x <- diff(which(x==0))  # Time to go back to zero
  r <- density(x[x<100])
  plot(r$x, log(r$y), 
       xlim = c(0,20),
       ylim = c(-6,0), 
       type = "l",
       xlab = "Hitting time",
       ylab = "log(density)",
       main = "(Discretized) Levy distribution")
  v <- r$x[ which.max(r$y) ]  # Mode
  abline(v = v, lty=3, lwd=3)
  curve( dnorm(x, mean = v, sd = 1, log = TRUE) -
         dnorm(0,sd=1,log=T) + log(max(r$y)) , 
         add = TRUE, col = "blue", lwd = 3, lty = 2 )
dev.off()

png(file="g608.png", width=600, height=600)
  curve(dexp(x), xlim=c(0,3), ylim=c(0,2))
  curve(dweibull(x,1), lty=3, lwd=3, add=T)
  curve(dweibull(x,2), col='red', add=T)
  curve(dweibull(x,.8), col='blue', add=T)
  title(main="Weibull Probability Distribution Function")
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('Exponential', 'Weibull, shape=1',
           'Weibull, shape=2', 'Weibull, shape=.8'),
         lwd=c(1,3,1,1),
         lty=c(1,3,1,1),
         col=c(par("fg"), par("fg"), 'red', 'blue'))
dev.off()

png(file="g609.png", width=600, height=600)
  curve( dgamma(x,1,1), xlim=c(0,5) )
  curve( dgamma(x,2,1), add=T, col='red' )
  curve( dgamma(x,3,1), add=T, col='green' )
  curve( dgamma(x,4,1), add=T, col='blue' )
  curve( dgamma(x,5,1), add=T, col='orange' )
  title(main="Gamma probability distribution function")
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('k=1 (Exponential distribution)', 'k=2', 'k=3', 'k=4', 'k=5'),
         lwd=1, lty=1,
         col=c(par('fg'), 'red', 'green', 'blue', 'orange') )
dev.off()

png(file="g610.png", width=600, height=600)
  n <- 500
  x1 <- rexp(n,17)
  x2 <- rexp(n,17)
  x3 <- rexp(n,17)
  x <- x1 + x2 + x3
  # Simpler, but less readable:
  # k <- 3
  # x <- drop(apply( matrix( rexp(n*k,17), nr=n, nc=k ), 1, sum))
  y <- qgamma(ppoints(n),3,17)
  plot( sort(x) ~ sort(y), log='xy' )
  abline(0,1, col='red')
  title("Comparision: gamma distribution and sum of exponential r.v.")
dev.off()

png(file="g611.png", width=600, height=600)
  curve( dbeta(x,1,1), xlim=c(0,1), ylim=c(0,4) )
  curve( dbeta(x,2,1), add=T, col='red' )
  curve( dbeta(x,3,1), add=T, col='green' )
  curve( dbeta(x,4,1), add=T, col='blue' )
  curve( dbeta(x,2,2), add=T, lty=2, lwd=2, col='red' )
  curve( dbeta(x,3,2), add=T, lty=2, lwd=2, col='green' )
  curve( dbeta(x,4,2), add=T, lty=2, lwd=2, col='blue' )
  curve( dbeta(x,2,3), add=T, lty=3, lwd=3, col='red' )
  curve( dbeta(x,3,3), add=T, lty=3, lwd=3, col='green' )
  curve( dbeta(x,4,3), add=T, lty=3, lwd=3, col='blue' )
  title(main="Beta distribution")
  legend(par('usr')[1], par('usr')[4], xjust=0,
         c('(1,1)', '(2,1)', '(3,1)', '(4,1)', 
           '(2,2)', '(3,2)', '(4,2)',
           '(2,3)', '(3,3)', '(4,3)' ),
         lwd=1, #c(1,1,1,1, 2,2,2, 3,3,3),
         lty=c(1,1,1,1, 2,2,2, 3,3,3),
         col=c(par('fg'), 'red', 'green', 'blue', 
               'red', 'green', 'blue', 
               'red', 'green', 'blue' ))
dev.off()

png(file="g612.png", width=600, height=600)
  N <- 500
  n <- 5
  y <- drop(apply( matrix( runif(n*N), nr=N, nc=n), 1, max ))
  x <- qbeta(ppoints(N), n, 1)
  plot( sort(y) ~ x )
  abline(0,1, col='red')
  title("Order statistic and Beta distribution")
dev.off()

png(file="g613.png", width=600, height=600)
  N <- 500
  n <- 5
  k <- 3
  y <- drop(apply( matrix( runif(n*N), nr=n, nc=N), 2, sort )[n-k,])
  x <- qbeta(ppoints(N), n-k, k+1) # Exercice: Where do those
                                   # coefficients come from?
  plot( sort(y) ~ x )
  abline(0,1, col='red')
  title("Order statistics and Beta distribution")
dev.off()

png(file="g614.png", width=600, height=600)
  # I admit it: I found the coefficients above by trial-and-error...
  op <- par(mfrow=c(5,5), mar=c(0,0,0,0) )
  for (i in 1:5) {
    for (j in 1:5) {
      plot( sort(y) ~ qbeta(ppoints(N), j, i), xlab='', ylab='', axes=F )
      abline(0,1, col='red')
      box()
      text( (par('usr')[1]+par('usr')[2])/2, 
            (par('usr')[3]+par('usr')[4])/2, 
            paste(j,i),
            cex=3, col='blue' )
    }
  }
  par(op)
dev.off()

png(file="g615.png", width=600, height=600)
  curve(dbeta(x,8,4),xlim=c(0,1))
  title(main="posterior distrobution of p")
dev.off()

png(file="g616.png", width=600, height=600)
  curve(dbeta(x,10,10), xlim=c(0,1), lwd=3)
  curve(dbeta(x,1,1), add=T, col='red', lwd=3)
  curve(dbeta(x,2,2), add=T, col='green', lwd=3)
  curve(dbeta(x,5,2), add=T, col='blue',lwd=3)
  curve(dbeta(x,.1,.5), add=T, col='orange')
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('B(10,10)', 'B(1,1)', 'B(2,2)', 'B(5,2)', 'B(.1,.5)'),
         lwd=c(3,3,3,3,1), lty=1,
         col=c(par('fg'),'red','green','blue','orange'))
  title("A few beta probability distributions")
dev.off()

png(file="g617.png", width=600, height=600)
  # Example from the manual
  library(nor1mix)
  ppos <- which("package:nor1mix" == search())
  nms <- ls(pat="^MW.nm", pos = ppos)
  nms <- nms[order(as.numeric(substring(nms,6)))]
  op <- par(mfrow=c(4,4), mgp = c(1.2, 0.5, 0), tcl = -0.2,
            mar = .1 + c(2,2,2,1), oma = c(0,0,3,0))
  for(n in nms) { plot(get(n, pos = ppos)) }
    mtext("The Marron-Wand Densities", outer = TRUE,
          font = 2, cex = 1.6)
  par(op)
dev.off()

png(file="g618.png", width=600, height=600)
  rmnorm <- function (R, C, mu=rep(0,dim(C)[1])) {
    A <- t(chol(C))
    n <- dim(C)[1]
    t(mu + A %*% matrix(rnorm(R*n),nr=n))
  }
  C <- matrix(c(2,.5,.5,1),nr=2)
  mu <- c(2,1)
  y <- rmnorm(1000,C, mu)
  cov(y)
  plot(y)

  abline(h=mu[2], lty=3, lwd=3, col='blue')  
  abline(v=mu[1], lty=3, lwd=3, col='blue')  
  e <- eigen(C)
  r <- sqrt(e$values)
  v <- e$vectors
  N <- 100
  theta <- seq(0,2*pi, length=N)
  x <- mu[1] + r[1]*v[1,1]*cos(theta) +
       r[2]*v[1,2]*sin(theta)
  y <- mu[2] + r[1]*v[2,1]*cos(theta) +
       r[2]*v[2,2]*sin(theta)
  lines(x,y, lwd=3, col='blue')
dev.off()

png(file="g619.png", width=600, height=600)
  library(MASS)
  N <- 200
  mu <- c(1, .5)
  S <- matrix(c(4,2,2,2), nc=2)
  x <- mvrnorm(N, mu, S)
  plot(x)
  # Cloud center
  points(mu[1], mu[2], pch=3, cex=3, lwd=5, col="blue")
  # Ellipse axes
  e <- eigen(S)
  a <- mu + sqrt(e$values[1]) * e$vectors[,1]
  b <- mu - sqrt(e$values[1]) * e$vectors[,1]
  segments(a[1], a[2], b[1], b[2], col="blue", lwd=3)
  a <- mu + sqrt(e$values[2]) * e$vectors[,2]
  b <- mu - sqrt(e$values[2]) * e$vectors[,2]
  segments(a[1], a[2], b[1], b[2], col="blue", lwd=3)
  # TODO: Draw the ellipse...
dev.off()

png(file="g620.png", width=600, height=600)
  x <- read.table("SUNW.csv", header=T, sep=",")
  x <- x$Close
  x <- diff(log(x))   # Compute the returns
  qqnorm(x, main="Normality of stock returns (15 years)")
  qqline(x, col='red')
dev.off()

png(file="g621.png", width=600, height=600)
  qqnorm(x[1:90], main="Normality of stock returns (3 months)")
  qqline(x[1:90], col='red')
dev.off()

png(file="g622.png", width=600, height=600)
  plot(density(x), main="Stock returns density vs normal density")
  m <- mean(x)
  s <- sd(x)
  curve( dnorm(x, m, s), col='red', add=T)
dev.off()

png(file="g623.png", width=600, height=600)
  plot(density(x), log='y', main="Stock returns density vs normal density")
  curve( dnorm(x, m, s), col='red', add=T)
dev.off()

png(file="g624.png", width=600, height=600)
  plot(density(x), log='y', ylim=c(1e-1,1.5e1), xlim=c(-.2,.2),
       main="Stock returns density vs normal density (detail)")
  curve( dnorm(x, m, s), col='red', add=T)
dev.off()

png(file="g625.png", width=600, height=600)
  library(evd)
  curve(dfrechet(x, shape=1), lwd=3, xlim=c(-1,2), ylim=c(0,1),
        ylab="", main="density of the Frechet distribution")
  for (s in seq(.1,2,by=.1)) {
    curve(dfrechet(x, shape=s), add=T)
  }
  curve(dfrechet(x, shape=2), lwd=3, add=T, col='red')
  curve(dfrechet(x, shape=.5), lwd=3, add=T, col='blue')
dev.off()

png(file="g626.png", width=600, height=600)
  curve(drweibull(x, shape=1), lwd=3, xlim=c(-2,1), ylim=c(0,1),
        ylab="", main="density of the (reverse) Weibull distribution")
  for (s in seq(.1,2,by=.1)) {
    curve(drweibull(x, shape=s), add=T)
  }
  curve(drweibull(x, shape=2), lwd=3, add=T, col='red')
  curve(drweibull(x, shape=.5), lwd=3, add=T, col='blue')
dev.off()

png(file="g627.png", width=600, height=600)
  curve(dgumbel(x), lwd=3, xlim=c(-2,2), ylim=c(0,1),
        ylab="", main="density of the Gumbel distribution")
dev.off()

png(file="g628.png", width=600, height=600)
  curve(dgev(x, shape=0), lwd=3, xlim=c(-2,2), ylim=c(0,1),
        ylab="", main="density of the GEV distribution")
  for (s in seq(-2,2,by=.1)) {
    curve(dgev(x, shape=s), add=T)
  }
  curve(dgev(x, shape=-1), lwd=3, add=T, col='red')
  curve(dgev(x, shape=1), lwd=3, add=T, col='blue')
dev.off()

png(file="g629.png", width=600, height=600)
  x <- read.table("SUNW.csv", header=T, sep=",")
  x <- x$Close
  x <- diff(log(x))   # Compute the returns
  x <- sort(x)
  n <- length(x)
  m <- floor(.9*n)
  y <- (1:n)/n
  plot(x,y, type='l', main="empirical cdf of stock returns")
  lines(x[m:n], y[m:n], col='red', lwd=3)
dev.off()

png(file="g630.png", width=600, height=600)
  plot(x[m:n], y[m:n], col='red', lwd=3,
       type='l', main="empirical conditionnal excess cdf of stock returns")
dev.off()

png(file="g631.png", width=600, height=600)
  plot(x[m:n], (y[m:n] - y[m])/(1-y[m]) , col='red', lwd=3,
       type='l', main="rescaled empirical conditionnal excess cdf of stock returns")
dev.off()

png(file="g632.png", width=600, height=600)
  curve(dgpd(x, shape=0), lwd=3, xlim=c(-.1,2), ylim=c(0,2),
        ylab="", main="density of the Generalized Pereto Distribution (GPD)")
  for (s in seq(-2,2,by=.1)) {
    curve(dgpd(x, shape=s), add=T)
  }
  curve(dgpd(x, shape=-1), lwd=3, add=T, col='red')
  curve(dgpd(x, shape=1), lwd=3, add=T, col='blue')
dev.off()

png(file="g633.png", width=600, height=600)
  curve(pgpd(x, shape=0), lwd=3, xlim=c(-.1,2), ylim=c(0,1),
        ylab="", main="cdf of the Generalized Pereto Distribution (GPD)")
  for (s in seq(-2,2,by=.1)) {
    curve(pgpd(x, shape=s), add=T)
  }
  curve(pgpd(x, shape=-1), lwd=3, add=T, col='red')
  curve(pgpd(x, shape=1), lwd=3, add=T, col='blue')
dev.off()

png(file="g634.png", width=600, height=600)
  x <- read.table("SUNW.csv", header=T, sep=",")
  x <- x$Close
  x <- diff(log(x))   # Compute the returns
  x <- sort(x)
  n <- length(x)
  m <- floor(.9*n)
  y <- (1:n)/n
  op <- par(mfrow=c(3,3), mar=c(2,2,2,2))
  for (s in seq(0,2,length=9)) {
    plot(qgpd(ppoints(n-m+1),shape=s), x[m:n],
         xlab='', ylab='')
  }
  par(op)
dev.off()

png(file="g635.png", width=600, height=600)
  qqnorm(x)
dev.off()

png(file="g636.png", width=600, height=600)
  x <- sort(x)
  e <- rep(NA, length(x))
  for (i in seq(along=x)) {
    u <- x[i]
    e[i] <- mean( (x-u)[x>u] )
  }
  plot(x, e, type='o', main="Sample mean excess plot")
dev.off()

png(file="g637.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (f in c("SUNW", "ADBE", "ADPT", 
              "PCAR", "COST", "INTC",
              "MSFT", "ADCT", "BMET")) {
    x <- read.table(paste(f,".csv", sep=''), header=T, sep=",")
    x <- x$Close
    x <- diff(log(x))   # Compute the returns
    x <- sort(x)
    e <- rep(NA, length(x))
    for (i in seq(along=x)) {
      u <- x[i]
      e[i] <- mean( (x-u)[x>u] )
    }
    plot(x, e, type='o', main="Sample mean excess plot")
  }
  par(op)
dev.off()

png(file="g638.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (f in c("SUNW", "ADBE", "ADPT", 
              "PCAR", "COST", "INTC",
              "MSFT", "ADCT", "BMET")) {
    x <- read.table(paste(f,".csv", sep=''), header=T, sep=",")
    x <- x$Close
    x <- -diff(log(x))   # Compute the returns
    x <- sort(x)
    e <- rep(NA, length(x))
    for (i in seq(along=x)) {
      u <- x[i]
      e[i] <- mean( (x-u)[x>u] )
    }
    plot(x, e, type='o', main="Sample mean excess plot")
  }
  par(op)
dev.off()
detach.everything()

png(file="g639.png", width=600, height=600)
  colorie <- function (x, y1, y2, N=1000, ...) {
    for (t in (0:N)/N) {
      lines(x, t*y1+(1-t)*y2, ...)
    }
  }
  # No, there is already a function to do this
  colorie <- function (x, y1, y2, ...) {
    polygon( c(x, x[length(x):1]), c(y1, y2[length(y2):1]), ... )
  }
  x <- seq(-6,6, length=100)
  y <- dnorm(x)
  plot(y~x, type='l')
  i = x<qnorm(.025)
  colorie(x[i],y[i],rep(0,sum(i)) ,col='red')
  i = x>qnorm(.975)
  colorie(x[i],y[i],rep(0,sum(i)) ,col='red')
  lines(y~x)
  title(main="Type I error")
dev.off()

png(file="g640.png", width=600, height=600)
  x <- seq(-6,6, length=1000)
  y <- dnorm(x)
  plot(y~x, type='l')
  y2 <- dnorm(x-.5)
  lines(y2~x)
  i <- x>qnorm(.025) & x<qnorm(.975)
  colorie(x[i],y2[i],rep(0,sum(i)), col='red')
  segments( qnorm(.025),0,qnorm(.025),dnorm(qnorm(.025)), col='red' )
  segments( qnorm(.975),0,qnorm(.975),dnorm(qnorm(.975)), col='red' )
  lines(y~x)
  lines(y2~x)
  title(main="High risk of type II error")
dev.off()

png(file="g641.png", width=600, height=600)
  x <- seq(-6,6, length=1000)
  y <- dnorm(x)
  plot(y~x, type='l')
  y2 <- dnorm(x-3.5)
  lines(y2~x)
  i <- x>qnorm(.025) & x<qnorm(.975)
  colorie(x[i],y2[i],rep(0,sum(i)), col='red')
  segments( qnorm(.025),0,qnorm(.025),dnorm(qnorm(.025)), col='red' )
  segments( qnorm(.975),0,qnorm(.975),dnorm(qnorm(.975)), col='red' )
  lines(y~x)
  lines(y2~x)
  title(main="Lower risk of type II error")
dev.off()

png(file="g642.png", width=600, height=600)
  delta <- seq(-1.5, 1.5, length=500)
  p <- NULL
  for (d in delta) {
    p <- append(p, 
                power.t.test(delta=abs(d), sd=1, sig.level=0.05, n=20,
                             type='one.sample')$power)
  }
  plot(1-p~delta, type='l',
       xlab='mean difference', ylab="propability of a type II error",
       main="type II error in a Student T test")
  abline(h=0,lty=3)
  abline(h=0.05,lty=3)
  abline(v=0,lty=3)
dev.off()

png(file="g643.png", width=600, height=600)
  delta <- seq(0, 1.5, length=100)
  p <- NULL
  for (d in delta) {
    p <- append(p, 
                power.t.test(delta=d, sd=1, sig.level=0.05, n=20,
                             type='one.sample')$power)
  }
  plot(p~delta, type='l',
       ylab='power', main='Power of a one-sample t-test')
dev.off()

png(file="g644.png", width=600, height=600)
  N <- seq(10,200, by=5)
  delta <- NULL
  for (n in N) {
    delta <- append(delta, 
                    power.t.test(n=n, sd=1, sig.level=.05, 
                                 power=.80, type='one.sample')$delta
                   )
  }
  plot(delta~N, type='l', xlab='sample size')
  delta <- NULL
  for (n in N) {
    delta <- append(delta, 
                    power.t.test(n=n, sd=1, sig.level=.01, 
                                 power=.80, type='one.sample')$delta
                   )
  }
  lines(delta~N, col='red')
  delta <- NULL
  for (n in N) {
    delta <- append(delta, 
                    power.t.test(n=n, sd=1, sig.level=.001, 
                                 power=.80, type='one.sample')$delta
                   )
  }
  lines(delta~N, col='blue')
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('significance level=.05', 'significance level=.01', 'significance level=.001'),
         col=c(par('fg'), 'red', 'blue'),
         lwd=1, lty=1)
  title(main='Sample size and difference detected for tests of pover 0.80')
dev.off()

png(file="g645.png", width=600, height=600)
  p <- c()
  for (i in 1:1000) {
    x <- rnorm(200)
    p <- append(p, t.test(x)$p.value)
  }
  hist(p, col='light blue')
dev.off()

png(file="g646.png", width=600, height=600)
  p <- sort(p)
  p[950]
  p[50]
  x <- 1:1000
  plot(p ~ x, main="p-value of a Student T test when H0 is true")
dev.off()

png(file="g647.png", width=600, height=600)
  # Sample size
  n <- 10
  # Number of simulations
  # (sufficiently large to get a good approximation
  # of the probability)
  m <- 1000
  # Number of points to draw the curve
  k <- 50
  # Maximum value for the actual mean
  M <- 2
  r <- vector()
  for (j in M*(0:k)/k) {
    res <- 0
    for (i in 1:m) {
      x <- rnorm(10, mean=j)
      if( t.test(x)$p.value > .05 ){
        res <- res + 1
      }
    }
    r <- append(r, res/m)
  }
  rr <- M*(0:k)/k
  plot(r~rr, type="l",
       xlab='difference in mean',
       ylab="type II error probability")

  # Comparison with the curve produced by "power.t.test"
  x <- seq(0,2,length=200)
  y <- NULL
  for (m in x) {
    y <- append(y, 1-power.t.test(delta=m, sd=1, n=10, sig.level=.05, 
                                  type='one.sample')$power)
  }
  lines(x,y,col='red')

  # Theoretical curve
  # (This is a Z test, not too different)
  r2 <- function (p,q,conf,x) {
    p(q(1-conf/2)-x) - p(q(conf/2)-x)
  }
  f <- function(x) {
    p <- function (t) { pnorm(t, sd=1/sqrt(10)) }
    q <- function (t) { qnorm(t, sd=1/sqrt(10)) }
    r2(p,q,.05,x)
  }
  curve( f(x) , add=T, col="blue" )

  # Theoretical curve (T test)
  f <- function(x) {
    p <- function (t) { pt(sqrt(10)*t, 10) } # Is this correct?
    q <- function (t) { qt(t, 10)/sqrt(10) }
    r2(p,q,.05,x)
  }
  curve( f(x) , add=T, col="green" )

  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('simulation', 'power.t.test', '"exact" value, Z test',
           'excat value'),
         col=c(par('fg'),'red','blue','green'),
         lwd=1,lty=1)
  title(main="Type II error risk in a Student T test")
dev.off()

png(file="g648.png", width=600, height=600)
  N <- 10000
  x <- 100*(1:N)/N
  plot( x~I(x/100), type='n', ylab="cumulated percents", xlab="p-value" )
  do.it <- function (m, col) {
    p <- c()
    for (i in 1:N) {
      x <- m+rnorm(200)
      p <- append(p, t.test(x)$p.value)
    }
    p <- sort(p)
    x <- 100*(1:N)/N
    lines(x ~ p, type='l', col=col)
  }
  do.it(0,   par('fg'))
  do.it(.05, 'red')
  do.it(.1,  'green')
  do.it(.15, 'blue')
  do.it(.2,  'orange')
  abline(v=.05)
  title(main='p-value distribution')
  legend(par('usr')[2],par('usr')[3],xjust=1,yjust=1,
         c('m=0', 'm=0.05', 'm=.01', 'm=.015', 'm=.02'),
         col=c(par('fg'), 'red', 'green', 'blue', 'orange'),
         lty=1,lwd=1)
dev.off()

png(file="g649.png", width=600, height=200)
  x <- 1:10
  y <- c(7:20, 200)
  boxplot(x,y, horizontal=T)
dev.off()

png(file="g650.png", width=600, height=200)
  boxplot(x,y, log="x", horizontal=T)
dev.off()

png(file="g651.png", width=600, height=600)
  curve(dnorm(x), from=-5, to=5, add=F, col="orange", lwd=3, lty=2)  
  curve(dt(x,100), from=-5, to=5, add=T, col=par('fg'))  
  curve(dt(x,5),  from=-5, to=5, add=T, col="red")  
  curve(dt(x,2),  from=-5, to=5, add=T, col="green")  
  curve(dt(x,1),  from=-5, to=5, add=T, col="blue")  
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('gaussian', 'df=100', 'df=5', 'df=2', 'df=1'),
         col=c('orange', par('fg'), 'red', 'green', 'blue'),
         lwd=c(3,1,1,1,1),
         lty=c(2,1,1,1,1))
  title(main="Student's T probability distribution function")
dev.off()

png(file="g652.png", width=600, height=600)
  N <- 50
  n <- 5
  v <- matrix(c(0,0),nrow=2)
  for (i in 1:N) {
    x <- rnorm(n)
    v <- cbind(v, t.test(x)$conf.int)
  }
  v <- v[,2:(N+1)]
  plot(apply(v,2,mean), ylim=c(min(v),max(v)),
       ylab='Confidence interval', xlab='')
  abline(0,0)
  c <- apply(v,2,min)>0 | apply(v,2,max)<0
  segments(1:N,v[1,],1:N,v[2,], col=c(par('fg'),'red')[c+1], lwd=3)
  title(main="The population mean need not be in the confidence interval")
dev.off()

png(file="g653.png", width=600, height=600)
  N <- 1000
  n <- 10
  delta <- .8
  v <- vector()
  w <- vector()
  for (i in 1:N) {
    x <- delta+runif(n, min=-1, max=1)
    v <- append(v, t.test(x)$p.value)
    w <- append(w, wilcox.test(x)$p.value)
  }
  plot(sort(v), type='l', lwd=3, lty=2, ylab="p-valeur")
  lines(sort(w), col='red')
  legend(par('usr')[1], par('usr')[4], xjust=0,
         c('Student', 'Wilcoxon'),
         lwd=c(2,1),
         lty=c(2,1),
         col=c(par("fg"), 'red'))
dev.off()

png(file="g654.png", width=600, height=600)
  N <- 1000
  n <- 100
  delta <- .1
  v <- vector()
  w <- vector()
  for (i in 1:N) {
    x <- delta+runif(n, min=-1, max=1)
    v <- append(v, t.test(x)$p.value)
    w <- append(w, wilcox.test(x)$p.value)
  }
  plot(sort(v), type='l', lwd=3, lty=2)
  lines(sort(w), col='red')
  legend(par('usr')[1], par('usr')[4], xjust=0,
         c('Student', 'Wilcoxon'),
         lwd=c(2,1),
         lty=c(2,1),
         col=c(par("fg"), 'red'))
dev.off()

png(file="g655.png", width=600, height=600)
  N <- 1000
  n <- 100
  delta <- .8
  v <- vector()
  w <- vector()
  for (i in 1:N) {
    x <- delta+runif(n, min=-1, max=1)
    v <- append(v, t.test(x)$p.value)
    w <- append(w, wilcox.test(x)$p.value)
  }
  plot(sort(v), type='l', lwd=3, lty=2)
  lines(sort(w), col='red')
  legend(par('usr')[1], par('usr')[4], xjust=0,
         c('Student', 'Wilcoxon'),
         lwd=c(2,1),
         lty=c(2,1),
         col=c(par("fg"), 'red'))
dev.off()

png(file="g656.png", width=600, height=600)
  N <- 1000
  n <- 10
  delta <- 1
  v <- vector()
  w <- vector()
  for (i in 1:N) {
    x <- delta+rcauchy(n)
    v <- append(v, t.test(x)$p.value)
    w <- append(w, wilcox.test(x)$p.value)
  }
  plot(sort(v), type='l', lwd=3, lty=2)
  lines(sort(w), col='red')
dev.off()

png(file="g657.png", width=600, height=600)
  N <- 1000
  n <- 100
  delta <- 1
  v <- vector()
  w <- vector()
  for (i in 1:N) {
    x <- delta+rcauchy(n)
    v <- append(v, t.test(x)$p.value)
    w <- append(w, wilcox.test(x)$p.value)
  }
  plot(sort(v), type='l', lwd=3, lty=2)
  lines(sort(w), col='red')
dev.off()

png(file="g658.png", width=600, height=600)
  N <- 1000
  n <- 10
  v <- 100
  a <- NULL
  b <- NULL
  for (i in 1:N) {
    x <- rnorm(n)
    y <- rnorm(n, 0, v)
    a <- append(a, t.test(x,y)$p.value)
    b <- append(b, t.test(x/var(x), y/var(y))$p.value)
  }
  plot(sort(a), type='l', col="green")
  points(sort(b), type="l", col="red")
  abline(0,1/N)
dev.off()

png(file="g659.png", width=600, height=600)
  N <- 1000
  n <- 10
  v <- 100
  a <- NULL
  b <- NULL
  c <- NULL
  d <- NULL
  for (i in 1:N) {
    x <- rnorm(n)
    y <- rnorm(n, 100, v)
    a <- append(a, t.test(x,y)$p.value)
    b <- append(b, t.test(x/var(x), y/var(y))$p.value)
    c <- append(c, t.test(x, y/10000)$p.value)
    d <- append(d, wilcox.test(x, y)$p.value)
  }
  plot(sort(a), type='l', col="green")
  points(sort(b), type="l", col="red")
  points(sort(c), type="l", col="blue")
  points(sort(d), type="l", col="orange")
  abline(0,1/N)
  legend(par('usr')[1], par('usr')[4], 
         c('Student Test', 'Renormalized Student Test',
           'Student Test renormalized with the sample variances',
           "(Non-Parametric) Wilcoxson's U Test"),
         col=c('green', 'blue', 'red', 'orange'),
         lwd=1,lty=1)
  title(main="Student's T test on non-equivariant samples")
dev.off()

png(file="g660.png", width=600, height=300)
  data(sleep)
  boxplot(extra ~ group, data=sleep,
          horizontal=T,
          xlab='extra', ylab='group')
dev.off()

png(file="g661.png", width=600, height=600)
  curve(dchisq(x,2),  from=0, to=5, add=F, col="red",
        ylab="dchisq(x,i)") 
  n <- 10
  col <- rainbow(n)
  for (i in 1:n) { 
    curve(dchisq(x,i),  from=0, to=5, add=T, col=col[i])
  }
  legend(par('usr')[2], par('usr')[4], xjust=1,
          paste('df =',1:n),
         lwd=1,
         lty=1,
         col=col)
  title(main="Chi^2 Probability Distribution Function")
dev.off()

png(file="g662.png", width=600, height=600)
  chisq.var.test <- function (x, var=1, conf.level=.95,
                              alternative='two.sided') {
    result <- list()
    alpha <- 1-conf.level
    n <- length(x)
    v <- var(x)
    result$var <- v
    result$sd <- sd(x)
    chi2 <- (n-1)*v/var
    result$chi2 <- chi2
    p <- pchisq(chi2,n-1)
    if( alternative == 'less' ) {
      stop("Not implemented yet")
    } else if (alternative == 'greater') {
      stop("Not implemented yet")
    } else if (alternative == 'two.sided') {
      if(p>.5) 
        p <- 1-p
      p <- 2*p
      result$p.value <- p
      result$conf.int.var <- c(
        (n-1)*v / qchisq(alpha/2, df=n-1, lower.tail=F),
        (n-1)*v / qchisq(alpha/2, df=n-1, lower.tail=T),
      )
    }
    result$conf.int.sd <- sqrt( result$conf.int.var )
    result
  }
  x <- rnorm(100)
  chisq.var.test(x)

  # We can check tha the results are correct by looking at
  # the distribution of the p-values: it should be uniform
  # in [0,1].
  v <- NULL
  for (i in 1:1000) { 
    v <- append(v, chisq.var.test(rnorm(100))$p.value) 
  } 
  plot(sort(v))
dev.off()

png(file="g663.png", width=600, height=600)
  p1 <- NULL
  p2 <- NULL
  for (i in 1:100) {
    x <- rnorm(10)
    p1 <- append(p1, chisq.var.test(x)$p.value)
    p2 <- append(p2, var.test(x, rnorm(10000))$p.value)
  }
  plot( p1 ~ p2 )
  abline(0,1,col='red')
dev.off()

png(file="g664.png", width=600, height=600)
  p <- .3
  col.values <- c(par('fg'),'red', 'blue', 'green', 'orange')
  n.values <- c(5,10,20,50,100)
  plot(0, type='n', xlim=c(0,1), ylim=c(0,1), xlab='exact', ylab='approximate')
  for (i in 1:length(n.values)) {
    n <- n.values[i]
    x <- NULL
    y <- NULL
    for (a in 0:n) {
      x <- append(x, binom.test(a,n,p)$p.value)
      y <- append(y, prop.test(a,n,p)$p.value)
    }
    o <- order(x)
    lines(x[o],y[o], col=col.values[i])
  }
  legend(par('usr')[1],par('usr')[4],
         as.character(n.values),
         col=col.values,
         lwd=1,lty=1)
  title(main="Comparing the binomial test and its approximation")
dev.off()

png(file="g665.png", width=600, height=600)
  p <- .3
  n <- 5
  N <- 1000
  e <- rbinom(N, n, p)
  x <- y <- NULL
  for (a in e) {
    x <- append(x, binom.test(a,n,p)$p.value)
    y <- append(y, prop.test(a,n,p)$p.value)
  }
  x <- sort(x)
  y <- sort(y)
  plot(x, type='l', lwd=3, ylab='p-value')
  lines(y, col='red')
  legend(par('usr')[2], par('usr')[3], xjust=1, yjust=0,
         c('exact', 'approximation'),
         lwd=c(3,1),
         lty=1,
         col=c(par("fg"),'red'))
  title(main="Binomial test (H0 true)")
dev.off()

png(file="g666.png", width=600, height=600)
  p1 <- .3
  p2 <- .5
  n <- 5
  N <- 1000
  e <- rbinom(N, n, p1)
  x <- y <- NULL
  for (a in e) {
    x <- append(x, binom.test(a,n,p2)$p.value)
    y <- append(y, prop.test(a,n,p2)$p.value)
  }
  x <- sort(x)
  y <- sort(y)
  plot(x, type='l', lwd=3, ylab='p-value')
  lines(y, col='red')
  legend(par('usr')[2], par('usr')[3], xjust=1, yjust=0,
         c('exact', 'approximation'),
         lwd=c(3,1),
         lty=1,
         col=c(par("fg"),'red'))
  title(main="Binomial test (H0 false)")
dev.off()

png(file="g667.png", width=600, height=600)
  p <- .3
  col.values <- c(par('fg'),'red', 'blue', 'green', 'orange')
  n.values <- c(5,10,20,50,100)
  plot(0, type='n', xlim=c(0,1), ylim=c(0,1), xlab='exact', ylab='approximate')
  for (i in 1:length(n.values)) {
    n <- n.values[i]
    x <- NULL
    y <- NULL
    z <- NULL
    for (a in 0:n) {
      x <- append(x, binom.test(a,n,p)$p.value)
      y <- append(y, chisq.test(c(a,n-a),p=c(p,1-p))$p.value)
      z <- append(z, prop.test(a,n,p)$p.value)
    }
    o <- order(x)
    lines(x[o],y[o], col=col.values[i])
    lines(x[o],z[o], col=col.values[i], lty=3)
  }
  legend(par('usr')[1],par('usr')[4],
         as.character(c(n.values, "prop.test", "chisq.test")),
         col=c(col.values, par('fg'), par('fg')),
         lwd=1,
         lty=c(rep(1,length(n.values)), 1,3) 
        )
  title(main="The binomial test and its approximations")
dev.off()

png(file="g668.png", width=600, height=600)
  # A Monte Carlo multinomial test
  multinom.test <- function (x, p, N=1000) {
    n <- sum(x)
    m <- length(x)
    chi2 <- sum( (x-n*p)^2/(n*p) )
    v <- NULL
    for (i in 1:N) {
      x <- table(factor(sample(1:m, n, replace=T, prob=p), levels=1:m))
      v <- append(v, sum( (x-n*p)^2/(n*p) ))
    }
    sum(v>=chi2)/N
  }
  multinom.test( c(25,40,25,25), p=c(.25,.25,.25,.25) )   # 0.13
  chisq.test(    c(25,40,25,25), p=c(.25,.25,.25,.25) )   # 0.12

  N <- 100
  m <- 4
  n <- 10
  p <- c(.25,.25,.1,.4)
  x <- NULL
  y <- NULL
  for (i in 1:N) {
    a <- table( factor(sample(1:m, n, replace=T, prob=p), levels=1:m) )
    x <- append(x, multinom.test(a,p))
    y <- append(y, chisq.test(a,p=p)$p.value)
  }
  plot(y~x)
  abline(0,1,col='red')
  title("Monte Carlo Multinomial Test and Chi^2 Test")
dev.off()

png(file="g669.png", width=600, height=600)
  # We sample 10 subjects in a 4-class population.
  # We repeat the experiment 100 times.
  N <- 1000
  m <- 4
  n <- 10
  p <- c(.24,.26,.1,.4)
  p.valeur.chi2 <- rep(NA,N)
  for (i in 1:N) {
    echantillon <- table(factor(sample(1:m, replace=T, prob=p), levels=1:m))
    p.valeur.chi2[i] <- chisq.test(echantillon,p=p)$p.value
  }
  plot( sort(p.valeur.chi2), type='l', lwd=3 )
  abline(0, 1/N, lty=3, col='red', lwd=3)
  title(main="p-values in a Chi^2 test")
dev.off()

png(file="g670.png", width=600, height=600)
  foo <- function (N) {
    population1 <- c(rep('A',10), rep('B',20), rep('A',60), rep('B',10))
    population1 <- factor(population1, levels=c('A','B'))
    population2 <- c(rep('C',10), rep('C',20), rep('D',60), rep('D',10))
    population2 <- factor(population2, levels=c('C','D'))
    o <- sample(1:100, N, replace=T)
    table( population2[o], population1[o] )
  }
  a <- foo(1000)
  op <- par(mfcol=c(1,2))
  plot( a, shade=T )
  plot( t(a), shade=T )
  par(op)
dev.off()

png(file="g671.png", width=600, height=600)
  n1 <- 10
  n2 <- 100
  N <- 1000
  x1 <- rep(NA,N)
  x2 <- rep(NA,N)
  for (i in 1:N) {
    x1[i] <- fisher.test(foo(n1))$p.value
    x2[i] <- fisher.test(foo(n2))$p.value
  }
  plot( sort(x1), type='l', lwd=3, ylab='p-valeur')
  lines( sort(x2), col='blue', lwd=3 )
  abline(0,1/N,col='red',lwd=3,lty=3)
  abline(h=c(0,.05),lty=3)
  abline(v=c(0,N*.05),lty=3)
  title(main="p-value of a Fisher test, H0 false")
  legend(par('usr')[1],par('usr')[4],
         c("n=10", "n=100"),
         col=c(par('fg'), 'blue'),
         lwd=3,
         lty=1)
dev.off()

png(file="g672.png", width=600, height=600)
  foo <- function (N) {
    population1 <- c(rep('A',2), rep('B',8), rep('A',18), rep('B',72))
    population1 <- factor(population1, levels=c('A','B'))
    population2 <- c(rep('C',2), rep('C',8), rep('D',18), rep('D',72))
    population2 <- factor(population2, levels=c('C','D'))
    o <- sample(1:100, N, replace=T)
    table( population2[o], population1[o] )
  }
  n1 <- 10
  n2 <- 100
  N <- 1000
  x1 <- rep(NA,N)
  x2 <- rep(NA,N)
  for (i in 1:N) {
    x1[i] <- fisher.test(foo(n1))$p.value
    x2[i] <- fisher.test(foo(n2))$p.value
  }
  plot( sort(x1), type='l', lwd=3, ylab='p-valeur', ylim=c(0,1))
  lines( sort(x2), col='blue', lwd=3 )
  abline(0,1/N,col='red',lwd=3,lty=3)
  abline(h=c(0,.05),lty=3)
  abline(v=c(0,N*.05),lty=3)
  title(main="p-valueof a Fisher test, H0 true")
  legend(par('usr')[2], .2, xjust=1, yjust=0,
         c("n=10", "n=100"),
         col=c(par('fg'), 'blue'),
         lwd=3,
         lty=1)
dev.off()

png(file="g673.png", width=600, height=600)
  sign.test <- function (x, mu=0) {
    n <- length(x)
    y <- sum(x<mu) # should warn about ties!
    p.value <- min(c( pbinom(y,n,.5), pbinom(y,n,.5,lower.tail=F) ))*2
    p.value
  }
  N <- 500
  n <- 200
  res <- rep(NA,N)
  for (i in 1:N) {
    res[i] <- sign.test(rlnorm(n),mu=1)
  }
  plot(sort(res))
  abline(0,1/N,lty=2)
dev.off()

png(file="g674.png", width=600, height=600)
  N <- 500
  n <- 10
  res <- rep(NA,N)
  for (i in 1:N) {
    res[i] <- sign.test(rlnorm(n),mu=2)
  }
  plot(sort(res), ylim=c(0,1))
  abline(0,1/N,lty=2)
dev.off()

png(file="g675.png", width=600, height=600)
  n <- 5
  x <- y <- z <- vector()
  for(i in 1:10000){
    t <- rnorm(n)
    z <- append(z, sum(t*t)/n)
    t <- t - mean(t)
    t <- t*t
    x <- append(x, sum(t)/n)
    y <- append(y, sum(t)/(n-1))
  }
  boxplot(x,y,z)
dev.off()

png(file="g676.png", width=600, height=600)
  plot(density(x))
  points(density(y), type="l", col="red")
  points(density(z), type="l", col="blue")
dev.off()

png(file="g677.png", width=600, height=600)
  op <- par( mfrow = c(3,1) )
  hist(x, xlim=c(0,5), breaks=20)
  hist(y, xlim=c(0,5), breaks=20)
  hist(z, xlim=c(0,5), breaks=20)
  par(op)
dev.off()

png(file="g678.png", width=600, height=600)
  x <- sqrt(x)
  y <- sqrt(x)
  z <- sqrt(z)
  boxplot(x,y,z)
dev.off()

png(file="g679.png", width=600, height=600)
  plot(density(x))
  points(density(y), type="l", col="red")
  points(density(z), type="l", col="blue")
dev.off()

png(file="g680.png", width=600, height=600)
  op <- par( mfrow = c(3,1) )
  hist(x, xlim=c(0,5), breaks=20)
  hist(y, xlim=c(0,5), breaks=20)
  hist(z, xlim=c(0,5), breaks=20)
  par(op)
dev.off()

png(file="g681.png", width=600, height=600)
  # The population mean
  m <- runif(1, min=-1, max=1)
  # The n-element sample
  n <- 5
  v <- rnorm(n, mean=m)
  # Likelihood
  N <- 1000
  l <- seq(-2,2, length=N)
  y <- vector()
  for (i in l) {
    y <- append(y, prod(dnorm(v,mean=i)))
  }
  plot(y~l, type='l')
  # Population mean
  points(m, prod(dnorm(v,mean=m)), lwd=3)
  # Sample mean
  points(mean(v), prod(dnorm(v,mean=mean(v))), col='red', lwd=3)
dev.off()

png(file="g682.png", width=600, height=600)
  f <- function (x, p, m1, s1, m2, s2) {
    p*dnorm(x,mean=m1,sd=s1) + (1-p)*dnorm(x,mean=m2,sd=s2)
  }
  data(faithful)
  fn <- function(arg) {
    prod(f(faithful$eruptions, arg[1], arg[2], arg[3], arg[4], arg[5]))
  }
  start <- c(.5, 
             min(faithful$eruptions), var(faithful$eruptions),
             max(faithful$eruptions), var(faithful$eruptions),
            )
  p <- optim(start, function(a){-fn(a)/fn(start)}, control=list(trace=1))$par
  hist(faithful$eruptions, breaks=20, probability=T, col='light blue')
  lines(density(faithful$eruptions,bw=.15), col='blue', lwd=3)
  curve(f(x, p[1], p[2], p[3], p[4], p[5]), add=T, col='red', lwd=3)
  #curve(dnorm(x, mean=p[2], sd=p[3]), add=T, col='red', lwd=3, lty=2)
  #curve(dnorm(x, mean=p[4], sd=p[5]), add=T, col='red', lwd=3, lty=2)
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c('sample density', 'theoretical density'),
         lwd=3, lty=1,
         col=c('blue', 'red'))
dev.off()

png(file="g683.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  curve(dnorm(x-5)+dnorm(x+5), xlim=c(-10,10))
  curve(dnorm(x-5)+.4*dnorm(x+5), xlim=c(-10,10))
  curve(dnorm(5*(x-5)) + .5*dnorm(x+5), xlim=c(-10,10))
  par(op)
dev.off()

png(file="g684.png", width=600, height=600)
  get.sample <- function (n=10, p=1/100) {
    ifelse( runif(n)>p, runif(n), 2 )
  }
  N <- 1000
  d <- rep(NA,N) 
  for (i in 1:N) {
    d[i] <- max(get.sample())
  }
  hist(d, probability=T, ylim=c(0, max(density(d)$y)), col='light blue')
  lines(density(d), type='l', col='red', lwd=3)
dev.off()

png(file="g685.png", width=600, height=600)
  get.sample <- function (n=100, p=1/100) {
    ifelse( runif(n)>p, runif(n), 2 )
  }
  N <- 1000
  d <- rep(NA,N) 
  for (i in 1:N) {
    d[i] <- max(get.sample())
  }
  #hist(d, breaks=seq(0,3,by=.02),
  #     probability=T, ylim=c(0, max(density(d)$y)), col='light blue')
  #lines(density(d), type='l', col='red', lwd=3)
  plot(density(d), type='l', col='red', lwd=3)
dev.off()

png(file="g686.png", width=600, height=600)
  get.sample <- function (n=100, p=1/100) {
    ifelse( runif(n)>p, runif(n), 2+2*runif(n) )
  }
  N <- 1000
  d <- rep(NA,N) 
  for (i in 1:N) {
    d[i] <- max(get.sample())
  }
  #hist(d, breaks=seq(0,4,by=.05),
  #     probability=T, ylim=c(0, max(density(d)$y)), col='light blue')
  #lines(density(d), type='l', col='red', lwd=3)
  plot(density(d), type='l', col='red', lwd=3)
dev.off()

png(file="g687.png", width=600, height=600)
  curve(ifelse(x < .3, 0, ifelse(x > .7, 0, 1)),
        lwd = 3, col = "blue", 
        xlim = c(0,1),
        main = "A classical membership function",
        xlab="", ylab="")
  abline(h = 0, lty = 3)
dev.off()

png(file="g688.png", width=600, height=600)
  curve( dnorm( x - .5, sd = .1 ) / dnorm(0, sd=.1),
         lwd = 3, col = "blue", 
         xlim = c(0, 1),
         main = "A fuzzy membership function",
         xlab="", ylab="" )
  abline(h=0, lty=3)
dev.off()
detach.everything()

png(file="g689.png", width=600, height=600)
  N <- 1000
  a <- 1
  b <- -1
  Z <- rnorm(N)
  epsilon <- rnorm(N)
  eta <- rnorm(N)
  aa <- runif(1)
  bb <- runif(1)
  X <- (aa + bb * Z + epsilon) + eta
  Y <- a + b * X + epsilon
  plot(X,Y)
  abline(a,b, lty=2, lwd=3)
  abline(lm(Y~X), col="red", lwd=3)
dev.off()

png(file="g690.png", width=600, height=600)
  library(sem)
  r <- tsls(Y ~ X, instruments = ~ Z)
  plot(X,Y)
  abline(a,b, lty=2, lwd=3)
  abline(lm(Y~X), col="red", lwd=3)
  abline(r$coef[1], r$coef[2], col="blue", lwd=3)
dev.off()

png(file="g691.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (n in c(10,1e2,1e3,1e4)) {
    x <- runif(n)
    y <- 1 - x + .2*rnorm(n)
    plot(y~x, main=paste("sample with", n, "observations"))
  }
  par(op)
dev.off()

png(file="g692.png", width=600, height=600)
  data(anscombe)
  ff <- y ~ x
  op <- par(mfrow = c(2,2), 
            mar = .1 + c(4,4,1,1), 
            oma = c(0,0,2,0))
  for(i in 1:4) {
    ff[2:3] <- lapply(paste(c("y","x"), i, sep=""), as.name)
    ## or   ff[[2]] <- as.name(paste("y", i, sep=""))
    ##      ff[[3]] <- as.name(paste("x", i, sep=""))
    assign(paste("lm.",i,sep=""), 
           lmi <- lm(ff, data= anscombe))
    #print(anova(lmi))
  }
  for(i in 1:4) {
    ff[2:3] <- lapply(paste(c("y","x"), i, sep=""), as.name)
    plot(ff, data = anscombe, 
         col = "red", pch = 21, bg = "orange", cex = 1.2,
         xlim = c(3,19), ylim = c(3,13))
    abline(get(paste("lm.",i,sep="")), col="blue")
  }
  mtext("Anscombe's 4 Regression data sets", 
        outer = TRUE, cex = 1.5)
  par(op)
dev.off()

png(file="g693.png", width=600, height=600)
  data(cars)
  x <- cars$speed
  y <- cars$dist
  plot(y~x)
  o <- order(x)
  n <- length(x)
  m <- floor(n/2)
  p1 <- c( median(x[o[1:m]]), median(y[o[1:m]]) )
  m <- ceiling(n/2)
  p2 <- c( median(x[o[m:n]]), median(y[o[m:n]]) )
  p <- rbind(p1,p2)
  points(p, pch='+', lwd=3, cex=5, col='red' )
  lines(p, col='red', lwd=3)
  title(main="Brown-Mood Line")
dev.off()

png(file="g694.png", width=600, height=600)
  three.group.resistant.line <- function (y, x) { 
    o <- order(x)
    n <- length(x)
    o1 <- o[1:floor(n/3)]
    o2 <- o[ceiling(n/3):floor(2*n/3)]
    o3 <- o[ceiling(2*n/3):n]
    p1 <- c( median(x[o1]), median(y[o1]) )
    p2 <- c( median(x[o2]), median(y[o2]) )
    p3 <- c( median(x[o3]), median(y[o3]) )
    p <- rbind(p1,p2,p3)
    g <- apply(p,2,mean)
    plot(y~x)
    points(p, pch='+', lwd=3, cex=3, col='red')
    polygon(p, border='red')
    a <- (p3[2] - p1[2])/(p3[1] - p1[1])
    b <- g[2]-a*g[1]
    abline(b,a,col='red')
  }
  three.group.resistant.line(cars$dist, cars$speed)
dev.off()

png(file="g695.png", width=600, height=600)
  n <- 100
  x <- runif(n,min=0,max=2)
  y <- x*(1-x) + rnorm(n)
  three.group.resistant.line(y,x)
dev.off()

png(file="g696.png", width=600, height=600)
  median.line <- function (y,x) {
    n <- length(x)
    b <- matrix(NA, nr=n, nc=n)
    # Exercise: Write this without a loop
    for (i in 1:n) {
      for (j in 1:n) {
        if(i!=j)
          b[i,j] <- ( y[i] - y[j] )/( x[i] - x[j] )
      }
    }
    b <- median(b, na.rm=T)
    a <- median(y-b*x)
    plot(y~x)
    abline(a,b, col='red')
    title(main="The median line")
  }
  median.line(cars$dist, cars$speed)
dev.off()

png(file="g697.png", width=600, height=600)
  other.median.line <- function (y,x) {
   n <- length(x)
    b <- matrix(NA, nr=n, nc=n)
    for (i in 1:n) {
      for (j in 1:n) {
        if(i!=j)
          b[i,j] <- ( y[i] - y[j] )/( x[i] - x[j] )
      }
    }
    b <- median( apply(b, 1, median, na.rm=T), na.rm=T ) # Only change
    a <- median(y-b*x)
    plot(y~x)
    abline(a,b, col='red')
    title(main="The median line")
  }
  other.median.line(cars$dist, cars$speed)
dev.off()

png(file="g698.png", width=600, height=600)
  do.it <- function (x, y) {
    plot(x,y, main=paste("cor =", round(cor(x,y), digits=2)))
    abline(lm(y~x), col='red', lwd=3)
  }

  n <- 100
  x <- runif(n)
  x <- x[order(x)]
  y <- x
  do.it(x,y)
  abline(0,1,lty=2)
dev.off()

png(file="g699.png", width=600, height=600)
  y <- rnorm(n,x,.1)
  do.it(x,y)
  abline(0,1,lty=2)
dev.off()

png(file="g700.png", width=600, height=600)
  y <- rnorm(n,x,1)
  do.it(x,y)
  abline(0,1,lty=2)
dev.off()

png(file="g701.png", width=600, height=600)
  y <- runif(n)
  do.it(x,y)
dev.off()

png(file="g702.png", width=600, height=600)
  x <- runif(n,-1,1)
  y <- x*x
  do.it(x,y)
dev.off()

png(file="g703.png", width=600, height=600)
  y <- rnorm(n, x*x, .1)
  do.it(x,y)
  x <- sort(x)
  lines(x,x*x,lty=2)
dev.off()

png(file="g704.png", width=600, height=600)
  y <- rnorm(n, x*x, 1)
  do.it(x,y)
  x <- sort(x)
  lines(x,x*x,lty=2)
dev.off()

png(file="g705.png", width=600, height=600)
  x <- runif(n)
  y <- rnorm(n,-x,1)
  do.it(x,y)
  abline(0,-1,lty=2)
dev.off()

png(file="g706.png", width=600, height=600)
  y <- rnorm(n,-x,.1)
  do.it(x,y)
  abline(0,-1,lty=2)
dev.off()

png(file="g707.png", width=600, height=600)
  y <- -x
  do.it(x,y)
  abline(0,-1,lty=2)
dev.off()

png(file="g708.png", width=600, height=600)
  library(mvtnorm)
  k <- 100                # Number of samples for each correlation
  N <- 20                 # Size of the samples
  r <- seq(-1, 1, by=.2)  # The true correlations
  n <- length(r)
  rr <- matrix(NA, nr=n, nc=k)
  for (i in 1:n) {
    for (j in 1:k) {
      x <- rmvnorm(N, sigma=matrix(c(1, r[i], r[i], 1), nr=2, nc=2))
      rr[i,j] <- cor( x[,1], x[,2] )      
    }
  }
  estimated.correlation <- as.vector(rr)
  true.correlation <- r[row(rr)] 
  boxplot(estimated.correlation ~ true.correlation, 
          col = "pink",
          xlab = "True correlation", 
          ylab = "Estimated correlation" )
dev.off()

png(file="g709.png", width=600, height=600)
  boxplot(estimated.correlation - true.correlation ~ true.correlation, 
          col = "pink",
          xlab = "True correlation", ylab = "Error" )
dev.off()

png(file="g710.png", width=600, height=600)
  N <- 20    # Sample size
  n <- 1000  # Number of samples
  true.correlation <- runif(n, -1, 1)
  estimated.correlation <- rep(NA, n)
  for (i in 1:n) {
    x <- rmvnorm(N, sigma=matrix(c(1, true.correlation[i], 
                                   true.correlation[i], 1), nr=2, nc=2))
    estimated.correlation[i] <- cor( x[,1], x[,2] )      
  }
  plot(estimated.correlation ~ true.correlation)
  abline(0,1, col="blue", lwd=3)
dev.off()

png(file="g711.png", width=600, height=600)
  plot(estimated.correlation - true.correlation ~ true.correlation, ylab="Error")
  abline(h=0, lty=3)
dev.off()

png(file="g712.png", width=600, height=600)
  plot(abs(estimated.correlation - true.correlation) ~ true.correlation, 
       ylab="abs(Error)")
  lines(lowess(abs(estimated.correlation - true.correlation) ~ true.correlation),
        col="red", lwd=3)
  abline(h=0, lty=3)
dev.off()

png(file="g713.png", width=600, height=600)
  N <- 1000
  n <- 100
  x <- matrix(nr=N, nc=3)
  colnames(x) <- c("Pearson", "Spearman (rank)", 
                   "Kendall's tau")
  y1 <- 1:n
  for (i in 1:N) {
    y2 <- sample(y1)
    x[i,] <- c( cor(y1, y2), 
                cor(y1, y2, method="spearman"),
                cor(y1, y2, method="kendall") )
  }
  plot(x[,2:3], 
       xlab="Spearman rank correlation",
       ylab = "Kendall's tau",
       main="Rank correlation and Kendall's tau contain the same information")
dev.off()

png(file="g714.png", width=600, height=600)
  N <- 1000
  n <- 100
  x <- matrix(nr=N, nc=3)
  colnames(x) <- c("Pearson", "Spearman (rank)", 
                   "Kendall's tau")
  y1 <- 1:n
  for (i in 1:N) {
    # We only shuffle k elements of the vector
    k <- sample(2:n, 1)  # At least two elements to shuffle
    k <- sample(1:n, k)
    y2 <- y1
    y2[k] <- sample(y2[k])
    # In order to have negative correlations, we also 
    # reverse the vector, from time to time
    if (sample(c(T,F),1)) {
      y2 <- n + 1 - y2
    }
    x[i,] <- c( cor(y1, y2), 
                cor(y1, y2, method="spearman"),
                cor(y1, y2, method="kendall") )
  }
  plot(x[,2:3], 
       # Colour: usual (Pearson) correlation
       col = rainbow(nrow(x))[rank(x[,1])],
       xlab="Spearman rank correlation",
       ylab = "Kendall's tau",
       main="Rank correlation and Kendall's tau contain the same information")
  abline(h=0, v=0, lty=3)
  abline(0, 1, lwd=3)
dev.off()

png(file="g715.png", width=600, height=600)
  gkt <- function (x, y, n=3, ...) {
    q <- quantile(x, c(1/n, 1-1/n), na.rm=T)
    i1 <- which( x <= q[1] & ! is.na(y) )
    i2 <- which( x >= q[2] & ! is.na(y) )
    n <- 0
    for (i in i1) {
      n <- n + sum( y[i] <= y[i2] )
    }
    n <- n / length(i1) / length(i2)
    2 * n - 1
  }
  gkt(1:100, sample(1:100))

  N <- 1000
  n <- 100
  x <- matrix(nr=N, nc=4)
  colnames(x) <- c("Pearson", 
                   "Spearman (rank)", 
                   "Kendall's tau", 
                   "Generalized Kendall's tau")
  y1 <- 1:n
  for (i in 1:N) {
    # We only shuffle k elements of the vector
    k <- sample(2:n, 1)  # At least two elements to shuffle
    k <- sample(1:n, k)
    y2 <- y1
    y2[k] <- sample(y2[k])
    # In order to have negative correlations, we also 
    # reverse the vector, from time to time
    if (sample(c(T,F),1)) {
      y2 <- n + 1 - y2
    }
    x[i,] <- c( cor(y1, y2), 
                cor(y1, y2, method="spearman"),
                cor(y1, y2, method="kendall"),
                gkt(y1, y2) ) 
  }
  pairs(x)
dev.off()

png(file="g716.png", width=600, height=600)
  library(mvtnorm)
  set.seed(1)
  V <- matrix(c(
    1, 0, 0, .5,
    0, 1, .3, 0,
    0, .3, 1, 0,
    .5, 0, 0,  1), nr=4, nc=4)
  stopifnot( eigen(V)$values > 0 )
  n <- 200
  x <- rmvnorm(n, sigma=V)
  colnames(x) <- LETTERS[1:4]

  library(sma)
  plot.cor( cor(x), labels=T )
dev.off()

png(file="g717.png", width=600, height=600)
  RMT <- function (x ,  # One variable per column
                   main="") {
    k <- dim(x)[2]      # Number of variables
    N <- 100            # Number of permutations
    r <- cov(x)
    res0 <- eigen(r)$values
    res <- matrix(NA, nr=N, nc=k)
    for (i in 1:N) {
      y <- apply(x, 2, sample)
      res[i,] <- eigen(cov(y))$values
    }
    if (k>10) {
      res <- res[,1:10]
      res0 <- res0[1:10]
    }
    boxplot(as.vector(res) ~ as.vector(col(res)),
            ylim=range(res0, res, 0),
            col="pink", ylab="eigen values",
            main=main)
    lines(res0, col="blue", lwd=3)
  }

  k <- 10  # Number of variables
  n <- 20  # Number of observations
  x <- matrix(rnorm(n*k), nr=n, nc=k)
  RMT(x, "Independant variables")
dev.off()

png(file="g718.png", width=600, height=600)
  k <- 10  # Number of variables
  n <- 20  # Number of observations for each variable
  m <- 3   # Number of factors
  # Building the variance-covariance matrix
  correlation.matrix <- function (x) {
    # x contains the correlations, one column per factor
    k <- dim(x)[2]  # Number of factors
    n <- dim(x)[1]  # Number of variables
    x %*% t(x)
  }
  covariance.matrix <- function (correlations.with.the.factors, variances) {
    r <- correlations.with.the.factors %*% t(correlations.with.the.factors)
    sqrt(variances %*% t(variances)) * r
  }
  V <- covariance.matrix(matrix(runif(k*m, -1,1), nr=k, nc=m),
                         runif(k, 1,2))
  library(mvtnorm)
  x <- rmvnorm(n, sigma=V)
  RMT(x)
dev.off()

png(file="g719.png", width=600, height=600)
  op <- par(mfrow=c(2,2), mar=c(2,2,3,1))
  for (v in c(.1, .25, .5, 1)) {
    RMT(x + v*rnorm(n*k), main=paste("noise sd =", v))
  }
  par(op)
dev.off()

png(file="g720.png", width=600, height=600)
  N <- 1000
  x <- rt(N, df=1)
  y <- ifelse(sample( c(T,F), N, replace=T ), x, -x) + rt(N,df=1)

  plot(x,y, main="The distributions are too dispersed")
dev.off()

png(file="g721.png", width=600, height=600)
  uniformize <- function (x) {
    x <- rank(x, na.last="keep", ties.method="random")
    x / max(x, na.rm=T)
  }
  x <- uniformize(x)
  y <- uniformize(y)
  plot(x, y, main="After uniformization")
dev.off()

png(file="g722.png", width=600, height=600)
  r <- kde2d(x, y)
  image(r)
  contour(kde2d(x,y), add=T, lwd=3)
dev.off()

png(file="g723.png", width=600, height=600)
  op <- par(mfrow=c(2,2), mar=c(2,3,4,2))
  for (r in c(-.9, -.5, 0, .5)) {
    N <- 1000
    x <- rnorm(N)
    y <- rnorm(N)
    y <- cbind(x,y) %*% chol( matrix(c(1,r,r,1), nr=2) )[,2]  
    cor(x,y)
    x <- uniformize(x)
    y <- uniformize(y)
    s <- kde2d(x, y)
    image(s, main=paste("Correlation =", r))
    contour(s, add=T, lwd=3)
  }
  par(op)  
dev.off()

png(file="g724.png", width=600, height=600)
  do.it <- function (seed, k=3, N=10000) {
    set.seed( seed )
    centers <- matrix(rnorm(2*k), nc=2)
    cluster <- sample(1:k, N, replace=T)
    x <- centers[cluster,1] + rnorm(N)
    y <- centers[cluster,2] + rnorm(N)
    x <- uniformize(x)
    y <- uniformize(y)
    s <- kde2d(x,y)
    image(s)
  }
  do.it(1)
dev.off()

png(file="g725.png", width=600, height=600)
  do.it(2)
dev.off()

png(file="g726.png", width=600, height=600)
  library(bayesm)
  do.it <- function (seed=2, k=3, N=10000) {
    set.seed( seed )
    r <- list()
    for (i in 1:k) {
      m <- matrix(rnorm(4), nr=2)
      r <- append(r, list(list( rnorm(2), solve(chol( t(m) %*% m)))))
    }
    p <- runif(k)
    p <- p / sum(p)
    s <- rmixture(N, p, r)   # Very long...
    op <- par(mfrow=c(2,2), mar=c(2,3,4,2))
    plot(s$x, col=s$z, main="Mixture of gaussians", xlab="", ylab="")
    image(kde2d(s$x[,1], s$x[,2]), main="Density of a mixture of gaussians",
          col=rev(heat.colors(100)))
    box()
    x <- uniformize(s$x[,1])
    y <- uniformize(s$x[,2])
    plot(x, y, col=s$z, main="Uniformized variables", xlab="", ylab="")
    r <- kde2d(x,y)
    image(r, main="Mixture-of-gaussians copula")
    box()
    contour(r, add=T, lwd=3)
    par(op)
  }
  do.it()
dev.off()

png(file="g727.png", width=600, height=600)
  f <- function (u, v, a) {
    exp( -( (-log(u))^a + (-log(v))^a )^(1/a) )
  }
  N <- 50
  v <- u <- ppoints(N)
  uu <- rep(u, N)
  vv <- rep(v, each=N)
  op <- par(mfrow=c(2,2), mar=c(2,2,4,1))
  for (a in c(-10, -2, 0, 5)) {
    w <- matrix( f(uu, vv, a), nr=N )
    image(w, main=paste("Gumbel copula, a =", a))
  }
  par(op)
dev.off()

png(file="g728.png", width=600, height=600)
  my.lss <- function (x, y, ...) {
    n <- length(x)
    sx <- sum(x)
    sy <- sum(y)
    sxx <- sum(x*x)
    sxy <- sum(x*y)
    d <- n*sxx-sx*sx
    a <- (sxx*sy - sx*sxy)/d
    b <- (-sx*sy + n*sxy)/d
    plot(x,y, ...)
    abline(a, b, col='red', ...)
    c(a,b)
  }

  n <- 10
  x <- runif(n)
  y <- 1 - 2*x + .3*rnorm(n)
  my.lss(x,y)
dev.off()

png(file="g729.png", width=600, height=600)
  n <- 10
  x <- runif(n)
  y <- 1 - x + .2*rnorm(n)
  res <- lm( y ~ x )
  plot(y~x, pch=16)
  abline(res, col='red')
dev.off()

png(file="g730.png", width=600, height=600)
  data(cars)
  plot(cars)
  abline(lm(cars$dist ~ cars$speed), col='red')
  title(main="dist ~ speed regression")
dev.off()

png(file="g731.png", width=600, height=600)
  plot(cars)
  r <- lm(cars$dist ~ cars$speed)
  abline(r, col='red')
  r <- lm(cars$speed ~ cars$dist)
  a <- r$coefficients[1] # Intercept
  b <- r$coefficients[2] # slope
  abline(-a/b , 1/b, col="blue")
  title(main="dist ~ speed and speed ~ dist regressions")
dev.off()

png(file="g732.png", width=600, height=600)
  plot(cars)
  r <- lm(cars$dist ~ cars$speed)
  abline(r, col='red')
  segments(cars$speed, cars$dist, cars$speed, r$fitted.values,col="red")
  title(main="dist ~ speed: distances measured vertically")
dev.off()

png(file="g733.png", width=600, height=600)
  plot(cars)
  r <- lm(cars$speed ~ cars$dist)
  a <- r$coefficients[1] # Intercept
  b <- r$coefficients[2] # slope
  abline(-a/b , 1/b, col="blue")
  segments(cars$speed, cars$dist, r$fitted.values, cars$dist, col="blue")
  title(main="speed ~ dist: distances measured horizontally")
dev.off()

png(file="g734.png", width=600, height=600)
  plot(cars)
  r <- lm(cars$dist ~ cars$speed)
  abline(r, col='red')
  r <- lm(cars$speed ~ cars$dist)
  a <- r$coefficients[1] # Intercept
  b <- r$coefficients[2] # slope
  abline(-a/b , 1/b, col="blue")
  r <- princomp(cars)
  b <- r$loadings[2,1] / r$loadings[1,1]
  a <- r$center[2] - b * r$center[1]
  abline(a,b)
  title(main='Comparing three "regressions"')
dev.off()

png(file="g735.png", width=600, height=600)
  set.seed(1)
  x <- rnorm(100)
  y <- x + rnorm(100)
  plot(y~x)
  r <- lm(y~x)
  abline(r, col='red')
  r <- lm(x ~ y)
  a <- r$coefficients[1] # Intercept
  b <- r$coefficients[2] # slope
  abline(-a/b , 1/b, col="blue")
  r <- princomp(cbind(x,y))
  b <- r$loadings[2,1] / r$loadings[1,1]
  a <- r$center[2] - b * r$center[1]
  abline(a,b)
  title(main='Comparing three "regressions"')
dev.off()

png(file="g736.png", width=600, height=600)
  plot(y~x, xlim=c(-4,4), ylim=c(-4,4) )
  abline(a,b)
  # Change-of-base matrix
  u <- r$loadings
  # Projection onto the first axis
  p <- matrix( c(1,0,0,0), nrow=2 )
  X <- rbind(x,y)
  X <- r$center + solve(u, p %*% u %*% (X - r$center))
  segments( x, y, X[1,], X[2,] )
  title(main="PCA: distances measured perpendicularly to the line")
dev.off()

png(file="g737.png", width=600, height=600)
  library(nlme)
  data(Orthodont)
  fm1 <- lme(distance ~ age, Orthodont, random = ~ age | Subject)
  # standardized residuals versus fitted values by gender
  plot(fm1, resid(., type = "p") ~ fitted(.) | Sex, abline = 0)
dev.off()

png(file="g738.png", width=600, height=600)
  # box-plots of residuals by Subject
  plot(fm1, Subject ~ resid(.))
dev.off()

png(file="g739.png", width=600, height=600)
  # observed versus fitted values by Subject
  plot(fm1, distance ~ fitted(.) | Subject, abline = c(0,1))
dev.off()

png(file="g740.png", width=600, height=600)
  data(crabs)
  n <- length(crabs$RW)
  r <- lm(FL~RW, data=crabs)
  op <- par(mfrow=c(2,2))
  plot(r)
  par(op)
dev.off()

png(file="g741.png", width=600, height=600)
  data(LifeCycleSavings)
  plot(LifeCycleSavings)
dev.off()

png(file="g742.png", width=600, height=600)
  op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  plot(lm.SR <- lm(sr ~ pop15 + pop75 + dpi + ddpi, data = LifeCycleSavings))
  par(op)
dev.off()

png(file="g743.png", width=600, height=600)
  ## 4 plots on 1 page; allow room for printing model formula in outer margin:
  op <- par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  #plot(lm.SR)
  #plot(lm.SR, id.n = NULL)               # no id's
  #plot(lm.SR, id.n = 5, labels.id = NULL)# 5 id numbers
  plot(lm.SR, panel = panel.smooth)
  ## Gives a smoother curve
  #plot(lm.SR, panel = function(x,y) panel.smooth(x, y, span = 1))
  par(op)
dev.off()

png(file="g744.png", width=600, height=600)
  n <- 100
  x <- abs(rnorm(n)) ^ runif(1,min=0,max=2)
  kk <- seq(.01,5,by=.01)
  N <- length(kk)
  pv <- rep(NA, N)
  for (i in 1:N) {
    pv[i] <- shapiro.test( x^kk[i] )$p.value
  }
  plot( pv ~ kk, type='l', xlab='exponent', ylab='p-value' )
  seuil <- .05
  abline( v = kk[ range( (1:N)[pv>seuil] ) ], lty=3 )  
  abline( v = kk[pv==max(pv)], lty=3 )
  abline( h = seuil, lty=3 )
  title(main="Shapiro's test")
dev.off()

png(file="g745.png", width=600, height=600)
  x <- exp( rnorm(100) )
  y <- 1 + 2*x + .3*rnorm(100)
  y <- y^1.3
  library(MASS)
  boxcox(y~x, plotit=T)
dev.off()

png(file="g746.png", width=600, height=600)
  bc <- boxcox(y ~ x, plotit=F)
  a <- bc$x[ order(bc$y, decreasing=T)[1] ]
  op <- par( mfcol=c(1,2) )
  plot( y~x )
  plot( y^a ~ x )
  par(op)
dev.off()

png(file="g747.png", width=600, height=600)
  data(trees)
  plot(trees)
dev.off()

png(file="g748.png", width=600, height=600)
  boxcox(Volume ~ Girth, data = trees)
dev.off()

png(file="g749.png", width=600, height=600)
  boxcox(Volume ~ log(Height) + log(Girth), data = trees )
dev.off()

png(file="g750.png", width=600, height=600)
  n <- 100
  v <- .1
  x1 <- rlnorm(n)
  x2 <- rlnorm(n)
  x3 <- rlnorm(n)
  x4 <- x1 + x2 + x3 + v*rlnorm(n)
  m1 <- cbind(x1,x2,x3,x4)
  pairs(m1, main="No missing value")
dev.off()

png(file="g751.png", width=600, height=600)
  remove.higher.values <- function (x) {
    n <- length(x)
    ifelse( rbinom(n,1,(x-min(x))/(max(x)+1))==1 , NA, x)
  }
  x1 <- remove.higher.values(x1)
  x2 <- remove.higher.values(x2)
  x3 <- remove.higher.values(x3)
  x4 <- remove.higher.values(x4)
  m2 <- cbind(x1,x2,x3,x4)
  pairs(m2, main="A few missing values")
dev.off()

png(file="g752.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  library(Design)
  r <- naclus(data.frame(m2))
  naplot(r)
  par(op)
dev.off()

png(file="g753.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for(m in c("ward","complete","median")) {
    plot(naclus(data.frame(m2), method=m))
    title(m)
  }
  plot(naclus(data.frame(m2)))
  title("Default")
  par(op)
dev.off()

png(file="g754.png", width=600, height=600)
  library(Hmisc)
  op <- par(mfrow=c(2,2))
  w <- transcan(~x1+x2+x3+x4, imputed=T, transformed=T, trantab=T, impcat='tree',
                data=data.frame(x1,x2,x3,x4), pl=TRUE)
  par(op)
dev.off()

png(file="g755.png", width=600, height=600)
  library(acepack)
  TWOPI <- 8*atan(1)
  x <- runif(200,0,TWOPI)
  y <- exp(sin(x)+rnorm(200)/2)
  a <- ace(x,y)
  op <- par(mfrow=c(3,1))
  plot(a$y,a$ty)  # view the response transformation
  plot(a$x,a$tx)  # view the carrier transformation
  plot(a$tx,a$ty) # examine the linearity of the fitted model
  par(op)
dev.off()

png(file="g756.png", width=600, height=600)
  library(acepack)
  TWOPI <- 8*atan(1)
  x <- runif(200,0,TWOPI)
  y <- exp(sin(x)+rnorm(200)/2)
  a <- avas(x,y)
  op <- par(mfrow=c(3,1))
  plot(a$y,a$ty)  # view the response transformation
  plot(a$x,a$tx)  # view the carrier transformation
  plot(a$tx,a$ty) # examine the linearity of the fitted model
  par(op)
dev.off()

png(file="g757.png", width=600, height=800)
  m <- read.table("faraway")
  op <- par(mfrow=(c(6,1)))
  for (i in 1:6) {
    boxplot(m[,i], horizontal=T, main=names(m)[i])
  }
  par(op)
dev.off()

png(file="g758.png", width=600, height=800)
  m <- read.table("faraway")
  op <- par(mfrow=(c(3,2)))
  for (i in 1:6) {
    hist(m[,i], col='light blue', probability=T, xlab=names(m)[i])
    lines(density(m[,i]), col='red', lwd=2)
  }
  par(op)
dev.off()

png(file="g759.png", width=600, height=800)
  m <- read.table("faraway")
  op <- par(mfrow=(c(3,2)))
  for (i in 1:6) {
    qqnorm(m[,i], main=names(m)[i])
    qqline(m[,i], col='red')
  }
  par(op)
dev.off()
detach.everything()

png(file="g760.png", width=600, height=600)
  y <- cars$dist
  x <- cars$speed
  o = order(x)
  plot( y~x )
  do.it <- function (model, col) {
    r <- lm( model )
    yp <- predict(r)
    lines( yp[o] ~ x[o], col=col, lwd=3 )
  }
  do.it(y~x, col="red")
  do.it(y~x+I(x^2), col="blue")
  do.it(y~-1+I(x^2), col="green")
  legend(par("usr")[1], par("usr")[4], 
         c("affine function", "degree-2 polynomial", "degree 2 monomial"),
         lwd=3,
         col=c("red", "blue", "green"),
        )
dev.off()

png(file="g761.png", width=600, height=600)
  n <- 5
  p <- matrix( nrow=n, ncol=n+1 )
  for (i in 1:n) {
    p[i,1:(i+1)] <- summary(lm( y ~ poly(x,i) ))$coefficients[,4]
  }
  matplot(p, type='l', lty=1, lwd=3)
  legend( par("usr")[1], par("usr")[4],
          as.character(1:n),
          lwd=3, lty=1, col=1:n
        )
  title(main="Evolution of the p-values (orthonormal polynomials)")
dev.off()

png(file="g762.png", width=600, height=600)
  p <- matrix( nrow=n, ncol=n+1 )
  p[1,1:2] <- summary(lm(y ~ x) )$coefficients[,4]
  p[2,1:3] <- summary(lm(y ~ x+I(x^2)) )$coefficients[,4]
  p[3,1:4] <- summary(lm(y ~ x+I(x^2)+I(x^3)) )$coefficients[,4]
  p[4,1:5] <- summary(lm(y ~ x+I(x^2)+I(x^3)+I(x^4)) )$coefficients[,4]
  p[5,1:6] <- summary(lm(y ~ x+I(x^2)+I(x^3)+I(x^4)+I(x^5)) )$coefficients[,4]
  matplot(p, type='l', lty=1, lwd=3)
  legend( par("usr")[1], par("usr")[4],
          as.character(1:n),
          lwd=3, lty=1, col=1:n
        )
  title(main="Evolution of the p-values (non orthonormal polynomials)")
dev.off()

png(file="g763.png", width=600, height=600)
  # The matrix
  n <- 5
  m <- matrix( ncol=n+1, nrow=length(x) )
  for (i in 2:(n+1)) {
    m[,i] <- x^(i-1)
  }
  m[,1] <- 1
  # Orthonormalization (Gram--Schmidt)
  for (i in 1:(n+1)) {
    if(i==1) m[,1] <- m[,1] / sqrt( t(m[,1]) %*% m[,1] )
    else {
      for (j in 1:(i-1)) {
        m[,i] <- m[,i] - (t(m[,i]) %*% m[,j])*m[,j]
      }
      m[,i] <- m[,i] / sqrt( t(m[,i]) %*% m[,i] )
    }
  }
  p <- matrix( nrow=n, ncol=n+1 )
  p[1,1:2] <- summary(lm(y~ -1 +m[,1:2] ))$coefficients[,4]
  p[2,1:3] <- summary(lm(y~ -1 +m[,1:3] ))$coefficients[,4]
  p[3,1:4] <- summary(lm(y~ -1 +m[,1:4] ))$coefficients[,4]
  p[4,1:5] <- summary(lm(y~ -1 +m[,1:5] ))$coefficients[,4]
  p[5,1:6] <- summary(lm(y~ -1 +m[,1:6] ))$coefficients[,4]
  matplot(p, type='l', lty=1, lwd=3)
  legend( par("usr")[1], par("usr")[4],
          as.character(1:n),
          lwd=3, lty=1, col=1:n
        )
  title(main="Idem, orthonormalisation by hand")
dev.off()

png(file="g764.png", width=600, height=600)
  library(ts)
  data(beavers)
  y <- beaver1$temp
  x <- 1:length(y)
  plot(y~x)
  for (i in 1:10) {
    r <- lm( y ~ poly(x,i) )
    lines( predict(r), type="l", col=i )
  }
  summary(r)
dev.off()

png(file="g765.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- exp(x)
  plot(y~x)
  title(main="Non-polynomial relation")
dev.off()

png(file="g766.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- exp(x) + .1*rnorm(n)
  plot(y~x)
  title(main="Non-polynomial relation")
dev.off()

png(file="g767.png", width=600, height=600)
  library(modreg)
  plot(cars$speed, cars$dist)
  lines( smooth.spline(cars$speed, cars$dist), col='red')
  abline(lm(dist~speed,data=cars), col='blue', lty=2)
dev.off()

png(file="g768.png", width=600, height=600)
  plot(quakes$long, quakes$lat)
  lines( smooth.spline(quakes$long, quakes$lat), col='red', lwd=3)
dev.off()

png(file="g769.png", width=600, height=600)
  library(Design)
  # 4-node spline
  r3 <- lm( quakes$lat ~ rcs(quakes$long) )
  plot( quakes$lat ~ quakes$long )
  o <- order(quakes$long)
  lines( quakes$long[o], predict(r)[o], col='red', lwd=3 )
  r <- lm( quakes$lat ~ rcs(quakes$long,10) )
  lines( quakes$long[o], predict(r)[o], col='blue', lwd=6, lty=3 )
  title(main="Regression with rcs")
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("4 knots", "10 knots"),
          lwd=c(3,3), lty=c(1,3), col=c("red", "blue") )
dev.off()

png(file="g770.png", width=600, height=600)
  library(splines)
  data(quakes)
  x <- quakes[,2]
  y <- quakes[,1]
  o <- order(x)
  x <- x[o]
  y <- y[o]
  r1 <- lm( y ~ bs(x,df=10) )
  r2 <- lm( y ~ ns(x,df=6) )
  plot(y~x)
  lines(predict(r1)~x, col='red', lwd=3)
  lines(predict(r2)~x, col='green', lwd=3)
dev.off()

png(file="g771.png", width=600, height=600)
  # The manual asks us to be cautious with predictions 
  # for ither values of x: I do not see any problem.
  plot(y~x)
  xp <- seq(min(x), max(x), length=200)
  lines(predict(r2) ~ x, col='green', lwd=3)
  lines(xp, predict(r2,data.frame(x=xp)), col='blue', lwd=7, lty=3)
dev.off()

png(file="g772.png", width=600, height=600)
  n <- 100
  x <- seq(-2,2,length=n)
  y <- exp(x)
  plot(y~x, type='l', lwd=3)
  title(main='Expenential growth')
dev.off()

png(file="g773.png", width=600, height=600)
  x <- seq(-2,2,length=n)
  y <- exp(-x)
  plot(y~x, type='l', lwd=3)
  title(main='Exponential Decrease')
dev.off()

png(file="g774.png", width=600, height=600)
  x <- seq(-2,4,length=n)
  y <- 1-exp(-x)
  plot(y~x, type='l', lwd=3)
  title(main='Negative Exponential')
dev.off()

png(file="g775.png", width=600, height=600)
  x <- seq(0,5,length=n)
  u <- 1
  v <- 2
  y <- u/(u-v) * (exp(-v*x) - exp(-u*x))
  plot(y~x, type='l', lwd=3)
  title(main='Double Exponential')
dev.off()

png(file="g776.png", width=600, height=600)
  x <- seq(-5,5,length=n)
  y <- 1/(1+exp(-x))
  plot(y~x, type='l', lwd=3)
  title(main='Sigmoid Growth')
dev.off()

png(file="g777.png", width=600, height=600)
  x <- seq(-2,5,length=n)
  y <- exp(-exp(-x))
  plot(y~x, type='l', lwd=3)
  title(main='Less Symetric Sigmoid')
dev.off()

png(file="g778.png", width=600, height=600)
  library(nls)
  f <- function (x,p) {
    u <- p[1]
    v <- p[2]
    u/(u-v) * (exp(-v*x) - exp(-u*x))    
  }
  n <- 100
  x <- runif(n,0,2)
  y <- f(x, c(3.14,2.71)) + .1*rnorm(n)
  r <- nls( y ~ f(x,c(a,b)), start=c(a=3, b=2.5) )
  plot(y~x)
  xx <- seq(0,2,length=200)
  lines(xx, f(xx,r$m$getAllPars()), col='red', lwd=3)
  lines(xx, f(xx,c(3.14,2.71)), lty=2)
dev.off()

png(file="g779.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  p <- profile(r)
  plot(p, conf = c(95, 90, 80, 50)/100)
  plot(p, conf = c(95, 90, 80, 50)/100, absVal = FALSE)
  par(op)
dev.off()

png(file="g780.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,0,2)
    y <- SSbiexp(x,1,1,-1,2) + .01*rnorm(n)
    try(  r <- nls(y ~ SSbiexp(x,a,u,b,v))  )
  }
  plot(y~x)
  xx <- seq(0,2,length=200)
  lines(xx, SSbiexp(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2],
                    r$m$getAllPars()[3],
                    r$m$getAllPars()[4]
                    ), col='red', lwd=3)
  title(main='biexponential')
dev.off()

png(file="g781.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  try(plot(profile(r)))
  par(op)
dev.off()

png(file="g782.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,-5,5)
    y <- SSlogis(x,1,0,1) + .1*rnorm(n)
    try( r <- nls(y ~ SSlogis(x,a,m,s)) )
  }
  plot(y~x)
  xx <- seq(-5,5,length=200)
  lines(xx, SSlogis(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2],
                    r$m$getAllPars()[3]
                    ), col='red', lwd=3)
  title(main='logistic')
dev.off()

png(file="g783.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,-5,5)
    y <- SSfpl(x,-1,1,0,1) + .1*rnorm(n)
    try( r <- nls(y ~ SSfpl(x,a,b,m,s)) )
  }
  plot(y~x)
  xx <- seq(-5,5,length=200)
  lines(xx, SSfpl(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2],
                    r$m$getAllPars()[3],
                    r$m$getAllPars()[4]
                    ), col='red', lwd=3)
  title(main='4-parameter logistic model')
dev.off()

png(file="g784.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,-5,5)
    y <- SSmicmen(x,1,1) + .01*rnorm(n)
    try( r <- nls(y ~ SSmicmen(x,m,h)) )
  }
  plot(y~x, ylim=c(-5,5))
  xx <- seq(-5,5,length=200)
  lines(xx, SSmicmen(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2]
                    ), col='red', lwd=3)
  title(main='Michaelis-Menten')
dev.off()

png(file="g785.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,0,5)
    y <- SSmicmen(x,1,1) + .1*rnorm(n)
    try( r <- nls(y ~ SSmicmen(x,m,h)) )
  }
  plot(y~x)
  xx <- seq(0,5,length=200)
  lines(xx, SSmicmen(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2]
                    ), col='red', lwd=3)
  title(main='Michaelis-Menten')
dev.off()

png(file="g786.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,0,1)
    y <- SSfol(1,x,1,2,1) + .1*rnorm(n)
    try( r <- nls(y ~ SSfol(1,x,a,b,c)) )
  }
  plot(y~x)
  xx <- seq(0,1,length=200)
  lines(xx, SSfol(1,
                     xx,
                     r$m$getAllPars()[1],
                     r$m$getAllPars()[2],
                     r$m$getAllPars()[3]
                    ), col='red', lwd=3)
  title(main='SSfol')
dev.off()

png(file="g787.png", width=600, height=600)
  rm(r)
  while(!exists("r")) {
    x <- runif(n,0,2)
    y <- SSasymp(x,1,.5,1) + .1*rnorm(n)
    try( r <- nls(y ~ SSasymp(x,a,b,c)) )
  }
  plot(y~x, xlim=c(-.5,2),ylim=c(0,1.3))
  xx <- seq(-1,2,length=200)
  lines(xx, SSasymp(xx,
                    r$m$getAllPars()[1],
                    r$m$getAllPars()[2],
                    r$m$getAllPars()[3]
                   ), col='red', lwd=3)
  title(main='SSasymp')
  # See also SSasympOff et SSasympOrig
dev.off()

png(file="g788.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(0, 5, len = 101)
  yy <- 5 - 4 * exp(-xx/(2*log(2)))
  par(mar = c(0, 0, 4.1, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSasymp model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(-0.4, 1, 0, 1, lty = 2, lwd = 0.75)
  arrows(-0.3, 0.25, -0.3, 0, length = 0.07, angle = 25)
  arrows(-0.3, 0.75, -0.3, 1, length = 0.07, angle = 25)
  text(-0.3, 0.5, expression(phi[2]), adj = c(0.5, 0.5))
  segments(1, 3.025, 1, 4, lty = 2, lwd = 0.75)
  arrows(0.2, 3.5, 0, 3.5, length = 0.08, angle = 25)
  arrows(0.8, 3.5, 1, 3.5, length = 0.08, angle = 25)
  text(0.5, 3.5, expression(t[0.5]), adj = c(0.5, 0.5))
dev.off()

png(file="g789.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(0.5, 5, len = 101)
  yy <- 5 * (1 -  exp(-(xx - 0.5)/(2*log(2))))
  par(mar = c(0, 0, 4.0, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSasympOff model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(0.5, 0, 0.5, 3, lty = 2, lwd = 0.75)
  text(0.5, 3.1, expression(phi[3]), adj = c(0.5, 0))
  segments(1.5, 2.525, 1.5, 3, lty = 2, lwd = 0.75)
  arrows(0.7, 2.65, 0.5, 2.65, length = 0.08, angle = 25)
  arrows(1.3, 2.65, 1.5, 2.65, length = 0.08, angle = 25)
  text(1.0, 2.65, expression(t[0.5]), adj = c(0.5, 0.5))
dev.off()

png(file="g790.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(0, 5, len = 101)
  yy <- 5 * (1- exp(-xx/(2*log(2))))
  par(mar = c(0, 0, 3.5, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSasympOrig model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(1, 2.525, 1, 3.5, lty = 2, lwd = 0.75)
  arrows(0.2, 3.0, 0, 3.0, length = 0.08, angle = 25)
  arrows(0.8, 3.0, 1, 3.0, length = 0.08, angle = 25)
  text(0.5, 3.0, expression(t[0.5]), adj = c(0.5, 0.5))
dev.off()

png(file="g791.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(-0.5, 5, len = 101)
  yy <- 1 + 4 / ( 1 + exp((2-xx)))
  par(mar = c(0, 0, 3.5, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSfpl model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  abline(h = 1, lty = 2, lwd = 0)
  arrows(-0.3, 0.25, -0.3, 0, length = 0.07, angle = 25)
  arrows(-0.3, 0.75, -0.3, 1, length = 0.07, angle = 25)
  text(-0.3, 0.5, expression(phi[2]), adj = c(0.5, 0.5))
  segments(2, 0, 2, 3.3, lty = 2, lwd = 0.75)
  text(2, 3.3, expression(phi[3]), adj = c(0.5, 0))
  segments(3, 1+4/(1+exp(-1)) - 0.025, 3, 2.5, lty = 2, lwd = 0.75)
  arrows(2.3, 2.7, 2.0, 2.7, length = 0.08, angle = 25)
  arrows(2.7, 2.7, 3.0, 2.7, length = 0.08, angle = 25)
  text(2.5, 2.7, expression(phi[4]), adj = c(0.5, 0.5))
dev.off()

png(file="g792.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(-0.5, 5, len = 101)
  yy <- 5 / ( 1 + exp((2-xx)))
  par(mar = c(0, 0, 3.5, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSlogis model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(2, 0, 2, 4.0, lty = 2, lwd = 0.75)
  text(2, 4.0, expression(phi[2]), adj = c(0.5, 0))
  segments(3, 5/(1+exp(-1)) + 0.025, 3, 4.0, lty = 2, lwd = 0.75)
  arrows(2.3, 3.8, 2.0, 3.8, length = 0.08, angle = 25)
  arrows(2.7, 3.8, 3.0, 3.8, length = 0.08, angle = 25)
  text(2.5, 3.8, expression(phi[3]), adj = c(0.5, 0.5))
dev.off()

png(file="g793.png", width=600, height=600)
  # Copied from the manual
  xx <- seq(0, 5, len = 101)
  yy <- 5 * xx/(1+xx)
  par(mar = c(0, 0, 3.5, 0))
  plot(xx, yy, type = "l", axes = FALSE, ylim = c(0,6), xlim = c(-1, 5),
       xlab = "", ylab = "", lwd = 2,
       main = "Parameters in the SSmicmen model")
  usr <- par("usr")
  arrows(usr[1], 0, usr[2], 0, length = 0.1, angle = 25)
  arrows(0, usr[3], 0, usr[4], length = 0.1, angle = 25)
  text(usr[2] - 0.2, 0.1, "x", adj = c(1, 0))
  text(-0.1, usr[4], "y", adj = c(1, 1))
  abline(h = 5, lty = 2, lwd = 0)
  arrows(-0.8, 2.1, -0.8, 0, length = 0.1, angle = 25)
  arrows(-0.8, 2.9, -0.8, 5, length = 0.1, angle = 25)
  text(-0.8, 2.5, expression(phi[1]), adj = c(0.5, 0.5))
  segments(1, 0, 1, 2.7, lty = 2, lwd = 0.75)
  text(1, 2.7, expression(phi[2]), adj = c(0.5, 0))
dev.off()

png(file="g794.png", width=600, height=600)
  broken.line <- function (x, y) { 
    n <- length(x)
    n2 = floor(n/2)
    o <- order(x)
    m <- mean(c(x[o[n2+0:1]]))
    x1 <- x[o[1:n2]]
    y1 <- y[o[1:n2]]
    n2 <- n2+1
    x2 <- x[o[n2:n]]
    y2 <- y[o[n2:n]]
    r1 <- lm(y1~x1)
    r2 <- lm(y2~x2)
    plot(y~x)
    segments(x[o[1]], r1$coef[1] + x[o[1]]*r1$coef[2],
             m, r1$coef[1] + m *r1$coef[2],
             col='red')
    segments(m, r2$coef[1] + m*r2$coef[2],
             x[o[n]], r2$coef[1] + x[o[n]] *r2$coef[2],
             col='blue')
    abline(v=m, lty=3)
  }
  set.seed(1)
  n <- 10
  x <- runif(n)
  y <- 1-x+.2*rnorm(n)
  broken.line(x,y)
dev.off()

png(file="g795.png", width=600, height=600)
  z <- x*(1-x)
  broken.line(x,z)
dev.off()

png(file="g796.png", width=600, height=600)
  broken.line <- function (x,y) {
    n <- length(x)
    o <- order(x)
    n1 <- floor((n+1)/2)
    n2 <- ceiling((n+1)/2)
    m <- mean(c( x[o[n1]], x[o[n2]] ))
    plot(y~x)
    B1 <- function (xx) {
      x <- xx
      x[xx<m] <- m - x[xx<m]
      x[xx>=m] <- 0
      x
    }
    B2 <- function (xx) {
      x <- xx
      x[xx>m] <- x[x>m] -m
      x[xx<=m] <- 0
      x
    }
    x1 <- B1(x)
    x2 <- B2(x)
    r <- lm(y~x1+x2)
    xx <- seq(x[o[1]],x[o[n]],length=100)
    yy <- predict(r, data.frame(x1=B1(xx), x2=B2(xx)))
    lines( xx, yy, col='red' )
    abline(v=m, lty=3)        
  }
  broken.line(x,y)
dev.off()

png(file="g797.png", width=600, height=600)
  broken.line(x,z)
dev.off()

png(file="g798.png", width=600, height=600)
  set.seed(5)
  n <- 200
  x <- rnorm(n)
  k <- rnorm(1)
  x1 <- ifelse( x < k, x - k, 0 )
  x2 <- ifelse( x < k, 0, x - k )
  a <- rnorm(1)
  b1 <- rnorm(1)
  b2 <- rnorm(1)
  y <- a + b1 * x1 + b2 * x2 + .2*rnorm(n)
  plot( y ~ x, col = ifelse(x < k, "blue", "green") )
  abline(a - k*b1, b1, col="blue", lwd=2, lty=2)
  abline(a - k*b2, b2, col="green", lwd=2, lty=2)
dev.off()

png(file="g799.png", width=600, height=600)
  library(segmented)
  r <- segmented( lm(y~x), 
                  x,        # Variable along which we are allowed to cut
                  0         # Initial values of the cut-points 
                )           # (there can be several)
  plot(x, y)
  o <- order(x)
  lines( x[o], r$fitted.values[o], col="red", lwd=3 )
dev.off()

png(file="g800.png", width=600, height=600)
  set.seed(1)
  n <- 10
  x <- runif(n)
  y <- 1-x+.2*rnorm(n)
  z <- x*(1-x)  # Same data as above
  plot(y~x)
  lines(lowess(x,y), col="red")
dev.off()

png(file="g801.png", width=600, height=600)
  plot(z~x)
  lines(lowess(x,z), col="red")
dev.off()

png(file="g802.png", width=600, height=600)
  data(quakes)
  plot(lat~long, data=quakes)
  lines(lowess(quakes$long, quakes$lat), col='red', lwd=3)
dev.off()

png(file="g803.png", width=600, height=600)
  n <- 1000
  x <- runif(n, min=-1,max=1)
  y <- (x-1)*x*(x+1) + .5*rnorm(n)
  # Do not use this code for large data sets.
  moyenne.mobile <- function (y,x, r=.1) {
    o <- order(x)
    n <- length(x)
    m <- floor((1-r)*n)
    p <- n-m
    x1 <- vector(mode='numeric',length=m)
    y1 <- vector(mode='numeric',length=m)
    y2 <- vector(mode='numeric',length=m)
    y3 <- vector(mode='numeric',length=m)
    for (i in 1:m) {
      xx <- x[ o[i:(i+p)] ]
      yy <- y[ o[i:(i+p)] ]
      x1[i] <- mean(xx)
      y1[i] <- quantile(yy,.25)
      y2[i] <- quantile(yy,.5)
      y3[i] <- quantile(yy,.75)
    }
    plot(y~x)
    lines(x1,y2,col='red', lwd=3)
    lines(x1,y1,col='blue', lwd=2)
    lines(x1,y3,col='blue', lwd=2)
  }
  moyenne.mobile(y,x)
dev.off()

png(file="g804.png", width=600, height=600)
  library(gregmisc)
  bandplot(x,y)
dev.off()

png(file="g805.png", width=600, height=600)
  library(modreg)
  plot(cars$speed, cars$dist)
  lines(ksmooth(cars$speed, cars$dist, "normal", bandwidth=2), col='red')
  lines(ksmooth(cars$speed, cars$dist, "normal", bandwidth=5), col='green')
  lines(ksmooth(cars$speed, cars$dist, "normal", bandwidth=10), col='blue')
dev.off()

png(file="g806.png", width=600, height=600)
  curve(dnorm(x), xlim=c(-3,3), ylim=c(0,1.1))
  x <- seq(-3,3,length=200)
  D.Epanechikov <- function (t) {
    ifelse(abs(t)<1, 3/4*(1-t^2), 0)
  }
  lines(D.Epanechikov(x) ~ x, col='red')
  D.tricube <- function (t) { # aka "triweight kernel"
    ifelse(abs(t)<1, (1-abs(t)^3)^3, 0)
  }
  lines(D.tricube(x) ~ x, col='blue')
  legend( par("usr")[1], par("usr")[4], yjust=1,
          c("noyau gaussien", "noyau d'Epanechikov", "noyau tricube"),
          lwd=1, lty=1,
          col=c(par('fg'),'red', 'blue'))
  title(main="Differents kernels")
dev.off()

png(file="g807.png", width=600, height=600)
  # With real data...
  library(KernSmooth)
  data(quakes)
  x <- quakes$long
  y <- quakes$lat
  plot(y~x)
  bw <- dpill(x,y) # .2
  lines( locpoly(x,y,degree=0, bandwidth=bw), col='red' )
  lines( locpoly(x,y,degree=1, bandwidth=bw), col='green' )
  lines( locpoly(x,y,degree=2, bandwidth=bw), col='blue' )
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="Local Polynomial Regression")
dev.off()

png(file="g808.png", width=600, height=600)
  plot(y~x)
  bw <- .5
  lines( locpoly(x,y,degree=0, bandwidth=bw), col='red' )
  lines( locpoly(x,y,degree=1, bandwidth=bw), col='green' )
  lines( locpoly(x,y,degree=2, bandwidth=bw), col='blue' )
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="Local Polynomial Regression (wider window)")
dev.off()

png(file="g809.png", width=600, height=600)
  n <- 50
  x <- runif(n)
  f <- function (x) { cos(3*x) + cos(5*x) }
  y <- f(x) + .2*rnorm(n)
  plot(y~x)
  curve(f(x), add=T, lty=2)
  bw <- dpill(x,y)
  lines( locpoly(x,y,degree=0, bandwidth=bw), col='red' )
  lines( locpoly(x,y,degree=1, bandwidth=bw), col='green' )
  lines( locpoly(x,y,degree=2, bandwidth=bw), col='blue' )
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
dev.off()

png(file="g810.png", width=600, height=600)
  a <- locpoly(x,y,degree=0, bandwidth=bw)
  b <- locpoly(x,y,degree=1, bandwidth=bw)
  c <- locpoly(x,y,degree=2, bandwidth=bw)
  matplot( cbind(a$x,b$x,c$x), abs(cbind(a$y-f(a$x), b$y-f(b$x), c$y-f(c$x)))^2, 
           xlab='', ylab='',
           type='l', lty=1, col=c('red', 'green', 'blue') )
  legend( .8*par("usr")[1]+.2*par("usr")[2], par("usr")[4], yjust=1,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="MSE (Mean Square Error)")
dev.off()

png(file="g811.png", width=600, height=600)
  f <- function (x) { sqrt(abs(x-.5)) }
  y <- f(x) + .1*rnorm(n)
  plot(y~x)
  curve(f(x), add=T, lty=2)
  bw <- dpill(x,y)
  lines( locpoly(x,y,degree=0, bandwidth=bw), col='red' )
  lines( locpoly(x,y,degree=1, bandwidth=bw), col='green' )
  lines( locpoly(x,y,degree=2, bandwidth=bw), col='blue' )
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
dev.off()

png(file="g812.png", width=600, height=600)
  a <- locpoly(x,y,degree=0, bandwidth=bw)
  b <- locpoly(x,y,degree=1, bandwidth=bw)
  c <- locpoly(x,y,degree=2, bandwidth=bw)
  matplot( cbind(a$x,b$x,c$x), abs(cbind(a$y-f(a$x), b$y-f(b$x), c$y-f(c$x)))^2, 
           xlab='', ylab='',
           type='l', lty=1, col=c('red', 'green', 'blue') )
  legend( .8*par("usr")[1]+.2*par("usr")[2], par("usr")[4], yjust=1,
          c("degree = 0", "degree = 1", "degree = 2"),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="MSE (Mean Square Error)")
dev.off()

png(file="g813.png", width=600, height=600)
  # It is linear
  library(modreg)
  n <- 10
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    x <- rnorm(n)
    y <- 1-2*x+.3*rnorm(n)
    plot(y~x)
    lo <- loess(y~x)
    xx <- seq(min(x),max(x),length=100)
    yy <- predict(lo, data.frame(x=xx))
    lines(xx,yy, col='red')
    lo <- loess(y~x, family='sym')
    xx <- seq(min(x),max(x),length=100)
    yy <- predict(lo, data.frame(x=xx))
    lines(xx,yy, col='red', lty=2)
    lines(lowess(x,y),col='blue',lty=2)
  #  abline(lm(y~x), col='green')
  #  abline(1,-2, col='green', lty=2)
  }
  par(op)
dev.off()

png(file="g814.png", width=600, height=600)
  # It is not linear
  n <- 10
  op <- par(mfrow=c(2,2))
  for (i in 1:4) {
    x <- rnorm(n)
    y <- x*(1-x)+.3*rnorm(n)
    plot(y~x)
    lo <- loess(y~x)
    xx <- seq(min(x),max(x),length=100)
    yy <- predict(lo, data.frame(x=xx))
    lines(xx,yy, col='red')
    lo <- loess(y~x, family='sym')
    xx <- seq(min(x),max(x),length=100)
    yy <- predict(lo, data.frame(x=xx))
    lines(xx,yy, col='red', lty=2)
    lines(lowess(x,y),col='blue',lty=2)
  #  curve(x*(1-x), add = TRUE, col = "green", lty=2)
  }
  par(op)
dev.off()

png(file="g815.png", width=600, height=600)
  x <- cars$speed
  y <- cars$dist
  my.lar <- function (y,x) {
    f <- function (arg) {
      a <- arg[1]
      b <- arg[2]
      sum(abs(y-a-b*x)) 
    }
    r <- optim( c(0,0), f )$par
    plot( y~x )
    abline(lm(y~x), col='red', lty=2)
    abline(r[1], r[2])
    legend( par("usr")[1], par("usr")[4], yjust=1,
            c("Least Squares", "Least Absolute Values"),
            lwd=1, lty=c(2,1),
            col=c(par('fg'),'red'))
  }
  my.lar(y,x)
dev.off()

png(file="g816.png", width=600, height=600)
  library(MASS)
  n <- 20
  x <- rnorm(n)
  y <- 1 - 2*x + rnorm(n)
  y[ sample(1:n, floor(n/4)) ] <- 10
  plot(y~x)
  abline(1,-2,lty=3)
  abline(lm(rlm(y~x)), col='red')
  abline(lm(y~x), lty=3, lwd=3)
dev.off()

png(file="g817.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- 1 - 2*x + rcauchy(n,1)
  plot(y~x)
  abline(1,-2,lty=3)
  abline(lm(rlm(y~x)), col='red')
  abline(lm(y~x), lty=3, lwd=3)
dev.off()

png(file="g818.png", width=600, height=600)
  #library(lqs)   # now merged into MASS
  x <- rnorm(20)
  y <- 1 - x + rnorm(20)
  x <- c(x,10)
  y <- c(y,1)
  plot(y~x)
  abline(1,-1, lty=3)
  abline(lm(y~x))
  abline(rlm(y~x, psi = psi.bisquare, init = "lts"), col='orange',lwd=3)
  abline(rlm(y~x), col='red')
  abline(rlm(y~x, psi = psi.hampel, init = "lts"), col='green')
  abline(rlm(y~x, psi = psi.bisquare), col='blue')
  title(main='Huber regression (rlm)')
dev.off()

png(file="g819.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- 1 - 2*x + rnorm(n)
  y[ sample(1:n, floor(n/4)) ] <- 7
  plot(y~x)
  r1 <- lm(y~x)
  r2 <- lqs(y~x, method='lts')
  abline(r1, col='red')
  abline(r2, col='green')
  abline(1,-2,lty=3)
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("Linear regression", "Trimmed regression"),
          lty=1, lwd=1,
          col=c("red", "green") )
  title("Least Trimmed Squares (LTS)")
dev.off()

png(file="g820.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  plot(r1)
  par(op)
dev.off()

png(file="g821.png", width=600, height=600)
  # Evil data...
  n <- 10
  x <- rnorm(n)
  y <- 1 - 2*x + rnorm(n)
  x[1] <- 5
  y[1] <- 0
  my.wls <- function (y,x) {
    # A first regression
    r <- lm(y~x)$residuals
    # The weights
    w <- compute.weights(r)
    # A new regression
    lm(y~x, weights=w)
  }
  compute.weights <- function (r) {
    # Compute the weights as you want, as long as they are positive,
    # sum up to 1 and the high-residuals have low weights.
    # My choice might be neither standard nor relevant.
    w <- r*r
    w <- w/mean(w)
    w <- 1/(1+w)
    w <- w/mean(w)
  }
  plot(y~x)
  abline(1,-2, lty=3)
  abline(lm(y~x))
  abline(my.wls(y,x), col='red')
  title(main="Weighted Regression")
dev.off()

png(file="g822.png", width=600, height=600)
  # Situation in which you would like to use weighted least squares
  N <- 500
  x <- runif(N)
  y <- 1 - 2 * x + (2 - 1.5 * x) * rnorm(N)
  op <- par(mar = c(1,1,3,1))
  plot(y ~ x, axes = FALSE,
       main = "Weighted least squares")
  box()
  for (u in seq(-3,3,by=.5)) {
    segments(0, 1 + 2 * u, 1, -1 + .5 * u,
             col = "blue")
  }
  abline(1, -2, col = "blue", lwd = 3)
  par(op)
dev.off()

png(file="g823.png", width=600, height=600)
  my.irls.plot <- function (y,x, n=10) {
    plot(y~x)
    abline(lm(y~x))
    r <- lm(y~x)$residuals
    for (i in 1:n) {
      w <- compute.weights(r)
      print(w)
      r <- lm(y~x, weights=w)
      abline(r, col=topo.colors(n)[i], lwd=ifelse(i==n,3,1))
      r <- r$residuals
    }
    lm(y~x, weights=w)
  }
  my.irls.plot(y,x)
  abline(1,-2, lty=3)
  abline(my.wls(y,x), col='blue', lty=3, lwd=3)
  title(main="Iteratively Reweighted Least Squares")
dev.off()

png(file="g824.png", width=600, height=600)
  op <- par(mfrow=c(2,2), mar=c(2,4,4,2))
  r <- function (tau, x) { ifelse(x<0, (tau-1)*x, tau*x) }
  curve(r(0,x), xlim=c(-1,1), ylim=c(0,1), lwd=3, main="Minimum", xlab="")
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=c(0,.25,.5,.75), lty=3)
  curve(r(.25,x), xlim=c(-1,1), ylim=c(0,1), lwd=3, main="First quartile")
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=c(0,.25,.5,.75), lty=3)
  curve(r(.5,x), xlim=c(-1,1), ylim=c(0,1), lwd=3, main="Median")
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=c(0,.25,.5,.75), lty=3)
  curve(r(.75,x), xlim=c(-1,1), ylim=c(0,1), lwd=3, main="Third quartile")
  abline(0,1,lty=2)
  abline(0,-1,lty=2)
  abline(h=c(0,.25,.5,.75), lty=3)
  par(op)
dev.off()

png(file="g825.png", width=600, height=600)
  N <- 2000
  x <- runif(N)
  y <- rnorm(N)
  y <- -1 + 2 * x + ifelse(y>0, y+5*x^2, y-x^2)
  plot(x,y)
  abline(lm(y~x), col="red")
dev.off()

png(file="g826.png", width=600, height=600)
  library(quantreg)
  plot(y~x)
  for (a in seq(.1,.9,by=.1)) {
    abline(rq(y~x, tau=a), col="blue", lwd=3)
  }  
dev.off()

png(file="g827.png", width=600, height=600)
  plot(y~x)
  for (a in seq(.1,.9,by=.1)) {
    r <- lprq(x,y, 
              h=bw.nrd0(x),  # See ?density
              tau=a)
    lines(r$xx, r$fv, col="blue", lwd=3)
  }  
dev.off()

png(file="g828.png", width=600, height=600)
  op <- par(mar=c(3,2,4,1))
  r <- rq(y~x, tau=1:49/50)
  plot(summary(r), nrow=1)
  par(op)
dev.off()

png(file="g829.png", width=600, height=600)
  y <- -1 + 2 * x + rnorm(N)
  op <- par(mar=c(3,2,4,1))
  r <- rq(y~x, tau=1:49/50)
  plot(summary(r), nrow=1)
  par(op)
dev.off()

png(file="g830.png", width=600, height=600)
  # library(lqs)   # now part of MASS
  n <- 100
  x <- rnorm(n)
  y <- 1 - 2*x + rnorm(n)
  y[ sample(1:n, floor(n/4)) ] <- 7
  plot(y~x)
  abline(1,-2,lty=3)
  r1 <- lm(y~x)
  r2 <- lqs(y~x, method='lts')
  r3 <- lqs(y~x, method='lqs')
  r4 <- lqs(y~x, method='lms')
  r5 <- lqs(y~x, method='S')
  abline(r1, col='red')
  abline(r2, col='green')
  abline(r3, col='blue')
  abline(r4, col='orange')
  abline(r5, col='purple')
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("Linear Regression", "LTS",
          "lqs", "lms", "S"),
          lty=1, lwd=1,
          col=c("red", "green", "blue", "orange", "purple") )
  title("LTS and variants")
dev.off()

png(file="g831.png", width=600, height=600)
  # Example
  x <- rnorm(20)
  y <- 1 - x + rnorm(20)
  x <- c(x,10)
  y <- c(y,1)
  plot(y~x)
  abline(lm(y~x), col='blue')
  title(main="Classical regression")
dev.off()

png(file="g832.png", width=600, height=600)
  # How bad usual regression is
  plot(lm(y~x), which=4, 
       main="Cook's distance")
dev.off()

png(file="g833.png", width=600, height=600)
  plot(y~x)
  for (i in 1:length(x))
    abline(lm(y~x, subset= (1:length(x))!=i ), col='red')
  title(main="Classical regression minus one point")
dev.off()

png(file="g834.png", width=600, height=600)
  # line (in the eda package)
  library(eda)
  plot(y~x)
  abline(1,-1, lty=3)
  abline(lm(y~x))
  abline(coef(line(x,y)), col='red')
  title(main='"line", in the "eda" package')
dev.off()

png(file="g835.png", width=600, height=600)
  # Trimmed regression
  #library(lqs)
  plot(y~x)
  abline(1,-1, lty=3)
  abline(lm(y~x))
  abline(lqs(y~x), col='red')
  title(main='Trimmed regression (lqs)')
dev.off()

png(file="g836.png", width=600, height=600)
  # glm (from the manual, it should be IWLS, but we get the same
  # result...) 
  plot(y~x)
  abline(1,-1, lty=3)
  abline(lm(y~x))
  abline(glm(y~x), col='red')
  title(main='"glm" regression')
dev.off()

png(file="g837.png", width=600, height=600)
  plot(y~x)
  abline(1,-1, lty=3)
  abline(lm(y~x))
  abline(rlm(y~x, psi = psi.bisquare, init = "lts"), col='orange',lwd=3)
  abline(rlm(y~x), col='red')
  abline(rlm(y~x, psi = psi.hampel, init = "lts"), col='green')
  abline(rlm(y~x, psi = psi.bisquare), col='blue')
  title(main='Huber regression (rlm)')
dev.off()

png(file="g838.png", width=600, height=600)
  my.ridge <- function (y,x,k=0) {
    xm <- apply(x,2,mean)
    ym <- mean(y)
    y <- y - ym
    x <- t( t(x) - xm )
    ss <- function (b) {
      t( y - x %*% b ) %*% ( y - x %*% b ) + k * t(b) %*% b
    }
    b <- nlm(ss, rep(0,dim(x)[2]))$estimate
    c(ym-t(b)%*%xm, b)
  }
  my.ridge.test <- function (n=20, s=.1) {
    x <- rnorm(n)
    x1 <- x + s*rnorm(n)
    x2 <- x + s*rnorm(n)
    x <- cbind(x1,x2)
    y <- x1 + x2 + 1 + rnorm(n)
    lambda <- c(0, .001, .01, .1, .2, .5, 1, 2, 5, 10)
    b <- matrix(nr=length(lambda), nc=1+dim(x)[2])
    for (i in 1:length(lambda)) {
      b[i,] <- my.ridge(y,x,lambda[i])
    }
    plot(b[,2], b[,3], 
         type="b", 
         xlim=range(c(b[,2],1)), ylim=range(c(b[,3],1)))
    text(b[,2], b[,3], lambda, adj=c(-.2,-.2), col="blue")
    points(1,1,pch="+", cex=3, lwd=3)
    points(b[8,2],b[8,3],pch=15)
  }
  my.ridge.test()
dev.off()

png(file="g839.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    my.ridge.test()
  }
  par(op)
dev.off()

png(file="g840.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    my.ridge.test(20,10)
  }
  par(op)
dev.off()

png(file="g841.png", width=600, height=600)
  my.sample <- function (n=20) {
    x <- rnorm(n)
    x1 <- x + .1*rnorm(n)
    x2 <- x + .1*rnorm(n)
    y <- 0 + x1 - x2 + rnorm(n)
    cbind(y, x1, x2)
  }

  n <- 500
  r <- matrix(NA, nr=n, nc=3)
  s <- matrix(NA, nr=n, nc=3)
  for (i in 1:n) {
    m <- my.sample()
    r[i,] <- lm(m[,1]~m[,-1])$coef
    s[i,2:3] <- lm.ridge(m[,1]~m[,-1], lambda=.1)$coef
    s[i,1] <- mean(m[,1])
  }
  plot( density(r[,1]), xlim=c(-3,3),
        main="Multicolinearity: high variance")
  abline(v=0, lty=3)
  lines( density(r[,2]), col='red' )
  lines( density(s[,2]), col='red', lty=2 )
  abline(v=1, col='red', lty=3)
  lines( density(r[,3]), col='blue' )
  lines( density(s[,3]), col='blue', lty=2 )
  abline(v=-1, col='blue', lty=3)
  # We give the mean, to show that it is biased
  evaluate.density <- function (d,x, eps=1e-6) {
    density(d, from=x-eps, to=x+2*eps, n=4)$y[2]
  }
  x<-mean(r[,2]); points( x, evaluate.density(r[,2],x) )
  x<-mean(s[,2]); points( x, evaluate.density(s[,2],x) )
  x<-mean(r[,3]); points( x, evaluate.density(r[,3],x) )
  x<-mean(s[,3]); points( x, evaluate.density(s[,3],x) )
  legend( par("usr")[1], par("usr")[4], yjust=1,
          c("intercept", "x1", "x2"),
          lwd=1, lty=1, 
          col=c(par('fg'), 'red', 'blue') )
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c("classical regression", "ridge regression"),
          lwd=1, lty=c(1,2), 
          col=par('fg') )
dev.off()

png(file="g842.png", width=600, height=600)
  n <- 500
  v <- matrix(c(0,1,-1), nr=n, nc=3, byrow=T)
  mse <- NULL
  kl <- c(1e-4, 2e-4, 5e-4, 
          1e-3, 2e-3, 5e-3, 
          1e-2, 2e-2, 5e-2,   
          .1,  .2,  .3,  .4,  .5,  .6,  .7,  .8,  .9, 
          1,  1.2, 1.4, 1.6, 1.8,   2)
  for (k in kl) {
    r <- matrix(NA, nr=n, nc=3)
    for (i in 1:n) {
      m <- my.sample()
      r[i,2:3] <- lm.ridge(m[,1]~m[,-1], lambda=k)$coef
      r[i,1] <- mean(m[,1])
    }
    mse <- append(mse, apply( (r-v)^2, 2, mean )[2])
  }
  plot( mse ~ kl, type='l' )  
  title(main="MSE evolution")
dev.off()

png(file="g843.png", width=600, height=600)
  plot( mse-min(mse)+.01 ~ kl, type='l', log='y' )
  title(main="MSE Evolution")
dev.off()

png(file="g844.png", width=600, height=600)
  plot( rank(mse) ~ kl, type='l' )
  title(main="MSE Evolution")
dev.off()

png(file="g845.png", width=600, height=600)
  m <- my.sample()
  b <- matrix(NA, nr=length(kl), nc=2)
  for (i in 1:length(kl)) {
    b[i,] <- lm.ridge(m[,1]~m[,-1], lambda=kl[i])$coef
  }
  matplot(kl, b, type="l")
  abline(h=c(0,-1,1), lty=3)
  # Heuristic estimation for k...
  k <- min( lm.ridge(m[,1]~m[,-1], lambda=kl)$GCV )
  abline(v=k, lty=3)
  title(main="Bias towards 0 in ridge regression")
dev.off()

png(file="g846.png", width=600, height=600)
  m <- my.sample()
  b <- matrix(NA, nr=length(kl), nc=2)
  for (i in 1:length(kl)) {
    b[i,] <- lm.ridge(m[,1]~m[,-1], lambda=kl[i])$coef
  }
  matplot(kl, b, type="l", log='x')
  abline(h=c(0,-1,1), lty=3)
  k <- min( lm.ridge(m[,1]~m[,-1], lambda=kl)$GCV )  
  abline(v=k, lty=3)
  title(main="Bias towards 0 in ridge regression")
dev.off()

png(file="g847.png", width=600, height=600)
  my.lm.ridge.diag <- function (y, x, k=.1) {
    my <- mean(y)
    y <- y - my
    mx <- apply(x,2,mean)
    x <- x - matrix(mx, nr=dim(x)[1], nc=dim(x)[2], byrow=T)
    sx <- apply(x,2,sd)
    x <- x/matrix(sx, nr=dim(x)[1], nc=dim(x)[2], byrow=T)
    b <- solve( t(x) %*% x + diag(k, dim(x)[2]), t(x) %*% y)
    v <- solve( t(x) %*% x + diag(k, dim(x)[2]), 
         t(x) %*% x %*% solve( t(x) %*% x + diag(k, dim(x)[2]), 
         diag( var(y), dim(x)[2] ) ))
    ss <- t(b) %*% t(x) %*% y
    list( b = b, varb = v, ss = ss )
  }

  m <- my.sample()
  b <- matrix(NA, nr=length(kl), nc=2)
  v <- matrix(NA, nr=length(kl), nc=1)
  ss <- matrix(NA, nr=length(kl), nc=1)
  for (i in 1:length(kl)) {
    r <- my.lm.ridge.diag(m[,1], m[,-1], k=kl[i])
    b[i,] <- r$b
    v[i,] <- sum(diag(r$v))
    ss[i,] <- r$ss
  }
  matplot(kl, b, 
          type="l", lty=1, col=par('fg'), axes=F, ylab='')
  axis(1)
  abline(h=c(0,-1,1), lty=3)
  par(new=T)
  matplot(kl, v, type="l", col='red', axes=F, ylab='')
  par(new=T)
  matplot(kl, ss, type="l", col='blue', axes=F, ylab='')
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c("parameters", "variance", "sum of squares"),
          lwd=1, lty=1, col=c(par('fg'), "red", "blue") )
dev.off()

png(file="g848.png", width=600, height=600)
  matplot(log(kl), b, 
          type="l", lty=1, col=par('fg'), axes=F, ylab='')
  axis(1)
  abline(h=c(0,-1,1), lty=3)
  par(new=T)
  matplot(log(kl), v, type="l", col='red', axes=F, ylab='')
  par(new=T)
  matplot(log(kl), ss, type="l", col='blue', axes=F, ylab='')
  # I cannot put the legend if the scale is logarithmic...
  legend( par("usr")[1], 
          .9*par("usr")[3] + .1*par("usr")[4], 
          yjust=0, xjust=0,
          c("parameters", "variance", "sum of squares"),
          lwd=1, lty=1, col=c(par('fg'), "red", "blue") )
dev.off()

png(file="g849.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (j in 1:9) {
    m <- my.sample()
    b <- matrix(NA, nr=length(kl), nc=2)
    for (i in 1:length(kl)) {
      b[i,] <- lm.ridge(m[,1]~m[,-1], lambda=kl[i])$coef
    }
    matplot(kl, b, type="l", log='x')
    abline(h=c(0,-1,1), lty=3)
    k <- min( lm.ridge(m[,1]~m[,-1], lambda=kl)$GCV )  
    abline(v=k, lty=3)
  }
  par(op)
dev.off()

png(file="g850.png", width=600, height=600)
  m <- my.sample()
  N <- 20
  err <- matrix(nr=length(kl), nc=N)
  for (j in 1:N) {
    s <- sample(dim(m)[1], floor(3*dim(m)[1]/4))
    mm <- m[s,]
    mv <- m[-s,]
    for (i in 1:length(kl)) {
      r <- lm.ridge(mm[,1]~mm[,-1], lambda=kl[i])
      # BUG...
      b <- r$coef / r$scales
      a <- r$ym - t(b) %*% r$xm
      p <- rep(a, dim(mv)[1]) + mv[,-1] %*% b
      e <- p-mv[,1]
      err[i,j] <- sum(e*e)
    }
  }
  err <- apply(err, 1, mean)
  plot(err ~ kl, type='l', log='x')
dev.off()

png(file="g851.png", width=600, height=600)
  my.lasso <- function (y,x,k=0) {
    xm <- apply(x,2,mean)
    ym <- mean(y)
    y <- y - ym
    x <- t( t(x) - xm )
    ss <- function (b) {
      t( y - x %*% b ) %*% ( y - x %*% b ) + k * sum(abs(b))
    }
    b <- nlm(ss, rep(0,dim(x)[2]))$estimate
    c(ym-t(b)%*%xm, b)
  }

  my.lasso.test <- function (n=20) {
    s <- .1
    x <- rnorm(n)
    x1 <- x + s*rnorm(n)
    x2 <- x + s*rnorm(n)
    x <- cbind(x1,x2)
    y <- x1 + x2 + 1 + rnorm(n)
    lambda <- c(0, .001, .01, .1, .2, .5, 1, 2, 5, 10)
    b <- matrix(nr=length(lambda), nc=1+dim(x)[2])
    for (i in 1:length(lambda)) {
      b[i,] <- my.lasso(y,x,lambda[i])
    }
    plot(b[,2], b[,3], 
         type = "b", 
         xlim = range(c(b[,2],1)), 
         ylim = range(c(b[,3],1)))
    text(b[,2], b[,3], lambda, adj=c(-.2,-.2), col="blue")
    points(1,1,pch="+", cex=3, lwd=3)
    points(b[8,2],b[8,3],pch=15)
  }
  my.lasso.test()
dev.off()

png(file="g852.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    my.lasso.test()
  }
  par(op)
dev.off()

png(file="g853.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    my.lasso.test(1000)
  }
  par(op)
dev.off()

png(file="g854.png", width=600, height=600)
  get.sample <- function (n=20,s=.1) {
    x <- rnorm(n)
    x1 <- x + s*rnorm(n)
    x2 <- x + s*rnorm(n)
    y <- x1 + x2 + 1 + rnorm(n)
    data.frame(y,x1,x2)
  }
  lambda <- c(0, .001, .01, .1, .2, .5, 1, 2, 5, 10)

  do.it <- function (n=20,s=.1) {
    d <- get.sample(n,s)
    y <- d$y
    x <- cbind(d$x1,d$x2)
    ridge <- matrix(nr=length(lambda), nc=1+dim(x)[2])
    for (i in 1:length(lambda)) {
      ridge[i,] <- my.ridge(y,x,lambda[i])
    }
    lasso <- matrix(nr=length(lambda), nc=1+dim(x)[2])
    for (i in 1:length(lambda)) {
      lasso[i,] <- my.lasso(y,x,lambda[i])
    }
    xlim <- range(c( 1, ridge[,2], lasso[,2] ))        
    ylim <- range(c( 1, ridge[,3], lasso[,3] ))        
    plot(ridge[,2], ridge[,3], 
         type = "b", col = 'red', 
         xlim = xlim, ylim = ylim)
    points(ridge[8,2],ridge[8,3],pch=15,col='red')
    lines(lasso[,2], lasso[,3], type="b")
    points(lasso[8,2],lasso[8,3],pch=15)
    points(1,1,pch="+", cex=3, lwd=3)
  }

  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    do.it()
  }
  par(op)
dev.off()

png(file="g855.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    do.it(100)
  }
  par(op)
dev.off()

png(file="g856.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    do.it(20,10)
  }
  par(op)
dev.off()

png(file="g857.png", width=600, height=600)
  my.pcr <- function (y,x,k) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    ym <- mean(y)
    xm <- apply(x,2,mean)
    # Ideally, we should also normalize x and y... 
    # (exercise left to the reader)
    y <- y - ym 
    x <- t( t(x) - xm )
    pc <- princomp(x)
    b <- lm(y~pc$scores[,1:k]-1)$coef
    b <- c(b, rep(0,p-k))
    b <- pc$loadings %*% b
    names(b) <- colnames(x)
    b <- c(ym-t(b)%*%xm, b)
    b
  }
  get.sample <- function (n=20, s=.1) {
    x <- rnorm(n)
    x1 <- x + s*rnorm(n)
    x2 <- x + s*rnorm(n)
    x3 <- x + s*rnorm(n)
    x4 <- x + s*rnorm(n)
    x <- cbind(x1,x2,x3,x4)
    y <- x1 + x2 - x3 - x4 + 1 + rnorm(n)
    list(x=x,y=y)
  }
  pcr.test <- function (n=20, s=.1) {
    s <- get.sample(n,s)
    x <- s$x
    y <- s$y
    pcr <- matrix(nr=4,nc=5)
    for (k in 1:4) {
      pcr[k,] <- my.pcr(y,x,k)
    }
    plot(pcr[,2], pcr[,3], 
         type = "b", 
         xlim = range(c(pcr[,2],1)), 
         ylim = range(c(pcr[,3],1)))
    points(pcr[4,2], pcr[4,3], lwd=2)
    points(1,1, pch="+", cex=3, lwd=3)
  }
  pcr.test()
dev.off()

png(file="g858.png", width=600, height=600)
  pcr.test <- function (n=20, s=.1) {
    s <- get.sample(n,s)
    x <- s$x
    y <- s$y

    lambda <- c(0, .001, .01, .1, .2, .5, 1, 2, 5, 10)
    ridge <- matrix(nr=length(lambda), nc=1+dim(x)[2])
    for (i in 1:length(lambda)) {
      ridge[i,] <- my.ridge(y,x,lambda[i])
    }

    pcr <- matrix(nr=4,nc=5)
    for (k in 1:4) {
      pcr[k,] <- my.pcr(y,x,k)
    }

    xlim <- range(c( 1, ridge[,2], pcr[,2] ))        
    ylim <- range(c( 1, ridge[,3], pcr[,3] ))        
    plot(ridge[,2], ridge[,3], 
         type = "b", col = 'red', 
         xlim = xlim, ylim = ylim)
    points(ridge[4,2], ridge[4,3],
           pch = 15, col = 'red')

    lines(pcr[,2], pcr[,3], type="b")
    points(pcr[4,2], pcr[4,3], lwd=2)
    points(1,1, pch="+", cex=3, lwd=3)    
  }
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    pcr.test()    
  }
  par(op)
dev.off()
detach.everything()

png(file="g859.png", width=600, height=600)
  x <- runif(20)
  y <- 1-2*x+.1*rnorm(20)
  res <- lm(y~x)
  plot(y~x)
  new <- data.frame( x=seq(0,1,length=21) )
  p <- predict(res, new)
  points( p ~ new$x, type='l' )
  p <- predict(res, new, interval='confidence')
  points( p[,2] ~ new$x, type='l', col="green" )
  points( p[,3] ~ new$x, type='l', col="green" )
  p <- predict(res, new, interval='prediction')
  points( p[,2] ~ new$x, type='l', col="red" )
  points( p[,3] ~ new$x, type='l', col="red" )
  title(main="Confidence and prediction bands")
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("Confidence band", "Prediction band"),
          lwd=1, lty=1, col=c("green", "red") )
dev.off()

png(file="g860.png", width=600, height=600)
  plot(y~x, xlim=c(-1,2), ylim=c(-3,3))
  new <- data.frame( x=seq(-2,3,length=200) )
  p <- predict(res, new)
  points( p ~ new$x, type='l' )
  p <- predict(res, new, interval='confidence')
  points( p[,2] ~ new$x, type='l', col="green" )
  points( p[,3] ~ new$x, type='l', col="green" )
  p <- predict(res, new, interval='prediction')
  points( p[,2] ~ new$x, type='l', col="red" )
  points( p[,3] ~ new$x, type='l', col="red" )
  title(main="Confidence and prediction bands")
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("Confidence band", "Prediction band"),
          lwd=1, lty=1, col=c("green", "red") )
dev.off()

png(file="g861.png", width=600, height=600)
  plot(y~x, xlim=c(-5,6), ylim=c(-11,11))
  new <- data.frame( x=seq(-5,6,length=200) )
  p <- predict(res, new)
  points( p ~ new$x, type='l' )
  p <- predict(res, new, interval='confidence')
  points( p[,2] ~ new$x, type='l', col="green" )
  points( p[,3] ~ new$x, type='l', col="green" )
  p <- predict(res, new, interval='prediction')
  points( p[,2] ~ new$x, type='l', col="red" )
  points( p[,3] ~ new$x, type='l', col="red" )
  title(main="Confidence and prediction bands")
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("Confidence band", "Prediction band"),
          lwd=1, lty=1, col=c("green", "red") )
dev.off()

png(file="g862.png", width=600, height=600)
  N <- 100
  n <- 20
  x <- runif(N, min=-1, max=1)
  y <- 1 - 2*x + rnorm(N, sd=abs(x))
  res <- lm(y~x)
  plot(y~x)
  x0 <- seq(-1,1,length=n)
  new <- data.frame( x=x0 )
  p <- predict(res, new)
  points( p ~ x0, type='l' )
  p <- predict(res, new, interval='prediction')
  segments( x0, p[,2], x0, p[,3], col='red')
  p <- predict(res, new, interval='confidence')
  segments( x0, p[,2], x0, p[,3], col='green', lwd=3 )
dev.off()

png(file="g863.png", width=600, height=600)
  mySegments <- function(a,b,c,d,...) {
    u <- par('usr')
    e <- (u[2]-u[1])/100
    segments(a,b,c,d,...)
    segments(a+e,b,a-e,b,...)
    segments(c+e,d,c-e,d,...)
  }
  plot(y~x)
  p <- predict(res, new)
  points( p ~ x0, type='l' )
  p <- predict(res, new, interval='prediction')
  mySegments( x0, p[,2], x0, p[,3], col='red')
  p <- predict(res, new, interval='confidence')
  mySegments( x0, p[,2], x0, p[,3], col='green', lwd=3 )
dev.off()

png(file="g864.png", width=600, height=600)
  library(ellipse)
  my.confidence.region <- function (g, a=2, b=3) {
    e <- ellipse(g,c(a,b))
    plot(e,
         type="l",
         xlim=c( min(c(0,e[,1])), max(c(0,e[,1])) ),
         ylim=c( min(c(0,e[,2])), max(c(0,e[,2])) ),
        )
    x <- g$coef[a]
    y <- g$coef[b]
    points(x,y,pch=18)
    cf <- summary(g)$coefficients
    ia <- cf[a,2]*qt(.975,g$df.residual)
    ib <- cf[b,2]*qt(.975,g$df.residual)
    abline(v=c(x+ia,x-ia),lty=2)
    abline(h=c(y+ib,y-ib),lty=2)
    points(0,0)
    abline(v=0,lty="F848")
    abline(h=0,lty="F848")
  }

  n <- 20
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- x1+x2+x3+rnorm(n)
  g <- lm(y~x1+x2+x3)
  my.confidence.region(g)
dev.off()

png(file="g865.png", width=600, height=600)
  n <- 20
  x <- rnorm(n)
  x1 <- x+.2*rnorm(n)
  x2 <- x+.2*rnorm(n)
  y <- x1+x2+rnorm(n)
  g <- lm(y~x1+x2)
  my.confidence.region(g)
dev.off()

png(file="g866.png", width=600, height=600)
  my.confidence.region <- function (g, a=2, b=3, which=0, col='pink') {
    e <- ellipse(g,c(a,b))
    x <- g$coef[a]
    y <- g$coef[b]
    cf <- summary(g)$coefficients
    ia <- cf[a,2]*qt(.975,g$df.residual)
    ib <- cf[b,2]*qt(.975,g$df.residual)
    xmin <- min(c(0,e[,1]))
    xmax <- max(c(0,e[,1]))
    ymin <- min(c(0,e[,2]))
    ymax <- max(c(0,e[,2]))
    plot(e,
         type="l",
         xlim=c(xmin,xmax),
         ylim=c(ymin,ymax),
        )
    if(which==1){ polygon(e,col=col) }
    else if(which==2){ rect(x-ia,par('usr')[3],x+ia,par('usr')[4],col=col,border=col) }
    else if(which==3){ rect(par('usr')[1],y-ib,par('usr')[2],y+ib,col=col,border=col) }
    lines(e)
    points(x,y,pch=18)
    abline(v=c(x+ia,x-ia),lty=2)
    abline(h=c(y+ib,y-ib),lty=2)
    points(0,0)
    abline(v=0,lty="F848")
    abline(h=0,lty="F848")
  }
  my.confidence.region(g, which=1)
dev.off()

png(file="g867.png", width=600, height=600)
  my.confidence.region(g, which=2)
dev.off()

png(file="g868.png", width=600, height=600)
  my.confidence.region(g, which=3)
dev.off()

png(file="g869.png", width=600, height=600)
  n <- 20000
  x <- runif(n)
  y <- 4 - 8*x + rnorm(n)
  plot(y~x, pch='.')
  abline(lm(y~x), col='red')
  arrows( .1, -6, .1, 6, code=3, lwd=3, col='blue' )
  arrows( .9, -3.2-2, .9, -3.2+2, code=3, lwd=3, col='blue' )
  text( .1, 6, "TSS", adj=c(0,0), cex=2, col='blue' )  
  text( .9, -3.2+2, "RSS", adj=c(1,0), cex=2, col='blue' )  
dev.off()

png(file="g870.png", width=600, height=600)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- x1^2+rnorm(n)
  x4 <- 1/(1+x2^2)+.2*rnorm(n)
  y <- 1+x1-x2+x3-x4+.1*rnorm(n)
  pairs(cbind(x1,x2,x3,x4,y))  
dev.off()

png(file="g871.png", width=600, height=200)
  r <- lm(y~x1+x2+x3+x4)  
  boxplot(r$res, horizontal=T)
dev.off()

png(file="g872.png", width=600, height=600)
  hist(r$res)
dev.off()

png(file="g873.png", width=600, height=600)
  plot(r$res, main='Residuals')
dev.off()

png(file="g874.png", width=600, height=600)
  plot(rstandard(r), main='Standardized residuals')
dev.off()

png(file="g875.png", width=600, height=600)
  plot(rstudent(r), main="Studentized residuals")
dev.off()

png(file="g876.png", width=600, height=600)
  plot(r$res ~ r$fitted.values, 
       main="Residuals and predicted values")    
  abline(h=0, lty=3)
dev.off()

png(file="g877.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  plot(r$res ~ x1)    
  abline(h=0, lty=3)
  plot(r$res ~ x2)    
  abline(h=0, lty=3)
  plot(r$res ~ x3)    
  abline(h=0, lty=3)
  plot(r$res ~ x4)    
  abline(h=0, lty=3)
  par(op)
dev.off()

png(file="g878.png", width=600, height=600)
  n <- 100
  x1 <- rnorm(n)
  x2 <- 1:n
  y <- rnorm(1)
  for (i in 2:n) {
    y <- c(y, y[i-1] + rnorm(1))
  }
  y <- x1 + y
  r <- lm(y~x1+x2) # Or simply: lm(y~x1)
  op <- par(mfrow=c(2,1))
  plot( r$res ~ x1 )
  plot( r$res ~ x2 )
  par(op)
dev.off()

png(file="g879.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- 1-x+rnorm(n)
  r <- lm(y~x)
  plot(r$res ~ y)    
  abline(h=0, lty=3)
  abline(lm(r$res~y),col='red')
  title(main='Not a good idea...')
dev.off()

png(file="g880.png", width=600, height=600)
  partial.regression.plot <- function (y, x, n, ...) {
    m <- as.matrix(x[,-n])
    y1 <- lm(y ~ m)$res
    x1 <- lm(x[,n] ~ m)$res
    plot( y1 ~ x1, ... )
    abline(lm(y1~x1), col='red')
  }

  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- x1+x2+rnorm(n)
  x <- cbind(x1,x2,x3)
  y <- x1+x2+x3+rnorm(n)
  op <- par(mfrow=c(2,2))
  partial.regression.plot(y, x, 1)
  partial.regression.plot(y, x, 2)
  partial.regression.plot(y, x, 3)
  par(op)
dev.off()

png(file="g881.png", width=600, height=600)
  library(car)
  av.plots(lm(y~x1+x2+x3),ask=F)
dev.off()

png(file="g882.png", width=600, height=600)
  my.partial.residual.plot <- function (y, x, i, ...) {
    r <- lm(y~x)
    xi <- x[,i]
    # Y, minus the linear effects of X_j
    yi <- r$residuals + r$coefficients[i] * x[,i]
    plot( yi ~ xi, ... )  
  }
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- x1+x2+rnorm(n)
  x <- cbind(x1,x2,x3)
  y <- x1+x2+x3+rnorm(n)
  op <- par(mfrow=c(2,2))
  my.partial.residual.plot(y, x, 1)
  my.partial.residual.plot(y, x, 2)
  my.partial.residual.plot(y, x, 3)
  par(op)
dev.off()

png(file="g883.png", width=600, height=600)
  n <- 10
  x <- seq(0,1,length=n)
  y <- 1-2*x+.3*rnorm(n)
  plot(spline(x, y, n = 10*n), col = 'red', type='l', lwd=3)
  points(y~x, pch=16, lwd=3, cex=2)
  abline(lm(y~x))
  title(main='Overfit')
dev.off()

png(file="g884.png", width=600, height=600)
  x <- runif(100, -1, 1)
  y <- 1-x+x^2+.3*rnorm(100)
  plot(y~x)
  abline(lm(y~x), col='red')
dev.off()

png(file="g885.png", width=600, height=600)
  plot(y~x)
  lines(smooth.spline(x,y), col='red', lwd=2)
  title(main="Splines can help you spot non-linear relations")
dev.off()

png(file="g886.png", width=600, height=600)
  plot(y~x)
  lines(lowess(x,y), col='red', lwd=2)
  title(main='Non-linear relations and "lowess"')
dev.off()

png(file="g887.png", width=600, height=600)
  plot(y~x)
  xx <- seq(min(x),max(x),length=100)
  yy <- predict( loess(y~x), data.frame(x=xx) )
  lines(xx,yy, col='red', lwd=3)
  title(main='Non-linear relation and "loess"')
dev.off()

png(file="g888.png", width=600, height=600)
  r <- lm(y~x)
  plot(r$residuals ~ r$fitted.values,
        xlab='predicted values', ylab='residuals',
        main='Residuals and predicted values')
  lines(lowess(r$fitted.values, r$residuals), col='red', lwd=2)
  abline(h=0, lty=3)
dev.off()

png(file="g889.png", width=600, height=600)
  plot(r$residuals ~ x,
        xlab='x', ylab='residuals',
        main='Residuals and the predictive variable')
  lines(lowess(x, r$residuals), col='red', lwd=2)
  abline(h=0, lty=3)
dev.off()

png(file="g890.png", width=600, height=600)
  n <- 20
  done.outer <- F
  while (!done.outer) {
    done <- F
    while(!done) {
      x <- rnorm(n)
      done <- max(x)>4.5
    }
    y <- 1 - 2*x + x*rnorm(n)
    r <- lm(y~x)
    done.outer <- max(cooks.distance(r))>5
  }
  plot(y~x)
  abline(1,-2,lty=2)
  abline(lm(y~x),col='red',lwd=3)
  lm(y~x)$coef
dev.off()

png(file="g891.png", width=600, height=200)
  boxplot(x, horizontal=T)
dev.off()

png(file="g892.png", width=600, height=200)
  stripchart(x, method='jitter')
dev.off()

png(file="g893.png", width=600, height=600)
  hist(x, col='light blue', probability=T)
  lines(density(x), col='red', lwd=3)
dev.off()

png(file="g894.png", width=600, height=200)
  boxplot(y, horizontal=T)
dev.off()

png(file="g895.png", width=600, height=200)
  stripchart(y, method='jitter')
dev.off()

png(file="g896.png", width=600, height=600)
  hist(y, col='light blue', probability=T)
  lines(density(y), col='red', lwd=3)
dev.off()

png(file="g897.png", width=600, height=600)
  plot(hat(x), type='h', lwd=5)
dev.off()

png(file="g898.png", width=600, height=600)
  plot(dffits(r),type='h',lwd=3)
dev.off()

png(file="g899.png", width=600, height=600)
  plot(dfbetas(r)[,1],type='h',lwd=3)
dev.off()

png(file="g900.png", width=600, height=600)
  plot(dfbetas(r)[,2],type='h',lwd=3)
dev.off()

png(file="g901.png", width=600, height=600)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  yy <- x1 - x2 + rnorm(n)
  yy[1] <- 10
  r <- lm(yy~x1+x2)
  pairs(dfbetas(r))
dev.off()

png(file="g902.png", width=600, height=600)
  cd <- cooks.distance(r)
  plot(cd,type='h',lwd=3)
dev.off()

png(file="g903.png", width=600, height=600)
  n <- 100
  xx <- rnorm(n)
  yy <- 1 - 2 * x + rnorm(n)
  rr <- lm(yy~xx)
  cd <- cooks.distance(rr)
  plot(cd,type='h',lwd=3)
dev.off()

png(file="g904.png", width=600, height=200)
  boxplot(cd, horizontal=T)
dev.off()

png(file="g905.png", width=600, height=200)
  stripchart(cd, method='jitter')
dev.off()

png(file="g906.png", width=600, height=600)
  hist(cd, probability=T, breaks=20, col='light blue')
dev.off()

png(file="g907.png", width=600, height=600)
  plot(density(cd), type='l', col='red', lwd=3)
dev.off()

png(file="g908.png", width=600, height=600)
  qqnorm(cd)
  qqline(cd, col='red')
dev.off()

png(file="g909.png", width=600, height=600)
  half.qqnorm <- function (x) {
    n <- length(x)
    qqplot(qnorm((1+ppoints(n))/2), x)
  }
  half.qqnorm(cd)
dev.off()

png(file="g910.png", width=600, height=600)
  m <- max(cooks.distance(r))
  plot(y~x, cex=1+5*cooks.distance(r)/m)
dev.off()

png(file="g911.png", width=600, height=600)
  cd <- cooks.distance(r)
  # rescaled Cook's distance
  rcd <- (99/4) * cd*(cd+1)^2
  rcd[rcd>100] <- 100
  plot(r$res~r$fitted.values, cex=1+.05*rcd)
  abline(h=0,lty=3)
dev.off()

png(file="g912.png", width=600, height=600)
  m <- max(cd)
  plot(r$res, 
       cex=1+5*cd/m,
       col=heat.colors(100)[ceiling(70*cd/m)],
       pch=16,
      )
  points(r$res, cex=1+5*cd/m)
  abline(h=0,lty=3)
dev.off()

png(file="g913.png", width=600, height=600)
  plot(r$res, 
       cex=1+.05*rcd,
       col=heat.colors(100)[ceiling(rcd)],
       pch=16,
      )
  points(r$res, cex=1+.05*rcd)
  abline(h=0,lty=3)
dev.off()

png(file="g914.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  y <- 1 - 2*x + rnorm(n)
  r <- lm(y~x)
  cd <- cooks.distance(r)
  m <- max(cd)
  plot(r$res ~ r$fitted.values,
       cex=1+5*cd/m,
       col=heat.colors(100)[ceiling(70*cd/m)],
       pch=16,
      )
  points(r$res ~ r$fitted.values, cex=1+5*cd/m)
  abline(h=0,lty=3)
dev.off()

png(file="g915.png", width=600, height=600)
  op <- par(fg='white', bg='black', 
            col='white', col.axis='white', 
            col.lab='white', col.main='white', 
            col.sub='white')
  plot(r$res ~ r$fitted.values,
       cex=1+5*cd/m,
       col=heat.colors(100)[ceiling(100*cd/m)],
       pch=16,
      )
  abline(h=0,lty=3)
  par(op)
dev.off()

png(file="g916.png", width=600, height=600)
  # With Cook's distance
  x <- rnorm(20)
  y <- 1 + x + rnorm(20)
  x <- c(x,10)
  y <- c(y,1)
  r <- lm(y~x)
  d <- cooks.distance(r)
  d <- (99/4)*d*(d+1)^2 + 1
  d[d>100] <- 100
  d[d<20] <- 20
  d <- d/20
  plot( y~x, cex=d )
  abline(r)
  abline(coef(line(x,y)), col='red')
  abline(lm(y[1:20]~x[1:20]),col='blue')
dev.off()

png(file="g917.png", width=600, height=600)
  n <- 200
  s <- .2
  x <- runif(n)
  y1 <- 1 - 2 * x + s*rnorm(n)
  y2 <- 2 * x - 1 + s*rnorm(n)
  y <- ifelse( sample(c(T,F),n,replace=T,prob=c(.25,.75)), y1, y2 )
  plot(y~x)
  abline(1,-2,lty=3)
  abline(-1,2,lty=3)
dev.off()

png(file="g918.png", width=600, height=200)
  x <- runif(100)
  y <- 1 - 2*x + .3*exp(rnorm(100)-1)
  r <- lm(y~x)
  boxplot(r$residuals, horizontal=T)
dev.off()

png(file="g919.png", width=600, height=600)
  hist(r$residuals, breaks=20, probability=T, col='light blue')
  lines(density(r$residuals), col='red', lwd=3)
  f <- function(x) {
    dnorm(x,
          mean=mean(r$residuals),
          sd=sd(r$residuals),
    )
  }
  curve(f, add=T, col="red", lwd=3, lty=2)
dev.off()

png(file="g920.png", width=600, height=600)
  qqnorm(r$residuals)
  qqline(r$residuals, col='red')
dev.off()

png(file="g921.png", width=600, height=600)
  rcauchy.with.hole <- function (n) { 
    x <- rcauchy(n)
    x[x>0] <- 10+x[x>0]
    x[x<0] <- -10+x[x<0]
    x
  }
  n <- 20
  x <- rcauchy(n)
  y <- 1 - 2*x + .5*rcauchy.with.hole(n)
  plot(y~x)
  abline(1,-2)
  r <- lm(y~x)
  abline(r, col='red')
dev.off()

png(file="g922.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  hist(r$residuals, breaks=20, probability=T, col='light blue')
  lines(density(r$residuals), col='red', lwd=3)
  f <- function(x) {
    dnorm(x,
          mean=mean(r$residuals),
          sd=sd(r$residuals),
    )
  }
  curve(f, add=T, col="red", lwd=3, lty=2)
  qqnorm(r$residuals)
  qqline(r$residuals, col='red')
  plot(r$residuals ~ r$fitted.values)
  plot(r$residuals ~ x)
  par(op)
dev.off()

png(file="g923.png", width=600, height=600)
  done <- F
  while(!done) {
    # A situation where the prediction interval is not too
    # large, so that it appears on the plot.
    n <- 5
    x <- rcauchy(n)
    y <- 1 - 2*x + .5*rcauchy.with.hole(n)
    r <- lm(y~x)
    n <- 100000
    xp <- sort(runif(n,-50,50))
    yp <- predict(r, data.frame(x=xp), interval="prediction")
    done <- ( yp[round(n/2),2] > -75 & yp[round(n/2),3] < 75 )
  }
  yr <- 1 - 2*xp + .5*rcauchy.with.hole(n)
  plot(yp[,1]~xp, type='l',
       xlim=c(-50,50), ylim=c(-100,100))
  points(yr~xp, pch='.')
  lines(xp, yp[,2], col='blue')  
  lines(xp, yp[,3], col='blue')
  abline(r, col='red')
  points(y~x, col='orange', pch=16, cex=1.5)
  points(y~x, cex=1.5)
dev.off()

png(file="g924.png", width=600, height=600)
  done <- F
  while(!done) {
    # Even worse: the sign of the slope is incorrect
    n <- 5
    x <- rcauchy(n)
    y <- 1 - 2*x + .5*rcauchy.with.hole(n)
    r <- lm(y~x)
    n <- 100000
    xp <- sort(runif(n,-50,50))
    yp <- predict(r, data.frame(x=xp), interval="prediction")
    print(r$coef[2])
    done <- ( yp[round(n/2),2] > -75 &  
              yp[round(n/2),3] <  75 &  
              r$coef[2]>0 )
  }
  yr <- 1 - 2*xp + .5*rcauchy.with.hole(n)
  plot(yp[,1]~xp, type='l',
       xlim=c(-50,50), ylim=c(-100,100))
  points(yr~xp, pch='.')
  lines(xp, yp[,2], col='blue')  
  lines(xp, yp[,3], col='blue')
  abline(r, col='red')
  points(y~x, col='orange', pch=16, cex=1.5)
  points(y~x, cex=1.5)
dev.off()

png(file="g925.png", width=600, height=600)
  done <- F
  while (!done) {
    n <- 5
    x <- rcauchy(n)
    y <- 1 - 2*x + .5*rcauchy.with.hole(n)
    r <- lm(y~x)
    done <- T
  }
  n <- 10000
  xp <- sort(runif(n,-2,2))
  yp <- predict(r, data.frame(x=xp), interval="prediction")
  yr <- 1 - 2*xp + .5*rcauchy.with.hole(n)
  plot(c(xp,x), c(yp[,1],y), pch='.',
       xlim=c(-2,2), ylim=c(-50,50) )
  lines(yp[,1]~xp)
  abline(r, col='red')
  lines(xp, yp[,2], col='blue')  
  lines(xp, yp[,3], col='blue')
  points(yr~xp, pch='.')
  points(y~x, col='orange', pch=16)
  points(y~x)
dev.off()

png(file="g926.png", width=600, height=600)
  done <- F
  essais <- 0
  while (!done) {
    n <- 5
    x <- rcauchy(n)
    y <- 1 - 2*x + .5*rcauchy.with.hole(n)
    r <- lm(y~x)
    yp <- predict(r, data.frame(x=2), interval='prediction')
    done <- yp[3]<0
    essais <- essais+1
  }
  print(essais) # Around 20 or 30
  n <- 10000
  xp <- sort(runif(n,-2,2))
  yp <- predict(r, data.frame(x=xp), interval="prediction")
  yr <- 1 - 2*xp + .5*rcauchy.with.hole(n)
  plot(c(xp,x), c(yp[,1],y), pch='.',
       xlim=c(-2,2), ylim=c(-50,50) )
  lines(yp[,1]~xp)
  points(yr~xp, pch='.')
  abline(r, col='red')
  lines(xp, yp[,2], col='blue')  
  lines(xp, yp[,3], col='blue')
  points(y~x, col='orange', pch=16)
  points(y~x)
dev.off()

png(file="g927.png", width=600, height=600)
  done <- F
  e <- NULL
  for (i in 1:100) {
    essais <- 0
    done <- F
    while (!done) {
      n <- 5
      x <- rcauchy(n)
      y <- 1 - 2*x + .5*rcauchy.with.hole(n)
      r <- lm(y~x)
      yp <- predict(r, data.frame(x=2), interval='prediction')
      done <- yp[3]<0
      essais <- essais+1
    }
    e <- append(e,essais)
  }
  hist(e, probability=T, col='light blue')
  lines(density(e), col='red', lwd=3)
  abline(v=median(e), lty=2, col='red', lwd=3)
dev.off()

png(file="g928.png", width=600, height=600)
  x <- runif(100)
  y <- 1 - 2*x + .3*x*rnorm(100)
  plot(y~x)
  r <- lm(y~x)
  abline(r, col='red')
  title(main="Heteroscedasticity")
dev.off()

png(file="g929.png", width=600, height=600)
  plot(r$residuals ~ r$fitted.values)
dev.off()

png(file="g930.png", width=600, height=600)
  plot(abs(r$residuals) ~ r$fitted.values)
  lines(lowess(r$fitted.values, abs(r$residuals)), col='red')
dev.off()

png(file="g931.png", width=600, height=600)
  plot(abs(r$residuals) ~ x)
  lines(lowess(x, abs(r$residuals)), col='red')
dev.off()

png(file="g932.png", width=600, height=600)
  data(crabs)
  plot(FL~RW, data=crabs)
dev.off()

png(file="g933.png", width=600, height=600)
  r <- lm(FL~RW, data=crabs)
  plot(r, which=1)
dev.off()

png(file="g934.png", width=600, height=600)
  plot(r, which=3, panel = panel.smooth)
dev.off()

png(file="g935.png", width=600, height=600)
  library(car)
  spread.level.plot(r)
dev.off()

png(file="g936.png", width=600, height=600)
  x <- runif(100)
  y <- 1 - 2*x + .3*x*rnorm(100)
  r <- lm(y~x)
  n <- 10000
  xp <- sort(runif(n,))
  yp <- predict(r, data.frame(x=xp), interval="prediction")
  yr <- 1 - 2*xp + .3*xp*rnorm(n)

  plot(c(xp,x), c(yp[,1],y), pch='.')
  lines(yp[,1]~xp)
  abline(r, col='red')
  lines(xp, yp[,2], col='blue')  
  lines(xp, yp[,3], col='blue')
  points(yr~xp, pch='.')
  points(y~x, col='orange', pch=16)
  points(y~x)
  title(main="Consequences of heteroscedasticity on prediction intervals")
dev.off()

png(file="g937.png", width=600, height=600)
  n <- 100
  x <- runif(n)
  y <- 1 - 2*x + x*rnorm(n)
  plot(y~x)
  r <- lm(y~x)
  abline(r, col='red')
  title(main="Classical linear regression")
dev.off()

png(file="g938.png", width=600, height=600)
  plot(abs(r$res) ~ x)
  r2 <- lm( abs(r$res) ~ x )
  abline(r2, col="red")
  title(main="Heteroscedasticity of the residuals")
dev.off()

png(file="g939.png", width=600, height=600)
  # We assume the the standard deviation of the residuals 
  # is of the form a*x
  a <- lm( I(r$res^2) ~ I(x^2) - 1 )$coefficients
  w <- (a*x)^-2
  r3 <- lm( y ~ x, weights=w )
  plot(y~x)
  abline(1,-2, lty=3)
  abline(lm(y~x), lty=3, lwd=3)
  abline(lm(y~x, weights=w), col='red') 
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("acutal model", "least squares", 
            "weighted least squares"),
          lwd=c(1,3,1),
          lty=c(3,3,1),
          col=c(par("fg"), par("fg"), 'red') )
  title("Weighted least squares and heteroscedasticity")
dev.off()

png(file="g940.png", width=600, height=600)
  # Prediction intervals
  N <- 10000
  xx <- runif(N,min=0,max=2)
  yy <- 1 - 2*xx + xx*rnorm(N)
  plot(y~x, xlim=c(0,2), ylim=c(-3,2))
  points(yy~xx, pch='.')
  abline(1,-2, col='red')
  xp <- seq(0,3,length=100)
  yp1 <- predict(r, new=data.frame(x=xp), interval='prediction')
  lines( xp, yp1[,2], col='red', lwd=3 )
  lines( xp, yp1[,3], col='red', lwd=3 )
  yp3 <- predict(r3, new=data.frame(x=xp), interval='prediction')
  lines( xp, yp3[,2], col='blue', lwd=3 )
  lines( xp, yp3[,3], col='blue', lwd=3 )
  legend( par("usr")[1], par("usr")[3], yjust=0,
          c("least squares", "weighted least squares"),
          lwd=3, lty=1,
          col=c('red', 'blue') )
  title(main="Prediction band")
dev.off()

png(file="g941.png", width=600, height=600)
  my.acf.plot <- function (x, n=10, ...) {
    y <- rep(NA,n)
    l <- length(x)
    for (i in 1:n) {
      y[i] <- cor( x[1:(l-i)], x[(i+1):l] )
    }
    plot(y, type='h', ylim=c(-1,1),...)
  }
  n <- 100
  x <- runif(n)
  b <- .1*rnorm(n+1)
  y <- 1-2*x+b[1:n]
  my.acf.plot(lm(y~x)$res, lwd=10)
  abline(h=0, lty=2)
dev.off()

png(file="g942.png", width=600, height=600)
  z <- 1-2*x+.5*(b[1:n]+b[1+1:n])
  my.acf.plot(lm(z~x)$res, lwd=10)
  abline(h=0, lty=2)
dev.off()

png(file="g943.png", width=600, height=600)
  n <- 500
  x <- runif(n)
  b <- rep(NA,n)
  b[1] <- 0
  for (i in 2:n) {
    b[i] <- b[i-1] + .1*rnorm(1)
  }
  y <- 1-2*x+b[1:n]
  my.acf.plot(lm(y~x)$res, n=100)
  abline(h=0, lty=2)
  title(main='Very autocorrelated example')
dev.off()

png(file="g944.png", width=600, height=600)
  r <- lm(y~x)
  plot(r$res ~ r$fitted.values)
  title(main="Residuals of the very correlated example")
dev.off()

png(file="g945.png", width=600, height=600)
  r <- lm(y~x)
  plot(r$res)
  title(main="Residuals of the very correlated example")
dev.off()

png(file="g946.png", width=600, height=600)
  n <- 500
  x <- runif(n)
  b <- rep(NA,n)
  b[1] <- 0
  for (i in 2:n) {
    b[i] <- b[i-1] + .1*rnorm(1)
  }
  y <- 1-2*x+b[1:n]
  r <- lm(y~x)$res
  plot( r[1:(n-1)], r[2:n],
        xlab='i-th residual', 
        ylab='(i+1)-th residual' )
dev.off()

png(file="g947.png", width=600, height=600)
  n <- 500
  x <- runif(n)
  b <- rep(NA,n)
  b[1] <- 0
  b[2] <- 0
  for (i in 3:n) {
    b[i] <- b[i-2] + .1*rnorm(1)
  }
  y <- 1-2*x+b[1:n]
  r <- lm(y~x)$res
  plot(data.frame(x=r[3:n-2], y=r[3:n-1], z=r[3:n]))
dev.off()

png(file="g948.png", width=600, height=600)
  plot(r)
dev.off()

png(file="g949.png", width=600, height=600)
  data(randu)
  plot(randu)
  # Nothing visible
dev.off()

png(file="g950.png", width=600, height=600)
  m <- matrix( c(0.0491788982891203, -0.998585856299176, 0.0201921658647648,  
                 0.983046639705112,  0.0448184901961194, -0.177793720645666, 
                 -0.176637312387723,  -0.028593540105802, -0.983860594462783),
               nr=3, nc=3)
  plot( t( m %*% t(randu) )[,1:2] )
dev.off()

png(file="g951.png", width=600, height=600)
  data(EuStockMarkets)
  plot(EuStockMarkets)
dev.off()

png(file="g952.png", width=600, height=600)
  x <- EuStockMarkets[,1]
  y <- EuStockMarkets[,2]
  r <- lm(y~x)
  plot(y~x)
  abline(r, col='red', lwd=3)
dev.off()

png(file="g953.png", width=600, height=600)
  plot(r, which=1)
dev.off()

png(file="g954.png", width=600, height=600)
  plot(r, which=3)
dev.off()

png(file="g955.png", width=600, height=600)
  plot(r, which=4)
dev.off()

png(file="g956.png", width=600, height=600)
  r <- r$res
  hist(r, probability=T, col='light blue')
  lines(density(r), col='red', lwd=3)
dev.off()

png(file="g957.png", width=600, height=600)
  plot(r)
dev.off()

png(file="g958.png", width=600, height=600)
  acf(r)
dev.off()

png(file="g959.png", width=600, height=600)
  pacf(r)
dev.off()

png(file="g960.png", width=600, height=600)
  r <- as.vector(r)
  x <- r[1:(length(r)-1)]
  y <- r[2:length(r)]
  plot(x,y, xlab='x[i]', ylab='x[i+1]')
dev.off()

png(file="g961.png", width=600, height=600)
  n <- 100
  x <- rnorm(n)
  e <- vector()
  e <- append(e, rnorm(1))
  for (i in 2:n) {
    e <- append(e, .6 * e[i-1] + rnorm(1) )
  }
  y <- 1 - 2*x + e
  i <- 1:n
  plot(y~x)
dev.off()

png(file="g962.png", width=600, height=600)
  r <- lm(y~x)$residuals
  plot(r)
dev.off()

png(file="g963.png", width=600, height=600)
  library(nlme)
  plot(y~x)
  abline(lm(y~x))
  abline(gls(y~x, correlation = corAR1(form= ~i)), col='red')
dev.off()

png(file="g964.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  e <- vector()
  e <- append(e, rnorm(1))
  for (i in 2:n) {
    e <- append(e, 1 * e[i-1] + rnorm(1) )
  }
  y <- 1 - 2*x + e
  i <- 1:n
  plot(lm(y~x)$residuals)
dev.off()

png(file="g965.png", width=600, height=600)
  plot(y~x)
  abline(lm(y~x))
  abline(gls(y~x, correlation = corAR1(form= ~i)), col='red')
  abline(1,-2, lty=2)
dev.off()

png(file="g966.png", width=600, height=600)
  check.multicolinearity <- function (M) {
    a <- NULL
    n <- dim(M)[2]
    for (i in 1:n) {
      m <- as.matrix(M[, 1:n!=i])
      y <- M[,i]
      a <- append(a, summary(lm(y~m))$adj.r.squared)
    }
    names(a) <- names(M)
    print(round(a,digits=2))
    invisible(a)
  }
  data(freeny)
  names(freeny) <- paste(
    names(freeny), 
    " (",
    round(check.multicolinearity(freeny), digits=2),
    ")",
    sep='')
  pairs(freeny,
    upper.panel=panel.smooth,
    lower.panel=panel.smooth)
dev.off()

png(file="g967.png", width=600, height=600)
  n <- 100
  v <- .1
  x <- rnorm(n)
  x1 <- x + v*rnorm(n)
  x2 <- rnorm(n)
  x3 <- x + v*rnorm(n)
  y <- x1+x2-x3 + rnorm(n)
  m <- summary(lm(y~x1+x2+x3), correlation=T)$correlation
  plot(col(m), row(m), cex=10*abs(m),
       xlim=c(0, dim(m)[2]+1), 
       ylim=c(0, dim(m)[1]+1),
       main="Correlation matrix of the coefficients of a regression")
dev.off()

png(file="g968.png", width=600, height=600)
  n <- 100
  x <- runif(n)
  z <- ifelse(x>.5,1,0)
  y <- 2*z -x + .1*rnorm(n)
  plot( y~x, col=c('red','blue')[1+z] )
dev.off()

png(file="g969.png", width=600, height=600)
  n <- 20
  x <- rnorm(n)
  y <- 1 - 2*x - .1*x^2 + rnorm(n)
  #summary(lm(y~poly(x,10)))
  plot(y~x, xlim=c(-20,20), ylim=c(-30,30))
  r <- lm(y~x)
  abline(r, col='red')
  xx <- seq(-20,20,length=100)
  p <- predict(r, data.frame(x=xx), interval='prediction')
  lines(xx,p[,2],col='blue')
  lines(xx,p[,3],col='blue')
  title(main="Widening of the prediction band")
dev.off()

png(file="g970.png", width=600, height=600)
  plot(y~x, xlim=c(-20,20), ylim=c(-30,30))
  r <- lm(y~x)
  abline(r, col='red')
  xx <- seq(-20,20,length=100)
  yy <- 1 - 2*xx - .1*xx^2 + rnorm(n)
  p <- predict(r, data.frame(x=xx), interval='prediction')
  points(yy~xx)
  lines(xx,p[,2],col='blue')
  lines(xx,p[,3],col='blue')
  title(main="Extrapolation problem: it is not linear...")
dev.off()

png(file="g971.png", width=600, height=600)
  data(cars)
  y <- cars$dist
  x <- cars$speed
  o <- x<quantile(x,.25)
  x1 <- x[o]
  y1 <- y[o]
  r <- lm(y1~x1)
  xx <- seq(min(x),max(x),length=100)
  p <- predict(r, data.frame(x1=xx), interval='prediction')
  plot(y~x)
  abline(r, col='red')
  lines(xx,p[,2],col='blue')
  lines(xx,p[,3],col='blue')
dev.off()

png(file="g972.png", width=600, height=600)
  n <- 100
  e <- .5
  x0 <- rnorm(n)
  x <- x0 + e*rnorm(n)
  y <- 1 - 2*x0 + e*rnorm(n)
  plot(y~x)
  points(y~x0, col='red')
  abline(lm(y~x))
  abline(lm(y~x0),col='red')
dev.off()

png(file="g973.png", width=600, height=600)
  n <- 100
  k <- 5
  x <- matrix(rnorm(n*k),nc=k)
  y <- x[,1] + x[,2] - sqrt(abs(x[,3]*x[,4]))
  y <- y-median(y)
  y <- factor(y>0)
  pairs(x, col=as.numeric(y)+1)
dev.off()

png(file="g974.png", width=600, height=600)
  library(nlme)
  n <- 20
  m <- 15
  d <- as.data.frame(matrix(rnorm(n*m),nr=n,nc=m))
  # i <- sample(1:m, 3)
  i <- 1:3
  d <- data.frame(y=apply(d[,i],1,sum)+rnorm(n), d)
  k <- m
  res <- matrix(nr=k, nc=5)
  for (j in 1:k) {
    r <- lm(d$y ~ as.matrix(d[,2:(j+1)]))
    res[j,] <- c( logLik(r), AIC(r), BIC(r), 
                  summary(r)$r.squared, 
                  summary(r)$adj.r.squared )
  }
  colnames(res) <- c('logLik', 'AIC', 'BIC', 
                     "R squared", "adjusted R squared")
  res <- t( t(res) - apply(res,2,mean) )
  res <- t( t(res) / apply(res,2,sd) )
  matplot(0:(k-1), res, 
          type = 'l', 
          col = c(par('fg'),'blue','green', 'orange', 'red'), 
          lty = 1,
          xlab = "Number of variables")
  legend(par('usr')[2], par('usr')[3], 
         xjust = 1, yjust = 0,
         c('log-vraissemblance', 'AIC', 'BIC', 
           "R^2", "adjusted R^2" ),
         lwd = 1, lty = 1,
         col = c(par('fg'), 'blue', 'green', "orange", "red") )
  abline(v=3, lty=3)
dev.off()

png(file="g975.png", width=600, height=600)
  get.sample <- function () {
    # Number of observations
    n <- 20 
    # Number of variables
    m <- 10
    # Number of the variables that actually appear in the model
    k <- sample(1:m, 5)
    print(k)    
    # Coefficients
    b <- rnorm(m); b <- round(sign(b)+b); b[-k] <- 0
    x <- matrix(nr=n, nc=m, rnorm(n*m))
    y <- x %*% b + rnorm(n)
    data.frame(y=y, x)
  }

  my.variable.selection <- function (y,x, p=.05) {
    nvar <- dim(x)[2]
    nobs <- dim(x)[1]
    v <- rep(FALSE, nvar)
    p.values <- matrix(NA, nr=nvar, nc=nvar)
    res1 <- list()
    res2 <- list()
    done <- FALSE
    while (!done) {
      print(paste("Iteration", sum(v)))
      done <- TRUE
      # Is there a p-value lower that 
      pmax <- 1
      imax <- NA
      for (i in 1:nvar) {
        if(!v[i]){
          # Compute the p-value
          m <- cbind(x[,v], x[,i])
          m <- as.matrix(m)
          pv <- 1
          try( pv <- summary(lm(y~m))$coefficients[ dim(m)[2]+1, 4 ] )
          if( is.nan(pv) ) pv <- 1
          if (pv<pmax) {
            pmax <- pv
            imax <- i
          }
          p.values[i,sum(v)+1] <- pv
        }
      }
      if (pmax<p) {
        print(paste("Adding variable", imax, "with p-value", pmax))
        m1 <- as.matrix(x[,v])
        res1[[ length(res1)+1 ]] <- NULL
        try( res1[[ length(res1)+1 ]] <- data.frame(res=lm(y~m1)$res,xi=x[,imax]) )
        v[imax] <- TRUE
        done <- FALSE
        m2 <- as.matrix(cbind(x[,v], x[,imax]))
        res2[[ length(res2)+1 ]] <- data.frame(res=lm(y~m2)$res,xi=x[,imax])
      }
    }
    list(variables=v, p.values=p.values[,1:sum(v)], res1=res1, res2=res2)
  }

  d <- get.sample()
  y <- d$y
  x <- d[,-1]
  res <- my.variable.selection(y,x)

  k <- ceiling(length(res$res1)/3)
  op <- par(mfrow=c(k,3))
  for (i in 1:length(res$res1)) {
    r1 <- res$res1[[i]]
    r2 <- res$res2[[i]]
    plot(r1[,1] ~ r1[,2], ylab="res", xlab=names(r1)[2])
    points(r2[,1] ~ r2[,2],  col='red')
  }
  par(op)
dev.off()

png(file="g976.png", width=600, height=600)
  matplot(t(res$p.values), type='l', lty=1, lwd=1+2*res$variables)
  abline(h=.05, lty=3)
dev.off()

png(file="g977.png", width=600, height=600)
  library(ade4)
  data(microsatt)
  x <- microsatt$tab   # 18 observations, 112 variables, a lot of zeroes...
  y <- x[,3]
  x <- x[,-3]
  yn <- y/sqrt(sum(y*y))
  xn <- t(t(x)/sqrt(apply(x*x, 2, sum)))
  plot( sort(as.vector(t(yn) %*% xn)), type='h')
dev.off()

png(file="g978.png", width=600, height=600)
  library(leaps)
  library(car)
  get.sample <- function () {
    # Number of observations
    n <- 20 
    # Number of variables
    m <- 10
    # Number of the variables that actually appear in the model
    k <- sample(1:m, 5)
    print(k)    
    # Coefficients
    b <- rnorm(m); b <- round(sign(b)+b); b[-k] <- 0
    x <- matrix(nr=n, nc=m, rnorm(n*m))
    y <- x %*% b + rnorm(n)
    list(y=y, x=x, k=k, b=b)
  }
  d <- get.sample()
  x <- d$x
  y <- d$y
  k <- d$k
  b <- d$b
  subsets(regsubsets(x,y), statistic='bic', legend=F)
  title(main=paste(sort(k),collapse=', '))
dev.off()

png(file="g979.png", width=600, height=600)
  set.seed(1)
  n <- 20
  x <- runif(n, -1, 1)
  y <- 1 - x^2 + .2*rnorm(n)
  X <- runif(10000, -1, 1)
  Y <- 1 - X^2 + .2*rnorm(1000)
  N <- n
  res <- matrix(NA, nc=N, nr=2)
  dimnames(res) <- list(
    c("In-sample error", "Out-of-sample error"), 
    "Model complexity" = as.character(1:N)
  )
  r <- lm(y~x)
  res[1,1] <- mean(abs(residuals(r)))
  res[2,1] <- mean(abs(predict(r, data.frame(x=X)) - Y))  
  for (i in 2:N) {
    r <- lm(y ~ poly(x,i-1))
    res[1,i] <- mean(abs(residuals(r)))
    res[2,i] <- mean(abs(predict(r, data.frame(x=X)) - Y))  
  }
  
  op <- par(mar=c(5,4,4,4))
  ylim <- c(0, 1.5*max(res[1,]))
  plot(res[1,], col="blue", type="l", lwd=3,
       ylim=ylim,
       axes=F, 
       xlab="Model complexity", 
       ylab="", 
       main="In- and out-of-sample error")
  axis(1)
  axis(2, col="blue")
  par(new=TRUE)
  plot(res[2,], col="red", type="b", lwd=3,
       ylim=ylim,
       axes=F, xlab="", ylab="", main="")
  axis(4, col="red")
  mtext("In-sample error",     2, line=2, col="blue", cex=1.2)
  mtext("Out-of-sample error", 4, line=2, col="red",  cex=1.2)
  par(op)
dev.off()

png(file="g980.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)
  u <- sqrt(x^2+y^2)
  u <- ifelse( u<.5*mean(u), 1, 2)
  plot(y~x, col=u, pch=15)
dev.off()

png(file="g981.png", width=600, height=600)
  library(e1071)
  u <- factor(u)
  r <- svm(u~x+y)
  {
    # The "plot.svm" calls "browser()": Why???
    # (I use it when I debug my code, it is probably the
    # same for them.)
    # And then it crashes...
    # (e1071 version 1.3-16, 1.3-16)
    browser <- function () {}
    try(  plot(r, data.frame(x,y,u), y~x)  )
  }
dev.off()

png(file="g982.png", width=600, height=600)
  n <- 200
  x <- runif(n, -1,1)
  y <- runif(n, -1,1)
  u <- abs(x-y)
  u <- ifelse( u<.5*mean(u), 1, 2)
  plot(y~x, col=u, pch=15)
dev.off()

png(file="g983.png", width=600, height=600)
  u <- factor(u)
  r <- svm(u~x+y)
  { 
    browser <- function () {} 
    try(  plot(r, data.frame(x,y,u), y~x)  )
  }
dev.off()

png(file="g984.png", width=600, height=600)
  n <- 200
  x1 <- runif(n,-3,3)
  x2 <- runif(n,-3,3)
  x3 <- runif(n,-3,3)
  f1 <- sin; f2 <- cos; f3 <- abs;
  y <- f1(x1) + f2(x2) + f3(x3) + rnorm(n)
  pairs(cbind(y,x1,x2,x3)) # Nothing really visible...
dev.off()

png(file="g985.png", width=600, height=600)
  library(mgcv)
  r <- gam(y~s(x1)+s(x2)+s(x3))
  x <- seq(-3,3,length=200)
  z <- rep(0,200)
  m.theoretical <- cbind(f1(x),f2(x),f3(x))
  m.experimental <- cbind( 
    predict(r, data.frame(x1=x,x2=z,x3=z)),
    predict(r, data.frame(x1=z,x2=x,x3=z)),
    predict(r, data.frame(x1=z,x2=z,x3=x))
  )
  matplot(m.theoretical, type='l', lty=1)
  matplot(m.experimental, type='l', lty=2, add=T)
dev.off()

png(file="g986.png", width=600, height=600)
  zero.mean <- function (m) {
    t(t(m)-apply(m,2,mean))
  }
  matplot(zero.mean(m.theoretical), 
          type='l', lty=1)
  matplot(zero.mean(m.experimental), 
          type='l', lty=2, add=T)
  title(main="GAM")
  legend(par('usr')[2], par('usr')[3], 
         xjust=1, yjust=0,
         c('theoretical curves', 'experimental curves'),
         lty=c(1,2))
dev.off()

png(file="g987.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (i in 1:3) {
    plot(r, select=i)
    lines(zero.mean(m.theoretical)[,i] ~ x, 
          lwd=3, lty=3, col='red')
  }
  par(op)
dev.off()

png(file="g988.png", width=600, height=600)
  res <- residuals(r)
  op <- par(mfrow=c(2,2))
  plot(res)
  plot(res ~ predict(r))
  hist(res, col='light blue', probability=T)
  lines(density(res), col='red', lwd=3)
  rug(res)
  qqnorm(res)
  qqline(res)
  par(op)
dev.off()

png(file="g989.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  plot(res ~ x1)
  plot(res ~ x2)
  plot(res ~ x3)
  par(op)
dev.off()

png(file="g990.png", width=600, height=600)
  library(rpart)
  data(kyphosis)
  r <- rpart(Kyphosis ~ ., data=kyphosis)
  plot(r)
  text(r)
dev.off()

png(file="g991.png", width=600, height=600)
  library(randomForest)
  library(mlbench)
  library(mva)
  data(Vowel)
  r <- randomForest(
    Class ~ ., 
    data = Vowel, 
    importance = TRUE
  ) # Trois minutes...
  m <- r$confusion
  # We have a confusion matrix instead of a distance matrix.
  # We try to tweak the matrix to get something that looks
  # like a distance matrix -- this might not be the best way
  # to proceed.
  m <- m[,-dim(m)[2]]
  m <- m+t(m)-diag(diag(m))
  n <- dim(m)[1]
  m <- m/( matrix(diag(m),nr=n,nc=n) + 
    matrix(diag(m),nr=n,nc=n, byrow=T) )
  m <- 1-m
  diag(m) <- 0
  mds <- cmdscale(m,2)
  plot(mds, type='n')
  text(mds, colnames(m))
  rm(Vowel)
dev.off()

png(file="g992.png", width=600, height=600)
  op <- par(mfrow = c(2,2))
  m <- r$importance
  n <- dim(m)[1]
  for (i in 1:4) {
    plot(m[,i],
         type = "h", lwd=3, col='blue', axes=F,
         ylab='', xlab='Variables',
         main = colnames(r$importance)[i] )
    axis(2)
    axis(1, at=1:n, labels=rownames(m))
    # The two highest values in red
    a <- order(m[,i], decreasing=T)
    m2 <- m[,i]
    m2[ -a[1:2] ] <- NA
    lines(m2, type='h', lwd=5, col='red')
  }
  par(op)
dev.off()

png(file="g993.png", width=600, height=600)
  library(mda)
  library(mlbench)
  data(BostonHousing)
  x <- BostonHousing
  x[,4] <- as.numeric(x[,4])
  pairs(x)
dev.off()

png(file="g994.png", width=600, height=600)
  op <- par(mfrow=c(4,4))
  for (i in 1:14) {
    hist(x[,i],probability=T,
         col='light blue', main=paste(i,names(x)[i]))
    lines(density(x[,i]),col='red',lwd=3)
    rug(jitter(x[,i]))
  }
    hist(log(x[,1]),probability=T,
         col='light blue', main="log(x1)")
    lines(density(log(x[,1])),col='red',lwd=3)
    rug(jitter(log(x[,1])))
  par(op)
dev.off()

png(file="g995.png", width=600, height=600)
  op <- par(mfrow=c(4,4))
  for (i in 1:14) {
    qqnorm(x[,i], main=paste(i,names(x)[i]))
    qqline(x[,i], col='red')
  }
    qqnorm(log(x[,1]), main="log(x1)")
    qqline(log(x[,1]), col='red')
  par(op)
dev.off()

png(file="g996.png", width=600, height=600)
  x[,1] <- log(x[,1])
  n <- dim(x)[1]
  k <- sample(1:n, 100)
  d1 <- x[k,]
  d2 <- x[-k,]
  r <- mars(d1[,-1],d1[,1])
  p <- predict(r, d2[,-1])
  res <- d2[,1] - p

  op <- par(mfrow=c(4,4))
  plot(res)
  plot(res~p)
  for (i in 2:14) {
    plot(res~d2[,i])
  }  
  par(op)
dev.off()

png(file="g997.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  qqnorm(r$fitted.values, main="fitted values")
  qqline(r$fitted.values)
  hist(r$fitted.values,probability=T,
       col='light blue', main='fitted values')
  lines(density(r$fitted.values), col='red', lwd=3)
  rug(jitter(r$fitted.values))
  qqnorm(res, main="residuals")
  qqline(res)
  hist(res, probability=T,
       col='light blue', main='residuals')
  lines(density(res),col='red',lwd=3)
  rug(jitter(res))
  par(op)
dev.off()

png(file="g998.png", width=600, height=600)
  library(superpc)

  # Some simulated data
  n <- 50    # Number of observations (e.g., patients)
  m <- 500   # Number of values to predict
  p <- 1000  # Number of variables or "features" (e.g., genes)
  q <- 20    # Number of useful variables
  x <- matrix( rnorm((n+m)*p), nr=n+m )
  z <- svd(x[,1:q])$u
  y <- 1 - 2 * z[,1]
  d <- list(
    x = t(x[1:n,]),
    y = y[1:n],
    featurenames = paste("V", as.character(1:p),sep="")
  )
  new.x <- list(
    x = t(x[-(1:n),]),
    y = y[-(1:n)],
    featurenames = paste("V", as.character(1:p),sep="")
  )

  # Compute the correlation (more precisely, a score) of
  # each variable with the outcome.
  train.obj <- superpc.train(d, type="regression")
  hist(train.obj$feature.score, 
       col = "light blue",
       main = "SPCA (Supervised Principal Component Analysis)")
dev.off()

png(file="g999.png", width=600, height=600)
  ## PROBLEM: This should not look like that...

  # Compute the threshold to be used on the scores to select
  # the variables to retain. 
  cv.obj <- superpc.cv(train.obj, d)
  superpc.plotcv(cv.obj)
  title("SPCA")
dev.off()

png(file="g1000.png", width=600, height=600)
  fit.cts <- superpc.predict(
    train.obj, 
    d, 
    new.x,
    threshold = 0.7, 
    n.components = 3, 
    prediction.type = "continuous"
  )
  r <- superpc.fit.to.outcome(train.obj, new.x, fit.cts$v.pred)
  plot(r$results$fitted.values, new.x$y,
       xlab="fitted values", ylab="actual values",
       main="SPCA forecasts")
  abline( lm(new.x$y ~ r$results$fitted.values),
          col="red", lwd=3 )
  abline(0, 1, lty=2, lwd=2)
dev.off()
detach.everything()

png(file="g1001.png", width=600, height=600)
  n <- 10
  N <- 10
  s <- .3
  m1 <- rnorm(n, c(0,0))
  a <- rnorm(2*N*n, m1, sd=s)
  m2 <- rnorm(n, c(1,1))
  b <- rnorm(2*N*n, m2, sd=s)
  x1 <- c( a[2*(1:(N*n))], b[2*(1:(N*n))] )
  x2 <- c( a[2*(1:(N*n))-1], b[2*(1:(N*n))-1] )
  y <- c(rep(0,N*n), rep(1,N*n))
  plot( x1, x2, col=c('red','blue')[1+y] )
dev.off()

png(file="g1002.png", width=600, height=600)
  plot( x1, x2, col=c('red','blue')[1+y] )
  r <- lm(y~x1+x2)
  abline( (.5-coef(r)[1])/coef(r)[3], -coef(r)[2]/coef(r)[3] )
dev.off()

png(file="g1003.png", width=600, height=600)
  # I need a function to draw conics...
  conic.plot <- function (a,b,c,d,e,f, xlim=c(-2,2), ylim=c(-2,2), n=20, ...) {
    x0 <- seq(xlim[1], xlim[2], length=n)
    y0 <- seq(ylim[1], ylim[2], length=n)
    x <- matrix( x0, nr=n, nc=n )
    y <- matrix( y0, nr=n, nc=n, byrow=T )
    z <- a*x^2 + b*x*y + c*y^2 + d*x + e*y + f
    contour(x0,y0,z, nlevels=1, levels=0, drawlabels=F, ...)
  }
  r <- lm(y~x1+x2+I(x1^2)+I(x1*x2)+I(x2^2))$coef
  plot( x1, x2, col=c('red','blue')[1+y] )
  conic.plot(r[4], r[5], r[6], r[2], r[3], r[1]-.5, 
             xlim=par('usr')[1:2], ylim=par('usr')[3:4], add=T)
dev.off()

png(file="g1004.png", width=600, height=600)
  M <- 100
  d <- function (a,b, N=10) {
    mean( y[ order( (x1-a)^2 + (x2-b)^2 )[1:N] ] )
  }
  myOuter <- function (x,y,f) {
    r <- matrix(nrow=length(x), ncol=length(y))
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        r[i,j] <- f(x[i],y[j])
      }
    }
    r
  }
  cx1 <- seq(from=min(x1), to=max(x1), length=M)
  cx2 <- seq(from=min(x2), to=max(x2), length=M)
  plot( x1, x2, col=c('red','blue')[1+y] )
  contour(cx1, cx2, myOuter(cx1,cx2,d), levels=.5, add=T)
dev.off()

png(file="g1005.png", width=600, height=600)
  n <- 20
  N <- 10
  s <- .1
  m1 <- rnorm(n, c(0,0))
  a <- rnorm(2*N*n, m1, sd=s)
  m2 <- rnorm(n, c(0,0))
  b <- rnorm(2*N*n, m2, sd=s)
  x1 <- c( a[2*(1:(N*n))], b[2*(1:(N*n))] )
  x2 <- c( a[2*(1:(N*n))-1], b[2*(1:(N*n))-1] )
  y <- c(rep(0,N*n), rep(1,N*n))
  plot( x1, x2, col=c('red','blue')[1+y] )
  M <- 100
  cx1 <- seq(from=min(x1), to=max(x1), length=M)
  cx2 <- seq(from=min(x2), to=max(x2), length=M)
  #text(outer(cx1,rep(1,length(cx2))),
  #     outer(rep(1,length(cx1)), cx2),
  #     as.character(myOuter(cx1,cx2,d)))
  contour(cx1, cx2, myOuter(cx1,cx2,d), levels=.5, add=T)
  # Color the various areas
  points(matrix(cx1, nr=M, nc=M),
         matrix(cx2, nr=M, nc=M, byrow=T),
         pch='.',
         col=c("red", "blue")[ as.numeric( myOuter(cx1,cx2,d) >.5 )+1])
dev.off()

png(file="g1006.png", width=600, height=600)
  pastel <- .9
  plot( x1, x2, col=c('red','blue')[1+y] )
  points(matrix(cx1, nr=M, nc=M),
         matrix(cx2, nr=M, nc=M, byrow=T),
         pch=16,
         col=c(rgb(1,pastel,pastel), rgb(pastel,pastel,1))
             [ as.numeric( myOuter(cx1,cx2,d) >.5 )+1])
  points(x1, x2, col=c('red','blue')[1+y] )
  contour(cx1, cx2, myOuter(cx1,cx2,d), levels=.5, add=T)
dev.off()

png(file="g1007.png", width=600, height=600)
  plot( x1, x2, col=c('red','blue')[1+y] )
  v <- c(3,10,50)
  for (i in 1:length(v)) {
    contour(cx1, cx2, 
            myOuter(cx1,cx2, function(a,b){d(a,b,v[i])}), 
            levels=.5, add=T, drawlabels=F, col=i+1)
  }
  legend(min(x1),max(x2),as.character(v),col=1+(1:length(v)), lty=1)
dev.off()

png(file="g1008.png", width=600, height=600)
  get.model <- function (n=10, m=2, s=.1) {
    list( n=n, m=m, x=runif(n), y=runif(n), z=sample(1:m,n,replace=T), s=s )
  }
  get.sample <- function (model, n=200) {
    i <- sample( 1:model$n, n, replace=T )
    data.frame( x=rnorm(n, model$x[i], model$s),
                y=rnorm(n, model$y[i], model$s),
                z=model$z[i] )
  }
  nearest.neighbour.predict <- function (x, y, d, k=10) {
    o <- order( (d$x-x)^2 + (d$y-y)^2 )[1:k]
    s <- apply(outer(d[o,]$z, 1:max(d$z), '=='), 2, sum)
    order(s, decreasing=T)[1]
  }
  m <- get.model()
  d <- get.sample(m)
  N <- 1000
  d.test <- get.sample(m,N)
  n <- 50
  r <- rep(0, n)
  # Very slow
  for (k in 1:n) {
    for(i in 1:N) {
      r[k] <- r[k] + 
        (nearest.neighbour.predict(d.test$x[i], d.test$y[i], d, k) != d.test$z[i] )
    }
  }
  plot(r/N, ylim=c(0,1), type='l', xlab="Error rate")
  abline(h=c(0,.5,1), lty=3)
  rm(d.test)
dev.off()

png(file="g1009.png", width=600, height=600)
  m <- get.model()
  d <- get.sample(m, 20)
  N <- 1000
  d.test <- get.sample(m,N)
  n <- 50
  r <- rep(0, n)
  # Very slow
  for (k in 1:n) {
    for(i in 1:N) {
      r[k] <- r[k] + 
        (nearest.neighbour.predict(d.test$x[i], d.test$y[i], d, k) != d.test$z[i] )
    }
  }
  plot(r/N, ylim=c(0,1), type='l', xlab="Error rate")
  abline(h=c(0,.5,1), lty=3)
  rm(d.test)
dev.off()

png(file="g1010.png", width=600, height=600)
  nearest.neighbour.plot <- function (d, k=10, model=NULL) {
    col <- rainbow(max(d$z))
    plot( d$x, d$y, col=col )
    cx <- seq(min(d$x), max(d$x), length=100)
    cy <- seq(min(d$y), max(d$y), length=100)
    pastel <- .8
    colp <- round(pastel*255 + (1-pastel)*col2rgb(col))
    colp <- rgb(colp[1,], colp[2,], colp[3,], max=255)
    points(matrix(cx,nr=100,nc=100),
           matrix(cy,nr=100,nc=100,byrow=T),
           col = colp[ 
             myOuter(cx,cy, function(a,b){
               nearest.neighbour.predict(a,b,d,k)
             })
           ],
           pch=16
          )
    points( d$x, d$y, col=col )
    if(!is.null(model)){
      points(model$x,model$y,pch='+',cex=3,lwd=3,col=col[model$z])
    }
  }
  m <- get.model(n=10, m=4)
  d <- get.sample(m)
  nearest.neighbour.plot(d, model=m)
dev.off()

png(file="g1011.png", width=600, height=600)
  library(MASS)
  n <- 100
  k <- 5
  x1 <- runif(k,-5,5) + rnorm(k*n)
  x2 <- runif(k,-5,5) + rnorm(k*n)
  x3 <- runif(k,-5,5) + rnorm(k*n)
  x4 <- runif(k,-5,5) + rnorm(k*n)
  x5 <- runif(k,-5,5) + rnorm(k*n)
  y <- factor(rep(1:5,n))
  plot(lda(y~x1+x2+x3+x4+x5))
dev.off()

png(file="g1012.png", width=600, height=600)
  n <- 100
  x <- c(rnorm(n), 1+rnorm(n))
  y <- c(rep(0,n), rep(1,n))
  plot(y~x)
  abline(lm(y~x), col='red')
dev.off()

png(file="g1013.png", width=600, height=600)
  x <- seq(0,1, length=100)
  x <- x[2:(length(x)-1)]
  logit <- function (t) {
    log( t / (1-t) )
  }
  plot(logit(x) ~ x, type='l')
dev.off()

png(file="g1014.png", width=600, height=600)
  curve(logit(x), col='blue', add=F)
  curve(qnorm(x), col='red', add=T)
  a <- par("usr")
  legend( a[1], a[4], c("logit","probit"), col=c("blue","red"), lty=1)
dev.off()

png(file="g1015.png", width=600, height=600)
  curve(logit(x), col='blue', add=F)
  curve(qnorm(x), col='red', add=T)
  curve(log(-log(1-x)), col='green', add=T)
  abline(h=0, lty=3)
  abline(v=0, lty=3)
  a <- par("usr")
  legend( a[1], a[4], 
          c("logit","probit", "log-log"), 
          col=c("blue","red","green"), 
          lty=1)
dev.off()

png(file="g1016.png", width=600, height=600)
  ilogit <- function (l) {
    exp(l) / ( 1 + exp(l) )
  }
  fakelogit <- function (l) {
    ifelse(l>.5, 1e6, -1e6)
  }
  n <- 100
  x <- c(rnorm(n), 1+rnorm(n))
  y <- c(rep(0,n), rep(1,n))
  yy <- fakelogit(y)
  xp <- seq(min(x),max(x),length=200)
  yp <- ilogit(predict(lm(yy~x), data.frame(x=xp)))
  yp[is.na(yp)] <- 1
  plot(y~x)
  lines(xp,yp, col='red', lwd=3)
dev.off()

png(file="g1017.png", width=600, height=600)
  n <- 100
  x <- c(rnorm(n), 1+rnorm(n))
  y <- c(rep(0,n), rep(1,n))
  f <- function (a) {
    -sum(log(ilogit(a[1]+a[2]*x[y==1]))) - sum(log(1-ilogit(a[1]+a[2]*x[y==0])))
  }
  r <- optim(c(0,1),f)
  a <- r$par[1]
  b <- r$par[2]
  plot(y~x)
  curve( dnorm(x,1,1)*.5/(dnorm(x,1,1)*.5+dnorm(x,0,1)*(1-.5)), add=T, col='red')
  curve( ilogit(a+b*x), add=T )
  legend( .95*par('usr')[1]+.05*par('usr')[2],
          .9,
          c('theoretical curve', 'MLE'),
          col=c('red', par('fg')),
          lty=1, lwd=1)
  title(main="Logistic regression, by hand")
dev.off()

png(file="g1018.png", width=600, height=600)
  #
  # BEWARE:
  # Do not forget the "family" argument -- otherwise, it would be a
  # linear regression -- the very thing we are trying to avoid. 
  #
  r <- glm(y~x, family=binomial)
  plot(y~x)
  abline(lm(y~x),col='red',lty=2)
  xx <- seq(min(x), max(x), length=100)
  yy <- predict(r, data.frame(x=xx), type='response')
  lines(xx,yy, col='blue', lwd=5, lty=2)
  lines(xx, ilogit(r$coef[1]+xx*r$coef[2]))
  legend( .95*par('usr')[1]+.05*par('usr')[2],
          .9,
          c('linear regression', 
            'prediction with "predict"',
            "prediction with the coefficients"),
          col=c('red', 'blue', par('fg')),
          lty=c(2,2,1), lwd=c(1,5,1))
  title(main='Logistic regression with the "glm" function')
dev.off()

png(file="g1019.png", width=600, height=600)
  n <- 100
  x <- c(rnorm(n), 1+rnorm(n))
  y <- c(rep(0,n), rep(1,n))
  plot(y~x)
  # Brutal prediction
  m1 <- mean(x[y==0])
  m2 <- mean(x[y==1])
  m <- mean(c(m1,m2))
  if(m1<m2) a <- 0
  if(m1>m2) a <- 1
  if(m1==m2) a <- .5
  lines( c(min(x),m,m,max(x)),
         c(a,a,1-a,1-a),
         col='blue')  
  # Linear regression
  abline(lm(y~x), col='red')
  # Logistic regression
  xp <- seq(min(x),max(x),length=200)
  r <- glm(y~x, family=binomial)
  yp <- predict(r, data.frame(x=xp), type='response')
  lines(xp,yp, col='orange') 
  # Theoretical curve
  curve( dnorm(x,1,1)*.5/(dnorm(x,1,1)*.5+dnorm(x,0,1)*(1-.5)), add=T, lty=3, lwd=3)
  legend( .95*par('usr')[1]+.05*par('usr')[2],
          .9, #par('usr')[4],
          c('Brutal prediction', "Linear regression", "Logistic regression",
            "Theoretical curve"),
          col=c('blue','red','orange', par('fg')),
          lty=c(1,1,1,3),lwd=c(1,1,1,3))
  title(main="Comparing linear and logistic regression")
dev.off()

png(file="g1020.png", width=600, height=600)
  n <- 100
  x <- c(rnorm(n), 1+rnorm(n))
  y <- c(rep(0,n), rep(1,n))
  r <- glm(y~x, family=binomial)
  plot(r, which=1)
dev.off()

png(file="g1021.png", width=600, height=600)
  n <- 1000
  a <- -2
  b <- 1
  x <- runif(n, -4, 5)
  y <- exp(a*x+b + rnorm(n))
  y <- y/(1+y)
  y <- rbinom(n,1,y)
  plot(y~x)
dev.off()

png(file="g1022.png", width=600, height=600)
  boxplot(x~y, horizontal=T)
dev.off()

png(file="g1023.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  hist(x[y==1], probability=T, col='light blue')
  lines(density(x[y==1]),col='red',lwd=3)
  hist(x[y==0], probability=T, col='light blue')
  lines(density(x[y==0]),col='red',lwd=3)
  par(op)
dev.off()

png(file="g1024.png", width=600, height=600)
  rt <- glm(y~x, family=binomial)
  plot(rt, which=1)
dev.off()

png(file="g1025.png", width=600, height=600)
  plot(rt, which=2)
dev.off()

png(file="g1026.png", width=600, height=600)
  hist(rt$residuals, breaks=seq(min(rt$residuals),max(rt$residuals)+1,by=.5),
       xlim=c(-10,10),
       probability=T, col='light blue')
  points(density(rt$residuals, bw=.5), type='l', lwd=3, col='red')
dev.off()

png(file="g1027.png", width=600, height=600)
  plot(rt, which=3)
dev.off()

png(file="g1028.png", width=600, height=600)
  plot(rt, which=4)
dev.off()

png(file="g1029.png", width=600, height=600)
  # It is supposed not to be the same as in the linear situation.
  # Here, it seems to be the same...
  plot(hat(x), type='h')
dev.off()

png(file="g1030.png", width=600, height=600)
  plot(rt, which=4)
dev.off()

png(file="g1031.png", width=600, height=600)
  n <- 1000
  y <- factor(sample(0:1,n,replace=T))
  x <- rnorm(n)
  r <- glm(y~x,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F)
  par(op)
dev.off()

png(file="g1032.png", width=600, height=600)
  library(Design)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1033.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  a <- 1
  b <- -2
  p <- exp(a+b*x)/(1+exp(a+b*x))
  y <- factor(ifelse( runif(n)<p, 1, 0 ), levels=0:1)
  r <- glm(y~x,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F)
  par(op)
dev.off()

png(file="g1034.png", width=600, height=600)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1035.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  a <- 1
  b <- -2
  p <- exp(a+b*x)/(1+exp(a+b*x))
  y <- ifelse( runif(n)<p, 1, 0 )
  i <- runif(n)<.1
  y <- ifelse(i, 1-y, y)
  y <- factor(y, levels=0:1)
  col=c(par('fg'),'red')[1+as.numeric(i)]
  r <- glm(y~x,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F, col=col)
  par(op)
dev.off()

png(file="g1036.png", width=600, height=600)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1037.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  a <- 1
  b <- -2
  p <- exp(a+b*x)/(1+exp(a+b*x))
  y <- ifelse( runif(n)<p, 1, 0 )
  i <- runif(n)<.5 & abs(x)>1
  y <- ifelse(i, 1-y, y)
  y <- factor(y, levels=0:1)
  col=c(par('fg'),'red')[1+as.numeric(i)]
  r <- glm(y~x,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F, col=col)
  par(op)
dev.off()

png(file="g1038.png", width=600, height=600)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1039.png", width=600, height=600)
  n <- 1000
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  a <- 1
  b <- -2
  p <- exp(a+b*x1)/(1+exp(a+b*x1))
  y <- factor(ifelse( runif(n)<p, 1, 0 ), levels=0:1)
  r <- glm(y~x1+x2,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F)
  par(op)
dev.off()

png(file="g1040.png", width=600, height=600)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1041.png", width=600, height=600)
  n <- 1000
  x <- rnorm(n)
  a <- 1
  b1 <- -1
  b2 <- -2
  p <- 1/(1+exp(-(a+b1*x+b2*x^2)))
  y <- factor(ifelse( runif(n)<p, 1, 0 ), levels=0:1)
  r <- glm(y~x,family=binomial)
  op <- par(mfrow=c(2,2))
  plot(r,ask=F)
  par(op)
dev.off()

png(file="g1042.png", width=600, height=600)
  r <- lrm(y~x,y=T,x=T)
  P <- resid(r,"gof")['P']
  resid(r,"partial",pl=T)
  title(signif(P))
dev.off()

png(file="g1043.png", width=600, height=600)
  # NO!
  n <- 100
  x <- c( rnorm(n), 1+rnorm(n), 2.5+rnorm(n) )
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n) ))
  r <- glm(y~x, family=binomial)
  plot(as.numeric(y)-1~x)
  xp <- seq(-5,5,length=200)
  yp <- predict(r,data.frame(x=xp), type='response')
  lines(xp,yp)
dev.off()

png(file="g1044.png", width=600, height=600)
  n <- 100
  x <- c( rnorm(n), 10+rnorm(n), 25+rnorm(n), -7 + rnorm(n) )
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  y1 <- factor(c( rep(0,n), rep(0,n), rep(1,n), rep(1,n) ))
  y2 <- factor(c( rep(0,n), rep(1,n), rep(0,n), rep(1,n) ))
  r1 <- glm(y1~x, family=binomial)
  r2 <- glm(y2~x, family=binomial)
  xp <- seq(-50,50,length=500)
  y1p <- predict(r1,data.frame(x=xp), type='response')
  y2p <- predict(r2,data.frame(x=xp), type='response')

  plot(as.numeric(y)-1~x)
  lines(xp,y1p+2*y2p)
  lines(xp,y1p, col='red')
  lines(xp,y2p, col='blue')
dev.off()

png(file="g1045.png", width=600, height=600)
  plot(as.numeric(y1)-1~x)
  lines(xp,y1p, col='red')
dev.off()

png(file="g1046.png", width=600, height=600)
  n <- 100
  x <- c( rnorm(n), 10+rnorm(n), 25+rnorm(n), -7 + rnorm(n) )
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  y1 <- factor(c( rep(0,n), rep(1,n), rep(1,n), rep(1,n) ))
  y2 <- factor(c( rep(0,n), rep(0,n), rep(1,n), rep(1,n) ))
  y3 <- factor(c( rep(0,n), rep(0,n), rep(0,n), rep(1,n) ))
  r1 <- glm(y1~x, family=binomial)
  r2 <- glm(y2~x, family=binomial)
  r3 <- glm(y3~x, family=binomial)
  xp <- seq(-50,50,length=500)
  y1p <- predict(r1,data.frame(x=xp), type='response')
  y2p <- predict(r2,data.frame(x=xp), type='response')
  y3p <- predict(r3,data.frame(x=xp), type='response')

  plot(as.numeric(y)-1~x)
  lines(xp,y1p+y2p+y3p)
dev.off()

png(file="g1047.png", width=600, height=600)
  plot(as.numeric(y1)-1~x)
  lines(xp,y1p, col='red')
dev.off()

png(file="g1048.png", width=600, height=600)
  n <- 100
  x <- c( -7+rnorm(n), rnorm(n), 10+rnorm(n), 25+rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  y1 <- factor(c( rep(0,n), rep(1,n), rep(1,n), rep(1,n) ))
  y2 <- factor(c( rep(0,n), rep(0,n), rep(1,n), rep(1,n) ))
  y3 <- factor(c( rep(0,n), rep(0,n), rep(0,n), rep(1,n) ))
  r1 <- glm(y1~x, family=binomial)
  r2 <- glm(y2~x, family=binomial)
  r3 <- glm(y3~x, family=binomial)
  xp <- seq(-50,50,length=500)
  y1p <- predict(r1,data.frame(x=xp), type='response')
  y2p <- predict(r2,data.frame(x=xp), type='response')
  y3p <- predict(r3,data.frame(x=xp), type='response')

  plot(as.numeric(y)-1~x)
  lines(xp,y1p+y2p+y3p)
dev.off()

png(file="g1049.png", width=600, height=600)
  n <- 100
  x <- c( -.7+rnorm(n), rnorm(n), 1+rnorm(n), 2.5+rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  y1 <- factor(c( rep(0,n), rep(1,n), rep(1,n), rep(1,n) ))
  y2 <- factor(c( rep(0,n), rep(0,n), rep(1,n), rep(1,n) ))
  y3 <- factor(c( rep(0,n), rep(0,n), rep(0,n), rep(1,n) ))
  r1 <- glm(y1~x, family=binomial)
  r2 <- glm(y2~x, family=binomial)
  r3 <- glm(y3~x, family=binomial)
  xp <- seq(-5,5,length=500)
  y1p <- predict(r1,data.frame(x=xp), type='response')
  y2p <- predict(r2,data.frame(x=xp), type='response')
  y3p <- predict(r3,data.frame(x=xp), type='response')

  plot(as.numeric(y)-1~x)
  lines(xp,y1p+y2p+y3p)
  lines(xp,round(y1p+y2p+y3p, digits=0), col='red')
dev.off()

png(file="g1050.png", width=600, height=600)
  n <- 100
  x <- c( -.7+rnorm(n), rnorm(n), 1+rnorm(n), 2.5+rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  ordinal.regression.one <- function (y,x) {
    xp <- seq(min(x),max(x), length=100)
    yi <- matrix(nc=length(levels(y)), nr=length(y))
    ri <- list();
    ypi <- matrix(nc=length(levels(y)), nr=100)
    for (i in 1:length(levels(y))) {
      yi[,i] <- as.numeric(y) >= i
      ri[[i]] <- glm(yi[,i] ~ x, family=binomial)
      ypi[,i] <- predict(ri[[i]], data.frame(x=xp), type='response')
    }
    plot(as.numeric(y) ~ x)
    lines(xp, apply(ypi,1,sum), col='red', lwd=3)
  }
  ordinal.regression.one(y,x)
dev.off()

png(file="g1051.png", width=600, height=600)
  n <- 100
  v <- .2
  x <- c( -.7+v*rnorm(n), v*rnorm(n), 1+v*rnorm(n), 2.5+v*rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  ordinal.regression.one(y,x)
dev.off()

png(file="g1052.png", width=600, height=600)
  n <- 100
  x <- c( -.7+rnorm(n), rnorm(n), 1+rnorm(n), 2.5+rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  ordinal.regression.two <- function (y,x) {
    xp <- seq(min(x),max(x), length=100)
    yi <- list();
    ri <- list();
    ypi <- matrix(nc=length(levels(y)), nr=100)
    for (i in 1:length(levels(y))) {
      ya <- as.numeric(y)
      o <- ya >= i
      ya <- ya[o]
      xa <- x[o]
      yi[[i]] <- ya == i
      ri[[i]] <- glm(yi[[i]] ~ xa, family=binomial)
      ypi[,i] <- predict(ri[[i]], data.frame(xa=xp), type='response')
    }

    # The plot is trickier to draw than earlier
    plot(as.numeric(y) ~ x)
    p <- matrix(0, nc=length(levels(y)), nr=100)
    for (i in 1:length(levels(y))) {
      p[,i] = ypi[,i] * (1 - apply(p,1,sum))
    }
    for (i in 1:length(levels(y))) {
      p[,i] = p[,i]*i
    }
    lines(xp, apply(p,1,sum), col='red', lwd=3)
  }
  ordinal.regression.two(y,x)
dev.off()

png(file="g1053.png", width=600, height=600)
  n <- 100
  v <- .1
  x <- c( -.7+v*rnorm(n), v*rnorm(n), 1+v*rnorm(n), 2.5+v*rnorm(n))
  y <- factor(c( rep(0,n), rep(1,n), rep(2,n), rep(3,n) ))
  ordinal.regression.two(y,x)
dev.off()
detach.everything()

png(file="g1054.png", width=600, height=300)
  n <- 200
  x <- sample(0:1, n, replace=T)
  y <- 2*(x==0) +rnorm(n)
  plot( y ~ factor(x), 
        horizontal = TRUE, 
        xlab = 'y', 
        ylab = 'x',
        col = "pink" )
dev.off()

png(file="g1055.png", width=600, height=600)
  plot(density(y[x==0]), 
       lwd = 3, 
       xlim = c(-3,5), 
       ylim = c(0,.5), 
       col = 'blue',
       main = "Density in each group",
       xlab = "x")
  lines(density(y[x==1]), 
        lwd = 3, 
        col = 'red')
dev.off()

png(file="g1056.png", width=600, height=600)
  plot(y ~ x,
       main = "lm(y~x) when x is qualitative")
  abline(lm(y ~ x), 
         col = 'red')
dev.off()

png(file="g1057.png", width=600, height=600)
  n <- 200
  x1 <- sample(0:1, n, replace=T)
  x2 <- rnorm(n)
  y <- 1 - 2*(x1==0) + x2 +rnorm(n)
  plot( y ~ x2, 
        col = c(par('fg'), 'red')[1+x1],
        main = "y ~ (1 | x1) + x2")
  rc <- lm( y ~ x1+x2 )$coef
  abline(rc[1], rc[3])
  abline(rc[1]+rc[2], rc[3], 
         col = 'red')
  legend(par('usr')[2], par('usr')[3], # Bottom right
         xjust = 1, yjust = 0,
         c("x1 = 0", "x1 = 1"),
         col = c("red", par("fg")),
         lty = 1,
         lwd = 3)
dev.off()

png(file="g1058.png", width=600, height=600)
  x1 <- ifelse(x1==1, "x1 = 1", "x1 = 0")
  library(lattice)
  xyplot( y ~ x2 | x1,
          panel = function (x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.lmline(x, y)
          },
          main = "y ~ (1 | x1) + x2"
         )
dev.off()

png(file="g1059.png", width=600, height=600)
  n <- 200
  x1 <- sample(0:1, n, replace=T)
  x2 <- rnorm(n)
  y <- 1 + x2 - (x1==1)*(2+3*x2) + rnorm(n)
  plot( y ~ x2, 
        col = c(par('fg'), 'red')[1+x1],
        main = "y ~ x2 + (x2 | x1)" )
  # With no interaction term (dotted lines)
  rc <- lm( y ~ x1 + x2 )$coef
  abline(rc[1], rc[3], 
         lty = 3)
  abline(rc[1]+rc[2], rc[3], 
         col = 'red', lty = 3)
  # with
  rc <- lm( y ~ x1 + x2 + x1:x2 )$coef
  abline(rc[1], rc[3])
  abline(rc[1]+rc[2], rc[3]+rc[4], 
         col='red')
dev.off()

png(file="g1060.png", width=600, height=600)
  x1 <- ifelse(x1==1, "x1 = 1", "x1 = 0")
  xyplot( y ~ x2 | x1,
          panel = function (x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.lmline(x, y)
          },
          main = "y ~ x2 + (x2 | x1)"
         )
dev.off()

png(file="g1061.png", width=600, height=600)
  # One might also want to compare the intercepts
  xyplot( y ~ x2 | x1,
          panel = function (x, y, ...) {
            panel.xyplot(x, y, ...)
            panel.lmline(x, y)
            panel.abline(v=0, h=0, lty=3)
          },
          main = "y ~ x2 + (x2 | x1)"
         )
dev.off()

png(file="g1062.png", width=600, height=300)
  n <- 1000
  l <- 0:1
  v <- .5
  x1 <- factor( sample(l, n, replace=T), levels=l )
  x2 <- factor( sample(l, n, replace=T), levels=l )
  y <- ifelse( x1==0, 
         ifelse( x2==0, 1+v*rnorm(n),   v*rnorm(n) ),
         ifelse( x2==0, v*rnorm(n),   1+v*rnorm(n) )
       )
  boxplot( y ~ x1, 
           horizontal = TRUE,
           col = "pink",
           main = "y does not seem to depend on x1",
           xlab = "y", 
           ylab = "x1" )  
dev.off()

png(file="g1063.png", width=600, height=300)
  boxplot( y ~ x2,
           horizontal = TRUE,
           col = "pink",
           main = "y does not seem to depend on x2",
           xlab = "y", 
           ylab = "x2" )  
dev.off()

png(file="g1064.png", width=600, height=600)
  plot( as.numeric(x1) - 1 + .2*rnorm(n), 
        y, 
        col = c("red", "blue")[as.numeric(x2)], 
        xlab = "x1 (jittered)",
        main = "Interactions between x1 and x2")
  legend(par('usr')[2], par('usr')[3], # Bottom right
         xjust = 1, yjust = 0,
         c("x2 = 0", "x2 = 1"),
         col = c("red", "blue"),
         lty = 1,
         lwd = 3)
dev.off()

png(file="g1065.png", width=600, height=400)
  n <- 200
  x <- sample(0:3, n, replace=T)
  y <- 2*(x==0) + 5*(x==2) -2*(x==3) + rnorm(n)
  plot( y ~ factor(x), 
        horizontal = TRUE, 
        col = "pink",
        xlab = 'y', 
        ylab = 'x',
        main = "x now has four values" )
dev.off()

png(file="g1066.png", width=600, height=600)
  n <- 1000
  x1 <- sample(1:4, n, replace=T)
  x2 <- runif(n, min=-1, max=1)
  y <- (x1==1) * (2 - 2*x2) +
       (x1==2) * (1 + x2) +
       (x1==3) * (-.5 + .5*x2) +
       (x1==4) * (-1) +
       rnorm(n)
  cols <- c(par('fg'), 'red', 'blue', 'orange')
  plot( y ~ x2, 
        col = cols[x1],
        main = "y ~ x2 + (x2 | x1)" )
  x1 <- factor(x1)
  r <- lm( y ~ x1*x2 )
  xx <- seq( min(x2), max(x2), length = 100 )
  for (i in levels(x1)) {
    yy <- predict(r, data.frame(
      x1 = rep(i, length(xx)), 
      x2 = xx
    ))
    lines(xx, yy, 
          col = cols[as.numeric(i)], 
          lwd = 3)
  }
dev.off()

png(file="g1067.png", width=600, height=600)
  N <- 10
  a <- rnorm(N)
  b <- rnorm(N)
  c <- rnorm(N)
  d <- rnorm(N)
  df <- data.frame( 
    y = c(a,b,c,d), 
    x = factor(c(rep(1,N), rep(2,N), rep(3,N), rep(4,N))) 
  )
  plot( y ~ x, 
        data = df,
        main = "Our 4 samples" )
dev.off()

png(file="g1068.png", width=600, height=600)
  plot(density(df$y,bw=1), 
       xlim = c(-6,6), 
       ylim = c(0,.5), 
       type = 'l',
       main = "our 4 samples",
       xlab = "y")
  points(density(a,bw=1), type='l', col='red')
  points(density(b,bw=1), type='l', col='green')
  points(density(c,bw=1), type='l', col='blue')
  points(density(d,bw=1), type='l', col='orange')
dev.off()

png(file="g1069.png", width=600, height=600)
  curve( dnorm(x-2), from=-6, to=6, col='red',
         xlab = "", ylab = "" )
  curve( dnorm(x+2),    add = T, col = 'green')
  curve( dnorm(x+2+.2), add = T, col = 'blue')
  curve( dnorm(x+2-.3), add = T, col = 'orange')
  x <- c(2, -2, -2.2, -1.7)
  segments(x, c(0,0,0,0), 
           x, rep(dnorm(0),4), 
           col = c('red', 'green', 'blue', 'orange') )
  title("The means are significantly different")
dev.off()

png(file="g1070.png", width=600, height=600)
  s <- 3
  curve( dnorm(x-2, sd=s), from=-6, to=6, col='red',
         xlab = "", ylab = "" )
  curve( dnorm(x+2, sd=s),    add=T, col='green')
  curve( dnorm(x+2+.2, sd=s), add=T, col='blue')
  curve( dnorm(x+2-.3, sd=s), add=T, col='orange')
  x <- c(2, -2, -2.2, -1.7)
  segments( x, c(0,0,0,0), x, rep(dnorm(0, sd=s),4), 
            col=c('red', 'green', 'blue', 'orange') )
  title("The means are not significantly different")
dev.off()

png(file="g1071.png", width=600, height=600)
  N <- 5000  # Iterations
  n <- 5   # Sample size
  k <- 3     # Number of groups
  r <- rep(NA, N)
  for (i in 1:N) {
    l <- matrix(rnorm(n*k), nc=k)
    r[i] <- n * var(apply(l,2,mean)) / 
            mean(apply(l,2,var))
  }
  plot(sort(r), qf(ppoints(N),k-1,(k-1)*n),
    main = "When I know the distribution but not its parameters")
  abline(0,1,col="red")
dev.off()

png(file="g1072.png", width=600, height=600)
  N <- 5000  # Iterations
  n <- 5     # Sample size
  k <- 3     # Number of groups
  r <- rep(NA, N)
  for (i in 1:N) {
    l <- matrix(rnorm(n*k), nc=k)
    r[i] <- n * var(apply(l,2,mean)) / 
            mean(apply(l,2,var))
  }
  plot(sort(r), qf(ppoints(N),k-1,k*(n-1)),
    main = "When I know the distribution but not its parameters")
  abline(0,1,col="red")
dev.off()

png(file="g1073.png", width=600, height=600)
  curve( df(x, 2, 2), from=0, to=4, ylim=c(0,1),
         xlab = "", ylab = "", main = "" )
  curve( df(x, 2, 10), add=T, col='red' )
  curve( df(x, 4, 2),  add=T, col='green' )
  curve( df(x, 4, 6),  add=T, col='green' )
  curve( df(x, 4, 10), add=T, col='green' )
  curve( df(x, 4, 20), add=T, col='green' )
  curve( df(x, 6, 2),  add=T, col='blue' )
  curve( df(x, 6, 6),  add=T, col='blue' )
  curve( df(x, 6, 10), add=T, col='blue' )
  curve( df(x, 6, 20), add=T, col='blue' )
dev.off()

png(file="g1074.png", width=600, height=600)
  n <- 30
  x <- sample(LETTERS[1:3],n,replace=T, p=c(3,2,1)/6)
  x <- factor(x)
  y <- rnorm(n)  
  plot(y ~ x, 
       col = 'pink',
       xlab = "", ylab = "",
       main = "Simple anova: y ~ x")
dev.off()

png(file="g1075.png", width=600, height=600)
  n <- 100
  x1 <- sample(LETTERS[1:3],n,replace=T,p=c(3,2,1)/6)
  x1 <- factor(x1)
  x2 <- sample(LETTERS[1:2],n,replace=T,p=c(3,1)/4)
  x2 <- factor(x2)
  y <- rnorm(n)
  i <- which(x1=='A' & x2=='B')
  y[i] <- rnorm(length(i),.5)
  library(lattice)
  bwplot( ~ y | x1 * x2, 
         layout = c(1,6),
         main = "Double anova: y ~ x1 + x2")
dev.off()

png(file="g1076.png", width=600, height=600)
  x1 <- gl(2,2,4,c(0,1))
  x2 <- gl(2,1,4,c(0,1))
  n <- function (x) {
    as.numeric(as.vector(x))
  }
  y1 <- n(x1) + n(x2)
  interaction.plot(
    x1, x2, y1, 
    main = "Interaction plot: No interaction"
  )
dev.off()

png(file="g1077.png", width=600, height=600)
  y2 <- n(x1) + n(x2) - 2*n(x1)*n(x2)
  interaction.plot(
    x1, x2, y2, 
    main = "Interaction plot: Interaction"
  )
dev.off()

png(file="g1078.png", width=600, height=600)
  n <- 2000    # Number of experiments
  k <- 20      # Number of subjects
  l <- 4       # Number of groups
  kl <- sample(1:l, k, replace=T)  # Group of each subject
  x1 <- sample(1:k, n, replace=T)  
  x2 <- kl[x1]
  A <- rnorm(1,sd=4)
  B <- rnorm(k,sd=4)
  C <- rnorm(l,sd=4)
  y <- A + B[x1] + C[x2] + rnorm(n)
  x1 <- factor(x1)
  x2 <- factor(x2)
  op <- par(mfrow=c(1,2))
  plot(y~x1, col='pink')
  plot(y~x2, col='pink')
  par(op)
  mtext("Hierarchical anova", line=1.5, font=2, cex=1.2)
dev.off()

png(file="g1079.png", width=600, height=600)
  # If the data were real, we wouldn't know kl.
  # One may recover it that way.
  kl <- tapply(x2, 
               x1, 
               function (x) { 
                 a <- table(x)
                 names(a)[which(a==max(a))[1]] 
               })
  kl <- factor(kl, levels=levels(x2))
  plot( y ~ x1, col = rainbow(l)[kl],
        main = "Hierarchical anova")
dev.off()

png(file="g1080.png", width=600, height=600)
  x1 <- factor(x1, levels=order(kl))
  kl <- sort(kl)
  plot( y ~ x1, col = rainbow(l)[kl],
        main = "Hierarchical anova")
  legend(par('usr')[2], par('usr')[3], # Bottom right
         xjust = 1, yjust = 0,
         1:l,   # Is this right?
         col = rainbow(l),
         lty = 1,
         lwd = 3)
dev.off()

png(file="g1081.png", width=600, height=600)
  n <- 1000
  k <- 9
  subject <- factor(sample(1:k,n,replace=T), levels=1:k)
  a <- rnorm(k,sd=4)
  b <- rnorm(1)
  x <- rnorm(n)
  y <- a[subject] + b * x + rnorm(n)
  plot(y~x, main="Fixed effects")
  abline(lm(y~x))
dev.off()

png(file="g1082.png", width=600, height=600)
  col <- rainbow(k)
  plot(y~x, main="Random effects")
  for (i in 1:k) {
    points(y[subject==i] ~ x[subject==i], col=col[i])
    abline( lm(y~x, subset=subject==i), col=col[i] )
  }
dev.off()

png(file="g1083.png", width=600, height=600)
  library(lattice)
  xyplot(y~x|subject, main = "Random effects")
dev.off()

png(file="g1084.png", width=600, height=600)
  n <- 30
  k <- 10
  subject <- gl(k,n/k,n)
  x <- gl(n/k,1,n)
  y <- rnorm(k)
  y <- rnorm(n, y[subject])
  y[ x==2 ] <- y[ x==2 ] + 1

  # I do not see anything
  plot(y ~ x, 
       col = 'pink',
       main = "Repeated measures: y ~ x | subject")
dev.off()

png(file="g1085.png", width=600, height=600)
  interaction.plot(
    x, subject, y,
    main = "Repeated measures: y ~ x | subject"
  )
dev.off()

png(file="g1086.png", width=600, height=600)
  n <- 50
  k <- 10
  subjet <- gl(k,n/k,n)
  group <- sample(1:3, k, replace=T)[subjet]
  group <- factor(group)
  x <- gl(n/k,1,n)
  y <- rnorm(n)
  y[x==1] <- y[x==1] + 1
  plot(y~x, col='pink',
       main = "y ~ x | subject/group")
dev.off()

png(file="g1087.png", width=600, height=600)
  n <- 10
  subject <- gl(n,2)
  traitement <- gl(2,1,2*n, 0:1)
  duree <- ifelse( traitement==1, rpois(2*n,2), rpois(2*n,3) )
  sans <- duree[2*(1:n)-1]
  avec <- duree[2*(1:n)]
  plot(sans+rnorm(n), avec+rnorm(n))
  abline(0,1)
dev.off()

png(file="g1088.png", width=600, height=600)
  # There are 12 cities
  n.cities <- 12

  # The area of those cities (more reasonably, the logarithm
  # of their areas) are gaussian, independant variables.
  area.moyenne <- 5
  area.sd <- 1
  area <- rnorm(n.cities, area.moyenne, area.sd)

  a <- rnorm(n.cities)
  b <- rnorm(n.cities)

  # 200 inhabitants sampled in each city
  n.inhabitants <- 20
  city <- rep(1:n.cities, each=n.inhabitants)

  # The age are independant gaussian variables, mean=40, sd=10
  # We could have chosen a different distribution for each city. 
  # (either randomly, or depending on their area or population).

  age <- rnorm(n.cities*n.inhabitants, 40, 10)

  # The income (the variable we try to explain) is a function of the
  # age, but the coefficients depend on the city
  # Here, the coefficients are taken at random, but they could 
  # depend on the city area or population.
  # Here, the coefficients are independant -- this is rarely the case
  a <- rnorm(n.cities, 20000, sd=2000)
  b <- rnorm(n.cities, sd=20)
  income <- 200*area[city] + a[city] + b[city]*age + 
            rnorm(n.cities*n.inhabitants, sd=200)

  plot(income ~ age, col=rainbow(n.cities)[city], pch=16)
dev.off()

png(file="g1089.png", width=600, height=600)
  library(lattice)
  xyplot(income ~ age | city)
dev.off()

png(file="g1090.png", width=600, height=600)
  a <- rnorm(n.cities, 20000, sd=2000)
  b <- rep(rnorm(1, sd=20), n.cities)
  income <- 200*area[city] + a[city] + b[city]*age + 
            rnorm(n.cities*n.inhabitants, sd=200)
  xyplot(income ~ age | city)
dev.off()

png(file="g1091.png", width=600, height=600)
  a <- rep( rnorm(1, 20000, sd=2000), n.cities )
  b <- rnorm(n.cities, sd=20)
  income <- 200*area[city] + a[city] + b[city]*age + 
            rnorm(n.cities*n.inhabitants, sd=200)
  xyplot(income ~ age | city)
dev.off()

png(file="g1092.png", width=600, height=600)
  a <- rep( rnorm(1, 20000, sd=2000), n.cities )
  b <- rep(rnorm(1, sd=20), n.cities)
  income <- 200*area[city] + a[city] + b[city]*age + 
            rnorm(n.cities*n.inhabitants, sd=200)
  xyplot(income ~ age | city)
dev.off()

png(file="g1093.png", width=600, height=600)
  n.cities <- 12
  area <- rnorm(n.cities, 5, 1)
  a <- rnorm(n.cities)
  b <- rnorm(n.cities)
  n.inhabitants <- 20
  city <- rep(1:n.cities, each=n.inhabitants)
  age <- rnorm(n.cities*n.inhabitants, 40, 10)
  a <- rnorm(n.cities, 20000, sd=2000)
  b <- rep(rnorm(1, sd=20), n.cities)
  income <- 200*area[city] + a[city] + b[city]*age +
  rnorm(n.cities*n.inhabitants, sd=200)
  xyplot(income ~ age | city)
dev.off()

png(file="g1094.png", width=600, height=600)
  library(nlme)
  d <- data.frame(income, age, city, area=area[city])
  r <- lmList(income ~ age | city, data=d)
  plot(intervals(r))
dev.off()
detach.everything()

png(file="g1095.png", width=600, height=600)
  # No random effects
  N <- 200
  k <- 4
  x <- rnorm(N)
  g <- sample(1:k, N, replace=TRUE)
  a <- rep(runif(1,-1,1), k)
  b <- rep(runif(1,-1,1), k)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  plot(x, y, col=rainbow(k)[g],
       main="(Non-)mixed model: no random effects")
  for (i in 1:k) {
    abline(lm(y[g==i] ~ x[g==i]), 
           col=rainbow(k)[i], lwd=2)
  }
dev.off()

png(file="g1096.png", width=600, height=600)
  # Random intercept
  a <- runif(k,-1,1)
  b <- rep(runif(1,-1,1), k)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  plot(x, y, col=rainbow(k)[g],
       main="Mixed model: random intercept")
  for (i in 1:k) {
    abline(lm(y[g==i] ~ x[g==i]), 
           col=rainbow(k)[i], lwd=2)
  }
dev.off()

png(file="g1097.png", width=600, height=600)
  # Random slope
  a <- rep(runif(1,-1,1), k)
  b <- runif(k,-1,1)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  plot(x, y, col=rainbow(k)[g],
       main="Mixed model: random slope")
  for (i in 1:k) {
    abline(lm(y[g==i] ~ x[g==i]), 
           col=rainbow(k)[i], lwd=2)
  }
dev.off()

png(file="g1098.png", width=600, height=600)
  # Random intercept and slope
  a <- runif(k,-1,1)
  b <- runif(k,-1,1)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  plot(x, y, col=rainbow(k)[g],
       main="Mixed model: random intercept and slope")
  for (i in 1:k) {
    abline(lm(y[g==i] ~ x[g==i]), 
           col=rainbow(k)[i], lwd=2)
  }
dev.off()

png(file="g1099.png", width=600, height=600)
  library(lme4)
  xyplot(Reaction ~ Days, groups = Subject, 
         data = sleepstudy, 
         type = "l")
dev.off()

png(file="g1100.png", width=600, height=600)
  xyplot(Reaction ~ Days | Subject, 
         data = sleepstudy, 
         type = "l")
dev.off()

png(file="g1101.png", width=600, height=600)
  # A regression in each group
  r <- lmList(Reaction ~ Days | Subject, data = sleepstudy)
  # The regression coefficients
  do.call("rbind", lapply(r, function (x) { x$coef }))
  # To check if they are different, let us compute 
  # their 95% confidence interval.
  lmListIntervals <- function (r, level=.95) {
    s <- array(numeric(0), 
               dim=c(length(r), 
               dim(confint(r[[1]]))))
    dimnames(s)[2:3] <- dimnames(confint(r[[1]]))
    dimnames(s)[[1]] <- names(r)
    names(dimnames(s)) <- c("Group", "Variable", "Interval")
    for (i in 1:length(r)) {
      s[i,,] <- confint(r[[i]], level=level)
    }
    s
  }
  s <- lmListIntervals(r)
  aperm(s, c(1,3,2))

  lmListIntervalsPlot <- function (s, i=1) {
    plot.new()
    plot.window( xlim = range(s[,i,]), ylim=c(1,length(r)) )
    segments( s[,i,1], 1:length(r), 
              s[,i,2], 1:length(r) )
    axis(1)
    axis(2)
    box()
    level <- diff(as.numeric(gsub(" .*", "", dimnames(s)[[3]])))
    title(xlab = dimnames(s)[[2]][i],
          ylab = "Subject",
          main = paste(level, "% confidence intervals", sep=""))
  }
  lmListIntervalsPlot(s,1)
dev.off()

png(file="g1102.png", width=600, height=600)
  lmListIntervalsPlots <- function (s) {
    k <- dim(s)[2]
    op <- par(mfrow=c(1,k))
    for (i in 1:k) {
      lmListIntervalsPlot(s,i)
    }
    par(op)
  }
  lmListIntervalsPlots(s)
dev.off()

png(file="g1103.png", width=600, height=600)
  # No random effects
  N <- 200
  k <- 4
  x <- rnorm(N)
  g <- sample(1:k, N, replace=TRUE)
  a <- rep(runif(1,-1,1), k)
  b <- rep(runif(1,-1,1), k)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  d <- data.frame(x=x, y=y, g=as.factor(g))
  lmListIntervalsPlots(lmListIntervals(
    lmList( y ~ x | g, data = d )
  ))
  mtext("No random effects", 
        side=3, line=3, font=2, cex=1.2)
dev.off()

png(file="g1104.png", width=600, height=600)
  # Random intercept
  a <- runif(k, -1, 1)
  b <- rep(runif(1,-1,1), k)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  d <- data.frame(x=x, y=y, g=as.factor(g))
  lmListIntervalsPlots(lmListIntervals(
    lmList( y ~ x | g, data = d )
  ))
  mtext("Random intercept", 
        side=3, line=3, font=2, cex=1.2)
dev.off()

png(file="g1105.png", width=600, height=600)
  # Random slope
  a <- rep(runif(1,-1,1), k)
  b <- runif(k, -1, 1)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  d <- data.frame(x=x, y=y, g=as.factor(g))
  lmListIntervalsPlots(lmListIntervals(
    lmList( y ~ x | g, data = d )
  ))
  mtext("Random slope", 
        side=3, line=3, font=2, cex=1.2)
dev.off()

png(file="g1106.png", width=600, height=600)
  # Random intercept and slope
  a <- runif(k, -1, 1)
  b <- runif(k, -1, 1)
  y <- a[g] + b[g] * x + .2 * rnorm(N)
  d <- data.frame(x=x, y=y, g=as.factor(g))
  lmListIntervalsPlots(lmListIntervals(
    lmList( y ~ x | g, data = d )
  ))
  mtext("Random intercept and slope", 
        side=3, line=3, font=2, cex=1.2)
dev.off()

png(file="g1107.png", width=600, height=600)
  data(Orthodont, package="nlme")
  xyplot( distance ~ age | Sex, group = Subject, 
          data = Orthodont, type="l" )
dev.off()

png(file="g1108.png", width=600, height=600)
  data(BodyWeight, package="nlme")
  xyplot(weight ~ Time | Rat, 
         data = BodyWeight, 
         type = "l",  aspect = 8)
dev.off()

png(file="g1109.png", width=600, height=600)
  xyplot(weight ~ Time | Diet, group = Rat, 
         data = BodyWeight, 
         type = "l",  aspect = 3)
dev.off()

png(file="g1110.png", width=600, height=600)
  data(Cefamandole, package="nlme")
  xyplot( conc ~ Time | Subject, 
          data = Cefamandole,
          type = "l" )
dev.off()

png(file="g1111.png", width=600, height=600)
  data(Dialyzer, package="nlme")
  xyplot( rate ~ pressure | QB, group = Subject, 
          data = Dialyzer, type = "l")
dev.off()

png(file="g1112.png", width=600, height=600)
  data(Earthquake, package="nlme")
  xyplot(accel ~ distance | Quake, data=Earthquake )
dev.off()

png(file="g1113.png", width=600, height=600)
  plus <- function (x) {
    ifelse( x >= 0, x, 0 )
  }
  op <- par(mfrow=c(2,2), mar=c(3,3,2,2), oma=c(0,0,2,0))
  plot.basis.function <- function (a) {
    curve( plus(x - a),  
           xlim = c(-1, 1), 
           ylim = c(-1, 1), 
           lwd = 3, 
           col = "blue" )
    abline(h=0, v=0, lty=3)
  }
  plot.basis.function(-1)
  plot.basis.function(-.5)
  plot.basis.function(0)
  plot.basis.function(.5)
  par(op)
  mtext("Linear splines: basis functions", 
        font = 2, line = 2, cex = 1.5)
dev.off()

png(file="g1114.png", width=600, height=600)
  N <- 100
  x <- runif(N, -1, 1)
  y <- sin(2*pi*x) + .5*rnorm(N)
  plot(x,y)
  z <- apply(as.matrix(seq(-1,1,by=.2)), 1, function (a)  { plus(x-a) }) 
  z <- z[,-1] # This is already captured by the intercept
  y.pred <- predict(lm(y ~ z))
  o <- order(x)
  lines( x[o], y.pred[o], col="red", lwd=3 )
  title(main="Broken line regression")
dev.off()

png(file="g1115.png", width=600, height=600)
  broken.line.regression <- function (
      x, y, 
      knots = max(2,ceiling(length(x)/10)), 
      method = c("equispaced", "quantiles")
  ) {
    method <- match.arg(method)
    if (length(knots) == 1) {
      if (method == "quantiles") {
        knots <- quantile(x, knots+2)[ - c(1, knots+2) ]
      } else {
        knots <- seq(min(x), max(x), length=knots+2) [ - c(1, knots+2) ]
      }
    }
    z <- apply(as.matrix(knots), 1, function (a)  { plus(x-a) })
    lm( y ~ z )
  }
  plot.broken.line.regression <- function (x, y, ...) {
    o <- order(x)
    plot(x, y)
    r <- broken.line.regression(x, y, ...)
    lines( x[o], predict(r)[o], col="red", lwd=3 )
    invisible(r)
  }
  plot.broken.line.regression(x,y, knots = 3)
dev.off()

png(file="g1116.png", width=600, height=600)
  plot.broken.line.regression(x,y, knots = 3)
  title(main="Smoothing with linear splines, 3 knots")
dev.off()

png(file="g1117.png", width=600, height=600)
  plot.broken.line.regression(x,y, knots = 10)
  title(main="Smoothing with linear splines, 3 knots")
dev.off()

png(file="g1118.png", width=600, height=600)
  plot.broken.line.regression(x,y, knots = 30)
  title(main="Smoothing with linear splines, 3 knots")
dev.off()

png(file="g1119.png", width=600, height=600)
  N <- 1000
  n <- 20   # Sample size
  k <- 4    # number of samples
  res <- matrix(NA, nr=N, nc=2)
  for (a in 1:N) {
    x <- factor(rep(1:k,1,each=n))
    y <- rnorm(k*n)  
    res[a,1] <- summary(aov(y~x))[[1]][1,5]
    p <- 1
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
      p <- min(c( p, t.test(y[x==i],y[x==j])$p.value ))
      }
    }
    res[a,2] <- p
  }
  plot(res, xlab="p-value (anova)", ylab="p-value (multiple Student's T tests)")
dev.off()

png(file="g1120.png", width=600, height=600)
  plot(sort(res[,1]), 
       main = "p-values (anova)")
  abline(0, 1/N, col = 'red')
dev.off()

png(file="g1121.png", width=600, height=600)
  plot(sort(res[,2]), 
       main = "p-values (multiple Student's T tests)")
  abline(0, 1/N, col = 'red')
dev.off()

png(file="g1122.png", width=600, height=600)
  plot(sort(k*(k-1)/2*res[,2]), 
       ylim = c(0,1), 
       main = "p-values (Bonferronni)")
  abline(0, 1/N, col = 'red')
dev.off()

png(file="g1123.png", width=600, height=600)
  plot(sort(1-(1-res[,2])^(k*(k-1)/2)), 
       main = "p-values (*********)")
  abline(0, 1/N, col = 'red')
dev.off()

png(file="g1124.png", width=600, height=600)
  # Very, very long...
  N <- 1000
  n <- 20   # Sample size
  k <- 20   # Number of samples
  res <- matrix(NA, nr=N, nc=2)
  for (a in 1:N) {
    x <- factor(rep(1:k,1,each=n))
    y <- rnorm(k*n)  
    res[a,1] <- summary(aov(y~x))[[1]][1,5]
    p <- 1
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
      p <- min(c( p, t.test(y[x==i],y[x==j])$p.value ))
      }
    }
    res[a,2] <- p
  }

  op <- par(mfrow=c(2,2))
  plot(res, 
       xlab = "p-value (anova)", 
       ylab = "p-value (multiple Student tests)")

  # plot(sort(res[,1]), main="p-values (anova)")
  # abline(0,1/N, col='red')

  plot(sort(res[,2]), 
       ylim = c(0,1), 
       main = "p-values (multiple Student tests)")
  abline(0, 1/N, col = 'red')
  abline(h = .05, lty = 3)

  plot(sort(k*(k-1)/2*res[,2]), 
       ylim = c(0,1), 
       main = "p-values (Bonferronni)")
  abline(0, 1/N, col = 'red')
  abline(h = .05, lty = 3)

  plot(sort(1-(1-res[,2])^(k*(k-1)/2)), 
       main = "p-values (********)")
  abline(0, 1/N, col = 'red')
  abline(h = .05, lty = 3)
  par(op)
dev.off()

png(file="g1125.png", width=600, height=600)
  library(nlme)
  library(lattice)
  data(Orthodont)
  formula(Orthodont)
  plot(Orthodont)
dev.off()

png(file="g1126.png", width=600, height=600)
  plot(Orthodont, outer=T)
dev.off()

png(file="g1127.png", width=600, height=600)
  data(Machines)
  bwplot(score~Machine|Worker, data=Machines)
dev.off()

png(file="g1128.png", width=600, height=600)
  data(CO2)
  plot(CO2)
dev.off()

png(file="g1129.png", width=600, height=600)
  plot(CO2, outer=T)
dev.off()

png(file="g1130.png", width=600, height=600)
  data(InsectSprays)
  y <- InsectSprays$count
  x <- InsectSprays$spray
  boxplot(y~x, col='pink')
dev.off()

png(file="g1131.png", width=600, height=600)
  library(lattice)
  histogram(~ y | x)
dev.off()

png(file="g1132.png", width=600, height=600)
  qqmath(~ y | x)
  # TODO: qqmathline
dev.off()

png(file="g1133.png", width=600, height=600)
  n <- length(levels(x))
  res <- matrix(NA, nr=n, nc=n)
  rownames(res) <- levels(x)
  colnames(res) <- levels(x)
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      res[i,j] <- t.test( y[ x == levels(x)[i] ], y[ x == levels(x)[j] ] )$p.value
    }
  }
  #res <- res * n*(n-1)/2; res <- ifelse(res>1,1,res)
  res <- 1 - (1-res)^n
  image(res, col=topo.colors(255))  # TODO: unreadable

                                    # Write a function that plots
                                    # distance matrices, correlation
                                    # matrices (with the corresponding
                                    # p-values).
                                    # With the name of the
                                    # rows/columns, with a legend for
                                    # the colors. 
  round(res,3)
dev.off()

png(file="g1134.png", width=600, height=600)
  d <- 1-res
  d <- ifelse(is.na(d), t(d), d)
  diag(d) <- 0
  p <- isoMDS(d)$points
  plot(p, pch=16)
  text(p, levels(x), pos=c(1,2,1,1,1,3))
dev.off()
detach.everything()

png(file="g1135.png", width=600, height=300)
  plot(LakeHuron,
       ylab = "",
       main = "Level of Lake Huron")
dev.off()

png(file="g1136.png", width=600, height=600)
  x <- window(sunspots, start=1750, end=1800)
  plot(x,
       ylab = "",
       main = "Sunspot numbers")
dev.off()

png(file="g1137.png", width=600, height=600)
  plot(x, 
       type = 'p', 
       ylab = "",
       main = "Sunspot numbers")
  k <- 20
  lines( filter(x, rep(1/k,k)), 
         col = 'red', 
         lwd = 3 )
dev.off()

png(file="g1138.png", width=600, height=300)
  data(UKgas)
  plot.band <- function (x, ...) {
    plot(x, ...)
    a <- time(x)
    i1 <- floor(min(a))
    i2 <- ceiling(max(a))
    y1 <- par('usr')[3]
    y2 <- par('usr')[4]
    if( par("ylog") ){
      y1 <- 10^y1
      y2 <- 10^y2
    }
    for (i in seq(from=i1, to=i2-1, by=2)) {
      polygon( c(i,i+1,i+1,i), 
               c(y1,y1,y2,y2), 
               col = 'grey', 
               border = NA )
    }
    par(new=T)
    plot(x, ...)
  }
  plot.band(UKgas, 
            log = 'y', 
            ylab = "",
            main = "UK gas consumption")
dev.off()

png(file="g1139.png", width=600, height=300)
  x <- LakeHuron
  op <- par(mfrow = c(1,2),
            mar = c(5,4,1,2)+.1,
            oma = c(0,0,2,0))
  hist(x, 
       col = "light blue",
       xlab = "",
       main = "")
  qqnorm(x,
         main = "")
  qqline(x, 
         col = 'red')
  par(op)
  mtext("Lake Huron levels", 
        line = 2.5, 
        font = 2, 
        cex = 1.2)
dev.off()

png(file="g1140.png", width=600, height=300)
  x <- diff(LakeHuron)
  op <- par(mfrow = c(1,2),
            mar = c(5,4,1,2)+.1,
            oma = c(0,0,2,0))
  hist(x, 
       col = "light blue",
       xlab = "",
       main = "")
  qqnorm(x,
         main = "")
  qqline(x, 
         col = 'red')
  par(op)
  mtext("Lake Huron level increments", 
        line = 2.5, 
        font = 2, 
        cex = 1.2)
dev.off()

png(file="g1141.png", width=600, height=200)
  boxplot(x,
          horizontal = TRUE, 
          col = "pink",
          main = "Lake Huron levels")
dev.off()

png(file="g1142.png", width=600, height=300)
  plot(x,
       ylab = "",
       main = "Lake Huron levels")
dev.off()

png(file="g1143.png", width=600, height=600)
  n <- length(x)
  k <- 5
  m <- matrix(nr=n+k-1, nc=k)
  colnames(m) <- c("x[i]", "x[i-1]", "x[i-2]",
                   "x[i-3]", "x[i-4]")
  for (i in 1:k) {
    m[,i] <- c(rep(NA,i-1), x, rep(NA, k-i))
  }
  pairs(m, 
        gap = 0,
        lower.panel = panel.smooth,
        upper.panel = function (x,y) {
          panel.smooth(x,y)
          par(usr = c(0, 1, 0, 1))
          a <- cor(x,y, use='pairwise.complete.obs')
          text(.1,.9, 
               adj=c(0,1),
               round(a, digits=2),
               col='blue',
               cex=2*a)
        })
  title("Lake Huron levels: autocorrelations", 
        line = 3)
dev.off()

png(file="g1144.png", width=600, height=600)
  op <- par(mfrow = c(3,1),
            mar = c(2,4,1,2)+.1)
  acf(x,      xlab = "")
  pacf(x,     xlab = "")
  spectrum(x, xlab = "", main = "")
  par(op)
dev.off()

png(file="g1145.png", width=600, height=600)
  op <- par(mfrow = c(3,3),
            mar = .1 + c(0,0,0,0))

  n <- 100
  k <- 5
  N <- k*n
  x <- (1:N)/n
  y1 <- rnorm(N)
  plot(ts(y1), 
       xlab="", ylab="", main="", axes=F)
  box()

  y2 <- cumsum(rnorm(N))
  plot(ts(y2),
       xlab="", ylab="", main="", axes=F)
  box()

  y3 <- cumsum(rnorm(N))+rnorm(N)
  plot(ts(y3),
       xlab="", ylab="", main="", axes=F)
  box()

  y4 <- cumsum(cumsum(rnorm(N)))
  plot(ts(y4),
       xlab="", ylab="", main="", axes=F)
  box()

  y5 <- cumsum(cumsum(rnorm(N))+rnorm(N))+rnorm(N)
  plot(ts(y5),
       xlab="", ylab="", main="", axes=F)
  box()

  # With a trend
  y6 <- 1 - x + cumsum(rnorm(N)) + .2 * rnorm(N)
  plot(ts(y6),
       xlab="", ylab="", main="", axes=F)
  box()

  y7 <- 1 - x - .2*x^2 + cumsum(rnorm(N)) + 
        .2 * rnorm(N)
  plot(ts(y7),
       xlab="", ylab="", main="", axes=F)
  box()


  # With a seasonnal component
  y8 <- .3 + .5*cos(2*pi*x) - 1.2*sin(2*pi*x) +
        .6*cos(2*2*pi*x) + .2*sin(2*2*pi*x) +
        -.5*cos(3*2*pi*x) + .8*sin(3*2*pi*x)
  plot(ts(y8+ .2*rnorm(N)),
       xlab="", ylab="", main="", axes=F)
  box()
  lines(y8, type='l', lty=3, lwd=3, col='red')

  y9 <- y8 + cumsum(rnorm(N)) + .2*rnorm(N)
  plot(ts(y9),
       xlab="", ylab="", main="", axes=F)
  box()

  par(op)
dev.off()

png(file="g1146.png", width=600, height=400)
  my.acf <- function (
      x, 
      lag.max = ceiling(5*log(length(x)))
  ) {
    m <- matrix( 
      c( NA, 
         rep( c(rep(NA, lag.max-1), x),
              lag.max ),
         rep(NA,, lag.max-1)
       ),
      byrow=T,
      nr=lag.max)
    x0 <- m[1,]
    apply(m,1,cor, x0, use="complete")
  }
  n <- 200
  x <- rnorm(n)
  plot(my.acf(x), 
       xlab = "Lag",
       type = 'h')
  abline(h=0)
dev.off()

png(file="g1147.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  acf(x, main="ACF of white noise")
  x <- LakeHuron
  acf(x, main="ACF of a time series (Lake Huron)")
  par(op)
dev.off()

png(file="g1148.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  set.seed(1)
  x <- rnorm(100)
  # Default plot
  acf(x, main = "ACF with a distracting horizontal line")
  # Without the axis, with larger bars
  r <- acf(x, plot = FALSE)
  plot(r$lag, r$acf, 
       type = "h", lwd = 20, col = "grey",
       xlab = "lag", ylab = "autocorrelation",
       main = "Autocorrelation without the y=0 line")
  ci <- .95 
  clim <- qnorm( (1+ci) / 2 ) / sqrt(r$n.used)
  abline(h = c(-1,1) * clim, 
         lty = 2, col = "blue", lwd = 2)
dev.off()

png(file="g1149.png", width=600, height=600)
  op <- par(mfrow=c(3,3), mar=c(0,0,0,0))
  for (y in sample(list(y1,y2,y3,y4,y5,y6,y7,y8,y9))) {
    acf(y, 
        xlab="", ylab="", main="", axes=F)
    box(lwd=2)
  }
  par(op)
dev.off()

png(file="g1150.png", width=600, height=300)
  my.plot.ts <- function (x, main="") {
    op <- par(mar=c(2,2,4,2)+.1)
    layout( matrix(c(1,2),nr=1,nc=2), widths=c(3,1) )
    plot(x, xlab="", ylab="")
    abline(h=0, lty=3)
    title(main=main)
    hist(x, col="light blue", main='', ylab="", xlab="")
    par(op)
  }
  n <- 100
  x <- ts(rnorm(n))
  my.plot.ts(x, "Gaussian iid noise")
dev.off()

png(file="g1151.png", width=600, height=300)
  n <- 100
  x <- ts(runif(n,-1,1))
  my.plot.ts(x, "Non gaussian iid noise")
dev.off()

png(file="g1152.png", width=600, height=300)
  x <- ts(rnorm(100)^3)
  my.plot.ts(x, "Non gaussian iid noise")
dev.off()

png(file="g1153.png", width=600, height=300)
  n <- 100
  x <- rep(0,n)
  z <- rnorm(n)
  for (i in 2:n) {
    x[i] <- z[i] * sqrt( 1 + .5 * x[i-1]^2 )
  }
  my.plot.ts(x, "Non iid noise")
dev.off()

png(file="g1154.png", width=600, height=300)
  n <- 100
  x <- rep(.7, n)
  for (i in 2:n) {
    x[i] <- 4 * x[i-1] * ( 1 - x[i-1] )
  }
  my.plot.ts(x, "A deterministic time series")
dev.off()

png(file="g1155.png", width=600, height=600)
  n <- 1000
  tn <- cumsum(rexp(n))
  # A C^infinity function defined as a sum 
  # of gaussian densities
  f <- function (x) { 
    # If x is a single number: sum(dnorm(x-tn))
    apply( dnorm( outer(x,rep(1,length(tn))) - 
                  outer(rep(1,length(x)),tn) ),
           1,
           sum )
  }
  op <- par(mfrow=c(2,1))
  curve( 
    f(x), 
    xlim = c(1,500),
    n = 1000,
    main = "From far away, it looks random..." 
  )
  curve( 
    f(x), 
    xlim = c(1,10),
    n = 1000,
    main="...but it is not: it is a C^infinity function" 
  )
  par(op)
dev.off()

png(file="g1156.png", width=600, height=600)
  z <- rnorm(200)
  op <- par(mfrow=c(2,1), mar=c(5,4,2,2)+.1)
  plot(ts(z))
  acf(z, main = "")
  par(op)
dev.off()

png(file="g1157.png", width=600, height=600)
  x <- diff(co2)
  y <- diff(x,lag=12)
  op <- par(mfrow=c(2,1), mar=c(5,4,2,2)+.1)
  plot(ts(y))
  acf(y, main="")
  par(op)
dev.off()

png(file="g1158.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  plot.box.ljung <- function (
      z, 
      k = 15, 
      main = "p-value of the Ljung-Box test", 
      ylab = "p-value"
  ) {
    p <- rep(NA, k)
    for (i in 1:k) {
      p[i] <- Box.test(z, i, 
                       type = "Ljung-Box")$p.value
    }
    plot(p, 
         type = 'h', 
         ylim = c(0,1), 
         lwd = 3, 
         main = main, 
         ylab = ylab)
    abline(h = c(0,.05),
           lty = 3)
  }
  plot.box.ljung(z, main="Random data")
  plot.box.ljung(y, main="diff(diff(co2),lag=12)")
  par(op)
dev.off()

png(file="g1159.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  library(lmtest)
  plot(LakeHuron, 
       main = "Lake Huron")
  acf(
    LakeHuron, 
    main = paste(
      "Durbin-Watson: p =", 
      signif( dwtest( LakeHuron ~ 1 ) $ p.value, 3 )
    )
  )
  par(op)
dev.off()

png(file="g1160.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  op <- par(mfrow=c(2,1))
  x <- ts(x)
  plot(x, main="White noise", ylab="")
  acf(
    x, 
    main = paste(
      "Durbin-Watson: p =", 
      signif( dwtest( x ~ 1 ) $ p.value, 3)
    )
  )
  par(op)
dev.off()

png(file="g1161.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  op <- par(mfrow=c(2,1))
  y <- filter(x,.8,method="recursive")
  plot(y, main="AR(1)", ylab="")
  acf(
    y, 
    main = paste(
      "p =", 
       signif( dwtest( y ~ 1 ) $ p.value, 3 )
    )
  )
  par(op)
dev.off()

png(file="g1162.png", width=600, height=600)
  set.seed(1)
  n <- 200
  x <- rnorm(n)
  y <- filter(x, c(0,1), method="recursive")
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(
    y, 
    main = paste(
      "one-sided DW test: p =",
      signif( dwtest ( y ~ 1 ) $ p.value, 3 )
    )
  )
  acf( y, main="")
  pacf(y, main="")
  par(op)
dev.off()

png(file="g1163.png", width=600, height=600)
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  res <- dwtest( y ~ 1, alternative="two.sided")
  plot(
    y, 
    main = paste(
      "two-sided p =",
      signif( res$p.value, 3 )
    )
  )
  acf(y, main="")
  pacf(y, main="")
  par(op)
dev.off()

png(file="g1164.png", width=600, height=600)
  data(co2)
  r <- arima(
    co2, 
    order = c(0, 1, 1), 
    seasonal = list(order = c(0, 1, 1), period = 12)
  )
  tsdiag(r)
dev.off()

png(file="g1165.png", width=600, height=300)
  data(co2)
  plot(co2)
dev.off()

png(file="g1166.png", width=600, height=300)
  y <- as.vector(co2)
  x <- as.vector(time(co2))
  r <- lm( y ~ poly(x,1) + cos(2*pi*x) + sin(2*pi*x) )
  plot(y~x, type='l', xlab="time", ylab="co2")
  lines(predict(r)~x, lty=3, col='red', lwd=3)
dev.off()

png(file="g1167.png", width=600, height=600)
  plot( y-predict(r),
        main = "The residuals are not random yet",
        xlab = "Time", 
        ylab = "Residuals" )
dev.off()

png(file="g1168.png", width=600, height=300)
  r <- lm( y ~ poly(x,2) + cos(2*pi*x) + sin(2*pi*x) )
  plot(y~x, type='l', xlab="Time", ylab="co2")
  lines(predict(r)~x, lty=3, col='red', lwd=3)
dev.off()

png(file="g1169.png", width=600, height=600)
  plot( y-predict(r),
        main = "Better residuals -- but still not random",
        xlab = "Time",
        ylab = "Residuals" )
dev.off()

png(file="g1170.png", width=600, height=300)
  r <- lm( y ~ poly(x,2) + cos(2*pi*x) + sin(2*pi*x) 
                         + cos(4*pi*x) + sin(4*pi*x) )
  plot(y~x, type='l', xlab="Time", ylab="co2")
  lines(predict(r)~x, lty=3, col='red', lwd=3)
dev.off()

png(file="g1171.png", width=600, height=600)
  plot( y-predict(r), 
        type = 'l',
        xlab = "Time",
        ylab = "Residuals",
        main = "Are those residuals any better?" )
dev.off()

png(file="g1172.png", width=600, height=600)
  r1 <- lm( y ~ poly(x,2) + 
                cos(2*pi*x) + 
                sin(2*pi*x) )
  r2 <- lm( y ~ poly(x,2) + 
                cos(2*pi*x) +
                sin(2*pi*x) + 
                cos(4*pi*x) + 
                sin(4*pi*x) )
  op <- par(mfrow=c(2,1))
  acf(y - predict(r1))
  acf(y - predict(r2))
  par(op)
dev.off()

png(file="g1173.png", width=600, height=600)
  m <- tapply(co2, gl(12,1,length(co2)), mean)
  m <- rep(m, ceiling(length(co2)/12)) [1:length(co2)]
  m <- ts(m, start=start(co2), frequency=frequency(co2))
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2))
  plot(co2)
  plot(m, ylab = "Periodic component")
  plot(co2-m, ylab = "Withou the periodic component")
  r <- lm(co2-m ~ poly(as.vector(time(m)),2))
  lines(predict(r) ~ as.vector(time(m)), col='red')
  par(op)
dev.off()

png(file="g1174.png", width=600, height=800)
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2)+.1)
  plot(r$res, type = "l")
  acf(r$res,  main="")
  pacf(r$res, main="")
  spectrum(r$res, col=par('fg'), main="")
  abline(v=1:6, lty=3)
  par(op)
dev.off()

png(file="g1175.png", width=600, height=600)
  k <- 12  
  m <- matrix( as.vector(diag(k)), 
               nr = length(co2), 
               nc = k, 
               byrow = TRUE )
  m <- cbind(m, poly(as.vector(time(co2)),2))
  r <- lm(co2~m-1)
  summary(r)
  b <- r$coef
  y1 <- m[,1:k] %*% b[1:k]
  y1 <- ts(y1, 
           start=start(co2), 
           frequency=frequency(co2))
  y2 <- m[,k+1:2] %*% b[k+1:2]
  y2 <- ts(y2, 
           start=start(co2), 
           frequency=frequency(co2))
  res <- co2 - y1 - y2
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(co2)
  lines(y2+mean(b[1:k]) ~ as.vector(time(co2)), 
        col='red')
  plot(y1)
  plot(res)
  par(op)
dev.off()

png(file="g1176.png", width=600, height=600)
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  acf(res,  main="")
  pacf(res, main="")
  spectrum(res, col=par('fg'), main="")
  abline(v=1:10, lty=3)
  par(op)
dev.off()

png(file="g1177.png", width=600, height=800)
  innovations <- arima(res, c(2,0,0))$residuals
  op <- par(mfrow=c(4,1), mar=c(3,4,2,2))
  plot(innovations)
  acf(innovations)
  pacf(innovations)
  spectrum(innovations)
  par(op)
dev.off()

png(file="g1178.png", width=600, height=300)
  x <- co2
  n <- length(x)
  k <- 12
  m <- matrix( c(x, rep(NA,k)), nr=n+k-1, nc=k )
  y <- apply(m, 1, mean, na.rm=T)
  y <- y[1:n + round(k/2)]
  y <- ts(y, start=start(x), frequency=frequency(x))
  y <- y[round(k/4):(n-round(k/4))]
  yt <- time(x)[ round(k/4):(n-round(k/4)) ]
  plot(x, ylab="co2")
  lines(y~yt, col='red')
dev.off()

png(file="g1179.png", width=600, height=300)
  x <- co2
  plot(x, ylab="co2")
  k <- 12
  lines( filter(x, rep(1/k,k), side=1), col='red')
dev.off()

png(file="g1180.png", width=600, height=300)
  x <- co2
  plot(window(x, 1990, max(time(x))), ylab="co2")
  k <- 12
  lines( filter(x, rep(1/k,k)), 
         col='red', lwd=3)
  lines( filter(x, rep(1/k,k), sides=1), 
         col='blue', lwd=3)
  legend(par('usr')[1], par('usr')[4], xjust=0,
         c('smoother', 'filter'),
         lwd=3, lty=1,
         col=c('red','blue'))
dev.off()

png(file="g1181.png", width=600, height=600)
  library(fBasics) # RMetrics
  x <- yahooImport("s=IBM&a=11&b=1&c=1999&d=0&q=31&f=2000&z=IBM&x=.csv")
  x <- as.numeric(as.character(x@data$Close))
  x20 <- filter(x, rep(1/20,20), sides=1)
  x50 <- filter(x, rep(1/50,50), sides=1)
  matplot(cbind(x,x20,x50), type="l", lty=1, ylab="price", log="y")
  segments(1:length(x), x[-length(x)], 1:length(x), x[-1], lwd=ifelse(x20>x50,1,5)[-1], col=ifelse(x20>x50,"black","blue")[-1])
dev.off()

png(file="g1182.png", width=600, height=800)
  exponential.moving.average <- 
  function (x, a) {
    m <- x
    for (i in 2:n) {
      # Definition 
      # Exercise: use the "filter" function instead, 
      # with its "recursive" argument (that should be
      # much, much faster)
      m[i] <- a * x[i] + (1-a)*m[i-1]
    }
    m <- ts(m, start=start(x), frequency=frequency(x))
    m
  }
  plot.exponential.moving.average <- 
  function (x, a=.9, ...) {
    plot(exponential.moving.average(x,a), ...)
    par(usr=c(0,1,0,1))
    text(.02,.9, a, adj=c(0,1), cex=3)
  }
  op <- par(mfrow=c(4,1), mar=c(2,2,2,2)+.1)
  plot(x, main="Exponential Moving Averages")
  plot.exponential.moving.average(x,     main="", ylab="")
  plot.exponential.moving.average(x,.1,  main="", ylab="")
  plot.exponential.moving.average(x,.02, main="", ylab="")
  par(op)
dev.off()

png(file="g1183.png", width=600, height=800)
  r <- x - exponential.moving.average(x,.02)
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2)+.1)
  plot(r, main="Residuals of an exponential moving average")
  acf(r,  main="")
  pacf(r, main="")
  spectrum(r, main="")
  abline(v=1:10,lty=3)
  par(op)
dev.off()

png(file="g1184.png", width=600, height=600)
  x <- co2
  m <- HoltWinters(x, alpha=.1, beta=0, gamma=0)
  p <- predict(m, n.ahead=240, prediction.interval=T)
  plot(m, predicted.values=p)
dev.off()

png(file="g1185.png", width=600, height=600)
  y <- x20 > x50
  op <- par(mfrow=c(5,1), mar=c(0,0,0,0), oma=.1+c(0,0,0,0))
  plot(y, type="l", axes=FALSE); box()
  plot(filter(y, 50/100, "recursive", sides=1), axes=FALSE); box()
  plot(filter(y, 90/100, "recursive", sides=1), axes=FALSE); box()
  plot(filter(y, 95/100, "recursive", sides=1), axes=FALSE); box()
  plot(filter(y, 99/100, "recursive", sides=1), axes=FALSE); box()
  par(op)
dev.off()

png(file="g1186.png", width=600, height=600)
  op <- par(mfrow=c(3,1), mar=c(3,4,0,1), oma=c(0,0,2,0))
  r <- stl(co2, s.window="periodic")$time.series
  plot(co2)
  lines(r[,2], col='blue')
  lines(r[,2]+r[,1], col='red')
  plot(r[,3])
  acf(r[,3], main="residuals")
  par(op)
  mtext("STL(co2)", line=3, font=2, cex=1.2)
dev.off()

png(file="g1187.png", width=600, height=600)
  r <- decompose(co2)
  plot(r)
dev.off()

png(file="g1188.png", width=600, height=600)
  op <- par(mfrow=c(2,1), mar=c(0,2,0,2), oma=c(2,0,2,0))
  acf(r$random, na.action=na.pass, axes=F, ylab="")
  box(lwd=3)
  mtext("PACF", side=2, line=.5)
  pacf(r$random, na.action=na.pass, axes=F, ylab="")
  box(lwd=3)
  axis(1)
  mtext("PACF", side=2, line=.5)
  par(op)
  mtext("stl(co2): residuals", line=2.5, font=2, cex=1.2)
dev.off()

png(file="g1189.png", width=600, height=600)
  data(LakeHuron)
  x <- LakeHuron
  before <- window(x, end=1935)
  after <- window(x, start=1935)
  a <- .2
  b <- 0
  g <- 0
  model <- HoltWinters(
    before, 
    alpha=a, beta=b, gamma=g)
  forecast <- predict(
    model, 
    n.ahead=37, 
    prediction.interval=T)
  plot(model, predicted.values=forecast,
       main="Holt-Winters filtering: constant model")
  lines(after)
dev.off()

png(file="g1190.png", width=600, height=600)
  data(LakeHuron)
  x <- LakeHuron
  before <- window(x, end=1935)
  after <- window(x, start=1935)
  a <- .2
  b <- .2
  g <- 0
  model <- HoltWinters(
    before, 
    alpha=a, beta=b, gamma=g)
  forecast <- predict(
    model, 
    n.ahead=37, 
    prediction.interval=T)
  plot(model, predicted.values=forecast,
       main="Holt-Winters filtering: trend model")
  lines(after)
dev.off()

png(file="g1191.png", width=600, height=600)
  data(LakeHuron)
  x <- LakeHuron
  op <- par(mfrow=c(2,2), 
            mar=c(0,0,0,0),
            oma=c(1,1,3,1))
  before <- window(x, end=1935)
  after <- window(x, start=1935)
  a <- .2
  b <- .5
  g <- 0
  for (b in c(.02, .04, .1, .5)) {
    model <- HoltWinters(
      before, 
      alpha=a, beta=b, gamma=g)
    forecast <- predict(
      model, 
      n.ahead=37, 
      prediction.interval=T)
    plot(model, 
         predicted.values=forecast, 
         axes=F, xlab='', ylab='', main='')
    box()
   text( (4*par('usr')[1]+par('usr')[2])/5,
          (par('usr')[3]+5*par('usr')[4])/6,
          paste("beta =",b),
          cex=2, col='blue' )
    lines(after)
  }
  par(op)
  mtext("Holt-Winters filtering: different values for beta",
        line=-1.5, font=2, cex=1.2)
dev.off()

png(file="g1192.png", width=600, height=300)
  n <- 200
  plot(ts(cumsum(rnorm(n)) + rnorm(n)),
       main="Noisy random walk",
       ylab="")
dev.off()

png(file="g1193.png", width=600, height=600)
  op <- par(mfrow=c(2,1), 
            mar=c(3,4,0,2)+.1,
            oma=c(0,0,3,0))
  plot(ts(cumsum(rnorm(n, sd=1))+rnorm(n,sd=.1)),
       ylab="")
  plot(ts(cumsum(rnorm(n, sd=.1))+rnorm(n,sd=1)),
       ylab="")
  par(op)
  mtext("Noisy random walk", line=2, font=2, cex=1.2)
dev.off()

png(file="g1194.png", width=600, height=300)
  n <- 200
  plot(ts( cumsum( cumsum(rnorm(n))+rnorm(n) ) + 
           rnorm(n) ),
       main = "Local trend model",
       ylab="")
dev.off()

png(file="g1195.png", width=600, height=600)
  n <- 200
  op <- par(mfrow=c(3,1), 
            mar=c(3,4,2,2)+.1,
            oma=c(0,0,2,0))
  plot(ts( cumsum( cumsum(rnorm(n,sd=1))+rnorm(n,sd=1) )
       + rnorm(n,sd=.1) ),
       ylab="")
  plot(ts( cumsum( cumsum(rnorm(n,sd=1))+rnorm(n,sd=.1) )
       + rnorm(n,sd=1) ),
       ylab="")
  plot(ts( cumsum( cumsum(rnorm(n,sd=.1))+rnorm(n,sd=1) )
       + rnorm(n,sd=1) ),
       ylab="")  
  par(op)
  mtext("Local level models", line=2, font=2, cex=1.2) 
dev.off()

png(file="g1196.png", width=600, height=600)
  structural.model <- function (
      n=200, 
      sd1=1, sd2=1, sd3=1, sd4=1, sd5=200
  ) {
    sd1 <- 1
    sd2 <- 2
    sd3 <- 3
    sd4 <- 4
    mu <- rep(rnorm(1),n)
    nu <- rep(rnorm(1),n)
    g <- rep(rnorm(1,sd=sd5),n)
    x <- mu + g + rnorm(1,sd=sd1)
    for (i in 2:n) {
      if (i>3) {    
        g[i] <- -(g[i-1]+g[i-2]+g[i-3]) + rnorm(1,sd=sd4)
      } else {
        g[i] <- rnorm(1,sd=sd5)
      }
      nu[i] <- nu[i-1] + rnorm(1,sd=sd3)
      mu[i] <- mu[i-1] + nu[i-1] + rnorm(1,sd=sd2)
      x[i] <- mu[i] + g[i] + rnorm(1,sd=sd1)
    }
    ts(x)
  }
  n <- 200
  op <- par(mfrow=c(3,1), 
            mar=c(2,2,2,2)+.1,
            oma=c(0,0,2,0))
  plot(structural.model(n))
  plot(structural.model(n))
  plot(structural.model(n))
  par(op)
  mtext("Structural models", line=2.5, font=2, cex=1.2)
dev.off()

png(file="g1197.png", width=600, height=600)
  n <- 200
  op <- par(mfrow=c(3,1), 
            mar=c(2,2,2,2)+.1,
            oma=c(0,0,2,0))
  plot(structural.model(n,sd4=0))
  plot(structural.model(n,sd4=0))
  plot(structural.model(n,sd4=0))
  par(op)
  mtext("Structural models", line=2.5, font=2, cex=1.2)
dev.off()

png(file="g1198.png", width=600, height=300)
  data(AirPassengers)
  plot(AirPassengers)
dev.off()

png(file="g1199.png", width=600, height=300)
  plot(log(AirPassengers))
dev.off()

png(file="g1200.png", width=600, height=600)
  x <- log(AirPassengers)
  r <- StructTS(x)
  plot(x, main="AirPassengers", ylab="")
  f <- apply(fitted(r), 1, sum)
  f <- ts(f, frequency=frequency(x), start=start(x))
  lines(f, col='red', lty=3, lwd=3)
dev.off()

png(file="g1201.png", width=600, height=600)
  matplot( 
    (StructTS(x-min(x)))$fitted, 
    type = 'l',
    ylab = "",
    main = "Structural model decomposition of a time series"
  )
dev.off()

png(file="g1202.png", width=600, height=600)
  l <- 1956
  x <- log(AirPassengers)
  x1 <- window(x, end=l)
  x2 <- window(x, start=l)
  r <- StructTS(x1)
  plot(x)
  f <- apply(fitted(r), 1, sum)
  f <- ts(f, frequency=frequency(x), start=start(x))
  lines(f, col='red')
  p <- predict(r, n.ahead=100)
  lines(p$pred, col='red')
  lines(p$pred + qnorm(.025) * p$se, 
        col='red', lty=2)
  lines(p$pred + qnorm(.975) * p$se, 
        col='red', lty=2)
  title(main="Forecasting with a structural model (StructTS)")
dev.off()

png(file="g1203.png", width=600, height=600)
  # A function to look at a time series
  eda.ts <- function (x, bands=FALSE) {
    op <- par(no.readonly = TRUE)
    par(mar=c(0,0,0,0), oma=c(1,4,2,1))
    # Compute the Ljung-Box p-values
    # (we only display them if needed, i.e., 
    # if we have any reason of
    # thinking it is white noise).
    p.min <- .05
    k <- 15
    p <- rep(NA, k)
    for (i in 1:k) {
      p[i] <- Box.test(
        x, i, type = "Ljung-Box"
      )$p.value
    }
    if( max(p)>p.min ) {
      par(mfrow=c(5,1))
    } else {
      par(mfrow=c(4,1))
    }
    if(!is.ts(x))
      x <- ts(x)
    plot(x, axes=FALSE);
    axis(2); axis(3); box(lwd=2)
    if(bands) {
      a <- time(x)
      i1 <- floor(min(a))
      i2 <- ceiling(max(a))
      y1 <- par('usr')[3]
      y2 <- par('usr')[4]
      if( par("ylog") ){
        y1 <- 10^y1
        y2 <- 10^y2
      }
      for (i in seq(from=i1, to=i2-1, by=2)) {
        polygon( c(i,i+1,i+1,i), c(y1,y1,y2,y2), 
                 col='grey', border=NA )
      }
      lines(x)
    }
    acf(x, axes=FALSE)
    axis(2, las=2)
    box(lwd=2)
    mtext("ACF", side=2, line=2.5)
    pacf(x, axes=FALSE)
    axis(2, las=2)
    box(lwd=2)
    mtext("ACF", side=2, line=2.5)
    spectrum(x, col=par('fg'), log="dB", 
             main="", axes=FALSE )
    axis(2, las=2)
    box(lwd=2)
    mtext("Spectrum", side=2, line=2.5)
    abline(v=1, lty=2, lwd=2)
    abline(v=2:10, lty=3)
    abline(v=1/2:5, lty=3)
    if( max(p)>p.min ) {
      main <- 
      plot(p, type='h', ylim=c(0,1), 
           lwd=3, main="", axes=F)
      axis(2, las=2)
      box(lwd=2)
      mtext("Ljung-Box p-value", side=2, line=2.5)
      abline(h=c(0,.05),lty=3)
    }
    par(op)
  }
 
  data(co2)
  eda.ts(co2, bands=T)
dev.off()

png(file="g1204.png", width=600, height=600)
  data(AirPassengers)
  x <- AirPassengers
  plot(x)
dev.off()

png(file="g1205.png", width=600, height=600)
  plot(x)
  abline(lm(x~time(x)), col='red')
dev.off()

png(file="g1206.png", width=600, height=600)
  plot(lm(x~time(x))$res)
dev.off()

png(file="g1207.png", width=600, height=600)
  plot(diff(x))
dev.off()

png(file="g1208.png", width=600, height=600)
  x <- log(x)
  plot(x)
  abline(lm(x~time(x)),col='red')
dev.off()

png(file="g1209.png", width=600, height=600)
  plot(lm(x~time(x))$res)
dev.off()

png(file="g1210.png", width=600, height=600)
  plot(diff(x))
dev.off()

png(file="g1211.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  plot(diff(x,1,2))
  plot(diff(x,1,3))
  plot(diff(x,1,4))
  par(op)
dev.off()

png(file="g1212.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  plot(x)
  abline(h=0, v=1950:1962, lty=3)
  y <- diff(x)
  plot(y)
  abline(h=0, v=1950:1962, lty=3)
  plot(diff(y, 12,1))
  abline(h=0, v=1950:1962, lty=3)
  par(op)
dev.off()

png(file="g1213.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  plot(x)
  abline(lm(x~time(x)), col='red', lty=2)
  abline(h=0, v=1950:1962, lty=3)
  y <- x - predict(lm(x~time(x)))
  plot(y)
  abline(h=0, v=1950:1962, lty=3)
  plot(diff(y, 12,1))
  abline(h=0, v=1950:1962, lty=3)
  par(op)
dev.off()

png(file="g1214.png", width=600, height=600)
  z <- diff(y,12,1)
  op <- par(mfrow=c(3,1))
  plot(z)
  abline(h=0,lty=3)
  plot(diff(z))
  abline(h=0,lty=3)
  plot(diff(z,1,2))
  abline(h=0,lty=3)
  par(op)  
dev.off()

png(file="g1215.png", width=600, height=600)
  k <- 3
  op <- par(mfrow=c(k,2))
  zz <- z
  for(i in 1:k) {
    acf(zz, main=i-1)
    pacf(zz, main=i-1)
    zz <- diff(zz)
  }
  par(op)
dev.off()

png(file="g1216.png", width=600, height=900)
  data(sunspots)
  op <- par(mfrow=c(5,1))
  for (i in 10+1:5) {
    plot(diff(sunspots,i))
  }
  par(op)
dev.off()

png(file="g1217.png", width=600, height=600)
  op <- par(mfrow=c(2,1), mar=c(2,4,3,2)+.1)
  x <- ts(rnorm(200))
  plot(x, main="gaussian iidrv", 
       xlab="", ylab="")
  acf(x, main="")
  par(op)
dev.off()

png(file="g1218.png", width=600, height=600)
  op <- par(mfrow=c(2,1), mar=c(2,4,3,2)+.1)
  data(BJsales)
  plot(BJsales, xlab="", ylab="", main="BJsales")
  acf(BJsales, main="")
  par(op)
dev.off()

png(file="g1219.png", width=600, height=600)
  f <- 24
  x <- seq(0,10, by=1/f)
  y <- sin(2*pi*x)
  y <- ts(y, start=0, frequency=f)
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2)+.1)
  plot(y, xlab="", ylab="")
  acf(y,  main="")
  pacf(y, main="")
  spectrum(y, main="", xlab="")
  par(op)
dev.off()

png(file="g1220.png", width=600, height=600)
  f <- 24
  x <- seq(0,10, by=1/f)
  y <- x + sin(2*pi*x) + rnorm(10*f)
  y <- ts(y, start=0, frequency=f)
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2)+.1)
  plot(y, xlab="", ylab="")
  acf(y,  main="")
  pacf(y, main="")
  spectrum(y, main="", xlab="")
  par(op)
dev.off()

png(file="g1221.png", width=600, height=600)
  n <- 200  
  x <- rnorm(n)
  y <- ( x[2:n] + x[2:n-1] ) / 2
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="white noise")
  plot(ts(y), xlab="", ylab="MA(1)")
  acf(y, main="")
  par(op)
dev.off()

png(file="g1222.png", width=600, height=600)
  n <- 200  
  x <- rnorm(n)
  y <- ( x[1:(n-3)] + x[2:(n-2)] + x[3:(n-1)] + x[4:n] )/4
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="white noise")
  plot(ts(y), xlab="", ylab="MA(3)")
  acf(y, main="")
  par(op)
dev.off()

png(file="g1223.png", width=600, height=600)
  n <- 200  
  x <- rnorm(n)
  y <- x[2:n] - x[1:(n-1)]
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="white noise")
  plot(ts(y), xlab="", ylab="momentum(1)")
  acf(y, main="")
  par(op)
dev.off()

png(file="g1224.png", width=600, height=600)
  n <- 200  
  x <- rnorm(n)
  y <- x[3:n] - 2 * x[2:(n-1)] + x[1:(n-2)]
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="white noise")
  plot(ts(y), xlab="", ylab="Momentum(2)")
  acf(y, main="")
  par(op)
dev.off()

png(file="g1225.png", width=600, height=600)
  n <- 200  
  x <- rnorm(n)
  y <- filter(x, c(1,-2,1))
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="White noise")
  plot(ts(y), xlab="", ylab="Momentum(2)")
  acf(y, na.action=na.pass, main="")
  par(op)
dev.off()

png(file="g1226.png", width=600, height=600)
  n <- 200
  x <- rnorm(n)
  y <- cumsum(x)
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="")
  plot(ts(y), xlab="", ylab="AR(1)")
  acf(y, main="")
  par(op)  
dev.off()

png(file="g1227.png", width=600, height=600)
  n <- 200
  a <- .7
  x <- rep(0,n)
  for (i in 2:n) {
    x[i] <- a*x[i-1] + rnorm(1)
  }
  y <- x[-1]
  x <- x[-n]
  r <- lm( y ~ x -1)
  plot(y~x)
  abline(r, col='red')
  abline(0, .7, lty=2)
dev.off()

png(file="g1228.png", width=600, height=600)
  n <- 200
  x <- rep(0,n)
  for (i in 4:n) {
    x[i] <- .3*x[i-1] -.7*x[i-2] + .5*x[i-3] + rnorm(1)
  }
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="AR(3)")
  acf(x,  main="", xlab="")
  pacf(x, main="", xlab="")
  par(op)
dev.off()

png(file="g1229.png", width=600, height=600)
  n <- 200
  x <- arima.sim(list(ar=c(.3,-.7,.5)), n)
  op <- par(mfrow=c(3,1), mar=c(2,4,2,2)+.1)
  plot(ts(x), xlab="", ylab="AR(3)")
  acf(x,  xlab="", main="")
  pacf(x, xlab="", main="")
  par(op)
dev.off()

png(file="g1230.png", width=600, height=600)
  n <- 200
  x <- seq(0,2,length=n)
  trend <- ts(sin(x))
  plot(trend, 
        ylim=c(-.5,1.5), 
       lty=2, lwd=3, col='red', 
       ylab='')
  r <- arima.sim(
    list(ar = c(0.5,-.3), ma = c(.7,.1)), 
    n, 
    sd=.1 
  )
  lines(trend+r)
dev.off()

png(file="g1231.png", width=600, height=600)
  n <- 200
  k <- 10
  x <- 1:n
  r <- matrix(nr=n,nc=k)
  for (i in 1:k) {
    r[,i] <- cumsum(rnorm(n))
  }
  matplot(x, r,
          type = 'l', 
          lty = 1, 
          col = par('fg'),
          main = "A random walk is not stationnary")
  abline(h=0,lty=3)
dev.off()

png(file="g1232.png", width=600, height=600)
  n <- 200
  ma <- 2
  mai <- 1/ma
  op <- par(mfrow=c(4,1), mar=c(2,4,1,2)+.1)
  x <- arima.sim(list(ma=ma),n)
  plot(x, xlab="", ylab="")
  acf(x, xlab="", main="")
  lines(0:n, 
        ARMAacf(ma=ma, lag.max=n), 
        lty=2, lwd=3, col='red')
  x <- arima.sim(list(ma=mai),n)
  plot(x, xlab="", ylab="")
  acf(x, main="", xlab="")
  lines(0:n, 
        ARMAacf(ma=mai, lag.max=n), 
        lty=2, lwd=3, col='red')
  par(op)
dev.off()

png(file="g1233.png", width=600, height=600)
  sym.poly <- function (z,k) {
    # Sum of the products of k 
    # distinct elements of the vector z
    if (k==0) {
      r <- 1
    } else if (k==1) {
      r <- sum(z)
    } else {
      r <- 0
      for (i in 1:length(z)) {
        r <- r + z[i]*sym.poly(z[-i],k-1)
      }
      r <- r/k   # Each term appeared k times
    }
    r
  }
  sym.poly( c(1,2,3), 1 )   #  6
  sym.poly( c(1,2,3), 2 )   # 11
  sym.poly( c(1,2,3), 3 )   #  6

  roots.to.poly <- function (z) {
    n <- length(z)
    p <- rep(1,n)
    for (k in 1:n) {
      p[n-k+1] <- (-1)^k * sym.poly(z,k)
    }
    p <- c(p,1)
    p
  }
  roots.to.poly(c(1,2))   # 2 -3 1
  round(
    Re(polyroot( roots.to.poly(c(1,2,3)) )), 
    digits=1
  )

  # After this interlude, we can finally 
  # construct an MA process and one of 
  # its inverses
  n <- 200
  k <- 3
  ma <- runif(k,-1,1)
  # The roots
  z <- polyroot(c(1,-ma))
  # The inverse of the roots
  zi <- 1/z
  # The polynomial
  p <- roots.to.poly(zi)
  # The result should be real, but because 
  # of rounding errors, it is not.
  p <- Re(p)
  # We want the constant term to be 1.
  p <- p/p[1]
  mai <- -p[-1]

  op <- par(mfrow=c(4,1), mar=c(2,4,1,2)+.1)
  x <- arima.sim(list(ma=ma),n)
  plot(x, xlab="")
  acf(x, main="", xlab="")
  lines(0:n, ARMAacf(ma=ma, lag.max=n), 
        lty=2, lwd=3, col='red')
  x <- arima.sim(list(ma=mai),n)
  plot(x, xlab="")
  acf(x, main="", xlab="")
  lines(0:n, ARMAacf(ma=mai, lag.max=n), 
        lty=2, lwd=3, col='red')
  par(op)
dev.off()

png(file="g1234.png", width=600, height=600)
  data(Nile)
  op <- par(mfrow=c(2,1), mar=c(2,4,3,2)+.1)
  plot(Nile, main="There is no trend", xlab="")
  acf(Nile, main="", xlab="")
  par(op)
dev.off()

png(file="g1235.png", width=600, height=600)
  data(BJsales)
  op <- par(mfrow=c(3,1), mar=c(2,4,3,2)+.1)
  plot(BJsales, xlab="",
      main="The trend disappears if we differentiate")
  acf(BJsales, xlab="", main="")
  acf(diff(BJsales), xlab="", main="", 
      ylab="ACF(diff(BJsales)")
  par(op)
dev.off()

png(file="g1236.png", width=600, height=600)
  n <- 2000
  x <- arima.sim(
    model = list(
      ar = c(.3,.6), 
      ma = c(.8,-.5,.2), 
      order = c(2,1,3)), 
    n
  )
  x <- ts(x)
  op <- par(mfrow=c(3,1), mar=c(2,4,3,2)+.1)
  plot(x, main="It suffices to differentiate once",
       xlab="", ylab="")
  acf(x, xlab="", main="")
  acf(diff(x), xlab="", main="",
      ylab="ACF(diff(x))")
  par(op)
dev.off()

png(file="g1237.png", width=600, height=800)
  n <- 10000
  x <- arima.sim(
    model = list(
      ar = c(.3,.6), 
      ma = c(.8,-.5,.2), 
      order = c(2,2,3)
    ), 
    n
  )
  x <- ts(x)
  op <- par(mfrow=c(4,1), mar=c(2,4,3,2)+.1)
  plot(x, main="One has to differentiate twice",
       xlab="", ylab="")
  acf(x, main="", xlab="")
  acf(diff(x), main="", xlab="",
      ylab="ACF(diff(x))")
  acf(diff(x,differences=2), main="", xlab="",
      ylab="ACF(diff(diff(x)))")
  par(op)
dev.off()

png(file="g1238.png", width=600, height=600)
  acf.exp <- function (x, lag.max=NULL, lag.max.reg=lag.max, ...) {
    a <- acf(x, lag.max=lag.max.reg, plot=F)
    b <- acf(x, lag.max=lag.max, ...)
    r <- lm( log(a$acf) ~ a$lag -1)
    lines( exp( b$lag * r$coef[1] ) ~ b$lag, lty=2 )
  }
  data(BJsales)
  acf.exp(BJsales,
          main="Exponential decay of the ACF")
dev.off()

png(file="g1239.png", width=600, height=600)
  acf.exp(BJsales, 
          lag.max=40,
          main="Exponential decay of the ACF")
dev.off()

png(file="g1240.png", width=600, height=600)
  data(Nile)
  acf.exp(Nile, lag.max.reg=10, main="Nile")
dev.off()

png(file="g1241.png", width=600, height=800)
  x <- diff(co2, lag=12)
  op <- par(mfrow=c(4,1), mar=c(2,4,3,2)+.1)
  plot(x, ylab="", xlab="")
  acf(x, xlab="", main="")
  pacf(x, xlab="", main="")
  spectrum(x, xlab="", main="",
           col=par('fg'))
  par(op)
dev.off()

png(file="g1242.png", width=600, height=800)
  y <- diff(x)
  op <- par(mfrow=c(4,1), mar=c(2,4,3,2)+.1)
  plot(y, xlab="", ylab="")
  acf(y,  xlab="", main="")
  pacf(y, xlab="", main="")
  spectrum(y, col=par('fg'),
           xlab="", main="")
  par(op)
dev.off()

png(file="g1243.png", width=600, height=600)
  n <- 200
  x <- arima.sim(
    list(
      order=c(2,1,2), 
      ar=c(.5,-.8), 
      ma=c(.9,.6)
    ), 
    n
  )
  op <- par(mfrow=c(3,1), mar=c(2,4,4,2)+.1)
  acf(x, main="You will have to defferentiate once")
  acf(diff(x), main="First derivative")
  acf(diff(x, differences=2), main="Second derivative")
  par(op)
dev.off()

png(file="g1244.png", width=600, height=600)
  n <- 200
  x <- arima.sim(
    list(
      order=c(2,2,2), 
      ar=c(.5,-.8), 
      ma=c(.9,.6)
    ), 
    n
  )
  op <- par(mfrow=c(3,1), mar=c(2,4,4,2)+.1)
  acf(x, main="You will have to differentiate twice")
  acf(diff(x), main="First derivative")
  acf(diff(x, differences=2), main="Second derivative")
  par(op)
dev.off()

png(file="g1245.png", width=600, height=600)
  data(sunspot)
  op <- par(mfrow=c(4,1), mar=c(2,4,3,2)+.1)
  plot(sunspot.month, xlab="", ylab="sunspot")
  acf(sunspot.month,  xlab="", main="")
  plot(diff(sunspot.month), 
       xlab="", ylab="diff(sunspot)")
  acf(diff(sunspot.month), xlab="", main="")
  par(op) 
dev.off()

png(file="g1246.png", width=600, height=600)
  data(JohnsonJohnson)
  x <- log(JohnsonJohnson)
  op <- par(mfrow=c(4,1), mar=c(2,4,3,2)+.1)
  plot(x, xlab="", ylab="JJ")
  acf(x, main="")
  plot(diff(x), ylab="diff(JJ)")
  acf(diff(x), main="")
  par(op) 
dev.off()

png(file="g1247.png", width=600, height=600)
  data(BJsales)
  x <- BJsales
  op <- par(mfrow=c(6,1), mar=c(2,4,0,2)+.1)
  plot(x)
  acf(x)
  plot(diff(x))
  acf(diff(x))
  plot(diff(x, difference=2))
  acf(diff(x, difference=2))
  par(op) 
dev.off()

png(file="g1248.png", width=600, height=600)
  data(austres)
  x <- austres
  op <- par(mfrow=c(6,1), mar=c(2,4,0,2)+.1)
  plot(x)
  acf(x)
  plot(diff(x))
  acf(diff(x))
  plot(diff(x, difference=2))
  acf(diff(x, difference=2))
  par(op) 
dev.off()

png(file="g1249.png", width=600, height=600)
  # In the preceding example, there was a linear trend:
  # let ut remove it.
  data(austres)
  x <- lm(austres ~ time(austres))$res
  op <- par(mfrow=c(6,1), mar=c(2,4,0,2)+.1)
  plot(x)
  acf(x)
  plot(diff(x))
  acf(diff(x))
  plot(diff(x, difference=2))
  acf(diff(x, difference=2))
  par(op) 
dev.off()

png(file="g1250.png", width=600, height=600)
  my.sarima.sim <- function (
      n = 20, 
      period = 12, 
      model, 
      seasonal
  ) {    
    x <- arima.sim( model, n*period )
    x <- x[1:(n*period)]
    for (i in 1:period) {
      xx <- arima.sim( seasonal, n )
      xx <- xx[1:n]
      x[i + period * 0:(n-1)] <- 
        x[i + period * 0:(n-1)] + xx
    }
    x <- ts(x, frequency=period)
    x
  }
  op <- par(mfrow=c(3,1))
  x <- my.sarima.sim(
    20, 
    12, 
    list(ar=.6, ma=.3, order=c(1,0,1)),
    list(ar=c(.5), ma=c(1,2), order=c(1,0,2))
  )
  eda.ts(x, bands=T)
dev.off()

png(file="g1251.png", width=600, height=600)
  x <- my.sarima.sim(
    20, 
    12, 
    list(ar=c(.5,-.3), ma=c(-.8,.5,-.3), order=c(2,1,3)),
    list(ar=c(.5), ma=c(1,2), order=c(1,0,2))
  )
  eda.ts(x, bands=T)
dev.off()

png(file="g1252.png", width=600, height=600)
  x <- my.sarima.sim(
    20, 
    12, 
    list(ar=c(.5,-.3), ma=c(-.8,.5,-.3), order=c(2,1,3)),
    list(ar=c(.5), ma=c(1,2), order=c(1,1,2))
  )
  eda.ts(x, bands=T)
dev.off()

png(file="g1253.png", width=600, height=600)
  x <- co2
  eda.ts(x)
dev.off()

png(file="g1254.png", width=600, height=600)
  eda.ts(diff(x))
dev.off()

png(file="g1255.png", width=600, height=600)
  eda.ts(diff(diff(x),lag=12))
dev.off()

png(file="g1256.png", width=600, height=600)
  eda.ts(diff(x,lag=12))
dev.off()

png(file="g1257.png", width=600, height=600)
  r3 <- arima(
    co2, 
    order = c(1, 1, 1), 
    seasonal = list(
      order = c(0, 1, 1), 
      period = 12
    )
  )
  eda.ts(r3$res)
dev.off()

png(file="g1258.png", width=600, height=600)
  x1 <- window(co2, end = 1990)
  r <- arima(
    x1, 
    order = c(1, 1, 1), 
    seasonal = list(
      order = c(0, 1, 1), 
      period = 12
    )
  )
  plot(co2)
  p <- predict(r, n.ahead=100)
  lines(p$pred, col='red')
  lines(p$pred+qnorm(.025)*p$se, col='red', lty=3)
  lines(p$pred+qnorm(.975)*p$se, col='red', lty=3)
dev.off()

png(file="g1259.png", width=600, height=600)
  # On the contrary, I do not know what to do with 
  # this plots (it looks like integrated noise).
  eda.ts(co2-p$pred)
dev.off()

png(file="g1260.png", width=600, height=600)
  r <- arima(
    co2, 
    order = c(1, 1, 1), 
    seasonal = list(
      order = c(0, 1, 1), 
      period = 12
    )
  )
  p <- predict(r, n.ahead=150)
  plot(co2, 
       xlim=c(1959,2010), 
       ylim=range(c(co2,p$pred)))
  lines(p$pred, col='red')
  lines(p$pred+qnorm(.025)*p$se, col='red', lty=3)
  lines(p$pred+qnorm(.975)*p$se, col='red', lty=3)
dev.off()

png(file="g1261.png", width=600, height=600)
  op <- par(mfrow=c(4,2), mar=c(2,4,4,2))
  n <- 200
  for (i in 1:4) {
    x <- NULL
    while(is.null(x)) {
      model <- list(ar=rnorm(1))
      try( x <- arima.sim(model, n) )
    }
    acf(x, 
        main = paste(
          "ARMA(1,0)",
          "AR:",
          round(model$ar, digits = 1)
        ))
    points(0:50, 
           ARMAacf(ar=model$ar, lag.max=50), 
           col='red')
    pacf(x, main="")
    points(1:50, 
           ARMAacf(ar=model$ar, lag.max=50, pacf=T), 
           col='red')
  } 
  par(op)  
dev.off()

png(file="g1262.png", width=600, height=600)
  op <- par(mfrow=c(4,2), mar=c(2,4,4,2))
  n <- 200
  for (i in 1:4) {
    x <- NULL
    while(is.null(x)) {
      model <- list(ar=rnorm(2))
      try( x <- arima.sim(model, n) )
    }
    acf(x, 
        main=paste("ARMA(2,0)","AR:",
                   round(model$ar[1],digits=1),
                   round(model$ar[2],digits=1)
                   ))
    points(0:50, 
           ARMAacf(ar=model$ar, lag.max=50), 
           col='red')
    pacf(x, main="")
    points(1:50, 
           ARMAacf(ar=model$ar, lag.max=50, pacf=T), 
           col='red')
  } 
  par(op)  
dev.off()

png(file="g1263.png", width=600, height=600)
  op <- par(mfrow=c(4,2), mar=c(2,4,4,2))
  n <- 200
  for (i in 1:4) {
    x <- NULL
    while(is.null(x)) {
      model <- list(ma=rnorm(1))
      try( x <- arima.sim(model, n) )
    }
    acf(x, 
        main = paste(
          "ARMA(0,1)",
          "MA:",
          round(model$ma, digits=1)
        ))
    points(0:50, 
           ARMAacf(ma=model$ma, lag.max=50), 
           col='red')
    pacf(x, main="")
    points(1:50, 
           ARMAacf(ma=model$ma, lag.max=50, pacf=T), 
           col='red')
  } 
  par(op)  
dev.off()

png(file="g1264.png", width=600, height=600)
  op <- par(mfrow=c(4,2), mar=c(2,4,4,2))
  n <- 200
  for (i in 1:4) {
    x <- NULL
    while(is.null(x)) {
      model <- list(ma=rnorm(2))
      try( x <- arima.sim(model, n) )
    }
    acf(x, main=paste("ARMA(0,2)","MA:",
                   round(model$ma[1],digits=1),
                   round(model$ma[2],digits=1)
                   ))
    points(0:50, 
           ARMAacf(ma=model$ma, lag.max=50), 
           col='red')
    pacf(x, main="")
    points(1:50, 
           ARMAacf(ma=model$ma, lag.max=50, pacf=T), 
           col='red')
  } 
  par(op)  
dev.off()

png(file="g1265.png", width=600, height=600)
  op <- par(mfrow=c(4,2), mar=c(2,4,4,2))
  n <- 200
  for (i in 1:4) {
    x <- NULL
    while(is.null(x)) {
      model <- list(ma=rnorm(1),ar=rnorm(1))
      try( x <- arima.sim(model, n) )
    }
    acf(x, main=paste("ARMA(1,1)",
                   "AR:", round(model$ar,digits=1),
                   "AR:", round(model$ma,digits=1)
                   ))
    points(0:50, 
           ARMAacf(ar=model$ar, ma=model$ma, lag.max=50), 
           col='red')
    pacf(x, main="")
    points(1:50, 
           ARMAacf(ar=model$ar, ma=model$ma, lag.max=50, pacf=T), 
           col='red')
  } 
  par(op)  
dev.off()

png(file="g1266.png", width=600, height=600)
  r <- arima(co2, order=c(0,1,1), list(order=c(0,1,1), period=12 ) ) 
  eda.ts(r$res)
dev.off()

png(file="g1267.png", width=600, height=600)
  x1 <- window(co2,end=1990)
  r <- arima(x1, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  plot(co2)
  p <- predict(r, n.ahead=100)
  lines(p$pred, col='red')
  lines(p$pred+qnorm(.025)*p$se, col='red', lty=3)
  lines(p$pred+qnorm(.975)*p$se, col='red', lty=3)
dev.off()

png(file="g1268.png", width=600, height=600)
  eda.ts(co2-p$pred)
dev.off()

png(file="g1269.png", width=600, height=600)
  r <- arima(co2, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  p <- predict(r, n.ahead=150)
  plot(co2, xlim=c(1959,2010), ylim=range(c(co2,p$pred)))
  lines(p$pred, col='red')
  lines(p$pred+qnorm(.025)*p$se, col='red', lty=3)
  lines(p$pred+qnorm(.975)*p$se, col='red', lty=3)
dev.off()

png(file="g1270.png", width=600, height=600)
  r1 <- arima(co2, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  r2 <- arima(co2, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12))
  p1 <- predict(r1, n.ahead=150)
  p2 <- predict(r2, n.ahead=150)
  plot(co2, xlim=c(1959,2010), ylim=range(c(co2,p$pred)))
  lines(p1$pred, col='red')
  lines(p2$pred, col='green', lty=3, lwd=4)
dev.off()

png(file="g1271.png", width=600, height=300)
  spectrum(co2)
  abline(v=1:10, lty=3)
dev.off()

png(file="g1272.png", width=600, height=600)
  signal.and.spectrum <- function (x, main="") {
    op <- par(mfrow=c(2,1), 
              mar=c(2,4,2,2)+.1,
              oma=c(0,0,2,0))
    plot(x, type="l", main="", ylab="Signal")
    spectrum(x, main="", xlab="")
    abline(v=.1*1:10, lty=3)
    par(op)
   mtext(main, line=1.5, font=2, cex=1.2)
  }
  N <- 100
  x <- 10 * (1:N / N)
  signal.and.spectrum(sin(2*pi*x), 
                      "Sine wave, period=10")
dev.off()

png(file="g1273.png", width=600, height=600)
  signal.and.spectrum(x - floor(x), 
                      "sawtooth, period=10")
dev.off()

png(file="g1274.png", width=600, height=600)
  signal.and.spectrum(abs(x - floor(x)-.5), 
                      "triangle, period=10")
dev.off()

png(file="g1275.png", width=600, height=600)
  signal.and.spectrum(x-floor(x)>.5, 
                      "square, period=10")
dev.off()

png(file="g1276.png", width=600, height=300)
  data(sunspots)
  plot(sunspots)
dev.off()

png(file="g1277.png", width=600, height=300)
  spectrum(sunspots)
dev.off()

png(file="g1278.png", width=600, height=300)
  spectrum(sunspots, spans=10)
dev.off()

png(file="g1279.png", width=600, height=300)
  spectrum(sunspots, spans=10, xlim=c(0,1))
  abline(v=1/1:12, lty=3)
dev.off()

png(file="g1280.png", width=600, height=300)
  spectrum(sunspots, spans=10, xlim=c(0,1), 
           main="10: Not quite")
  abline(v=1:3/10, lty=3)
dev.off()

png(file="g1281.png", width=600, height=300)
  spectrum(sunspots, spans=10, xlim=c(0,1), 
           main="11: Better")
  abline(v=1:3/11, lty=3)
dev.off()

png(file="g1282.png", width=600, height=300)
  spectrum(sunspots, spans=10, xlim=c(0,1), 
           main="12: Not quite")
  abline(v=1:3/12, lty=3)
dev.off()

png(file="g1283.png", width=600, height=300)
  library(car)
  box.cox.powers(3+sunspots)
  y <- box.cox(sunspots,.3)
  spectrum(y, spans=10, xlim=c(0,1))
  abline(v=1:3/11, lty=3) # A single harmonic
dev.off()

png(file="g1284.png", width=600, height=300)
  x <- as.vector(time(y))
  y <- as.vector(y)
  r1 <- lm( y ~ sin(2*pi*x/11) + cos(2*pi*x/11) )
  r2 <- lm( y ~ sin(2*pi*x/11) + cos(2*pi*x/11) 
                + sin(2* 2*pi*x/11) + cos(2* 2*pi*x/11) 
                + sin(3* 2*pi*x/11) + cos(3* 2*pi*x/11) 
                + sin(4* 2*pi*x/11) + cos(4* 2*pi*x/11) 
          )
  plot(y~x, type='l')
  lines(predict(r1)~x, col='red')
  lines(predict(r2)~x, col='green', lty=3, lwd=3)
dev.off()

png(file="g1285.png", width=600, height=600)
  z <- y - predict(r1)
  op <- par(mfrow=c(2,2), mar=c(2,4,3,2)+.1)
  hist(z, probability=T, col="light blue")
  lines(density(z), lwd=3, col="red")
  qqnorm(z)
  acf(z,  main="")
  pacf(z, main="")
  par(op)
dev.off()

png(file="g1286.png", width=600, height=600)
  n <- 8
  x <- 1:n/n*2*pi
  k <- 0
  y <- cbind(rep(1,n),
             1 + cos(x),
             1 + cos(x) + sin(x),
             1 + cos(x) + sin(x) + cos(2*x),
             1 + cos(x) + sin(x) + cos(2*x) + sin(2*x)
             )
  z <- mvfft(y)  
  op <- par(mfrow=c(2,1))
  barplot(Re(t(z)), 
          beside = TRUE, 
          col = rainbow(dim(y)[2]), 
          main = "Real part of the FFT")
  barplot(Im(t(z)), 
          beside = TRUE, 
          col = rainbow(dim(y)[2]), 
          main = "Imaginary part of the FFT")
  par(op)
dev.off()

png(file="g1287.png", width=600, height=600)
  x <- rep(0,256)
  x[129:256] <- 1
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal",
       xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', lwd=3, col="blue",
       xlab="", ylab="Mod(fft(x))", 
       main="DFT (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1288.png", width=600, height=600)
  x <- rep(0,256)
  x[33:64] <- 1
  x[64+33:64] <- 1
  x[128+33:64] <- 1
  x[128+64+33:64] <- 1
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1289.png", width=600, height=600)
  x <- rep( 1:32/32, 8 )
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1290.png", width=600, height=600)
  x <- 1:256/256
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1291.png", width=600, height=600)
  x <- abs(1:256-128)
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1292.png", width=600, height=600)
  x <- 1:256/256
  x <- sin(2*pi*x) + cos(3*pi*x) + sin(4*pi*x+pi/3)
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1293.png", width=600, height=600)
  x <- 1:256/256
  x <- sin(16*pi*x)
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1294.png", width=600, height=600)
  x <- 1:256/256
  x <- sin(16*pi*x) + .3*cos(56*pi*x)
  op <- par(mfrow=c(2,1), mar=c(3,4,4,2)+.1)
  plot(x, type='l', lwd=3,
       main="Signal", xlab="", ylab="")
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', col="blue", lwd=3,
       ylab="Mof(fft(x))", xlab="",
       main="DTF (Discrete Fourrier Transform)")
  par(op)
dev.off()

png(file="g1295.png", width=600, height=600)
  n <- 1000
  x <- cumsum(rnorm(n))+rnorm(n)
  y <- fft(x)
  y[20:(length(y)-19)] <- 0
  y <- Re(fft(y, inverse=T)/length(y))
  op <- par(mfrow=c(3,1), mar=c(3,4,2,2)+.1)
  plot(x, type='l', 
       main="FFT: Removing the high frequencies from a signal")
  lines(y, col='red', lwd=3)
  plot(Mod(fft(x)[1: ceiling((length(x)+1)/2) ]), 
       type='l', ylab="FFT")
  plot(Mod(fft(y)[1: ceiling((length(y)+1)/2) ]), 
       type='l', ylab="Truncated FFT")
  par(op)
dev.off()

png(file="g1296.png", width=600, height=600)
  n <- 1000
  x <- cumsum(rnorm(n))+rnorm(n)
  plot(x, type='l', ylab="", 
       main="FFT: Removing more and more high frequencies")
  for (i in 1:10) {
    y <- fft(x)
    y[(1+i):(length(y)-i)] <- 0
    y <- Re(fft(y, inverse=T)/length(y))
    lines(y, col=rainbow(10)[i])
  }
dev.off()

png(file="g1297.png", width=600, height=600)
  n <- 1000
  x <- c(129:256,1:128)/256
  y <- fft(x)
  y[20:(length(y)-19)] <- 0
  y <- Re(fft(y, inverse=T)/length(y))
  plot(x, type='l', 
       ylim=c(-.2,1.2), 
       xlab="", ylab="", 
       main="Gibbs phenomenon")
  lines(y, col='red', lwd=3)
dev.off()

png(file="g1298.png", width=600, height=600)
  x <- co2
  y <- fft(x)
  y[20:(length(y)-19)] <- 0
  y <- Re(fft(y, inverse=T)/length(y))
  plot(x, type='l')
  lines(y, col='red', lwd=3)
dev.off()

png(file="g1299.png", width=600, height=600)
  p <- predict(lm(co2~time(co2)))
  x <- co2 - p
  y <- fft(x)
  y[20:(length(y)-19)] <- 0
  y <- Re(fft(y, inverse=T)/length(y))
  plot(x+p, type='l')
  lines(y+p, col='red', lwd=3)
  lines(x-y+mean(x+p), col='green', lwd=3)
dev.off()

png(file="g1300.png", width=600, height=600)
  x <- co2
  y <- fft(x)
  k <- 20
  n <- length(y)
  y[1:(k+1)] <- 0
  y[(n-k):n] <- 0
  y <- Re(fft(y, inverse=T)/length(y))
  plot(x, type='l')
  lines(y+mean(x), col='red', lwd=3)
dev.off()

png(file="g1301.png", width=600, height=600)
  library(sound)
  x <- loadSample("sample.wav")
  plot(x)
dev.off()

png(file="g1302.png", width=600, height=600)
  n <- length(x$sound)
  n <- round(n/3)
  y <- x$sound[ n:(n+2000) ]
  n <- length(y)
  op <- par(mfrow=c(2,1), mar=c(2,4,2,2)+.1)
  plot(y, type='l')
  #plot(Mod(fft(y)[1: ceiling((length(y)+1)/2) ]), type='l')
  plot(Mod(fft(y)[1:100]), type='l')
dev.off()

png(file="g1303.png", width=600, height=600)
  N <- 1024
  k <- 6
  x <- ( (1:N) - N/2 ) * 2 * pi * k / N
  y <- ifelse( x>0, sin(x), sin(3*x) )
  plot(y, type='l',
       xlab="", ylab="")
dev.off()

png(file="g1304.png", width=600, height=600)
  # With the "wd" function in the "wavethresh" package
  library(wavethresh)
  y <- ifelse( x>0, sin(x), sin(10*x) )
  plot(wd(y), main="")
dev.off()

png(file="g1305.png", width=600, height=600)
  father <- function (x) { 
    ifelse( x<0 | x>1, 0, 1 ) 
  }
  mother <- function (x) { 
    ifelse( x<0 | x>1, 
            0, 
            ifelse(x<.5, 1, -1)) 
  }
  jk <- function (f,x,j,k) { 
    f(2^j * x - k) 
  }

  op <- par(mfrow=c(4,3), mar=c(0,0,0,0))
  for (j in 1:3) {
    curve(jk(father,x,0,j),
          xlim=c(-.2,3.2), ylim=c(-1,1), 
          xlab="", ylab="", 
          axes=F, 
          lwd=3, col="blue"
    )
    abline(h=0,lty=3)
    abline(v=0:4,lty=3)
  }
  for (i in 1:3) {
    for (j in 1:3) {
      curve(jk(mother,x,i,j),
            xlim=c(-.2,3.2), ylim=c(-1,1), 
            xlab="", ylab="", 
            axes=F, 
            lwd=3, col="blue"
      )
      abline(h=0,lty=3)
      abline(v=0:4,lty=3)
      box()
    }
  }
dev.off()

png(file="g1306.png", width=400, height=400)
  op <- par(mar=c(0,0,3,0))
  curve(jk(father,x,0,1),
        xlim=c(-.2,3.2),
        ylim=c(-1,1),
        xlab="", ylab="", axes=F, 
        lwd=3, col="blue",
        main="phi, father wavelet")
  abline(h=0,lty=3)
  #abline(v=0:4,lty=3)
  par(op)
dev.off()

png(file="g1307.png", width=400, height=400)
  op <- par(mar=c(0,0,3,0))
  curve(jk(mother,x,0,1),
        xlim=c(-.2,3.2),
        ylim=c(-1,1),
        xlab="", ylab="", axes=F, lwd=3, col="blue",
        main="psi, mother wavelet")
  abline(h=0,lty=3)
  #abline(v=0:4,lty=3)
  par(op)
dev.off()

png(file="g1308.png", width=600, height=600)
  N <- 1024
  k <- 6
  x <- ( (1:N) - N/2 ) * 2 * pi * k / N
  y <- ifelse( x>0, sin(x), sin(3*x) )
  r <- wd(y)
  draw(r, col="blue", lwd=3, main="")
  abline(h=0, lty=3)
dev.off()

png(file="g1309.png", width=600, height=600)
  number.of.runs <- function (x) {
    1+sum(abs(diff(as.numeric(x))))
  }
  my.runs.test <- function (x, R=999) {
    if( is.numeric(x) )
      x <- factor(sign(x), levels=c(-1,1))
    if( is.logical(x) )
      x <- factor(x, levels=c(FALSE,TRUE))
    if(!is.factor(x))
      stop("x should be a factor")

    # Non-parametric (permutation) test
    n <- length(x)
    res <- rep(NA, R)
    for (i in 1:R) {
      res[i] <- number.of.runs(
        x[sample(1:n,n,replace=F)]
      )
    }
    t0 <- number.of.runs(x)
    n1 <- 1+sum(t0<=res)
    n2 <- 1+sum(t0>=res)
    p <- min( n1/R, n2/R )*2
    p <- min(1,p) # If more than half the values are identical

    # Parametric test, based on a formula found on the
    # internet...
    # People believe that Z follows a gaussian distribution
    # (this is completely wrong if the events are rare -- I
    # had first used it with mutations on a DNA sequence...)
    n1 <- sum(x==levels(x)[1])
    n2 <- sum(x==levels(x)[2])
    r <- number.of.runs(x)
    mr <- 2*n1*n2/(n1+n2) + 1
    sr <- sqrt( 2*n1*n2*(2*n1*n2-n1-n2)/
            (n1+n2)^2/(n1+n2-1) )
    z <- (r-mr)/sr
    pp <- 2*min(pnorm(z), 1-pnorm(z))

    r <- list(t0=t0, t=res, R=R,
              p.value.boot=p,
              n1=n1, n2=n2, r=r, mr=mr, sr=sr, z=z,
              p.value.formula=pp)
    class(r) <- "nstest"
    r
  }

  print.nstest <- function (d) {
    cat("Runs test\n");
    cat("  NS = ")
    cat(d$t0)
    cat("\n  p-value (")
    cat(d$R)
    cat(" samples) = ")
    cat(round(d$p.value.boot,digits=3))
    cat("\n")
    cat("  theoretical p-value = ")
    cat(d$p.value.formula)
    cat("\n")
  }

  plot.statistic <- function (t0, t, ...) {
    xlim <- range(c(t,t0))
    hist(t, col='light blue', xlim=xlim, ...)
    points(t0, par("usr")[4]*.8,
           type='h', col='red', lwd=5)
    text(t0, par("usr")[4]*.85, signif(t0,3))
  }

  plot.nstest <- function (
      d, main="Runs test",
      ylab="effectif", ...
    ) {
    plot.statistic(d$t0, d$t, main=main,
                   xlab=paste("Runs, p =",signif(d$p.value.boot,3)), 
                   ...)
  }

  # Example
  data(EuStockMarkets)
  x <- EuStockMarkets[,1]
  x <- diff(log(x))
  i <- abs(x)>median(abs(x))
  plot(my.runs.test(i))
dev.off()

png(file="g1310.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (k in 1:4) {
    x <- EuStockMarkets[,k]
    x <- diff(log(x))
    i <- abs(x)>median(abs(x))
    plot(my.runs.test(i))
  }
  par(op)
dev.off()

png(file="g1311.png", width=600, height=600)
  sde.sim <- function (t, f, ...) {
    n <- length(t)
    S <- rep(NA,n)
    S[1] <- 1
    for (i in 2:n) {
      S[i] <- S[i-1] + f(S[i-1], t[i]-t[i-1], sqrt(t[i]-t[i-1])*rnorm(1), ...)
    }
    S
  } 
  a <- 0
  b <- 10
  N <- 1000
  x <- sde.sim(seq(a,b,length=N),
               function (S,dt,dX,m=1,s=1) { m * S * dt + s * S * dX })
  x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2))
  plot(x)
  plot(log(x))
  plot(diff(log(x)))
  acf(diff(log(x)))
  par(op)
dev.off()

png(file="g1312.png", width=600, height=600)
  N <- 1000
  a <- 0
  b <- 3

  op <- par(mfrow=c(3,1))
  for (i in 1:3) {
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=1,s=1) { s * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Random walk")
  }
  par(op)
dev.off()

png(file="g1313.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  for (i in 1:3) {
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=1,s=1) { m * dt + s * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Random walk with a trend")
  }
  par(op)
dev.off()

png(file="g1314.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  for (i in 1:3) {
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=1,s=1) { m * S * dt + s * S * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Lognormal random walk")
  }
  par(op)
dev.off()

png(file="g1315.png", width=600, height=600)
  b <- 10
  op <- par(mfrow=c(3,1))
  for (i in 1:3) {
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=3,s=1,n=3) { (n - m*S) * dt + s * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Mean-reverting random walk")
    abline(h=1,lty=3)
  }
  par(op)
dev.off()

png(file="g1316.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  for (i in 1:3) {
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=3,s=1,n=3) { (n - m*S) * dt + s * sqrt(S) * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Mean-reverting random walk")
    abline(h=1,lty=3)
  }
  par(op)
dev.off()

png(file="g1317.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
    x <- sde.sim(seq(a,b,length=N),
                 function (S,dt,dX,m=1,s=1) { m * S * dt + s * S * dX })
    x <- ts(x, start=a, end=b, freq=(N-1)/(b-a))
    plot(x, main="Lognormal random walk")
    return <- diff(x) / x[ -length(x) ]
    y <- return^2
    plot(y)
    lines(predict(loess(y~time(y))) ~ as.vector(time(y)), col='red', lwd=3)
  par(op)
dev.off()

png(file="g1318.png", width=600, height=600)
  n <- 200
  a <- c(1,.2,.8,.5)

  a0 <- a[1]
  a <- a[-1]
  k <- length(a)
  u <- rep(NA,n)
  u[1:k] <- a0 * rnorm(k)
  for (i in (k+1):n) {
    u[i] <- sqrt( a0 + sum(a * u[(i-1):(i-k) ]^2 )) * rnorm(1)
  }
  u <- ts(u)
  eda.ts(u)
dev.off()

png(file="g1319.png", width=600, height=600)
  h <- rnorm(1000)^2
  x <- filter(h, rep(1,50))
  x <- x[!is.na(x)]
  eda.ts(x)
dev.off()

png(file="g1320.png", width=600, height=600)
  y <- rnorm(length(x),0,sqrt(x))
  eda.ts(y)
dev.off()

png(file="g1321.png", width=600, height=600)
  h <- c(rep(1,100), rep(2,100))
  y <- ts(rnorm(length(h), 0, sd=h))
  eda.ts(y)
dev.off()

png(file="g1322.png", width=600, height=600)
  plot(abs(y))
  lines(predict(loess(abs(y)~time(y))), col='red', lwd=3)
  k <- 20
  lines(filter(abs(y), rep(1/k,k)), col='blue', lwd=3, lty=2)
dev.off()

png(file="g1323.png", width=600, height=600)
  plot.ma <- function (x, k=20, ...) {
    plot(abs(x), ...)
    a <- time(x)
    b <- predict(loess(abs(x) ~ a))
    lines(b ~ as.vector(a), col='red', lwd=3)
    k <- 20
    lines(filter(abs(x), rep(1/k,k)), col='blue', lwd=3)
  }
  plot.ma(u)  
dev.off()

png(file="g1324.png", width=600, height=600)
  n <- 1000
  op <- par(mfrow=c(4,1), mar=c(2,4,2,2))
  plot.ma(ts(rnorm(n)))
  plot.ma(arima.sim(list(ar = c(.8,.1)), n))
  plot.ma(arima.sim(list(ma = c(.8,.1)), n))
  plot.ma(arima.sim(list(ma = c(.8,.4,.1), ar = c(.8,.1)), n))
  par(op)
dev.off()

png(file="g1325.png", width=600, height=600)
  op <- par(mfrow=c(1,4))
  data(EuStockMarkets)
  for (a in 1:4) {
    x <- EuStockMarkets[,a]
    x <- diff(log(x))

    n <- length(x)
    s <- rep(NA, n+1)
    s[ which(x>0) + 1 ] <- "+"
    s[ which(x<0) + 1 ] <- "-"
    i <- which( !is.na(s) )
    s <- factor(s[i])
    x <- x[i]
    boxplot(log(abs(x))~s, col='pink')
  }
  par(op)
dev.off()

png(file="g1326.png", width=600, height=600)
  data(EuStockMarkets)
  x <- EuStockMarkets[,1]
  op <- par(mfrow=c(3,1))
  plot(x, main="An index")
  y <- diff(log(x))
  plot(abs(y), main="Volatility")
  k <- 30
  z <- filter(abs(y),rep(1,k)/k)
  plot(z, ylim=c(0,max(z,na.rm=T)), col='red', type='l', 
       main="smoothed volatility (30 days)")
  par(op)
dev.off()

png(file="g1327.png", width=600, height=600)
  n <- 200
  v <- function (t) { .1*(.5 + sin(t/50)^2) }
  x <- 1:n
  y <- rnorm(n) * v(x)
  y <- ts(y)
  op <- par(mfrow=c(3,1))
  plot(ts(cumsum(y)), main="A simulation")
  plot(abs(y), main="volatility")
  plot(v(x)~x, 
       ylim=c(0,.3),
       type='l', lty=2, main="smoothed volatility")
  k <- c(5,10,20,50)
  col <- rainbow(length(k))
  for (i in 1:length(k)) {
    z <- filter(abs(y),rep(1,k[i])/k[i])
    lines(z, col=col[i])
  }
  par(op)
dev.off()

png(file="g1328.png", width=600, height=600)
  library(tseries)
  x <- EuStockMarkets[,1]
  y <- diff(x)/x
  r <- garch(y)
  # plot(r)   The plot function is only for interactive use...
  op <- par(mfrow=c(2,1))
  plot(y, main = r$series, ylab = "Series")
  plot(r$residuals, main = "Residuals", ylab = "Series")
  par(op)
dev.off()

png(file="g1329.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  hist(y, 
       main = paste("Histogram of", r$series), 
       xlab = "Series")
  hist(r$residuals, 
       main = "Histogram of Residuals", 
       xlab = "Series")
  par(op)
dev.off()

png(file="g1330.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  qqnorm(y, 
         main = paste("Q-Q Plot of", r$series), 
         xlab = "Gaussian Quantiles")
  qqnorm(r$residuals, 
         main = "Q-Q Plot of Residuals", 
         xlab = "Normal Quantiles")
  par(op)
dev.off()

png(file="g1331.png", width=600, height=600)
  op <- par(mfrow=c(2,1))
  acf(y^2, 
      main = paste("ACF of Squared", r$series))
  acf(r$residuals^2, 
      main = "ACF of Squared Residuals", 
      na.action = na.remove)
  par(op)
dev.off()

png(file="g1332.png", width=600, height=600)
  data(EuStockMarkets)
  x <- EuStockMarkets[,1]
  y <- EuStockMarkets[,2]
  acf(ts.union(x,y),lag.max=100)
dev.off()

png(file="g1333.png", width=600, height=600)
  x <- diff(x)
  y <- diff(y)
  acf(ts.union(x,y),lag.max=100)
dev.off()

png(file="g1334.png", width=600, height=600)
  f <- function (N, alpha = 2, reverting = FALSE) {
    if (N <= 1) {
      res <- 0
    } else {
      a <- rnorm(1) * (N/2)^(1/alpha)
      res <- c( seq(0, a, length = N/2) + f(N/2, reverting = TRUE),
                seq(a, 0, length = N/2) + f(N/2, reverting = TRUE) )
    }
    if (!reverting) {
      final <- rnorm(1) * N^(1/alpha)
      res <- res + seq(0, final, length = length(res))
    }  
    res
  }
  N <- 1024
  plot( f(N), 
        type = "l", 
        xlab = "time", ylab = "x", 
        main = "power law random walk, alpha = 2" )
dev.off()

png(file="g1335.png", width=600, height=600)
  do.it <- function (alpha = 2, N=1024) {
    op <- par()
    layout(matrix(c(1,2,3, 1,2,4), nc=2))
    par(mar=c(2,4,4,2)+.1)
    x <- f(N, alpha = alpha)
    plot(x, 
         type="l", 
         ylab="", 
         main = paste("Power law random walk, alpha =", alpha))
    plot(diff(x), type="l")
    par(mar=c(5,4,4,2)+.1)
    hist(diff(x), 
         col = "light blue", 
         probability = TRUE)
    lines(density(diff(x)), 
          col="red", lwd=3 )
    qqnorm(diff(x))
    qqline(diff(x), col="red", lwd=3)
    par(op)
  }
  do.it()
dev.off()

png(file="g1336.png", width=600, height=600)
  do.it(alpha = 1.5)
dev.off()

png(file="g1337.png", width=600, height=600)
  do.it(alpha = 1)
dev.off()

png(file="g1338.png", width=600, height=600)
  do.it(alpha = 2.5)
dev.off()

png(file="g1339.png", width=600, height=600)
  do.it(alpha = 3)
dev.off()

png(file="g1340.png", width=600, height=600)
  N <- 10000
  x <- cumsum(rnorm(N))  # Random walk
  y <- spectrum(x)
dev.off()

png(file="g1341.png", width=600, height=600)
  plot(y$spec ~ y$freq, 
       xlab = "Frequency", 
       ylab = "Spectral density",
       log = "xy")
dev.off()

png(file="g1342.png", width=600, height=600)
  N <- 1000
  n <- 100
  res <- rep(NA, n)
  for (i in 1:n) {
    x <- cumsum(rnorm(N))
    y <- spectrum(x, plot = FALSE)
    res[i] <- lm( log(y$spec) ~ log(y$freq) )$coef[2]
  }
  summary(res)
  hist(res, col="light blue",
       main = "Exponent of the spectral density decay")
dev.off()

png(file="g1343.png", width=600, height=600)
  N <- 10000
  x <- cumsum(rnorm(N)) + rnorm(N)
  y <- acf(x, lag.max=N/4)
dev.off()

png(file="g1344.png", width=600, height=600)
  LAG <- function (x, k = 1) {
    stopifnot( is.vector(x) )
    n <- length(x)
    stopifnot( abs(k) < n )
    if (k > 0) {
      x <- c( x[ -(1:k) ], rep(NA, k) )
    } else if ( k < 0 ) {
      k <- -k
      x <- c(rep(NA,k), x[ 1:(k-n) ])
    }
    x
  }

  x <- as.vector(sunspots)

  # Delay plots
  op <- par(mfrow=c(3,3))
  for (k in 1:9) {
    plot( LAG(x, k), x )
  }
  par(op)
dev.off()

png(file="g1345.png", width=600, height=600)
  # Principal Components plots
  N <- 20
  d <- x
  for (k in 1:N) {
    d <- cbind(d, LAG(x,k))
  }
  d <- d[ 1:(dim(d)[1]-N), ]
  r <- prcomp(d)
  plot(r$x[,1:2])
dev.off()

png(file="g1346.png", width=600, height=600)
  pairs(r$x[,1:4])
dev.off()

png(file="g1347.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1 + 1:9) {
    plot(r$x[,c(1,k)])
  }
  par(op)
dev.off()

png(file="g1348.png", width=600, height=600)
  recurrence_plot <- function (x, ...) {
    image(outer(x, x, function (a, b) abs(a-b)), ...)
    box()
  }
  N <- 500
  recurrence_plot( sin(seq(0, 10*pi, length=N)),
                   main = "Recurrence plot: sine")
dev.off()

png(file="g1349.png", width=600, height=600)
  recurrence_plot( rnorm(100), 
                   main = "Recurrence plot: noise" )
dev.off()

png(file="g1350.png", width=600, height=600)
  recurrence_plot( cumsum(rnorm(200)),
                   main = "Recurrence plot: random walk")
dev.off()

png(file="g1351.png", width=600, height=600)
  library(tseriesChaos)
  recurrence_plot(lorenz.ts[100:200],
                  main = "Recurrence plot: Lorentz attractor")
dev.off()

png(file="g1352.png", width=600, height=600)
  # Thresholded recurrence plot
  thresholded_recurrence_plot <- function (
    x, 
    threshold = 0, 
    FUN = function (x) x, 
    ...
  ) {
    image(-outer(
        x, x, 
        function (a, b) 
          ifelse(FUN(a-b)>threshold,1,0)
      ), 
      ...)
    box()
  }
  thresholded_recurrence_plot(
    lorenz.ts[1:100],
    0,
    main = "Thresholded recurrence plot"
  )
dev.off()

png(file="g1353.png", width=600, height=600)
  thresholded_recurrence_plot(
    lorenz.ts[1:100],
    5,
    abs,
    main = "Thresholded recurrence plot"
  )
dev.off()

png(file="g1354.png", width=600, height=600)
  # Recurrence plot
  recurrence_plot <- function (x, m, ...) {
    stopifnot (m >= 1, m == floor(m))
    res <- outer(x, x, function (a,b) (a-b)^2)
    i <- 2
    LAG <- function (x, lag) {
      stopifnot(lag > 0)
      if (lag >= length(x)) {
        rep(NA, length(x))
      } else {
        c(rep(NA, lag), x[1:(length(x)-lag)])
      }
    }
    while (i <= m) {
      res <- res + outer(LAG(x,i-1), LAG(x,i-1), function (a,b) (a-b)^2)
      i <- i + 1
    }
    res <- sqrt(res)
    if (m>1) {
      res <- res[ - (1:(m-1)), ] [ , - (1:(m-1)) ]
    }
    image(res, ...)
    box()
  }
  library(tseriesChaos)
  recurrence_plot(lorenz.ts[1:200], m=1,
                  main = "1-dimensional recurrence plot")
dev.off()

png(file="g1355.png", width=600, height=600)
  recurrence_plot(lorenz.ts[1:200], m=2)
dev.off()

png(file="g1356.png", width=600, height=600)
  recurrence_plot(lorenz.ts[1:200], m=10)
dev.off()

png(file="g1357.png", width=600, height=600)
  # A more complete function
  recurrence_plot <- function (
    x, 
    m=1,           # Dimension of the embedding
    t=1,           # Lag used to define this embedding
    epsilon=NULL,  # If non-NULL, threshold
    box=TRUE, ...
  ) {
    stopifnot( length(m) == 1, m >= 1, m == floor(m), 
               length(t) == 1, t >= 1, t == floor(t),
               is.null(epsilon) || (
                 length(epsilon) == 1 && epsilon > 0 ) )
    stopifnot( length(x) > m * t )
    res <- outer(x, x, function (a,b) (a-b)^2)
    i <- 2
    LAG <- function (x, lag) {
      stopifnot(lag > 0)
      if (lag >= length(x)) {
        rep(NA, length(x))
      } else {
        c(rep(NA, lag), x[1:(length(x)-lag)])
      }
    }
    while (i <= m) {
      y <- LAG(x,t*(i-1))
      res <- res + outer(y, y, function (a,b) (a-b)^2)
      i <- i + 1
    }
    res <- sqrt(res)
    if (!is.null(epsilon)) {
      res <- res > epsilon
    }
    if (m>1) {
      # TODO: Check this...
      res <- res[ - (1:(t*(m-1))), ] [ , - (1:(t*(m-1))) ]
    }
    image(res, ...)        
    if (box) {
      box()
    }
  }
  library(tseriesChaos)
  recurrence_plot(lorenz.ts[1:200], m=10)
  title("Recurrence plot")
dev.off()

png(file="g1358.png", width=600, height=600)
  op <- par(mfrow=c(5,5), mar=c(0,0,0,0), oma=c(0,0,2,0))
  for (i in 1:5) {
    for (j in 1:5) {
      recurrence_plot(lorenz.ts[1:200], m=i, t=j, 
                      axes=FALSE)
    }
  }
  par(op)
  mtext("Recurrence plots", line=3, font=2, cex=1.2)
dev.off()

png(file="g1359.png", width=600, height=600)
  phase_plane_plot <- function (
    x, 
    col=rainbow(length(x)-1), 
    xlab = "x", ylab = "dx/dt", 
    ...) {
    plot( x[-1], diff(x), col = col, 
          xlab = xlab, ylab = ylab, ... )
  }
  phase_plane_plot( sin(seq(0, 20*pi, length=200)), pch=16 )
dev.off()

png(file="g1360.png", width=600, height=600)
  x <- c(sin(seq(0, 5*pi, length=500)), 
         sin(seq(0, 5*pi, length=1000)) + .2*rnorm(1000),
         sin(seq(0, 2*pi, length=500)^2))
  phase_plane_plot(x)
dev.off()

png(file="g1361.png", width=600, height=600)
  library(tseriesChaos)
  phase_plane_plot(lorenz.ts, pch=16,
                   main = "Phase plot")
dev.off()

png(file="g1362.png", width=600, height=600)
  phase_plane_plot(lorenz.ts, type="l",
                   main = "Phase plot")
dev.off()

png(file="g1363.png", width=600, height=600)
  library(pixmap)
  z <- read.pnm("2006-08-27_Hello.pgm") # Created (by hand) with The Gimp
  z <- z@grey
  d <- cbind(
    x = col(z)[ ! z ],
    y = -row(z)[ ! z ]
  )
  N <- 10000
  d <- d[sample(1:nrow(d), N, replace = TRUE),]
  d <- d + rnorm(2*N)
  #plot(d)
  d <- apply(d, 2, scale)
  explode <- function (
    d,
    FUN = function (x) { rank(x, na.last="keep") / length(x) }
  ) {
    # Convert to polar coordinates
    d <- data.frame( 
      rho = sqrt(d[,1]^2 + d[,2]^2),
      theta = atan2(d[,2], d[,1])
    )
    d$rho <- FUN(d$rho)
    # Convert back to cartesian coordinates
    d <- cbind(
      x = d$rho * cos(d$theta),
      y = d$rho * sin(d$theta)
    )
    d
  }
  #plot(explode(d))
  d <- explode(d, FUN = function (x) x^4)
  d <- apply(d, 2, function (x) (x - min(x))/diff(range(x)))
  d <- rbind(d, matrix(rnorm(2*N), nc=2))

  # Exercice: Find the word in the following cloud of points...
  op <- par(mfrow=c(2,2), mar=c(.1,.1,.1,.1))
  plot(d, axes = FALSE)
  box()
  plot( rank(d[,1]), rank(d[,2]), axes = FALSE )
  box()
  plot(explode(d), axes = FALSE)
  box()
  plot(explode(d, atan), axes = FALSE)
  box()
  par(op)
dev.off()

png(file="g1364.png", width=600, height=600)
  one_dimensional_fish_eye <- function (x1, x2, y, method="natural") {
    n <- length(y)
    x <- seq(min(x1), max(x1), length=n)
    x3 <- splinefun(x1, x2, method = method)(x)
    if (! all(x3 == sort(x3))) {
      warning("Non monotonic transformation!")
    }
    d <- cbind(x=x3, y=y)
    op1 <- par(mar=c(.1,.1,.1,.1))
    plot(d, type="l", lwd=3, axes = FALSE)
    box()
    abline(v=d[seq(0,length(y),by=ceiling(length(y)/50)),1])
    op2 <- par(fig=c(.02,.2,.8,.98), mar=c(0,0,0,0), new=TRUE)
    plot(x, x3, type = "l", lwd = 3, axes = FALSE)
    polygon(rep(par("usr")[1:2], 2)[c(1,2,4,3)], 
            rep(par("usr")[3:4], each=2), 
            border = NA, col = "white")
    lines(x, x3, type = "l", lwd = 3, col="blue")
    box(lwd=3, col="blue")
    par(op2)
    par(op1)
  }
  library(Ecdat) # Some econometric data
  data(DM)
  y <- DM$s
  # More details in the middle
  one_dimensional_fish_eye(
    seq(0, 1, length = 4),
    c(0, .2, .8, 1),
    y
  )
dev.off()

png(file="g1365.png", width=600, height=600)
  # More details on the left 
  one_dimensional_fish_eye(
    c(0, .33, .67, 1),
    c(0, .6, .9, 1),
    y
  )
dev.off()

png(file="g1366.png", width=600, height=600)
  # More details on the right
  one_dimensional_fish_eye(
    seq(0, 1, length=4),
    c(0, .1, .4, 1),
    y
  )
dev.off()

png(file="g1367.png", width=600, height=600)
  two_time_scales <- function (x1, x2, y) {
    stopifnot( length(x1) == length(y),
               length(x2) == length(y),
               x1 == sort(x1),
               x2 == sort(x2) )
    op <- par(mfrow=c(2,1))    
    plot(x1, y, type="l", xlab="", axes=FALSE)
    box()
    axis(2, lwd=2)
    mtext(side=1, "x1", line=2.5, col="blue", font=2)
    mtext(side=3, "x2", line=2.5, col="red", font=2)
    x2lab <- pretty(x2, 10)
    axis(1, col="blue", lwd=2)
    axis(3, at = approx(x2, x1, x2lab)$y, labels=x2lab,
         col="red", lwd=2)
    plot(x2, y, type="l", axes=FALSE,
         xlab="")
    box()
    mtext(side=1, "x1", line=2.5, col="blue", font=2)
    mtext(side=3, "x2", line=2.5, col="red", font=2)
    x1lab <- pretty(x1, 10)
    axis(1, at = approx(x1,x2,x1lab)$y, labels=x1lab,
         col="blue", lwd=2)
    axis(2, lwd=2)
    axis(3, col="red", lwd=2)
    par(op)
  }
  N <- 100
  x1 <- seq(0, 2, length=N)
  x2 <- x1^2
  y <- x2
  two_time_scales(x1,x2,y)
dev.off()

png(file="g1368.png", width=600, height=600)
  x <- read.csv(gzfile("2006-08-27_tick_data.csv.gz"))
  op <- par(mfrow=c(3,1))
  x$DateTime <- as.POSIXct(paste(as.character(x$Date),
                                 as.character(x$Time)))
  x <- x[!is.na(x$TradePrice),]
  plot(TradePrice ~ DateTime,
       data = x,
       type = "l",
       xlab = "Clock time",
       main = "Is time a good choice for the X axis?")
  plot(x$TradePrice,
       type = "l",
       ylab = "TradePrice",
       xlab = "Transaction time")
  coalesce <- function (x,y) ifelse(!is.na(x),x,y)
  plot(cumsum(coalesce(x$TradeSize,0)),
       x$TradePrice,
       type = "l",
       xlab = "Volume time",
       ylab = "TradePrice")
  par(op)
dev.off()

png(file="g1369.png", width=600, height=600)
  two_time_scales(
    cumsum(coalesce(x$TradeSize,0)),
    as.numeric(x$DateTime - x$DateTime[1]) / 3600,
    x$TradePrice
  )
dev.off()

png(file="g1370.png", width=600, height=600)
  data_driven_time_warp <- function (y) {
    cbind(
      x = cumsum(c(0, abs(diff(y)))),
      y = y
    )
  }
  library(Ecdat) # Some econometric data
  data(DM)
  y <- DM$s
  i <- seq(1,length(y),by=10)
  op <- par(mfrow=c(2,1), mar=c(.1,.1,.1,.1))
  plot(y, type="l", axes = FALSE)
  abline(v=i, col="grey")
  lines(y, lwd=3)
  box()
  d <- data_driven_time_warp(y)
  plot(d, type="l", axes=FALSE)
  abline(v=d[i,1], col="grey")
  lines(d, lwd=3)
  box()
  par(op)
dev.off()

png(file="g1371.png", width=600, height=600)
  # Building a Markov chain
  markov <- function (x) {
    x <- strsplit(x,'')[[1]]
    x <- factor(x)
    aa <- strsplit("ACDEFGHIKLMNPSQRTVWY01",'')[[1]]
    n <- length(aa)
    m <- matrix(0, nr=n, nc=n)
    colnames(m) <- aa
    rownames(m) <- aa
    m["1","1"] <- 1
    m["0", x[1]] <- m["0", x[1]] +1
    for (i in 1:(length(x)-1)) {
      m[ x[i], x[i+1] ] <- m[ x[i], x[i+1] ] +1
    }
    m[ x[length(x)], "1" ] <- m[ x[length(x)], "1" ] +1
    # This is a contingency matrix, we want a probability
    # matrix, where the sum of each row is 1
    m <- m +.001   # Pas n
    m <- m/apply(m,1,sum)
    print(round(m, digits=2))
    invisible(m)
  }
  x <-"MAKGVAVLNSSEGVTGTIFFTQEGDGVTTVSGTVSGLKPGLHGFHVHALGDTTNGCMSTGPHFNPDGKTHGAPEDANRHAGDLGNITVGDDGTATFTITDCQIPLTGPNSIVGRAVVVHADPDDLGKGGHELSLATGNAGGRVACGIIGLQG"
  m <- markov(x)
  image(m)
dev.off()

png(file="g1372.png", width=600, height=600)
  n <- 2000
  x <- 1:n
  y <- sin(x/10) + cos(pi*x/20) + rnorm(n)
  op <- par(mfrow=c(4,1))
  plot(ts(y))
  acf(y)
  pacf(y)
  spectrum(y, col=par('fg'))
  abline(v=c(1/40,1/20/pi), lty=3)
  par(op)
dev.off()

png(file="g1373.png", width=600, height=800)
  see.ts <- function (name, ma=NULL, ar=NULL, d=0, n=2000) {
    order=c(length(ar), d, length(ma))
    x <- arima.sim(list(ma=ma, ar=ar, order=order), n)
    op <- par(mfrow=c(4,1))
    plot(x, main=name)
    acf(x)
    pacf(x)
    spectrum(x, spans=10, col=par('fg'))
    par(op)
  }
  n <- 200
  see.ts("MA(1) theta_1=.9", .9)
dev.off()

png(file="g1374.png", width=600, height=800)
  see.ts("MA(1) theta_1=.5", .5)
dev.off()

png(file="g1375.png", width=600, height=800)
  see.ts("MA(1) theta_1=.1", .1)
dev.off()

png(file="g1376.png", width=600, height=800)
  see.ts("MA(1) theta_1=-.9", -.9)
dev.off()

png(file="g1377.png", width=600, height=800)
  see.ts("AR(1) phi_1=.9", 0, .9)
dev.off()

png(file="g1378.png", width=600, height=800)
  see.ts("AR(1) phi_1=.8", 0, .9)
dev.off()

png(file="g1379.png", width=600, height=800)
  see.ts("AR(1) phi_1=.5", 0, .5)
dev.off()

png(file="g1380.png", width=600, height=800)
  see.ts("AR(1) phi_1=.1", 0, .1)
dev.off()

png(file="g1381.png", width=600, height=800)
  see.ts("AR(1) phi_1=-.9", 0, -.9)
dev.off()

png(file="g1382.png", width=600, height=800)
  see.ts("AR(1) phi_1=-.8", 0, -.8)
dev.off()

png(file="g1383.png", width=600, height=800)
  see.ts("ARMA(1,1) theta_1=.9 phi_1=.9", .9, .9)
dev.off()

png(file="g1384.png", width=600, height=800)
  see.ts("ARMA(1,1) theta_1=.9 phi_1=-.9", .9, -.9)
dev.off()

png(file="g1385.png", width=600, height=800)
  see.ts("ARMA(1,1) theta_1=-.9 phi_1=.9", -.9, .9)
dev.off()

png(file="g1386.png", width=600, height=800)
  see.ts("ARMA(1,1) theta_1=-.9 phi_1=-.9", -.9, -.9)
dev.off()

png(file="g1387.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(AirPassengers)
  plot(AirPassengers)
  acf(AirPassengers)
  pacf(AirPassengers)
  spectrum(AirPassengers);
  par(op)
dev.off()

png(file="g1388.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(austres)
  plot(austres)
  acf(austres)
  pacf(austres)
  spectrum(austres);
  par(op)
dev.off()

png(file="g1389.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(beavers)
  plot(ts(beaver1$temp))
  acf(beaver1$temp)
  pacf(beaver1$temp)
  spectrum(beaver1$temp);
  par(op)
dev.off()

png(file="g1390.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(BJsales)
  plot(BJsales)
  acf(BJsales)
  pacf(BJsales)
  spectrum(BJsales);
  par(op)
dev.off()

png(file="g1391.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(EuStockMarkets)
  plot(EuStockMarkets[,1])
  acf(EuStockMarkets[,1])
  pacf(EuStockMarkets[,1])
  spectrum(EuStockMarkets[,1]);
  par(op)
dev.off()

png(file="g1392.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(JohnsonJohnson)
  plot(JohnsonJohnson)
  acf(JohnsonJohnson)
  pacf(JohnsonJohnson)
  spectrum(JohnsonJohnson);
  par(op)
dev.off()

png(file="g1393.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(LakeHuron)
  plot(LakeHuron)
  acf(LakeHuron)
  pacf(LakeHuron)
  spectrum(LakeHuron);
  par(op)
dev.off()

png(file="g1394.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(lh)
  plot(lh)
  acf(lh)
  pacf(lh)
  spectrum(lh);
  par(op)
dev.off()

png(file="g1395.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(lynx)
  plot(lynx)
  acf(lynx)
  pacf(lynx)
  spectrum(lynx);
  par(op)
dev.off()

png(file="g1396.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(Nile)
  plot(Nile)
  acf(Nile)
  pacf(Nile)
  spectrum(Nile);
  par(op)
dev.off()

png(file="g1397.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(nottem)
  plot(nottem)
  acf(nottem)
  pacf(nottem)
  spectrum(nottem);
  par(op)
dev.off()

png(file="g1398.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  #data(sunspot, package=ts) # Il y en a aussi dans "boot"...
  data(sunspot)
  plot(sunspot.month)
  acf(sunspot.month)
  pacf(sunspot.month)
  spectrum(sunspot.month);
  par(op)
dev.off()

png(file="g1399.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(treering)
  plot(treering)
  acf(treering)
  pacf(treering)
  spectrum(treering);
  par(op)
dev.off()

png(file="g1400.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(UKDriverDeaths)
  plot(UKDriverDeaths)
  acf(UKDriverDeaths)
  pacf(UKDriverDeaths)
  spectrum(UKDriverDeaths);
  par(op)
dev.off()

png(file="g1401.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(UKLungDeaths)
  plot(ldeaths)
  acf(ldeaths)
  pacf(ldeaths)
  spectrum(ldeaths);
  par(op)
dev.off()

png(file="g1402.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(UKgas)
  plot(UKgas)
  acf(UKgas)
  pacf(UKgas)
  spectrum(UKgas);
  par(op)
dev.off()

png(file="g1403.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(USAccDeaths)
  plot(USAccDeaths)
  acf(USAccDeaths)
  pacf(USAccDeaths)
  spectrum(USAccDeaths);
  par(op)
dev.off()

png(file="g1404.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  data(WWWusage)
  plot(WWWusage)
  acf(WWWusage)
  pacf(WWWusage)
  spectrum(WWWusage);
  par(op)
dev.off()
detach.everything()

png(file="g1405.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  n <- 20
  lambda <- rep(2,n)
  x <- seq(0,2,length=n)
  plot(lambda ~ x, type='l', main=expression(lambda))
  plot(lambda*x ~ x, type='l', main=expression(Lambda))
  plot(exp(-lambda*x) ~ x, type='l', main="S")
  par(op)
dev.off()

png(file="g1406.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  n <- 20
  alpha <- 1
  g <- rep(2,n)
  x <- seq(0,2,length=n)
  plot(g * alpha * x^(g-1) ~ x, type='l', main=expression(lambda (gamma==2)))
  plot(alpha * x^g ~ x, type='l', main=expression(Lambda))
  plot(exp(-alpha*x^g) ~ x, type='l', main="S")
  par(op)
dev.off()

png(file="g1407.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  n <- 20
  alpha <- 1
  g <- rep(.5,n)
  x <- seq(0,2,length=n)
  plot(g * alpha * x^(g-1) ~ x, type='l', main=expression(lambda (gamma==.5)))
  plot(alpha * x^g ~ x, type='l', main=expression(Lambda))
  plot(exp(-alpha*x^g) ~ x, type='l', main="S")
  par(op)
dev.off()

png(file="g1408.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  n <- 200
  x <- seq(0,5,length=n)
  plot(1-plnorm(x) ~ x, type='l', main="S")
  L <- function (x) { -log(1-plnorm(x)) }
  plot(L(x) ~ x, type='l', main=expression(Lambda))
  h <- .001
  plot( (L(x+h)-L(x))/h ~ x, type='l', main=expression(lambda))
  par(op)
dev.off()

png(file="g1409.png", width=600, height=600)
  set.seed(87638)
  library(survival)
  # survfit <- survival:::survfit # Incompatibility Design/survival?
  try( detach("package:Design") )
  n <- 200
  x <- rweibull(n,.5)
  y <- rexp(n,1/mean(x))
  s <- Surv(ifelse(x<y,x,y), x<=y)
  plot(s) # not insightful
dev.off()

png(file="g1410.png", width=600, height=600)
  plot(survfit(s))
  lines(survfit(s, type='fleming-harrington'), col='red')
  r <- survfit(s)
  lines( 1-pweibull( r$time, .5 ) ~ r$time, lty=3, lwd=3, col='blue' )
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c("Kaplan-Meier", "Fleming-Harrington", "Theoretical survival function"),
          lwd=c(1,1,3), lty=c(1,1,3),
          col=c(par("fg"), 'red', 'blue'))
  title(main="Survival function (Kaplan-Meier estimator)")
dev.off()

png(file="g1411.png", width=600, height=600)
  op <- par(mfrow=c(3,1))
  r <- survfit(s)
  plot(r$surv ~ r$time, type='l', main="S")
  curve( 1-pweibull(x,.5,1), col='red', lty=2, add=T )
  plot(-log(r$surv) ~ r$time, type='l', main=expression(Lambda))
  curve( -log(1-pweibull(x,.5,1)), col='red', lty=2, add=T )
  # Before derivating, we smooth Lambda
  a <- -log(r$surv)
  b <- r$time
  # library(modreg) # Merged into stats...
  l <- loess(a~b)
  bb <- seq(min(b),max(b),length=200)
  aa <- predict(l, data.frame(b=bb))
  plot( diff(aa) ~ bb[-1], type='l', main=expression(lambda) )
  aa <- -log(1-pweibull(bb,.5,1))
  lines( diff(aa) ~ bb[-1], col='red', lty=2 )
  par(op)
dev.off()

png(file="g1412.png", width=600, height=600)
  plot(r, fun="log", main="log-survival curve")
dev.off()

png(file="g1413.png", width=600, height=600)
  plot(r, fun="event", main="cumulative events: f(y)=1-y")
dev.off()

png(file="g1414.png", width=600, height=600)
  plot(r, fun="cumhaz", main="cumulative hazard: f(y)=-log(y)=Lambda")
dev.off()

png(file="g1415.png", width=600, height=600)
  try(
    plot(r, fun="cloglog", main="complementary log-log plot")
  )
  # f(y)=log(-log(y)), log-scale on the x-axis
dev.off()

png(file="g1416.png", width=600, height=600)
  n <- 200
  x1 <- rweibull(n,.5)
  x2 <- rweibull(n,1.2)
  f <- factor( sample(1:2, n, replace=T), levels=1:2 )
  x <- ifelse(f==1,x1,x2)
  y <- rexp(n,1/mean(x))
  s <- Surv(ifelse(x<y,x,y), x<=y)
  plot(s, col=as.numeric(f))
dev.off()

png(file="g1417.png", width=600, height=600)
  plot( density(s[,1][ s[,2] == 1 & f == 1]), lwd=3,
        main="Survival analysis, one factor" )
  lines( density(s[,1][ s[,2] == 0 & f == 1]), lty=2, lwd=3 )
  lines( density(s[,1][ s[,2] == 1 & f == 2]), col='red', lwd=3 )
  lines( density(s[,1][ s[,2] == 0 & f == 2]), lty=2, col='red', lwd=3 )
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c("non censored, f = 1", "censored, f = 1",
            "non censored, f = 2", "censored, f = 2"),
          lty=c(1,2,1,2),
          lwd=1,
          col=c(par('fg'), par('fg'), 'red', 'red') )
dev.off()

png(file="g1418.png", width=600, height=600)
  plot(survfit(s ~ f), col=as.numeric(levels(f)))
dev.off()

png(file="g1419.png", width=600, height=600)
  data(lung)
  x <- Surv(lung$time, lung$status)
  plot(x)
dev.off()

png(file="g1420.png", width=600, height=600)
  plot(survfit(x))
dev.off()

png(file="g1421.png", width=600, height=600)
  data(lung)
  x <- Surv(lung$time, lung$status)
  f <- function (p,t) { dweibull(t,p[1],p[2]) }
  S <- function (p,t) { 1-pweibull(t,p[1],p[2]) }
  ll <- function (p) {
    time <- x[,1]
    status <- x[,2]
    censored <- 0
    dead <- 1
    # cat(p); cat("\n"); str(time); cat("\n");
    -2*( sum(log(f(p,time[status==dead]))) + sum(log(S(p,time[status==censored]))) )
  }

  # Estimations of the second parameter
  m <- 1     # Does not work
  m <- 100
  s <- survfit(x)
  m <- mean(s$time)
  m <- max(s$time[s$surv>.5])

  r <- optim( c(1,m), ll )

  # Plot the log-likelohood
  myOuter <- function (x,y,f) {
    r <- matrix(nrow=length(x), ncol=length(y))
    for (i in 1:length(x)) {
      for (j in 1:length(y)) {
        r[i,j] <- f(x[i],y[j])
      }
    }
    r
  }
  lll <- function (u,v) { 
    r <- ll(c(u,v)) 
    if(r==Inf) r <- NA
    r
  }
  a <- seq(1,1.6,length=50)
  b <- seq(100,700,length=50)
  ab <- myOuter(a,b,lll)
  persp(a,b,ab)
dev.off()

png(file="g1422.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in seq(0,360,length=10)[-10]) {
    persp(a,b,ab,theta=i)
  }
  par(op)
dev.off()

png(file="g1423.png", width=600, height=600)
  image(a,b,ab)
  points(r$par[1],r$par[2],lwd=3,cex=3)
dev.off()

png(file="g1424.png", width=600, height=600)
  n <- 255
  image(a,b,ab, col=topo.colors(n), breaks=quantile(ab,(0:n)/n, na.rm=T))
  points(r$par[1],r$par[2],lwd=3,cex=3)
dev.off()

png(file="g1425.png", width=600, height=600)
  image(a,b,ab, col=topo.colors(n), breaks=quantile(ab,((0:n)/n)^2,na.rm=T))
  points(r$par[1],r$par[2],lwd=3,cex=3)
dev.off()

png(file="g1426.png", width=600, height=600)
  plot(survfit(x))
  curve( 1-pweibull(x,r$par[1],r$par[2]), add=T, col='red', lwd=3, lty=2 )
dev.off()

png(file="g1427.png", width=600, height=600)
  ph.mle.weibull <- function (x) {
    f <- function (p,t) { dweibull(t,p[1],p[2]) }
    S <- function (p,t) { 1-pweibull(t,p[1],p[2]) }
    m <- mean(survfit(x)$time)
    ph.mle.optim(x,f,S,c(1,m))
  }
  ph.mle.exp <- function (x) {
    f <- function (p,t) { dexp(t,p) }
    S <- function (p,t) { 1-pexp(t,p) }
    m <- mean(survfit(x)$time)
    ph.mle.optim(x,f,S,1/m)
  }
  ph.mle.gamma <- function (x) {
    f <- function (p,t) { dgamma(t,p[1],p[2]) }
    S <- function (p,t) { 1-pgamma(t,p[1],p[2]) }
    m <- mean(survfit(x)$time)
    ph.mle.optim(x,f,S,c(1,1/m))
  }
  ph.mle.optim <- function (x,f,S,m) {
    ll <- function (p) {
      time <- x[,1]
      status <- x[,2]
      censored <- 0
      dead <- 1
      -2*( sum(log(f(p,time[status==dead]))) + sum(log(S(p,time[status==censored]))) )
    }
    optim(m,ll)
  }
  eda.surv <- function (x) {
    r <- survfit(x)
    plot(r)
    a1 <- ph.mle.exp(x)$par
    lines( 1-pexp(r$time,a1) ~ r$time, col='red' )
    a2 <- ph.mle.weibull(x)$par
    lines( 1-pweibull(r$time,a2[1],a2[2]) ~ r$time, col='green' )
    a3 <- ph.mle.gamma(x)$par
    lines( 1-pgamma(r$time,a3[1],a3[2]) ~ r$time, col='blue' )
    legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
            c(paste("Exponential(", signif(a1,2), ")", sep=''), 
              paste("Weibull(", signif(a2[1],2), ", ", signif(a2[2],2), ")", sep=''), 
              paste("Gamma(", signif(a3[1],2), ", ", signif(a3[2],2), ")", sep='')
            ),
            lwd=1, lty=1,
            col=c('red', 'green', 'blue'))
    title(main="Parametric estimation of PH models")
  }

  data(lung)
  x <- Surv(lung$time, lung$status)
  eda.surv(x)
dev.off()

png(file="g1428.png", width=600, height=600)
  x <- Surv(lung$time, lung$status)
  r <- survfit(x)
  a1 <- ph.mle.exp(x)$par
  t1 = 1-pexp(r$time,a1)
  a2 <- ph.mle.weibull(x)$par
  t2 <- 1-pweibull(r$time,a2[1],a2[2])
  a3 <- ph.mle.gamma(x)$par
  t3 <- 1-pgamma(r$time,a3[1],a3[2])
  plot( t1 ~ r$surv, col='red', xlab='sample', ylab='model')
  points( t2 ~ r$surv, col='green')
  points( t3 ~ r$surv, col='blue' )
  abline(0,1)
  legend( par("usr")[1], par("usr")[4], yjust=1, xjust=0,
          c(paste("Exponental(", signif(a1,2), ")", sep=''), 
            paste("Weibull(", signif(a2[1],2), ", ", signif(a2[2],2), ")", sep=''), 
            paste("Gamma(", signif(a3[1],2), ", ", signif(a3[2],2), ")", sep='')
          ),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="Parametric estimation of PH models")
dev.off()

png(file="g1429.png", width=600, height=600)
  plot( t1 - r$surv ~ t1, col='red', xlab='predicted values', ylab='residuals')
  points( t2 - r$surv ~ t2, col='green')
  points( t3 - r$surv ~ t3, col='blue' )
  abline(h=0, lty=3)
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c(paste("Exponential(", signif(a1,2), ")", sep=''), 
            paste("Weibull(", signif(a2[1],2), ", ", signif(a2[2],2), ")", sep=''), 
            paste("Gamma(", signif(a3[1],2), ", ", signif(a3[2],2), ")", sep='')
          ),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="Parametric estimation of PH models")  
dev.off()

png(file="g1430.png", width=600, height=600)
  plot( abs(t1 - r$surv) ~ t1, col='red', 
        xlab='predicted values', ylab=expression( abs(residuals) ))
  points( abs(t2 - r$surv) ~ t2, col='green')
  points( abs(t3 - r$surv) ~ t3, col='blue' )
  abline(h=0, lty=3)
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=1,
          c(paste("Exponental(", signif(a1,2), ")", sep=''), 
            paste("Weibull(", signif(a2[1],2), ", ", signif(a2[2],2), ")", sep=''), 
            paste("Gamma(", signif(a3[1],2), ", ", signif(a3[2],2), ")", sep='')
          ),
          lwd=1, lty=1,
          col=c('red', 'green', 'blue'))
  title(main="Parametric estimation of PH models")  
dev.off()

png(file="g1431.png", width=600, height=600)
  x <- runif(100)
  y <- runif(100)
  nearest_neighbour <- function (x, y, d=dist(cbind(x,y)), ...) {
    n <- length(x)
    stopifnot(length(x) == length(y))
    d <- as.matrix(d)
    stopifnot( dim(d)[1] == dim(d)[2] )
    stopifnot( length(x) == dim(d)[1] )
    i <- 1:n 
    j <- apply(d, 2, function (a) order(a)[2])
    segments(x[i], y[i], x[j], y[j], ...)
  }
  plot(x, y,
       main="Nearest neighbour graph",
       xlab = "", ylab = "")
  nearest_neighbour(x, y)
dev.off()

png(file="g1432.png", width=600, height=600)
  plot(x, y,
       main = "Minimum spanning tree",
       xlab = "", ylab = "")
  nearest_neighbour(x, y, lwd=10, col="grey")
  points(x,y)
  library(ape)
  r <- mst(dist(cbind(x, y)))
  i <- which(r==1, arr.ind=TRUE )
  segments(x[i[,1]], y[i[,1]], x[i[,2]], y[i[,2]],
           lwd = 2, col = "blue")
dev.off()

png(file="g1433.png", width=600, height=600)
  # Voronoi diagram
  library(tripack)
  plot(voronoi.mosaic(x, y))
  segments(x[i[,1]], y[i[,1]], x[i[,2]], y[i[,2]],
           lwd=3, col="grey")  
  points(x, y, pch=3, cex=2, lwd=3)
  box()
dev.off()

png(file="g1434.png", width=600, height=600)
  # Delaunay triangulation
  # See also the "deldir" package
  plot(tri.mesh(x,y))
  plot(voronoi.mosaic(x, y), add=T, col="grey")
  points(x, y, pch=3, cex=2, lwd=3)
dev.off()

png(file="g1435.png", width=600, height=600)
  library(tripack)
  set.seed(1)
  n <- 20
  x <- runif(n)
  y <- runif(n)
  v <- voronoi.mosaic(x, y)
  plot(v, main="Voronoi mosaic and tile centers")
  points(x,y, pch=3, cex=1.5, lwd=2)
  
  # Center of gravity of a convex polygon, given by 
  # the coordinates of its vertices.
  voronoi.center <- function (x,y) {
    stopifnot( length(x) == length(y) )
    n <- length(x)
    # A point inside the polygon
    x0 <- mean(x)
    y0 <- mean(y)
    # Reorder the vertices
    o <- order(atan2(y-y0, x-x0))
    x <- x[o]
    y <- y[o]
    # Duplicate the first point at the end, for the loop
    x <- c(x, x[1])
    y <- c(y, x[1])
    # Compute the center of gravity and the area of each triangle
    gx <- gy <- rep(NA, n)
    a <- rep(NA, n) 
    for (i in 1:n) {
      xx <- c( x0, x[i], x[i+1] )
      yy <- c( y0, y[i], y[i+1] )
      gx[i] <- mean(xx)
      gy[i] <- mean(yy)
      a[i] <- voronoi.polyarea(xx, yy)
    }
    # Compute the barycenter of those centers of gravity, with
    # weights proportionnal to the triangle areas.
    res <- c( 
       x = weighted.mean(gx, w=a),
       y = weighted.mean(gy, w=a)
     )
     attr(res, "x") <- x[1:n]
     attr(res, "y") <- y[1:n]
     attr(res, "G") <- c(x0,y0)
     attr(res, "gx") <- gx
     attr(res, "gy") <- gy
     res
  }
  voronoi.centers <- function (v) {
    ntiles <- length(v$tri$x)
    res <- matrix(NA, nc=2, nr=ntiles)
    for (i in 1:ntiles) {
      vs <- voronoi.findvertices(i, v)
      if (length(vs) > 0) {
        res[i,] <- voronoi.center(v$x[vs], v$y[vs])
      }
      else {
        res[i,] <- NA
      }
    }
    res
  }

  points( voronoi.centers(v), pch=16, col="red")
dev.off()

png(file="g1436.png", width=600, height=600)
  for (i in 1:10) {
    z <- voronoi.centers(v)
    x <- ifelse( is.na(z[,1]), x, z[,1] )
    y <- ifelse( is.na(z[,2]), y, z[,2] )
    v <- voronoi.mosaic(x, y)
  }
  plot(v, main="Centroidal Voronoi tessallation?")
  points(x,y, pch=3, cex=1.5, lwd=2)
  points( voronoi.centers(v), pch=16, col="red")
dev.off()

png(file="g1437.png", width=600, height=600)
  for (i in 1:10) {
    z <- voronoi.centers(v)
    x <- ifelse( is.na(z[,1]), x, z[,1] )
    y <- ifelse( is.na(z[,2]), y, z[,2] )
    v <- voronoi.mosaic(x, y)
  }
  plot(v, main="Centroidal Voronoi tessallation?")
  points(x,y, pch=3, cex=1.5, lwd=2)
  points( voronoi.centers(v), pch=16, col="red")
dev.off()

png(file="g1438.png", width=600, height=600)
  for (i in 1:10) {
    z <- voronoi.centers(v)
    x <- ifelse( is.na(z[,1]), x, z[,1] )
    y <- ifelse( is.na(z[,2]), y, z[,2] )
    v <- voronoi.mosaic(x, y)
  }
  plot(v, main="Centroidal Voronoi tessallation?")
  points(x,y, pch=3, cex=1.5, lwd=2)
  points( voronoi.centers(v), pch=16, col="red")
dev.off()

png(file="g1439.png", width=600, height=600)
  for (i in 1:10) {
    z <- voronoi.centers(v)
    x <- ifelse( is.na(z[,1]), x, z[,1] )
    y <- ifelse( is.na(z[,2]), y, z[,2] )
    v <- voronoi.mosaic(x, y)
  }
  plot(v, main="Centroidal Voronoi tessallation?")
  points(x,y, pch=3, cex=1.5, lwd=2)
  points( voronoi.centers(v), pch=16, col="red")
dev.off()

png(file="g1440.png", width=600, height=600)
  library(tripack)
  n <- 100    # Number of cells
  k <- 100    # Number of points in each cell
  x <- runif(k*n)
  y <- runif(k*n)
  z <- kmeans(cbind(x,y), n)
  v <- voronoi.mosaic( z$centers[,1], z$centers[,2] )
  plot(v, main="Centroidal Voronoi tessallation (k-means)")
dev.off()

png(file="g1441.png", width=600, height=600)
  plot(voronoi.mosaic( rnorm(n), rnorm(n) ),
       main = "Non-centroidal Voronoi tessallation")
dev.off()

png(file="g1442.png", width=600, height=600)
  pi.monte.carlo.plot <- function (n=10000, pch='.', ...) {
    x <- runif(n)
    y <- runif(n)
    interieur <- x^2 + y^2 <= 1
    p <- 4*sum(interieur)/n
    xc <- seq(0,1,length=200)
    yc <- sqrt(1-xc^2)
    plot( xc, yc, type='l' )
    lines( c(0,1,1,0,0), c(0,0,1,1,0) )
    abline(h=0, lty=3)
    abline(v=0, lty=3)
    points(x[interieur], y[interieur], col='red', pch=pch, ...)
    points(x[!interieur], y[!interieur], pch=pch, ...)
    title(main=paste("Monte Carlo Simulation: pi=",p,sep=''))
  }
  pi.monte.carlo.plot(100, pch='+', cex=3)
dev.off()

png(file="g1443.png", width=600, height=600)
  pi.monte.carlo.plot()
dev.off()

png(file="g1444.png", width=600, height=600)
  n <- 10
  N <- 1000
  r <- matrix(NA, nr=N, nc=6)
  for (i in 1:N) {
    x <- ifelse( runif(n)>1/6, rnorm(n,1,2), rnorm(n,-1,1) )
    r[i,] <- c( mean(x), quantile(x) )
  }
  colnames(r) <- c("mean", "0% (min)", "25%", "50% (median)", "75%", "100% (max)")
  rr <- apply(r,2,density)
  xlim <- range( sapply(rr,function(a){range(a$x)}) )
  ylim <- range( sapply(rr,function(a){range(a$y)}) )
  plot(NA, xlim=xlim, ylim=ylim, ylab='density')
  for (i in 1:6) {
    lines(rr[[i]], col=i)
  }
  legend(par('usr')[2], par('usr')[4], xjust=1, yjust=1,
         c('mean', "min", "1st quartile", "median", "3rd quartile", "max"),
         lwd=1, lty=1,
         col=1:6 )
dev.off()

png(file="g1445.png", width=600, height=600)
  op <- par(mfrow=c(2,3))
  for (i in 1:6) {
    qqnorm(r[,i], main=colnames(r)[i])
    qqline(r[,i], col='red')
    text( par("usr")[1], par("usr")[4], adj=c(-.2,2),
          round(shapiro.test(r[,i])$p.value, digits=4) )
  }
  par(op)
dev.off()

png(file="g1446.png", width=600, height=600)
  my.simulation <- function (get.sample, statistic, R) {
    res <- statistic(get.sample())
    r <- matrix(NA, nr=R, nc=length(res))
    r[1,] <- res
    for (i in 2:R) {
      r[i,] <- statistic(get.sample())
    }
    list(t=r, t0=apply(r,2,mean))
  }

  r <- my.simulation(
    function () { 
      n <- 200
      x1 <- rnorm(n)
      x2 <- rnorm(n)
      y <-  1 - x1 + 2 * x2 + rnorm(n)
     data.frame(y,x1,x2)
    },
    function (d) {
      lm(y~x1+x2, data=d)$coef
    },
    R=999
  )

  matdensityplot <- function (r, ylab='density', ...) {
    rr <- apply(r,2,density)
    n <- length(rr)
    xlim <- range( sapply(rr,function(a){range(a$x)}) )
    ylim <- range( sapply(rr,function(a){range(a$y)}) )
    plot(NA, xlim=xlim, ylim=ylim, ylab=ylab, ...)
    for (i in 1:n) {
      lines(rr[[i]], col=i)
    }
  }

  matdensityplot(r$t)
dev.off()

png(file="g1447.png", width=600, height=600)
  q <- runif(1,0,10)
  m1 <- runif(1,0,10)
  m2 <- q*m1
  r <- my.simulation(
    function () { 
      n1 <- 200
      n2 <- 100
      x1 <- m1*rlnorm(n1)
      x2 <- m2*rlnorm(n2)
      data.frame( x=c(x1,x2), c=factor(c(rep(1,n1),rep(2,n2))))
    },
    function (d) { 
      a <- tapply(d[,1],d[,2],mean) 
      a[2]/a[1]
    }, 
    R=999
  )
  hist(r$t, probability=T, col='light blue',
       main="Distribution of the ratio of two means")
  lines(density(r$t), col='red', lwd=3)
  abline(v=c(q,r$t0), lty=3, lwd=3, col=c('blue','red'))
  legend( par("usr")[2], par("usr")[4], xjust=1, yjust=1,
          c("Theoretical mean", "Experimental mean"),
          lwd=1, lty=3, col=c('blue','red') )
dev.off()

png(file="g1448.png", width=600, height=600)
  # 5% confidence interval of the preceding example
  hist(r$t, probability=T, col='light blue',
       main="Distribution of the ratio of two means")
  qt <- quantile(r$t, c(.025,.975))
  d <- density(r$t)
  o <- d$x>qt[1] & d$x<qt[2]
  lines(d$x[o], d$y[o], col='red', lwd=3)
  lines(d, col='red')
  abline(v=c(q,r$t0), lty=3, lwd=3, col=c('blue','red'))
  legend( par("usr")[2], par("usr")[4], xjust=1, yjust=1,
          c("Theoretical mean", "Experimental mean",
            "5% confidence interval"),
          lwd=c(1,1,3), lty=c(3,3,1), col=c('blue','red', 'red') )
dev.off()

png(file="g1449.png", width=600, height=600)
  # We have two samples, of different sizes, and we study the
  # quotient of their means
  n1 <- 200
  n2 <- 100
  q <- runif(1,0,10)
  m1 <- runif(1,0,10)
  m2 <- q*m1
  x1 <- m1*rlnorm(n1)
  x2 <- m2*rlnorm(n2)
  d <- data.frame( x=c(x1,x2), c=factor(c(rep(1,n1),rep(2,n2))))
  R <- 999
  library(boot)
  r <- boot(d, 
            function (d,i) { 
              a <- tapply(d[,1][i],d[,2][i],mean) 
              a[2]/a[1]
            }, 
            R=R)
  hist(r$t, probability=T, col='light blue',
       main="Distribution of the ratio of two means")
  qt <- quantile(r$t, c(.025,.975))
  d <- density(r$t)
  o <- d$x>qt[1] & d$x<qt[2]
  lines(d$x[o], d$y[o], col='red', lwd=3)
  lines(d, col='red')
  abline(v=c(q,r$t0), lty=3, lwd=3, col=c('blue','red'))
  legend( par("usr")[2], par("usr")[4], xjust=1, yjust=1,
          c("Theoretical mean", "Experimental mean",
            "5% confidence interval"),
          lwd=c(1,1,3), lty=c(3,3,1), col=c('blue','red', 'red') )
dev.off()

png(file="g1450.png", width=600, height=600)
  # The variables to predict are the first of the data.frame
  my.simulation.predict <- function (
    get.training.sample,
    get.test.sample=get.training.sample,
    get.model,
    get.predictions=predict,
    R=999) {
    r <- matrix(NA, nr=R, nc=3)
    colnames(r) <- c("biais", "variance", "quadratic error")
    for (i in 1:R) {
      d.train <- get.training.sample()
      d.test <-  get.test.sample()
      m <- get.model(d.train)
      p <- get.predictions(m, d.test)
      if(is.vector(p)){ p <- data.frame(p) }
      d.test <- d.test[,1:(dim(p)[2])]
      b <- apply(d.test-p, 2, mean)
      v <- apply(d.test-p, 2, var)
      e <- apply((d.test-p)^2, 2, mean)
      r[i,] <- c(b,v,e)
    }
    list(t=r, t0=apply(r,2,mean))
  }

  get.sample <- function () {
    n <- 200
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    y <- sin(x1) + cos(x2) - x1*x2 + rnorm(n)
    data.frame(y,x1,x2)
  }
  r <- my.simulation.predict(
    get.sample,
    get.model = function (d) {
      lm(y~x1+x2,data=d)
    }
  )
  hist(r$t[,1], probability=T, col='light blue', 
       main="Bias of a linear regression")
  lines(density(r$t), col='red', lwd=3)
  abline(v=r$t0[1], lty=3)
dev.off()

png(file="g1451.png", width=600, height=600)
  library(MASS)
  get.sample <- function () {
    n <- 20
    x <- rnorm(n)
    x1 <- x + .2*rnorm(n)
    x2 <- x + .2*rnorm(n)
    x3 <- x + .2*rnorm(n)
    x4 <- x + .2*rnorm(n)
    y <- x + (x1+x2+x3+x4)/4 + rnorm(n)
    data.frame(y,x1,x2,x3,x4)
  }
  get.model <- function (d) {
    r <- lm.ridge(y~., data=d)
  }
  get.predictions <- function (r,d) {
    m <- t(as.matrix(r$coef/r$scales))
    inter <- r$ym - m %*% r$xm
    drop(inter) + as.matrix(d[,-1]) %*% drop(m)
  }

  lambda <- c(0,.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2,5,10,20,50,100)
  n <- length(lambda)
  res <- rep(NA,n)
  for (i in 1:n) {  
    res[i] <- my.simulation.predict(
      get.sample,
      get.model = function (d) {
        lm.ridge(y~., data=d, lambda=lambda[i])
      },
      get.predictions = get.predictions,
      R=99
    )$t0[3]
  }
  plot(res~lambda, type='l', log='x')
  abline(h=res[1], lty=3)
  i <- (1:n)[ res == min(res) ]
  x <- lambda[i]
  y <- res[i]
  points(x,y,col='red',pch=15)
  title(main="Bootstrap simulations to choose the parameter in a ridge regression")
dev.off()

png(file="g1452.png", width=600, height=600)
  library(MASS)
  get.sample <- function () {
    n <- 20
    x <- rnorm(n)
    x1 <- x + .2*rnorm(n)
    x2 <- x + .2*rnorm(n)
    x3 <- x + .2*rnorm(n)
    x4 <- x + .2*rnorm(n)
    y <- 1+x1+2*x2+x3+2*x4 + rnorm(n)
    data.frame(y,x1,x2,x3,x4)
  }
  get.model <- function (d) {
    r <- lm.ridge(y~., data=d)
  }
  get.predictions <- function (r,d) {
    m <- t(as.matrix(r$coef/r$scales))
    inter <- r$ym - m %*% r$xm
    drop(inter) + as.matrix(d[,-1]) %*% drop(m)
  }

  lambda <- c(0,.001,.002,.005,.01,.02,.05,.1,.2,.5,1,2,5,10,20,50,100)
  lambda <- exp(seq(-7,7,length=100))
  n <- length(lambda)
  res <- rep(NA,n)
  for (i in 1:n) {  
    res[i] <- my.simulation.predict(
      get.sample,
      get.model = function (d) {
        lm.ridge(y~., data=d, lambda=lambda[i])
      },
      get.predictions = get.predictions,
      R=20
    )$t0[1]
  }
  plot(res~lambda, type='l', log='x')
  abline(h=0,lty=3)
  title("Bias in ridge regression forecasts")
dev.off()

png(file="g1453.png", width=600, height=600)
  get.sample <- function () {
    n <- 20
    x <- rnorm(n)
    x1 <- x + .2*rnorm(n)
    x2 <- x + .2*rnorm(n)
    x3 <- x + .2*rnorm(n)
    x4 <- x + .2*rnorm(n)
    y <- 1+2*x1+3*x2+4*x3+5*x4 + rnorm(n)
    data.frame(y,x1,x2,x3,x4)
  }
  s <- function(d,l) {
    r <- lm.ridge(y~., data=d, lambda=l)
    m <- t(as.matrix(r$coef/r$scales))
    inter <- r$ym - m %*% r$xm
    c(inter,m)
  }
  lambda <- exp(seq(-7,7,length=100))
  n <- length(lambda)
  res <- matrix(NA,nr=n,nc=5)
  for (i in 1:n) {
    res[i,] <- my.simulation(
      get.sample,
      function (d) { s(d,lambda[i]) },
      R=20
    )$t0
  }
  matplot(lambda,res,type='l',lty=1,log='x')
  abline(h=1:5,lty=3,col=1:5)
  title("Bias in ridge regression coefficients")
dev.off()

png(file="g1454.png", width=600, height=600)
  get.sample <- function (a,b) {
    n <- 20
    x <- rnorm(n)
    x1 <- x + .2* rnorm(n)
    x2 <- x + .2* rnorm(n)
    y <- a*x1 + b*x2 + rnorm(n)
    data.frame(y,x1,x2)
  }
  get.parameters <- function (d) {
    y <- d[,1]
    x <- d[,-1]
    p <- princomp(x)
    r <- lm(y ~ p$scores[,1] -1)
    drop(p$loadings %*% c(r$coef,0))
  }
  N <- 5
  a.min <- 0
  a.max <- 10
  res <- matrix(NA,nr=N,nc=N)
  a <- seq(a.min,a.max,length=N)
  b <- seq(a.min,a.max,length=N)
  for (i in 1:N) {
    for(j in 1:N) {
      res[i,j] <- my.simulation(
        function () { get.sample(a[i],b[j]) },
        get.parameters,
        R=20
      )$t0[1] - a[i]
    }
  }

  persp(res)
dev.off()

png(file="g1455.png", width=600, height=600)
  image(a,b,res,col=topo.colors(255))
  text( matrix(a,nr=N,nc=N),
        matrix(b,nr=N,nc=N, byrow=T),
        round(res) )
dev.off()

png(file="g1456.png", width=600, height=600)
  get.parameters <- function (d) {
    y <- d[,1]
    x <- d[,-1]
    p <- princomp(x)
    r <- lm(y ~ p$scores[,1] -1)
    res <- drop(p$loadings %*% c(r$coef,0))
    a <- res[1]
    b <- res[2]
    c( a - (b-a)/2, b - (a-b)/2 )
  }
  N <- 5
  a.min <- 0
  a.max <- 10
  res <- matrix(NA,nr=N,nc=N)
  a <- seq(a.min,a.max,length=N)
  b <- seq(a.min,a.max,length=N)
  for (i in 1:N) {
    for(j in 1:N) {
      res[i,j] <- my.simulation(
        function () { get.sample(a[i],b[j]) },
        get.parameters,
        R=20
      )$t0[1] - a[i]
    }
  }
  persp(res)
dev.off()

png(file="g1457.png", width=600, height=600)
  image(a,b,res,col=topo.colors(255))
  text( matrix(a,nr=N,nc=N),
        matrix(b,nr=N,nc=N, byrow=T),
        round(res) )
dev.off()

png(file="g1458.png", width=600, height=600)
  get.parameters <- function (d) {
    y <- d[,1]
    x <- d[,-1]
    p <- princomp(x)
    r <- lm(y ~ p$scores[,1] -1)
    drop(p$loadings %*% c(r$coef,0))
  }
  N <- 5
  a.min <- 0
  a.max <- 10
  res <- matrix(NA,nr=N,nc=N)
  res.a <- matrix(NA,nr=N,nc=N)
  res.b <- matrix(NA,nr=N,nc=N)
  a <- seq(a.min,a.max,length=N)
  b <- seq(a.min,a.max,length=N)
  for (i in 1:N) {
    for(j in 1:N) {
      r <- my.simulation(
        function () { get.sample(a[i],b[j]) },
        get.parameters,
        R=20
      )$t0
      res.a[i,j] <- r[1]
      res.b[i,j] <- r[2]
      res[i,j] <- r[1] - a[i]
    }
  }
  plot(res.a, res.b, type='n')
  text(res.a, res.b, round(res))
dev.off()

png(file="g1459.png", width=600, height=600)
  plot(as.vector(res) ~ as.vector(res.a-res.b))
dev.off()

png(file="g1460.png", width=600, height=600)
  my.simulation.cross.validation <- function (
    s, # sample
    k, # number of observations to fit the the model
    get.model,
    get.predictions=predict,
    R = 999 )
  {
    n <- dim(s)[1]
    r <- matrix(NA, nr=R, nc=3)
    colnames(r) <- c("bias", "variance", "MSE")
    for (i in 1:R) {
      j <- sample(1:n, k)
      d.train <- s[j,]
      d.test <-  s[-j,]
      m <- get.model(d.train)
      p <- get.predictions(m, d.test)
      if(is.vector(p)){ p <- data.frame(p) }
      d.test <- d.test[,1:(dim(p)[2])]
      b <- apply(d.test-p, 2, mean)
      v <- apply(d.test-p, 2, var)
      e <- apply((d.test-p)^2, 2, mean)
      r[i,] <- c(b,v,e)
    }
    list(t=r, t0=apply(r,2,mean))
  }

  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 2 - x1 - x2 - x3 + rnorm(n)
  s <- data.frame(y,x1,x2,x3)
  r <- my.simulation.cross.validation(s, 150, function(s){lm(y~.,data=s)}, R=99)
  hist( r$t[,3], probability=T, col='light blue', 
        main="MSE estimation by cross-validation" )
  lines(density(r$t[,3]),col='red',lwd=3)
  abline(v=mean(r$t[,3]),lty=2,lwd=3,col='red')
dev.off()

png(file="g1461.png", width=600, height=600)
  my.simulation.jackknife <- function (
    s,
    get.model,
    get.predictions=predict,
    R=999)  
  {
    n <- dim(s)[1]
    if(R<n) {
      j <- sample(1:n, R)
    } else {
      R <- n
      j <- 1:n
    }
    p <- get.predictions(get.model(s), s)
    if( is.vector(p) ) {
      p <- as.data.frame(p)
    }
    r <- matrix(NA, nr=R, nc=dim(s)[2]+dim(p)[2])
    colnames(r) <- c(colnames(s), colnames(p))
    for (i in j) {
      d.train <- s[-i,]
      d.test <-  s[i,]
      m <- get.model(d.train)
      p <- get.predictions(m, d.test)
      r[i,] <- as.matrix(cbind(s[i,], p))
    }
    r
  }

  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 2 - x1 - x2 - x3 + rnorm(n)
  s <- data.frame(y,x1,x2,x3)
  r <- my.simulation.jackknife(s, function(s){lm(y~.,data=s)})
  hist( r[,5]-r[,1], col='light blue' )
dev.off()

png(file="g1462.png", width=600, height=600)
  my.simulation.bootstrap <- function (
    s,
    get.model,
    get.predictions=predict,
    R = 999 )
  {
    n <- dim(s)[1]
    r <- matrix(NA, nr=R, nc=3)
    colnames(r) <- c("bias", "variance", "MSE")
    for (i in 1:R) {
      j <- sample(1:n, n, replace=T)
      d.train <- s[j,]
      d.test <-  s
      m <- get.model(d.train)
      p <- get.predictions(m, d.test)
      if(is.vector(p)){ p <- data.frame(p) }
      d.test <- d.test[,1:(dim(p)[2])]
      b <- apply(d.test-p, 2, mean)
      v <- apply(d.test-p, 2, var)
      e <- apply((d.test-p)^2, 2, mean)
      r[i,] <- c(b,v,e)
    }
    list(t=r, t0=apply(r,2,mean))
  }

  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 2 - x1 - x2 - x3 + rnorm(n)
  s <- data.frame(y,x1,x2,x3)
  r <- my.simulation.bootstrap(s, function(s){lm(y~.,data=s)}, R=99)
  hist( r$t[,3], col='light blue', 
        main="Using the bootstrap to estimate the MSE" )
  lines(density(r$t[,3]),col='red',lwd=3)
  abline(v=mean(r$t[,3]),lty=2,lwd=3,col='red')
dev.off()

png(file="g1463.png", width=600, height=600)
  library(boot)
  n <- 200
  x <- rlnorm(n)
  r <- boot(x, function(s,i){mean(s[i])}, R=99)
  hist(r$t, probability=T, col='light blue')
  lines(density(r$t),col='red',lwd=3)
  abline(v=mean(r$t),lty=2,lwd=3,col='red')
  title("Bootstrap to estimate the distribution of the mean of a sample")
dev.off()

png(file="g1464.png", width=600, height=600)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- -2 - x1 +0* x2 + x3 + rnorm(n)
  s <- data.frame(y,x1,x2,x3)
  r <- boot(s, function(s,i){lm(y~.,data=s[i,])$coef}, R=99)
  plot(NA,xlim=range(r$t), ylim=c(0,6.5))
  for (i in 1:4) {
    lines(density(r$t[,i]), lwd=2, col=i)
  }
dev.off()

png(file="g1465.png", width=600, height=600)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 1+x1+x2+x3+rnorm(n)
  d <- data.frame(y,x1,x2,x3)
  r <- boot(d, 
            function (d,i) {
              r <- lm(y~.,data=d[i,])
              p <- predict(r, d[-i,])
              mean((p-d[,1][-i])^2)
            },
            R=99)
  hist(r$t,probability=T,col='light blue',
       main="Mean Square Error (MSE)")
  lines(density(r$t),lwd=3,col='red')
  abline(v=mean(r$t),lty=3)
dev.off()

png(file="g1466.png", width=600, height=600)
  get.sample <- function (n=50) { rnorm(n) }
  get.statistic <- function (x) { mean(x) }
  # Plot the density of the estimator
  R <- 200
  x <- rep(NA,R)
  for (i in 1:R) {
    x[i] <- get.statistic(get.sample())
  }
  plot(density(x), ylim=c(0,4), lty=1,lwd=3,col='red',
       main="Relevance of the bootstrap")
  # Let us plot a few estimations of this density, 
  # obtained from the bootstrap
  for (i in 1:5) {
    x <- get.sample()
    r <- boot(x, function(s,i){mean(s[i])}, R=R)
    lines(density(r$t))
  }
  curve(dnorm(x,sd=1/sqrt(50)),add=T, lty=2,lwd=3,col='blue')
  legend(par('usr')[2], par('usr')[4], xjust=1, yjust=1,
         c('parametric bootstrap', 'non-parametric bootstrap',
           "theoretical curve"),
         lwd=c(3,1,2),
         lty=c(1,1,2),
         col=c("red", par("fg"), "blue"))
dev.off()

png(file="g1467.png", width=600, height=600)
  n <- 50
  x <- rnorm(n)
  r <- boot(x, function(s,i){max(s[i])}, R=999)
  hist(r$t, probability=T, col='light blue',
       main="A statistic for which the bootstrap is not very relevant...")
  lines(density(r$t,bw=.05), col='red', lwd=3)
dev.off()

png(file="g1468.png", width=600, height=600)
  n <- 20
  N <- 100
  r <- matrix(NA,nr=N,nc=2)
  for (i in 1:N) {
    x <- rnorm(n)
    r[i,] <- boot.ci(boot(x,function(x,i){mean(x[i])},R=99))$basic[c(4,5)]
  }
  plot(NA,xlim=c(-1.5,1.5), ylim=c(0,2),
       main="Confidence interval boundaries (bootstrap)")
  lines(density(r[,1]), lwd=3, col='red')
  lines(density(r[,2]), lwd=3, col='blue')
  abline(v=0,lty=3)
dev.off()

png(file="g1469.png", width=600, height=600)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 1+x1+x2+x3+rnorm(n)
  d <- data.frame(y,x1,x2,x3)
  r1 <- boot(d, function (d,i) {
               r <- lm(y~.,data=d[i,])
               p <- predict(r, d)
               mean(p-d[,1])
             },
             R=99)$t
  r2 <- boot(d, function (d,i) {
               r <- lm(y~.,data=d[i,])
               p <- predict(r, d[-i,])
               mean(p-d[,1][-i])
             },
             R=99)$t
  r3 <- .368*r1 + .632*r2
  hist(r1,probability=T,col='light blue',
       main="Bias estimation: bootstrap, oobb, .632")
  lines(density(r1),lwd=3,col='red')
  lines(density(r2),lwd=3,col='green')
  lines(density(r3),lwd=3,col='blue')
  legend(par('usr')[2], par('usr')[4], xjust=1, yjust=1,
         c("raw bootstrap", "out-of-the-bag bootstrap", ".632"),
         lwd=3,
         col=c("red", "green", "blue")
        )
dev.off()

png(file="g1470.png", width=600, height=600)
  library(boot)
  my.max <- function (x, i) { max(x[i]) }
  r <- boot(faithful$eruptions, my.max, R=100)
  plot(r)
dev.off()

png(file="g1471.png", width=600, height=600)
  jack.after.boot(r)
dev.off()

png(file="g1472.png", width=600, height=600)
  get.sample <- function (n=50) {
    x <- rnorm(n)
    y <- 1-x+rnorm(n)
    data.frame(x,y)
  }
  s <- get.sample()

  R <- 100
  n <- dim(s)[1]
  x <- s$x
  y <- s$y
  res <- rep(NA,R)
  res.oob <- rep(NA,R)
  for (k in 1:R) {
    i <- sample(n,n,replace=T)
    xx <- x[i]    
    yy <- y[i]
    r <- lm(yy~xx)
    ip <- rep(T,n)
    ip[i] <- F
    xp <- x[ip]
    yp <- predict(r,data.frame(xx=xp))
    try( res.oob[k] <- sd(yp - s$y[ip]) )
    res[k] <- sd( yy - predict(r) )
  }
  plot(density(res), main="Basic and out-of-the-box bootstrap")
  lines(density(res.oob),lty=2)
  abline(v=1,lty=3)
  legend( par("usr")[2], par("usr")[4], xjust=1, yjust=1,
          c("Basic bootstrap", "Out-of-the-box bootstrap"),
          lwd=1, lty=c(1,2) )
dev.off()

png(file="g1473.png", width=600, height=600)
  echantillon.multicolineaire <- function (n=100) {
    x0 <- rep(1,n)
    x1 <- x0 + .1*rnorm(n)
    x2 <- x0 + .1*rnorm(n)
    x3 <- x0 + .1*rnorm(n)
    y <- x0 + x1 + x2 + x3 + rnorm(n)
    x <- cbind(x0,x1,x2,x3)
    m <- cbind(y,x)
    m
  }
  m <- echantillon.multicolineaire()
  f <- function (m, k) { 
    x <- m[,-1]
    y <- m[,1]
    n <- dim(x)[2]
    solve( t(x) %*% x + k*diag(rep(1,n)), t(x) %*% y )
  }
  N <- 200
  variance <- rep(NA, N)
  biais <- rep(NA, N)
  exact <- f(m,0)
  k <- exp(seq(.01,4, length=N))-1
  for (j in 1:N) {
    g <- function (m, i) { f(m[i,], k[j]) }
    s <- boot(m, g, R=1000)
    variance[j] <- max(var(s$t))
    biais[j] <- max(exact - apply(s$t,2,mean))
  }
  plot(variance ~ k, log='x', type='l')
  par(new=T)
  plot(biais ~ k, type='l', log='x', col='red', xlab='', ylab='', axes=F)
  axis(4, col="red")
  legend( par("usr")[2], par("usr")[4], yjust=1, xjust=0,
          c("variance", "optimism"),
          lwd=1, lty=1,
          col=c(par('fg'), 'red') )
dev.off()

png(file="g1474.png", width=600, height=600)
  my.bootstrap <- function (data, estimator, size, ...) {
    values <- NULL
    n <- dim(data)[1]
    for (i in 1:size) {
      bootstrap.sample <- data[ sample(1:n, n, replace=T), ]
      values <- rbind(values, t(estimator(bootstrap.sample, ...)))
    }
    values
  }

  n <- 100
  x0 <- rep(1,n)
  x1 <- x0 + .1*rnorm(n)
  x2 <- x0 + .1*rnorm(n)
  x3 <- x0 + .1*rnorm(n)
  y <- x0 + x1 + x2 + x3 + rnorm(n)
  x <- cbind(x0,x1,x2,x3)
  f <- function (m, k) { 
    x <- m[,-1]
    y <- m[,1]
    n <- dim(x)[2]
    solve( t(x) %*% x + k*diag(rep(1,n)), t(x) %*% y )
  }
  N <- 200
  variance <- rep(NA, N)
  biais <- rep(NA, N)
  k <- exp(seq(0,6, length=N))-1
  for (i in 1:N) {
    s <- my.bootstrap(cbind(y,x), f, 100, k[i])
    variance[i] <- max(apply(s, 2, var))
    biais[i] <- max(abs( s - 
      matrix( as.vector(f(cbind(y,x),0)), nr=dim(s)[1], nc=dim(s)[2], byrow=T )
    ))
  }
  plot(variance ~ k, log='x', type='l')
  par(new=T)
  plot(biais ~ k, type='l', log='x', col='red', xlab='', ylab='', axes=F)
  axis(4)
dev.off()

png(file="g1475.png", width=600, height=600)
  s <- my.bootstrap(cbind(y,x), f, 100, 10)
  b <- apply(s, 2, mean)
  b0 <- as.vector(f(cbind(y,x),0))

  # Another sample from the same population
  # TODO: put this in a function
  n <- 100
  x0 <- rep(1,n)
  x1 <- x0 + .1*rnorm(n)
  x2 <- x0 + .1*rnorm(n)
  x3 <- x0 + .1*rnorm(n)
  y <- x0 + x1 + x2 + x3 + rnorm(n)
  x <- cbind(x0,x1,x2,x3)

  a1 <- x %*% b  
  a2 <- x %*% b0
  boxplot( c(as.vector(a1), as.vector(a2)) ~ 
           gl(2,n,2*n, c("ridge", "classic")),
           horizontal=T )
  title(main="Residues of a ridge regression")
dev.off()

png(file="g1476.png", width=600, height=600)
  n <- 50
  k <- 10 
  N <- 50
  get.sample <- function (n) {
    x <- runif(n)
    y <- sin(2*pi*x)+.3*rnorm(n)
    m <- data.frame(x=x, y=y)
    m
  }
  error <- rep(NA,k)
  for (i in 1:k) {
    e <- rep(NA,N)
    for (j in 1:N) {
      training.sample <- get.sample(n)
      x <- training.sample$x
      y <- training.sample$y
      r <- lm(y~poly(x,i))
      test.sample <- get.sample(1000)
      e[j] <- mean( (predict(r, data.frame(x=test.sample$x)) - test.sample$y)^2 )
    }
    error[i] <- mean(e)
  }
  plot(error, type='l', log='y',
       xlab="Model complexity",
       ylab="Errorr")
dev.off()

png(file="g1477.png", width=600, height=600)
  a <- 2
  curve(dnorm(x-a), xlim=c(-3,5), ylim=c(0,.5))
  abline(v=0, lty=3, col='red')
  abline(v=a, lty=3, col='blue')
  arrows(0,.45,a,.45, code=3, col='green')
  arrows(a-pnorm(1), dnorm(pnorm(1)), a+pnorm(1), dnorm(pnorm(1)), code=3, col='red')
  legend(par('usr')[2], par('usr')[4], xjust=1,
         c("true value", "estimator expectation", "bias", "variance"),
         lwd=1, lty=c(3,3,1,1),
         col=c('red', 'blue', 'green', 'red') )
  title(main="Bias and variance of a biased estimator")
dev.off()

png(file="g1478.png", width=600, height=600)
  my.gibbs <- function (x0=c(0,0), N=20, r=.5, plot.it=T) {
    res <- matrix(NA, nr=N, nc=2)
    res[1,] <- c(0,0)
    for (i in seq(2,N-1,by=2)) {
      x1 <- res[i-1,1]
      x2 <- res[i-1,2]
      res[i,2] <- x2
      # WE are supposed to take the new point at random 
      # and accept or reject it depending of the conditionnal
      # probability.
      # Here, I directly sample from the conditionnal probability.
      x1 <- rnorm(1, r*x2, sd=sqrt(1-r^2))
      res[i,1] <- x1
      res[i+1,1] <- x1
      x2 <- rnorm(1, r*x1, sd=sqrt(1-r^2))
      res[i+1,2] <- x2
    }
    if (plot.it) {
      plot(res, type='l', xlab="x1", ylab="x2", main="Gibbs sampler")
      points(res, pch=16)
      invisible(res)
    } else {
      res
    }
  }
  my.gibbs()
dev.off()

png(file="g1479.png", width=600, height=600)
  my.gibbs(N=1000)
dev.off()

png(file="g1480.png", width=600, height=600)
  r <- my.gibbs(N=40,plot.it=F)[2*(1:20),]
  plot(r, type='l', xlab="x1", ylab="x2",
       main="Path of a Metropolis simulation (almost)")
  points(r, pch=16)
dev.off()

png(file="g1481.png", width=600, height=600)
  r <- my.gibbs(N=2000,plot.it=F)[2*(1:1000),]
  plot(r, type='l', xlab="x1", ylab="x2",
       main="Path of a Metropolis simulation (almost)")
  points(r, pch=16)
dev.off()

png(file="g1482.png", width=600, height=600)
  N <- 200
  K <- 5
  res <- matrix(NA, nr=N, nc=K)
  for (i in 1:K) {
    r <- my.gibbs(x0=runif(2,-1,1), N=N, plot.it=F)
    res[,i] <- apply(r,1,cumsum)[1,] / 1:N
  }
  matplot(res, type='l', lty=1,
          ylab="Estimator", xlab="Simulation length",
          main="Gibbs sampler")
dev.off()

png(file="g1483.png", width=600, height=600)
  library(ts)
  MH <- function (N = 1000, 
                  voisin = function (x) { rnorm(1,x,1) },
                  p = dnorm,                 # Probability distribution to simulate
                  q = function (x,y) { 1 }   # Hastings Correction
                 ) {
    res <- rep(NA,N)
    x <- 0
    for (i in 1:N) {
      y <- voisin(x)
      u <- runif(1)
      if ( u < q(x,y)/q(y,x) * p(y)/p(x) ) {
        x <- y
      }
      res[i] <- x
    }
    ts(res)
  }

  x1 <- ts(rnorm(1000))
  x2 <- MH()
  x3 <- MH(voisin = function (x) { rnorm(1,x,.01) })
  x4 <- MH(voisin = function (x) { rnorm(1,x,50) })

  op <- par(mfrow=c(4,1))
  plot(x1, main="Random")
  plot(x2, main="MCMC")
  plot(x3, main="MCMC, neighbours too near")
  plot(x4, main="MCMC, neighbours too far")
  par(op)
dev.off()

png(file="g1484.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  qqnorm(x1, main="Random"); qqline(x1); 
  qqnorm(x2, main="MCMC"); qqline(x2); 
  qqnorm(x3, main="MCMC, neighbours too near"); qqline(x3); 
  qqnorm(x4, main="MCMC, neighbours too far away"); qqline(x4)
  par(op)
dev.off()

png(file="g1485.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  hist(x1, col='light blue', probability=T, main="Random")
  curve(dnorm(x),col='red',lwd=3, add=T)
  hist(x2, col='light blue', probability=T, main="MCMC")
  curve(dnorm(x),col='red',lwd=3, add=T)
  hist(x3, col='light blue', probability=T, main="MCMC, neighbours too near")
  curve(dnorm(x),col='red',lwd=3, add=T)
  hist(x4, col='light blue', probability=T, main="MCMC, neighnours too far away")
  curve(dnorm(x),col='red',lwd=3, add=T)
  par(op)
dev.off()

png(file="g1486.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  acf(x1, main="Random")
  acf(x2, main="MCMC")
  acf(x3, main="MCMC, neighbours too near")
  acf(x4, main="MCMC, neighnours too far away")
  par(op)
dev.off()

png(file="g1487.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  pacf(x1, main="Random")
  pacf(x2, main="MCMC")
  pacf(x3, main="MCMC, neighbours too near")
  pacf(x4, main="MCMC, neighnours too far away")
  par(op)
dev.off()

png(file="g1488.png", width=600, height=600)
  N <- 3000
  x <- MH(N=N)
  plot(x[2:N]~x[2:N-1],xlab='x(i)',ylab='x(i+1)')
  abline(0,1,col='red')
dev.off()

png(file="g1489.png", width=600, height=600)
  N <- 100000
  x <- MH(N=N)
  x <- x[100:length(x)]  # burn-in
  x <- x[seq(1,length(x),by=50)]
  N <- length(x)
  plot(x[2:N]~x[2:N-1],xlab='x(i)',ylab='x(i+1)')
  abline(0,1,col='red')
dev.off()

png(file="g1490.png", width=600, height=600)
  op <- par(mfrow=c(4,1))
  plot(ts(x))
  acf(x)
  pacf(x)
  spectrum(x)
  par(op)
dev.off()

png(file="g1491.png", width=600, height=600)
  system("jags jags_regression_test.jags")
  library(coda)
  x <- read.jags()
  plot(x)
dev.off()

png(file="g1492.png", width=600, height=600)
  autocorr.plot(x)
dev.off()

png(file="g1493.png", width=600, height=600)
  N <- 3
  op <- par(mfrow=c(N,1))
  for (i in 1:N) {
    pacf(x[,i], main=colnames(x)[i])
  }
  par(op)
dev.off()

png(file="g1494.png", width=600, height=800)
  f <- function (x, by=2) {
    i <- seq(1, length(x), by=by)
    x[i]
  }
  op <- par(mfrow=c(4,1), mar=c(2,4,4,2))
  for (i in 2:5) {
    pacf(f(x[,3],i), main=paste("sigma, thinning=", i, sep=""), xlab="")
  }
  par(op)
dev.off()

png(file="g1495.png", width=600, height=600)
  # One of the graphs is too large to be displayed...
  library(Rgraphviz)
  data(graphExamples)
  op <- par(mfrow=c(6,3), mar=c(0,0,0,0))
  for (g in graphExamples) {
    cat(object.size(g), "\n")
    if (object.size(g) < 100000) {
      cat("Drawing\n")
      try( plot(g) )
    } else {
      cat("Skipping\n")
    }
  }
  par(op)
dev.off()

png(file="g1496.png", width=600, height=600)
  library(Rgraphviz)
  data(graphExamples)
  op <- par(mfrow=c(3,3))
  set.seed(2)
  for (i in 1:9) {
    try(   plot(randomGraph(LETTERS[1:10], 1, .3))   )
  }
  par(op)
dev.off()

png(file="g1497.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    try(   plot(randomGraph(LETTERS[1:10], 1, .5))   )
  }
  par(op)
dev.off()

png(file="g1498.png", width=600, height=600)
  op <- par(mfrow=c(3,3))
  for (i in 1:9) {
    try(   plot(randomGraph(LETTERS[1:10], 1, .8))   )
  }
  par(op)
dev.off()

png(file="g1499.png", width=600, height=600)
  op <- par(mfrow=c(2,2))
  for (a in 1:4) {

    n <- 10
    V <- LETTERS[1:n]
    m <- matrix(sample(0:1,n*n,replace=T,prob=c(.8,.2)), nr=n, nc=n, dimnames=list(V,V))
    m <- m + t(m)
    m <- 1 - (m-1)*(m-2)/2
    m <- m - diag(diag(m))

    e <- vector("list", length=n)
    names(e) <- V
    for (i in 1:n) {
      e[[i]] <- list( edges = which(m[i,]==1) )  # Edge numbers, not names
    } 
    g <- new("graphNEL", nodes=V, edgeL=e)
    plot(g)

  }
  par(op)
dev.off()
detach.everything()
detach.everything()

png(file="g1500.png", width=600, height=600)
  x <- seq(0,1, length=100)
  x <- x[2:(length(x)-1)]
  logit <- function (t) {
    log( t / (1-t) )
  }
  plot(logit(x) ~ x, type='l')
dev.off()

png(file="g1501.png", width=600, height=600)
  curve(logit(x), col='blue', add=F)
  curve(qnorm(x), col='red', add=T)
  a <- par("usr")
  legend( a[1], a[4], c("logit","probit"), col=c("blue","red"), lty=1)
dev.off()

png(file="g1502.png", width=600, height=600)
  N <- 10000
  f <- NULL
  c <- NULL
  for (i in 1:N) {
    x <- sample( c(rep(0,3),rep(1,7)), 10, replace=T )
    y <- sample( c(0,1), 10, replace=T )
    x <- factor(x, levels=c(0,1))
    y <- factor(y, levels=c(0,1))
    t <- table(x,y)
    f <- append(f, fisher.test(t)$p.value)
    c <- append(c, chisq.test(t)$p.value)
  }
  plot(sort(f), type='l')
  points(sort(c), type='l', col='red')
dev.off()

png(file="g1503.png", width=600, height=600)
  x <- sample( c(rep(0,3),rep(1,7)), 10, replace=T )
  y <- sample( c(0,1), 10, replace=T )
  x <- factor(x, levels=c(0,1))
  y <- factor(y, levels=c(0,1))
  pValeur <- function (x,y, N=1000) {
    x <- factor(x)
    y <- factor(y)
    s <- numeric(N)
    for (i in 1:N) {
      xx <- sample( x, length(x), replace=T )
      yy <- sample( y, length(y), replace=T )
      s[i] <- as.numeric(chisq.test(table(xx,yy))$statistic)
    }
    t <- table(x,y)
    p <- sum( chisq.test(t)$statistic > s )/length(s)
    if( p>.5 ) p <- 1-p
    p
  }
  t <- table(x,y)
  c( ChiSq = chisq.test(t)$p.value, Fisher = fisher.test(t)$p.value, MonteCarlo = pValeur(x,y) )

  #Représentons maintenant nos trois p-valeurs sur un même graphique.

  N <- 100
  f <- NULL
  c <- NULL
  s <- NULL
  for (i in 1:N) {
    x <- sample( c(rep(0,3),rep(1,7)), 10, replace=T )
    y <- sample( c(0,1), 10, replace=T )
    x <- factor(x, levels=c(0,1))
    y <- factor(y, levels=c(0,1))
    t <- table(x,y)
    f <- append(f, fisher.test(t)$p.value)
    c <- append(c, chisq.test(t)$p.value)
    s <- append(s, pValeur(x,y,100))
  }
  plot(sort(f), type='l')
  points(sort(c), type='l', col='green')
  points(sort(s), type='l', col='red')
  
  # J'ai de TRÈS gros doutes sur la pertinence de mes calculs...
  # A FAIRE
dev.off()

png(file="g1504.png", width=600, height=600)
  n <- 100
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- x1 - x2*x2 + 2*x3 + rnorm(n)
  pairs(cbind(y,x1,x2))
dev.off()

png(file="g1505.png", width=600, height=600)
  y3  <- lm(y~x3)$residuals
  x13 <- lm(x1~x3)$residuals
  x23 <- lm(x2~x3)$residuals
  plot( y3 ~ x13,
        xlab="x1 sans les effets linéaires de x3",
        ylab="y sans les effets linéaires de x3" )
dev.off()

png(file="g1506.png", width=600, height=600)
  library(lattice)
  bwplot( ~ y3 | equal.count(x13), layout=c(1,6) )    
dev.off()

png(file="g1507.png", width=1200, height=600)
  bwplot( ~ y3 | equal.count(x13) + equal.count(x23) )    
dev.off()

png(file="g1508.png", width=600, height=600)
  n <- 1000
  x <- runif(n)
  x <- x>.5
  x <- cumsum(x)/(1:n)
  plot(x, ylim=c(0,1), type="l")
dev.off()

png(file="g1509.png", width=600, height=600)
  plot(x, ylim=c(0,1), type="l", log="x")
dev.off()

png(file="g1510.png", width=600, height=600)
  N <- 10
  n <- 20
  x <- rnorm(n*N)
  y <- rep(1:n, N)
  boxplot(x ~ y)
dev.off()

png(file="g1511.png", width=600, height=600)
  par( mfrow = c(4,5) )
  for (i in 1:n) {
    hist(x[ y==i ])
  }
dev.off()

png(file="g1512.png", width=600, height=600)
  N <- 100
  n <- 20
  x <- rnorm(n*N)
  y <- rep(1:n, N)
  boxplot(x ~ y)
dev.off()

png(file="g1513.png", width=600, height=600)
  par( mfrow = c(4,5) )
  for (i in 1:n) {
    hist(x[ y==i ])
  }
dev.off()

png(file="g1514.png", width=600, height=600)
  N <- 3
  n <- 20
  x <- rnorm(n*N)
  y <- rep(1:n, N)
  m <- matrix(x, nrow=N, byrow=T)
  c <- apply(m,2,min)>0 | apply(m,2,max)<0
  boxplot(x ~ y, col=c(0,6)[1+c])
dev.off()

png(file="g1515.png", width=600, height=600)
  x <- rnorm(10+100+1000+10000+100000)
  y <- c( rep(5,10), rep(4,100), rep(3,1000), rep(2,10000), rep(1,100000))
  boxplot(x~y, horizontal=T, axes=F)
  axis(1)
  axis(2, 1:5, c(10,100,1000,10000,100000) )
dev.off()

png(file="g1516.png", width=600, height=600)
  x <- 1:9
  y <- log(1+1/x)/log(10)
  plot(y~x, type="h")
dev.off()

png(file="g1517.png", width=600, height=600)
  a <- read.table("Cours.txt")
  a <- a[a>0]
  a <- as.vector(a)
  a <- floor(a*10^-floor(log(a)/log(10)))
  hist(a, probability=T)
  points((x+.5), y, type='h', col='red', lwd=3)
  chisq.test(table(factor(a)), p=y) # p-value = 0.06
dev.off()

png(file="g1518.png", width=600, height=600)
  a <- read.table("Volume.txt")
  a <- a[a>0]
  a <- as.vector(a)
  a <- floor(a*10^-floor(log(a)/log(10)))
  hist(a, probability=T)
  points((x+.5), y, type='h', col='red', lwd=3)
  chisq.test(table(factor(a)), p=y) # p-value = 0.85
dev.off()

png(file="g1519.png", width=600, height=600)
  plot(sort(runif(100)))
  for (i in 1:1000) {
    lines( sort(runif(100)) )
  }
dev.off()

png(file="g1520.png", width=600, height=600)
  n <- 1000
  m <- 200
  a <- matrix( runif(n*m), c(n,m) )
  b <- apply(a, 1, sort)
  c <- apply(b, 1, range)
  plot( c[1,], type="l" )
  lines( c[2,] )

  n <- 10
  m <- 200
  a <- matrix( runif(n*m), c(n,m) )
  b <- apply(a, 1, sort)
  c <- apply(b, 1, range)
  lines( c[1,] )
  lines( c[2,] )
dev.off()
detach.everything()
detach.everything()
detach.everything()
detach.everything()
