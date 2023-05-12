###############################################################################
###############################################################################
# Chapter 3
###############################################################################
###############################################################################

########################################
# 3.1 Discrete distributions
########################################

########################################
# 3.1.1 Bernoulli distribution
########################################

########################################
# 3.1.2 Binomial distribution
########################################

# Example 1. Load the TeachingDemos package...
# install.packages("TeachingDemos")
# Note that TeachingDemos requires XQuartz from www.xquartz.org on the Mac.
library(TeachingDemos)
vis.binom()

# Example 2. If we cross two heterozygous carriers...
choose(3,1)* 0.25^1* 0.75^2

for (k in 0:3) {
    print(dbinom(k,3,0.25))
}

pbinom(2,3,0.25)

# Example 3. After four coin tosses, we can....
pbinom(2,4,0.50)  # binomial cumulative distribution function (CDF)


# Example 4. microRNA consists of a sequence composed of nucleotides A, G, U, and C,...
dbinom(14, 22, 0.7)
pbinom(13, 22, 0.7)

sum(dbinom(11:22, 22, 0.7))
1 - pbinom(10, 22, 0.7)

x <- 0:22
plot(x,           # X-values
     dbinom(x,size=22,prob=.7), # binomial mass function
     type="h",      # histogram
     col="blue",
     lwd=4,         # make lines thicker
     xlab="x",
     ylab="Binomial probability mass function f(x)",
     cex.lab=1.5)    # make axis labels big

binomialCDF = stepfun(x, c(0,pbinom(x,size=22,prob=.7))) # create a step function from binomial CDF
plot(binomialCDF, # binomial cumulative function
     col="red",
     vertical=FALSE,
     lwd=4,         # make lines thicker
     xlab="x",
     ylab="Binomial cumulative distribution function F(x)",
     main=NULL,
     cex.lab=1.5)    # make axis labels big)

rbinom(1000, 22, 0.7)

########################################
# 3.1.3 Poisson distribution
########################################

# Poisson approximating the binomial
lambda=4
plot(x,              # X-values
     dpois(x,lambda),# Poisson density function
     type="h",       # histogram
     col=1,
     lwd=8,          # make lines thicker
     ylim=c(0,0.25), # range of the y-axis
     xlab="x",
     ylab="Density function f(x)",
     cex.lab=1.5)    # make axis labels big
for (counter in 1:10) {
    lines(x, dbinom(x, size=counter*10, prob=lambda/(counter*10)), col=counter+1, lwd=2)
}
legend("topright", lty=rep(1,11), lwd=rep(2,11), col=1:11,
       legend=c(expression(Poisson(lambda == 4)),
           as.expression(sapply(seq(10,100, 10), function(x) bquote(B(n ==.(x), p == 4/.(x)))))))

# Example 1. Teaching demo
vis.binom()

# Example 2. Daily Lottery
exp(-2)*2^0/factorial(0)
dpois(0,2)

exp(-2)*2^1/factorial(1)
dpois(1,2)

x <- 0:20
plot(x,           # X-values
     dpois(x,5),  # Poisson mass function
     type="h",    # histogram
     col="blue",
     lwd=4,        # make lines thicker
     xlab="x",
     ylab="Poisson probability mass function f(x)",
     cex.lab=1.5)    # make axis labels big

poissonCDF = stepfun(x, c(0,ppois(x,5))) # create a step function from poisson CDF
plot(poissonCDF, # poisson cumulative function
     col="red",
     vertical=FALSE,
     lwd=4,         # make lines thicker
     xlab="x",
     ylab="Poisson cumulative distribution function F(x)",
     main=NULL,
     cex.lab=1.5)    # make axis labels big)

# Example 3. Finding k-mers
exp(-3.9)*2^0/factorial(0)
dpois(0,3.9)

1 - dpois(0,3.9)

# Example 4. Gaps in an open reading frame
exp(-5)*5^2/factorial(2)
dpois(2,5)

ppois(3,5)

y <- rpois(10000,5)   # 1000 random samples with lambda=5
mean(y)
var(y)

qpois(c(0.025,0.975), lambda=5, lower.tail = TRUE)

sum(dpois(1:10,5))

########################################
# 3.2 Continuous distributions
########################################

########################################
# 3.2.1 Exponential distribution
########################################

# Example 1. Suppose that lambda=5
pexp(.5,5)

f <-function(x) { dexp(x,5) }
plot(f,  # function to plot
     0,  # first x-value
     2,  # last x-value
     xlab="x",
     ylab="Expontential mass function f(x)",
     cex.lab=1.5)  # make axis labels bigger
x1 <- seq(0,0.5,.1) # set of magenta x-values
x2 <- seq(0.5,2,.1) # set of blue x-values
polygon(c(0,x1,.5), c(0,f(x1),0), col="magenta")
polygon(c(.5,x2,2), c(0,f(x2),0), col="lightblue")
arrows(0.75,2.5,0.25,1.6, lwd=3, col="green")
text(0.86, 2.7 - c(0,0.7), cex = 1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 0.5)),
       expression(paste("  = ", integral(f(x) * dx, 0, 0.5)))))

F <- function(x) { pexp(x,5) }
plot(F,             # function
     0,             # start x
     2,             # end x
     cex.lab=1.5,   # make axis labels big
     col="red",
     lwd=6,         # make line thicker
     xlab="x",
     ylab="Exponential cumulative distribution function F(x)")
mtext("0.92",side=2,at=0.92, col="red")
arrows(0.8,0.6,0.5,0.9, lwd=3, col="green")
text(0.65, 0.5 - c(0,.11,.22,.33,.44), cex=1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 0.5)),
       expression(paste("  = ", integral(f(x) * dx, 0, 0.5))),
       expression(paste("  = ", bgroup("[", F(x) ,"]")[0]^0.5)),
       expression(paste("  = ", CDF(0.5))),
       expression(paste("  = ", 0.92))))

1-pexp(2,5)

pexp(2,5)-pexp(.25,5)

qexp(c(0.025,0.975), rate = 5, lower.tail = TRUE, log.p = FALSE)
pexp(0.737775891,5)-pexp(0.005063562,5)

########################################
# 3.2.2 Normal distribution
########################################

# Example 1 . ...
library(TeachingDemos)
vis.normal()

# Example 2. When population X is distributed as N(1.9, .25)...
x <- rnorm(10000,1.9,0.5)
mean(x)
sd(x)

# Example 3. Suppose that the CCND3 (Cyclin D3) gene expression...
pnorm(1.4, 1.9, 0.5)   # left-side tail of the Normal cumulative density function (CDF)

dnormFun <- function(x) { dnorm(x,1.9,0.5) }
x1<- seq(0,1.4,0.1)
x2<- seq(1.4,4,0.1)
plot(dnormFun,      # function
     0,             # start
     4,             # end
     cex.lab=1.5,   # make axis labels big
     xlab="x",
     ylab="Normal probability density function f(x)")
polygon(c(0.0,x1,1.4), c(0,dnormFun(x1),0), col="magenta")
polygon(c(1.4,x2,4.0), c(0,dnormFun(x2),0), col="lightblue")
mtext("1.4",side=1,at=1.4, col="red")
arrows(0.6,0.43,1.0,0.1, lwd=3, col="green")
text(0.3, 0.55 - c(0,.10), cex = 1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 1.4)),
       expression(paste("  = ", integral(f(x) * dx, -infinity, 1.4)))))


pnormFun <- function(x) { pnorm(x,1.9,0.5) }
plot(pnormFun,      # function
     0,             # start
     4,             # end
     cex.lab=1.5,   # make axis labels big
     col="red",
     lwd=6,         # make line thicker
     xlab="x",
     ylab="Normal cumulative distribution function F(x)")
mtext("1.4",side=1,at=1.4, col="red")
mtext("0.16",side=2,at=0.16, col="red")
arrows(0.9,0.45,1.4,0.16, lwd=3, col="green")
text(0.5, 0.9 - c(0,.1,.2,.3,.4), cex=1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 1.4)),
       expression(paste("  = ", integral(f(x) * dx, -infinity, 1.4))),
       expression(paste("  = ", bgroup("[", F(x) ,"]")[-infinity]^1.4)),
       expression(paste("  = ", CDF(1.4))),
       expression(paste("  = ", 0.16))))

1 - pnorm(2.4, 1.9, 0.5)   # right-side tail of normal cumulative density function (CDF)
pnorm(2.4, 1.9, 0.5) - pnorm(1.4, 1.9, 0.5)  # central area of the normal cumulative density function (CDF)

qnorm(0.025,1.9,0.5)   # Normal quantile function
pnorm(0.920018,1.9,0.5)    # left-side tail of the Normal cumulative density function (CDF)

# Normal approximation to the binomial.
x <- 0:20
np=11
xPolygon<- seq(0,20,0.1)
dnormFun <- function(x) { dnorm(x,np,sqrt(np)) }
plot(dnormFun,      # normal density function
     0,             # start
     20,            # end
     col="lightblue",
     ylim=c(0,0.4), # range of the y-axis
     xlab="x",
     ylab="Density function f(x)",
     cex.lab=1.5)   # make axis labels big
polygon(c(0.0,xPolygon,20), c(0,dnormFun(xPolygon),0), col="lightblue")
for (counter in 1:9) {
    lines(x, dbinom(x, size=counter*12, prob=np/(counter*12)), col=counter, lwd=2)
}
legend("topleft", lty=rep(1,10), lwd=rep(2,10), col=c("lightblue",1:9),
       legend=c(expression(Normal(mu == 11, sigma^2 == 11)),
           as.expression(sapply(seq(12,108, 12), function(x) bquote(B(n ==.(x), p == 11/.(x)))))))

# Normal approximation to the Poisson.
x <- 0:20
xPolygon<- seq(0,20,0.1)
dnormFun <- function(x) { dnorm(x,11,sqrt(11)) }
plot(dnormFun,      # normal density function
     0,             # start
     20,            # end
     col="lightblue",
     ylim=c(0,0.4), # range of the y-axis
     xlab="x",
     ylab="Density function f(x)",
     cex.lab=1.5)   # make axis labels big
polygon(c(0.0,xPolygon,20), c(0,dnormFun(xPolygon),0), col="lightblue")
for (counter in 1:11) {
    lines(x, dpois(x,counter), col=counter, lwd=2)
}
legend("topright", lty=rep(1,12), lwd=rep(2,12), col=c("lightblue",1:11),
       legend=c(expression(Normal(mu == 11, sigma^2 == 11)),
           as.expression(sapply(1:11, function(x) bquote(Poisson(lambda ==.(x)))))))

########################################
# 3.2.3 Chi-squared distribution
########################################

# Example 1. Load the TeachingDemos package to view...
library(TeachingDemos)
vis.gamma()

# Example 2. Let's consider the chi-squared variable with 5 degrees of freedom...
pchisq(8, 5)   # left-side tail of the Chi-squared cumulative density function (CDF)

dchisqFun<-function(x) { dchisq(x,5) }
plot(dchisqFun,     # function
     0,             # start
     25,            # end
     cex.lab=1.5,   # make axis labels big
     xlab="x",
     ylab="Chi-Squared probability density function f(x)")
x1 <- seq(0,8,0.1)
x2 <- seq(8,25,0.1)
polygon(c(0,x1,8),  c(0,dchisqFun(x1),0), col="magenta")
polygon(c(8,x2,25), c(0,dchisqFun(x2),0), col="lightblue")
mtext("8",side=1,at=8)
arrows(11,0.07,5,0.07, lwd=3, col="green")
text(13, 0.075 - c(0,.018), cex = 1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 8)),
       expression(paste("  = ", integral(f(x) * dx, 0, 8)))))

pchisqFun<-function(x){pchisq(x,5)}
plot(pchisqFun,     # function
     0,             # start
     20,            # end
     cex.lab=1.5,   # make axis labels big
     col="red",
     lwd=6,         # make line thicker
     xlab="x",
     ylab="Chi-Squared cumulative distribution function F(x)")
mtext("8",   side=1,at=8, col="red")
mtext("0.84",side=2,at=0.84, col="red")
arrows(11,0.74,8,0.84, lwd=3, col="green")
text(12, 0.67 - c(0,.11,.22,.33,.44), cex=1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 8)),
       expression(paste("  = ", integral(f(x) * dx, 0, 8))),
       expression(paste("  = ", bgroup("[", F(x) ,"]")[0]^8)),
       expression(paste("  = ", CDF(8))),
       expression(paste("  = ", 0.84))))

qchisq(0.025, 5, lower.tail=TRUE)

# Example 3. The chi-squared distribution is frequently used as a so-called goodness of fit measure...
library(multtest); data(golub)
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
ccnd3 = grep("CCND3",golub.gnames[ ,2], ignore.case = TRUE)
ccnd3
x <- golub[ccnd3,golubFactor=="ALL"]
z <- (x-1.90)/0.50
sum(z^2)
pchisq(sum(z^2),26, lower.tail=FALSE)

########################################
# 3.2.4 t-Distribution
########################################

# Example 1. Load the TeachingDemos package and execute...
library(TeachingDemos)
vis.t()

# Example 2. A quick NCBI scan makes it reasonable to assume that ...
n <- 11
gdf5 = grep("GDF5",golub.gnames[ ,2], ignore.case = TRUE)
gdf5
x <- golub[gdf5, golubFactor=="AML"]
t.value <- sqrt(n)*(mean(x)-0)/sd(x)
t.value

1 - pt(1.236324, 10)   # right-side tail of the t-distribution cumulative density function (CDF)

f<-function(x) { dt(x,10) }
plot(f,             # function
     -5,            # start
     5,             # end
     cex.lab=1.5,   # make axis labels big
     xlab="x",
     ylab="t-Distribution probability density function f(x)")
x1 <- seq(-5,1.24,0.01)
x2 <- seq(1.24,5,0.01)
polygon(c(-5,x1,1.24), c(0,f(x1),0), col="lightblue")
polygon(c(1.24,x2,5),  c(0,f(x2),0), col="magenta")
mtext("1.24",side=1,at=1.24, col="red")
arrows(2.7,0.20,2.7,0.01, lwd=3, col="green")
text(1.8, 0.3 - c(0,.043, .086), cex = 1.2, adj=c(0,0), col="blue",
     c(expression(P(X >= 1.24)),
       expression(paste("  = ", integral(f(x) * dx, 1.24, infinity))),
       expression(paste("  = ", 1 - integral(f(x) * dx, -infinity, 1.24)))))

F<-function(x) { pt(x,10) }
plot(F,             # function
     -5,            # start
     5,             # end
     cex.lab=1.5,   # make axis labels big
     col="red",
     lwd=6,         # make line thicker
     xlab="x",
     ylab="t-Distribution cumulative distribution function F(x)")
mtext("1.24",side=1,at=1.24, col="red")
mtext("0.88",side=2,at=0.88, col="red")
arrows(-1,0.88,1.24,0.88, lwd=3, col="green")
text(-3.7, 0.85 - c(0,.11,.22,.33,.44,.55), cex=1.2, adj=c(0,0), col="blue",
     c(expression(P(X >= 1.24)),
       expression(paste("  = ", 1 - integral(f(x) * dx, -infinity, 1.24))),
       expression(paste("  = ", 1 - bgroup("[", F(x) ,"]")[-infinity]^1.24)),
       expression(paste("  = ", 1 - CDF(1.24))),
       expression(paste("  = ", 1 - 0.88)),
       expression(paste("  = ", 0.12))))

pt(2, 10) - pt(-2, 10)
qt(0.025, n-1)

########################################
# 3.2.5 F-Distribution
########################################

# Example 1. When the two population variances are in fact equal...
var(golub[1042,golubFactor=="ALL"])/var(golub[1042,golubFactor=="AML"])

pf(0.7116441, 26, 10)   # left-side tail of the F-distribution cumulative density function (CDF)

f<-function(x) { df(x,26,10) }
plot(f,             # function
     0,             # start
     10,            # end
     cex.lab=1.5,   # make axis labels big
     xlab="x",
     ylab="F-Distribution probability density function f(x)")
mtext("0.71",side=1,at=.7,cex=1, col="red")
x1 <- seq(0,0.71,0.01)
x2 <- seq(0.71,10,0.01)
polygon(c(0,x1,.71),  c(0,f(x1),0), col="magenta")
polygon(c(.71,x2,10), c(0,f(x2),0), col="lightblue")
arrows(4.5,.50,0.55,.3, lwd=3, col="green")
text(5.0, 0.53 - c(0,.11), cex = 1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 0.71)),
       expression(paste("  = ", integral(f(x) * dx, 0, 0.71)))))

f<-function(x) { pf(x,26,10) }
plot(f,            # function
     0,            # start
     10,           # end
     cex.lab=1.5,  # make axis labels big
     col="red",
     lwd=6,         # make line thicker
     xlab="x",
     ylab="F-Distribution cumulative distribution function F(x)")
mtext("0.71",side=1,at=.7, cex=1, col="red")
mtext("0.23",side=2,at=.23,cex=1, col="red")
arrows(3,0.3,0.71,0.23, lwd=3, col="green")
text(4.0, 0.5 - c(0,.12,.24,.36,.48), cex=1.2, adj=c(0,0), col="blue",
     c(expression(P(X <= 0.71)),
       expression(paste("  = ", integral(f(x) * dx, 0, 0.71))),
       expression(paste("  = ", bgroup("[", F(x) ,"]")[0]^0.71)),
       expression(paste("  = ", CDF(0.71))),
       expression(paste("  = ", 0.23))))

qf(0.025, 26, 10)

