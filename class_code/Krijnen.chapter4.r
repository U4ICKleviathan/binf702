###############################################################################
###############################################################################
# Chapter 4
###############################################################################
###############################################################################

########################################
# 4.1 Statistical hypothesis testing
########################################

########################################
# 4.1.1 The Z-test
########################################

# Example 1. To illustrate the Z-test we shall concentrate on the Gdf5 ...
data(golub, package = "multtest")
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
sigma <- 0.25; n <- 27; mu0 <- 0
x <- golub[2058,golubFactor=="ALL"]
z.value <- sqrt(n)*(mean(x) - mu0)/sigma

2*pnorm(-abs(z.value),0,1)

f<-function(x){dnorm(x,0,1)}
x1 <- seq(-4,-1.96,0.01)
y1 <- dnorm(x1,0,1)
x2 <- seq(-1.96,1.96,0.01)
y2 <- dnorm(x2,0,1)
x3 <- seq(1.96,4,0.01)
y3 <- dnorm(x3,0,1)
plot(f,    # function
     -4,   # begin x-value
     4,    # end x-value
     cex.lab=1.5,  # make axis labels big
     xlab="x",
     ylab="Normal probability density function f(x)")
polygon(c(-4,x1,-1.96), c(0,y1,0), col="red")
polygon(c(-1.96,x2,1.96), c(0,y2,0), col="lightblue")
polygon(c(1.96,x3,4), c(0,y3,0), col="red")
arrows(-3,0.15,-3,0.03)
text(-3,0.23,"Rejection")
text(-3,0.20,"Region")
text(-3,0.17,expression(alpha/2))
arrows(3,0.15,3,0.03)
text(3,0.23,"Rejection")
text(3,0.20,"Region")
text(3,0.17,expression(alpha/2))
text(0,0.23,"Acceptance")
text(0,0.20,"Region")

# Example 2. Using the data from Example 1, the 95% confidence interval...
mean(x)+qnorm(c(0.025),0,1)*sigma/sqrt(n)
mean(x)+qnorm(c(0.975),0,1)*sigma/sqrt(n)

library(TeachingDemos)
z.test(x,mu=0,sd=0.25)

# Example 3. To develop intuition with respect to confidence intervals...
ci.examp(mean.sim =0, sd = 1, n = 25, reps = 100,
         method = "z", lower.conf=0.025, upper.conf=0.975)

########################################
# 4.1.2 One Sample t-Test
########################################

f<-function(x){dt(x,5)}
x1 <- seq(-4,qt(0.025,5),0.01)
y1 <- f(x1)
x2 <- seq(qt(0.025,5),qt(0.975,5),0.01)
y2 <- f(x2)
x3 <- seq(qt(0.975,5),4,0.01)
y3 <- f(x3)
plot(f,    # function
     -4,   # begin x-value
     4,    # end x-value
     xlab="x",
     ylab="t-Distribution probability density function f(x)")
polygon(c(-4,x1,qt(0.025,5)), c(0,y1,0), col="red")
polygon(c(qt(0.025,5),x2,qt(0.975,5)), c(0,y2,0), col="lightblue")
polygon(c(qt(0.975,5),x3,4), c(0,y3,0), col="red")
arrows(-3,0.15,-3,0.03)
text(-3,0.23,"Rejection")
text(-3,0.20,"Region")
text(-3,0.17,expression(alpha/2))
arrows(3,0.15,3,0.03)
text(3,0.23,"Rejection")
text(3,0.20,"Region")
text(3,0.17,expression(alpha/2))
text(0,0.23,"Acceptance")
text(0,0.20,"Region")
mtext(expression(t[0.025]),side=1,at=qt(0.025,5), col="red")
mtext(expression(t[0.975]),side=1,at=qt(0.975,5), col="red")

# Example 1. Let's test H0 : ¹ = 0 against H1 : ¹ 6= 0 for the ALL...
x <- golub[2058,golubFactor=="ALL"]; mu0 <- 0; n <- 27
t.value <- sqrt(n)*(mean(x) - mu0)/sd(x)
t.value

2 * pt(-0.0010, n-1)
c(qt(0.025, n-1), qt(0.975, n-1))

mean(x)+qt(0.025,26)*sd(x)/sqrt(n)
t.test(x,mu=0)

# Example 2. In Chapter 2 a box-and-whiskers plot revealed that the...
ccnd3 <- grep("CCND3",golub.gnames[,2], ignore.case = TRUE)
ccnd3
t.test(golub[ccnd3,golubFactor=="ALL"],mu=0, alternative = c("greater"))

########################################
# 4.1.3 Two-sample t-test with unequal variances
########################################

# Example 1. Golub et al. (1999) argue that gene CCND3 Cyclin D3 plays...
t.test(golub[ccnd3,] ~ golubFactor, var.equal=FALSE)

########################################
# 4.1.4 two-sample t-test with equal variances
########################################

# Example 1. The null hypothesis for gene CCND3 Cyclin D3 that the...
t.test(golub[ccnd3,] ~ golubFactor, var.equal = TRUE)


########################################
# 4.1.5 Two-sample paired t-test
########################################
before = c(112.9, 133.5, 122.8, 145.6, 167.2, 169.2, 174.6, 175.3, 184.4, 121.3)
after  = c(120.7, 123.6, 134.0, 141.2, 156.8, 170.0, 182.0, 185.9, 196.0, 115.1)
t.test(before, after, paired=TRUE)

########################################
# 4.1.6 F-test on equal variances
########################################

# Example 1. The null hypothesis for gene CCND3 Cyclin D3 that the...
var.test(golub[ccnd3,] ~ golubFactor)

########################################
# 4.1.7 Binomial test
########################################

# Example 1. A microRNA of length 22 contains 18 purines. The null...
1 - pbinom(17, 22, 0.7)

binom.test(18, 22, p = 0.7, alternative = c("greater"),
           conf.level = 0.95)

########################################
# 4.1.8 Chi-squared test
########################################

# Example 1. Suppose we want to test the hypothesis that the nucleotides...
library(ape)
zyxinfreq <- table(read.GenBank(c("X94991.1"),as.character=TRUE))
chisq.test(zyxinfreq)

f<-function(x){dchisq(x,3)}
plot(f,    # function
     0,    # begin x-value
     25,   # end x-value
     cex.lab=1.5,  # make axis labels big
     xlab="q",
     ylab="Chi-squared probability density function f(q)")
x1<- seq(0,qchisq(0.95,3),0.01)
x2<- seq(qchisq(0.95,3),20,0.01)
polygon(c(0,x1,qchisq(0.95,3)), c(0,f(x1),0), col="lightblue")
polygon(c(qchisq(0.95,3),x2,20), c(0,f(x2),0), col="red")
arrows(10,0.045,10,0.02)
text(10,0.06,"Rejection")
text(10,0.05,"Region")
text(3,0.04,"Acceptance")
text(3,0.03,"Region")
mtext("7.8",side=1,at=7.86)

# Example 2. In the year 1866 Mendel observed in large number of exper-...
pi <- c(0.75,0.25)
x <-c(5474, 1850)
chisq.test(x, p=pi)
qchisq(0.95, 3)

# Example 3. Given certain expression values for a healthy control group...
data <- matrix(c(5,5,5,5),2,byrow=TRUE)
chisq.test(data)
data <- matrix(c(8,2,2,8),2,byrow=TRUE)
chisq.test(data)

########################################
# 4.1.9 Fisher’s exact test
########################################

# Example 1. Example 1: Oncogenes on chromosome 1....
data <- matrix(c(300,500,3000,7000),2,byrow=TRUE)
fisher.test(data)

########################################
# 4.1.10 Normality tests
########################################

# Example 1. To test the hypothesis that the ALL gene expression values...
shapiro.test(golub[ccnd3, golubFactor=="ALL"])

library(nortest)
ad.test(golub[ccnd3,golubFactor=="ALL"])

########################################
# 4.1.11 Outliers test
########################################

# Example 1. From Figure 2.4 we have observed that expression values...
library(outliers)
grubbs.test(golub[ccnd3, golubFactor=="ALL"])

########################################
# 4.1.12 Non-Parametric Tests
########################################

# One-sample sign test
# Example 1: Small sample size.
# install.packages(BSDA)
library(BSDA)
x<-c(6003, 6304, 6478, 6245, 6134, 6204, 6150, 6345, 6298, 6194)
SIGN.test(x, md = 6000)

# One-sample Wilcoxon signed-rank test
# Example 1: Small sample size.
x <- c(6003, 6304, 6478, 6245, 6134, 6204, 6150)
wilcox.test(x,mu=6000)

# Two-sample Wilcoxon rank-sum test
# Example 1. CCND3 gene expression. The null hypothesis that the expression values for gene...
wilcox.test(golub[ccnd3,] ~ golubFactor)

# Two-sample paired sign test
# Example 1: pathway gene expression.
before = c(112.9, 133.5, 122.8, 145.6, 167.2, 169.2, 174.6, 175.3, 184.4, 121.3)
after  = c(120.7, 123.6, 134.0, 141.2, 156.8, 170.0, 182.0, 185.9, 196.0, 115.1)
SIGN.test(before, after, alternative = "two.sided")

# Two-sample paired Wilcoxon rank-sum test
# Example 1: pathway gene expression.
before = c(112.9, 133.5, 122.8, 145.6, 167.2, 169.2, 174.6, 175.3, 184.4, 121.3)
after  = c(120.7, 123.6, 134.0, 141.2, 156.8, 170.0, 182.0, 185.9, 196.0, 115.1)
wilcox.test(before, after, paired=TRUE)

# Non-Parametric Bootstrapping
# Example 1: Small sample size. Lets use the bootstrap...
library(boot)
x <- c(6003, 6304, 6478, 6245, 6134)
boot.mean <- boot(x, function(x,i){mean(x[i])}, R=9999)
boot.mean

mean(boot.mean$t)
sd(boot.mean$t)

boot.ci(boot.mean, conf = c(0.95), type = c("perc"))

pvec <- c(0.025,0.975)
quantile(boot.mean$t,pvec)

########################################
# 4.1.13 Robust Estimation
########################################

# Example 1: Normal distribution with an outlier. Lets construct...
x <- rnorm(20,10,3)
mean(x)
x[21] <- 1000
mean(x)

library(MASS)
unlist(huber(x))

library(boot)
boot.mean <- boot(x, function(x,i){mean(x[i])}, R=9999)
jack.after.boot(boot.mean)


########################################
# 4.1.14 Type I and type II errors
########################################

########################################
# 4.1.15 Power of a statistical test
########################################

# Example 1. With respect to the normal distribution it is...
qnorm(0.95,0,1)
pnorm(1.644854,3.5,1)

# Example 2: Teaching demonstration. To further visualize...
library(TeachingDemos)
run.power.examp()

# Example 3. If in case of a one sample t-Test an estimate...
power.t.test(n=4, delta=300, sd=175, type = c("one.sample"))

########################################
# 4.2 Maximum likelihood Estimation
########################################

# Example 5. Various parameters can be estimated...
x <- rpois(100,4)
fitdistr(x,"Poisson")

x <- c(6003, 6304, 6478, 6245, 6134)
fitdistr(x,"Normal")

x <- rexp(40000,5)
fitdistr(x,"exponential")


# Example 6. Frequently it is assumed that the data are...
n <- 100
data <- matrix(0,nrow=n,ncol = 1)
for (i in 1:n) {
  p <- sample(c(0,1), 1, replace = TRUE, prob=c(0.4,0.6))
  data[i] <- p * rnorm(1,3000,300) +  (1-p) * rnorm(1,6000,300)
}

dmixnor <- function(x,p,m1,s1,m2,s2){p * dnorm(x,m1,s1) +  (1-p) * dnorm(x,m2,s2)}
fitdistr(data,dmixnor,list(p=0.6,m1=3000,s1=300,m2=6000,s2=300))


########################################
# 4.3 Application of tests to a whole gene set
########################################

# Example 1. Having a data matrix with gene expression values, a ques-...
data(golub,package="multtest")
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
sh <- apply(golub[,golubFactor=="ALL"], 1, function(x) shapiro.test(x)$p.value)
sum(sh > 0.05)/nrow(golub) * 100

# Example 2. In case the gene expression data are non-normally dis-...
data(golub, package = "multtest")
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
pt <- apply(golub, 1, function(x) t.test(x ~ golubFactor)$p.value)
pw <- apply(golub, 1, function(x) wilcox.test(x ~ golubFactor)$p.value)
result <- data.frame(cbind(pw,pt))
result[pw<0.05 & abs(pt-pw)>0.2,]



