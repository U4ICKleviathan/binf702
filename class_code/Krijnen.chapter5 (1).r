###############################################################################
###############################################################################
# Chapter 5
###############################################################################
###############################################################################

########################################
# 5.1 Definition of linear models
########################################

# Example 1. A common manner to introduce the linear model is by writing...
library(TeachingDemos)
put.points.demo()

zyxin = grep("Zyxin",golub.gnames[,2], ignore.case = TRUE)
cmyb = grep("c-myb",golub.gnames[,2], ignore.case = TRUE)
x <- golub[zyxin,]
y <- golub[cmyb,]
leastSquares = lm(y  ~ x)  # linear model with least squares regression with an intercept
plot(x,
     y,
     pch=19,       # plot solid circles
     cex.lab=1.5,  # make axis labels big
     col="blue",
     xlab="Relative Zyxin gene expression",
     ylab="Relative c-MYB gene expression")
predicted <- predict(leastSquares)   # get the predicted values
for (i in 1:length(predicted)) {     # add the residuals
    segments(x[i],predicted[i], x[i],y[i], lwd=2, lty=2, col="green")
}
abline(leastSquares$coef, lwd=3, lty=2, col="red")    # add regression line

lmSummary = summary(leastSquares)
lmSummary

f.value1 = (coef(lmSummary)[2, "t value"])^2   # get the t-value of the 1 predictor variable
f.value1

# Example 2. Suppose we have the following artificial gene expressing values...
y <- c(2,3,1,2, 8,7,9,8, 11,12,13,12)
factor <- gl(3,4)
factor
model.matrix(y ~ factor - 1)
summary(lm(y ~ factor - 1))

########################################
# 5.2 One-way analysis of variance
########################################

# Example 1. Let's continue with the data from the previous example....
y <- c(2,3,1,2, 8,7,9,8, 11,12,13,12)

factor <- gl(3,4)
factor

model.matrix(y ~ factor - 1)

summary(lm(y ~ factor - 1))

y <- c(2,3,1,2, 8,7,9,8, 11,12,13,12)
factor <- gl(3,4)
groupMeans <- as.numeric(tapply(y, factor, mean))
groupMeans

groupMeans <- as.numeric(tapply(y, factor, mean))
g <- 3; n <- 4; N <-12; ssb <- 0
for (j in 1:g) {
    ssb  <-  ssb + (groupMeans[j] - mean(y))^2
}
SSB <- n*ssb
SSB

SSW <- 0
for (j in 1:g) {
    SSW  <-  SSW + sum((y[factor==j] - groupMeans[j])^2)
}
SSW
f.value <- (SSB/(g-1))/(SSW/(N-g))
f.value

1 - pf(f.value,g-1,N-g)

anova(lm(y ~ factor))

# Example 2. By the previous analysis of variance it is concluded that...
summary(lm(y ~ factor))

# Example 3. Let's sample data from the normal distribution with mean...
y <- rnorm(12,1.9,0.5)
round(y,2)
factor <- gl(3,4)
anova(lm(y ~ factor))$Pr[1]

# Example 4. B-cell ALL: SKI-like oncogene. To illustrate analysis of variance...
library(ALL); data(ALL)
data(ALL,package="ALL")
samplesB1toB3 <- ALL$BT %in% c("B1","B2","B3")
x <- as.numeric(exprs(ALL)[row.names(exprs(ALL))=="1866_g_at",samplesB1toB3])
factor <- factor(ALL$BT[samplesB1toB3],labels=c("B1","B2","B3"))
stripchart(x ~ factor,
           method="jitter",  # add random horizontal jitter
           cex.lab=1.5,      # make axis labels big
           vertical = TRUE,  # boxplots vertical
           col=c("red", "darkgreen", "blue"),
           xlab="B-cell ALL stage",
           ylab="SKI-like oncogene expression")

ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1866_g_at",]
summary(lm(y ~ ALLB123$BT))

# Example 5. B-cell ALL: Ets2 repressor. To illustrate a case where the means...
library(ALL); data(ALL)
data(ALL,package="ALL")
x <- as.numeric(exprs(ALL)[row.names(exprs(ALL))=="1242_at",samplesB1toB3])
factor <- factor(ALL$BT[samplesB1toB3],labels=c("B1","B2","B3"))
stripchart(x ~ factor,
           method="jitter",  # add random horizontal jitter
           cex.lab=1.5,      # make axis labels big
           vertical = TRUE,  # boxplots vertical
           col=c("red", "darkgreen", "blue"),
           xlab="B-cell ALL stage",
           ylab="Ets2 expression")

ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1242_at",]
summary(lm(y ~ ALLB123$BT))

anova(lm(y ~ ALLB123$BT))

# Example 6. An interesting question is ...
anova.pValues <- apply(exprs(ALLB123),1,function(x) anova(lm(x~ALLB123$BT))$Pr[1])
sum(anova.pValues<0.05)

########################################
# 5.3 Two-way analysis of variance
########################################

# Example 1. A two-way approach. Looking at the ALL data from Chiaretti...
library("ALL"); data(ALL)
ALLBm <- ALL[,which(ALL$BT %in% c("B","B1","B2","B3","B4") & ALL$mol.biol %in% c("BCR/ABL","NEG"))]
factorMolbio <- factor(ALLBm$mol.biol)
factorB <- factor(ceiling(as.integer(ALLBm$BT)/3),levels=1:2,labels=c("B012","B34"))

anova(lm(exprs(ALLBm)["32069_at",] ~ factorB * factorMolbio))
summary(lm(exprs(ALLBm)["32069_at",] ~ factorB * factorMolbio))
anova.pValues <- apply(exprs(ALLBm), 1, function(x) anova(lm(x ~ factorB * factorMolbio))$Pr[1:3])
anova.pValues.t <- data.frame(t(anova.pValues))
colnames(anova.pValues.t) <- c("mainEffectB","mainEffectMolbio","interaction")
sum(anova.pValues.t$mainEffectB < 0.05 & anova.pValues.t$mainEffectMolbio < 0.05 & anova.pValues.t$interaction < 0.05)

########################################
# 5.4 Checking assumptions
########################################

# Example 1. Testing normality of the residuals. From Figure 5.1 it can...
data(ALL,package="ALL"); library(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1866_g_at",]
shapiro.test(residuals(lm(y ~ ALLB123$BT)))

# Example 2. Testing homoscedasticity of the residuals. From Figure...
library(ALL); data(ALL);
library(lmtest)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1866_g_at",]
bptest(lm(y ~ ALLB123$BT),studentize = FALSE)

########################################
# 5.5 Robust & nonparametric tests
########################################

############################################################
# 5.5.1 One independent (predictor) variable (one-way)
############################################################

# Example 1. In Example 2 of the previous section the hypothesis of...
data(ALL,package="ALL"); library(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1866_g_at",]
oneway.test(y ~ ALLB123$BT)

# Example 2. In Example 1 of the previous section we rejected the hypoth-...
data(ALL,package="ALL"); library(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
y <- exprs(ALLB123)["1866_g_at",]
kruskal.test(y ~ ALLB123$BT)

kruskal.pValues <- apply(exprs(ALLB123),1,function(x) kruskal.test(x~ALLB123$BT)$p.value)
sum(kruskal.pValues<0.05)


###############################################################################
# 5.5.2 Repeated measures one-way and
#       two independent (predictor) variables (two-way)
###############################################################################

#Example 1: Gene expression by gender and cancer stage.
geneExpr = c(.78,.83,.38, 9.73,8.62,9.55,  2.3,3.5,4.2, 13.4,12.9,14.3)
gender = c("F","F","F","M","M","M",  "F","F","F","M","M","M")
stage  = c("1","2","3","1","2","3",  "4","5","6","4","5","6")
gender.df = factor(gender)
stage.df = factor(stage)
table(gender.df, stage.df)
friedman.test(geneExpr ~ gender.df | stage.df)
friedman.test(geneExpr ~ stage.df | gender.df)
