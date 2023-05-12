###############################################################################
###############################################################################
# Chapter 8
###############################################################################
###############################################################################

########################################
# 8.1 Classification of microRNA
########################################

# Example 1. Classification of Micro RNA. MicroRNA are small RNA molecules...

########################################
# 8.2 ROC types of curves
########################################

# Example 1. For the sake of illustration we consider the prediction of ALL...
data(golub, package = "multtest")
golubLabels <- factor(golub.cl,levels=0:1,labels= c("ALL","not ALL"))
ccnd3 <- grep("CCND3",golub.gnames[,2], ignore.case = TRUE)
golubPredictor <- factor(golub[ccnd3,]>1.27,levels=c("TRUE","FALSE"), labels=c("ALL","notALL"))
table(golubPredictor,golubLabels)

# Example 2. The expression values for gene CCND3 Cyclin D3 from the...
library(ROCR)
golubLabels <- factor(golub.cl,levels=0:1,labels= c("TRUE","FALSE"))
pred <- prediction(golub[ccnd3,], golubLabels)
perf <- performance(pred, "tpr", "fpr" )
plot(perf,
     lwd=4,
     col="magenta")
## dev.copy2eps(device=x11,file="ROCgolub1042.eps")
performance(pred,"auc")

index <- c(1,1)
gdf5 = grep("GDF5",golub.gnames[ ,2], ignore.case = TRUE)
gdf5
data <- c(sort(golub[gdf5,1:27],decreasing = TRUE),sort(golub[gdf5,28:38],decreasing = TRUE))
for (i in 1:38) {
    if (golub.cl[i] == 0) {
        index[i] <- 2 }
    else {
        index[i] <-golub.cl[i]
    }
}
pred <- prediction(data, index)
perf <- performance(pred, "tpr", "fpr" )
plot(perf,
     lwd=4,
     col="magenta")
## dev.copy2eps(device=x11,file="ROCgolub2058.eps")
performance(pred,"auc" )


########################################
# 8.3 Classification trees
########################################

# Example 1. Optimal gene expressions. Suppose microarray expres-...
library(rpart)
library(rpart.plot)
set.seed(123); n<-10 ; sigma <- 0.5
factor <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(factor) <- c("ALL1","ALL2","AML")
geneA <- c(rnorm(10,0,sigma),rnorm(10,2,sigma),rnorm(10,4,sigma))
data <- data.frame(factor,geneA)
rpartFit <- rpart(factor ~ geneA, method="class",data=data)
boxplot(geneA ~ factor, col=c("blue", "red", "darkgreen"))
## dev.copy2eps(device=x11,file="rpartGeneAboxplot.eps")
prp(rpartFit,
    branch.lwd=4, # wide branches
    branch.col="darkgreen",
    extra=101)    # plot label numbers and the percentage of observations
## dev.copy2eps(device=x11,file="rpartGeneATree.eps")

# Example 2. Gene selection. Another situation is where Gene A discrim-...
set.seed(123)
n<-10 ; sigma <- 0.5
factor <- factor(c(rep(1,n),rep(2,n),rep(3,n)))
levels(factor) <- c("ALL1","ALL2","AML")
geneA <- c(rnorm(20,0,sigma),rnorm(10,2,sigma))
geneB <- c(rnorm(10,0,sigma),rnorm(20,2,sigma))
geneC <- c(rnorm(30,1,sigma))
data <- data.frame(factor,geneA,geneB,geneC)
rpartFit <- rpart(factor ~ geneA + geneB + geneC, method="class",data=data)
boxplot(geneA ~ factor, col=c("blue", "blue", "red"))
## dev.copy2eps(device=x11,file="rpartGeneABCboxplot.eps")
prp(rpartFit,
    branch.lwd=4, # wide branches
    branch.col="blue",
    extra=101)    # plot label numbers and the percentage of observations
## dev.copy2eps(device=x11,file="rpartGeneABCTree.eps")

# Example 3. Prediction by CCND3 Cyclin D3 gene expression values....
library(rpart); library(multtest); data(golub)
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
ccnd3 <- grep("CCND3",golub.gnames[,2], ignore.case = TRUE)
golubRpart <- rpart(golubFactor ~ golub[ccnd3,] , method="class")
predictedclass <- predict(golubRpart, type="class")
table(predictedclass, golubFactor)

summary(golubRpart)
predict(golubRpart, type="class")
predict(golubRpart, type="prob")
boxplot(golub[1024,] ~ golubFactor, col=c("blue", "red"))
## dev.copy2eps(device=x11,file="rpartgolubBoxPlot1042.eps")
prp(golubRpart,
     branch.lwd=4, # wide, thick branches
     branch.col="blue",
     extra=101)    # plot label numbers and the percentage of observations
## dev.copy2eps(device=x11,file="rpartgolub1042.eps")

# Example 4. Gene selection of the Golub (1999) data. By recursive...
library(rpart); library(multtest); data(golub)
row.names(golub)<- paste("gene", 1:3051, sep = "")

golubData <- data.frame(t(golub[1:3051,]))
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
golubRpart <- rpart(golubFactor~., data=golubData, method="class", cp=0.001)
prp(golubRpart,
     branch.lwd=4, # wide, thick branches
     branch.col="blue",
     extra=101)    # plot label numbers and the percentage of observations
golub.gnames[896,]

# Example 5. Application to the Chiaretti (2004) data. With respect to...
library("hgu95av2.db"); library(ALL); data(ALL)
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
anova.pValue <- apply(exprs(ALLB123), 1, function(x) anova(lm(x ~ ALLB123$BT))$Pr[1])
names <- featureNames(ALL)[anova.pValue<0.00001]
symb <- mget(names, env = hgu95av2SYMBOL)
ALLBTnames <- ALLB123[names, ]
test<- exprs(ALLBTnames)
probeData <- as.matrix(exprs(ALLBTnames))
row.names(probeData)<-unlist(symb)

diagnosed <- factor(ALLBTnames$BT)
rpartFit <- rpart(diagnosed ~ ., data = data.frame(t(probeData)))
prp(rpartFit,
    branch.lwd=4,
    branch.col="blue",
    extra=101)
dev.copy2eps(device=x11,file="rpartOnALLB123.eps")
rpartPredict <- predict(rpartFit, type="class")
table(rpartPredict,diagnosed)

predicted.class <- predict(rpartFit, type="class")
predicted.probabilities <- predict(rpartFit, type="prob")
out <- data.frame(predicted.probabilities,predicted.class,
                  diagnosis=factor(ALLBTnames$BT))
print(out,digits=2)
?sample
# Example 6. Training and validation. In the setting of B-cell ALL data...
training <- sample(1:78, 39, replace = FALSE)
testing <- setdiff(1:78,training)
df <- data.frame(t(probeData))
rpartFit <- rpart(diagnosed ~ ., data = df, subset=training)

rpart.predictor.t <- predict(rpartFit, df[training,], type="class")
table(rpart.predictor.t,factor(ALLBTnames$BT[training]))

rpart.predictor.v <- predict(rpartFit,df[testing,], type="class")
table(rpart.predictor.v,factor(ALLBTnames$BT[testing]))

########################################
# 8.4 Random Forest
########################################
library(randomForest)
Y <- factor(ALLBTnames$BT);
X <- t(probeData)
rf1 <- randomForest(X, Y, ntree = 1000, importance=TRUE, proximity=TRUE)
rf1

varImpPlot(rf1,
           n.var = 15,
           pch=19,
           main=NULL,
           col="red",
           gcolor="blue",
           lcolor="darkgreen")
dev.copy2eps(device=x11,file="varImpPlotRandomForestOnALLB123.eps")

rft <- randomForest(X[training,], Y[training], ntree = 1000, importance = TRUE, proximity=TRUE)
table(rfpredt=rft$pred,Yt=Y[training])

pv <- predict(rft, newData=X[testing,], type="class", proximity=TRUE)
table(rfpredv=pv$pred,Yv=Y[testing])


########################################
# 8.5 Support Vector Machine
########################################

# Example 1. Application to the Chiaretti (2004) data. The parameters...
library(e1071)
df <- data.frame(Y = factor(ALLBTnames$BT), X =t(probeData))
Y <- factor(ALLBTnames$BT); X <- t(probeData)
Yt <- factor(ALLBTnames$BT)[training]; Yv <- factor(ALLBTnames$BT)[testing]
X <- t(probeData); Xt <- X[training,]; Xv <- X[testing,]
svmEst <- svm(X, Y, data=df, type = "C-classification", kernel = "linear")
svmPredict <- predict(svmEst, X, probability=TRUE)
table(svmPredict, Y)
summary(svmEst)
dim(svmEst$SV)
dim(svmEst$coefs)

# Example 2. Training and validation. A generally applied manner to...
Yt <- factor(ALLBTnames$BT)[training]; Yv <- factor(ALLBTnames$BT)[testing]
X <- t(probeData); Xt <- X[training,]; Xv <- X[testing,]
svmEst <- svm(Xt, Yt, type = "C-classification", kernel = "linear")
svmPredict.t <- predict(svmEst, Xt, probability=TRUE)
table(svmPredict.t, Yt)
svmPredict.v <- predict(svmEst, Xv, probability=TRUE)
table(svmPredict.v, Yv)

########################################
# 8.6 Neural Networks
########################################

# Example 1. Application to the Chiaretti (2004) data. The models can...
Y <- factor(ALLBTnames$BT); X <- t(probeData)
library(nnet)
df <- data.frame(Y = Y, X = X[, sample(ncol(X), 20)])
nnest <- nnet(Y ~ .,data = df, size = 5, maxit = 500, decay = 0.01, MaxNWts = 5000)
pred <- predict(nnest, type = "class")
table(pred, Y) # prints confusion ma

# Example 2. Training and validation. The results from cross validation...
nnest.t <- nnet(Y ~ ., data = df,subset=training, size = 5,decay = 0.01,maxit=500)
prednnt <- predict(nnest.t, df[training,], type = "class")
table(prednnt,Ytrain=Y[training])
prednnv <- predict(nnest.t, df[testing,], type = "class")
table(prednnv, Yval= Y[testing])

########################################
# 8.7 Generalized Linear Models
########################################

# Example 1. CCND3 Cyclin D3. In the Golub et al. (1999) data we...
library(faraway)
ccnd3 <- grep("CCND3",golub.gnames[,2], ignore.case = TRUE)
logitmod <- glm((-golub.cl + 1) ~ golub[ccnd3,],
                family=binomial(link = "logit"))
pchisq(deviance(logitmod),df.residual(logitmod),lower=FALSE)
plot((-golub.cl + 1) ~ golub[ccnd3,], xlim=c(-2,5), ylim = c(0,1),
     xlab="CCND3 expression values ", ylab="Probability of ALL")
x <- seq(-2,5,.1)
lines(x,ilogit(-4.844124 + 4.439953*x))
pchisq(deviance(logitmod),df.residual(logitmod),lower=FALSE)
summary(logitmod)

predictor <- predict(logitmod,type="response") > 0.5
predictor.fac <- factor(predictor,levels=c(TRUE,FALSE),labels=c("ALL","not ALL"))
table(predictor.fac,golubFactor)

# Example 2. Application to the Chiaretti (2004) data. With respect to...
library(nnet); library("hgu95av2.db"); library(ALL); data(ALL)
probe.names <- c("1389_at","35991_at","40440_at")
ALLB123 <- ALL[,ALL$BT %in% c("B1","B2","B3")]
probeData <- exprs(ALLB123)[probe.names,]
row.names(probeData) <- unlist(mget(probe.names, env = hgu95av2SYMBOL))
factor <- factor(ALLB123$BT,levels=c("B1","B2","B3"))
data <- data.frame(factor,t(probeData))
mnmod <- multinom(factor ~ ., family=binomial(link = "logit"),data=data)
summary(mnmod)

predict.mn <- predict(mnmod,type="class")
table(predict.mn,factor)


