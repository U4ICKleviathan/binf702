###############################################################################
# Answers to exercises in Chapter 7: Cluster Analysis and Trees.
###############################################################################

# 1. Cluster analysis on the "Zyxin" expression values of the Golub et al. (1999) data.
data(golub, package="multtest")
library("cluster")

zyxin <- grep("Zyxin",golub.gnames[,2])
zyxin
data <- data.frame(golub[zyxin,])
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
stripchart(golub[zyxin,]~golubFactor, pch=as.numeric(golubFactor), method="jitter")
plot(hclust(dist(data,method="euclidian"),method="single"))

cl <- kmeans(data, centers=2)
table(cl$cluster,golubFactor)
stripchart(golub[zyxin,]~factor(cl$cluster), method="jitter", col=1:2)
points(x=cl$centers[1], y=1, col = 1, pch = 8, cex=2)
points(x=cl$centers[2], y=2, col = 2, pch = 8, cex=2)

initial <- as.matrix(tapply(golub[zyxin,],golubFactor,mean), nrow = 2, ncol=1)
cl <- kmeans(data, centers=initial)
stripchart(golub[zyxin,]~factor(cl$cluster), method="jitter", col=1:2)
points(x=cl$centers[1], y=1, col = 1, pch = 8, cex=2)
points(x=cl$centers[2], y=2, col = 2, pch = 8, cex=2)
table(cl$cluster,golubFactor)

n <- nrow(data); nboot<-1000
boot.cl <- matrix(0,nrow=nboot,ncol = 2)

for (i in 1:nboot) {
    dat.star <- data[sample(1:n,replace=TRUE),]
    cl <- kmeans(dat.star, initial, nstart = 10)
    boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
}

quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))

# 2. Close to CCND3 Cyclin D3.
library("genefilter"); data(golub, package = "multtest")
ccnd3 <- grep("CCND3",golub.gnames[,2])
closeg <- genefinder(golub, ccnd3, 10, method = "euc", scale = "none")
closeg
golub.gnames[closeg[[1]][[1]],2:3]
boxplot(golub[closeg[[1]][[1]][1],] ~golubFactor)
boxplot(golub[closeg[[1]][[1]][2],] ~golubFactor)
boxplot(golub[closeg[[1]][[1]][3],] ~golubFactor)
boxplot(golub[closeg[[1]][[1]][4],] ~golubFactor)

cyclins <- grep("Cyclin",golub.gnames[,2])
golub.gnames[cyclins,2]
dist.cyclin <- dist(golub[cyclins,],method="euclidian")
distanceMatrix <- as.matrix(dist.cyclin)
rownames(distanceMatrix) <- colnames(distanceMatrix) <- golub.gnames[cyclins,3]
distanceMatrix[1:5,1:5]
# Most of the distances given by genefinder are smaller than those between the Cyclin genes.

# 3. MCM3.
data(golub, package = "multtest")
mcm3 <- grep("MCM3",golub.gnames[,2])
x <- golub[mcm3[1],]; y <- golub[mcm3[2],]
plot(x,y)
which.min(y) # the plot suggests the smallest y as the outlier
cor.test(x,y)
cor.test(x[-21],y[-21])

nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1)
data <- matrix(c(x,y),ncol=2,byrow=FALSE)
for (i in 1:nboot) {
    dat.star <- data[sample(1:nrow(data),replace=TRUE),]
    boot.cor[i,] <- cor(dat.star)[2,1]
}
mean(boot.cor)
quantile(boot.cor[,1],c(0.025,0.975))

nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1)
data <- matrix(c(x[-21],y[-21]),ncol=2,byrow=FALSE)
for (i in 1:nboot) {
    dat.star <- data[sample(1:nrow(data),replace=TRUE),]
    boot.cor[i,] <- cor(dat.star)[2,1]
}
mean(boot.cor)
quantile(boot.cor[,1],c(0.025,0.975))

# 4. Cluster analysis on part of Golub data.
library(multtest); data(golub)
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
o1 <- grep("oncogene",golub.gnames[,2])
plot(hclust(dist(t(golub[o1,]),method="euclidian"),method="single"), hang=-1)

o2 <- grep("antigen",golub.gnames[,2])
plot(hclust(dist(t(golub[o2,]),method="euclidian"),method="single"), hang=-1)

o3 <- grep("receptor",golub.gnames[,2])
plot(hclust(dist(t(golub[o3,]),method="euclidian"),method="single"), hang=-1)

# 5. Principal Components Analysis on part of the ALL data.
library(ALL); data(ALL)
ALLB <- ALL[,ALL$BT %in% c("B1","B2","B3")]
anova.pValues <- apply(exprs(ALLB), 1, function(x) { anova(lm(x ~ ALLB$BT))$Pr[1] } )
ALLBsp <- ALLB[anova.pValues<0.001,]
dim(exprs(ALLBsp))
min(cor(exprs(ALLBsp)))
eigen(cor(exprs(ALLBsp)))$values[1:5]
data <- exprs(ALLBsp); p <- ncol(data); n <- nrow(data) ; nboot<-1000



eigenvalues <- array(dim=c(nboot,p))

for (i in 1:nboot) {
    dat.star <- data[sample(1:n,replace=TRUE),]
    eigenvalues[i,] <- eigen(cor(dat.star))$values
}

for (j in 1:p) {
    print(quantile(eigenvalues[,j],c(0.025,0.975)))
}

library(ggplot2)
library(ggfortify)
library(devtools)
library(ggfortify);
library(ggplot2)
library(devtools)
install_github("vqv/ggbiplot", force=TRUE)


ALLBsp_exprs = t(exprs(ALLBsp))
ALLBsp_exprs_bt <- data.frame(ALLBsp_exprs, BT = ALLBsp$BT)
dim(ALLBsp_exprs_bt)
allb.exp.matrix <- matrix(data = NA, nrow = nrow(ALLBsp_exprs), ncol = ncol(ALLBsp_exprs))
for(i in 1:ncol(allb.exp.matrix)) {allb.exp.matrix[,i] = ALLBsp_exprs[,i]}

?data.frame

table(ALLBsp_exprs_bt$BT)
pca <- prcomp(allb.exp.matrix)
screeplot(pca)
summary(pca)

colorPalette <- c("#0072B2", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#000000", "#D55E00", "#CC79A7")

?autoplot
ggplot2::autoplot(pca,
         data = ALLBsp_exprs_bt,
         colour = 'BT') +
  scale_colour_manual(values = colorPalette) +
  ggtitle("PCA of Significantly Differentially Expressed Genes")

ggbiplot(pca)





autoplot(kmeans(ALLBsp_exprs, 3), data = ALLBsp_exprs)
library(ggpubr)
library(factoextra)
res.km <- kmeans(scale(ALLBsp_exprs), 3, nstart = 25)

fviz_cluster(res.km, data = ALLBsp_exprs,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
# https://www.datanovia.com/en/blog/k-means-clustering-visualization-in-r-step-by-step-guide/
# http://rstudio-pubs-static.s3.amazonaws.com/53162_cd16ee63c24747459ccd180f69f07810.html
# 6. Some correlation matrices.
eigen(matrix(c(1,-0.8,-0.8,1),nrow=2))
eigen(matrix(c(1,0.8,0.8,0.8,1,0.8,0.8,0.8,1),nrow=3))
eigen(matrix(c(1,-0.5,-0.5,-0.5,1,-0.5,-0.5,-0.5,1),nrow=3))
2.6/3 * 100
eigen(matrix(c(1,0.8,0.8,0.8,1,0.8,0.8,0.8,1),nrow=3))$vectors
