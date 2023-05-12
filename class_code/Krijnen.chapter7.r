###############################################################################
###############################################################################
# Chapter 7
###############################################################################
###############################################################################

########################################
# 7.1 Distance
########################################

# Example 1. Let a = 1 and b = 3. Then, obviously, the distance d(1, 3) = 2....

# Example 2. Suppose that a = (a1, a2) = (1, 1) and b = (b1, b2) = (4, 5)....

# Example 3. To compute the Euclidian distance between two vectors one...
a <- c(1,1); b <- c(4,5)
sqrt(sum((a-b)^2))

# Example 4. Distances between Cyclin gene expressions. By the build-in-...
library(multtest); data(golub)
cyclins <- grep("Cyclin",golub.gnames[,2])
golub.gnames[cyclins,2]
dist.cyclin <- dist(golub[cyclins,],method="euclidian")
distanceMatrix <- as.matrix(dist.cyclin)
rownames(distanceMatrix) <- colnames(distanceMatrix) <- golub.gnames[cyclins,3]
distanceMatrix[1:5,1:5]

# Example 5. Finding the ten closest genes to a given one. After selecting...
library("genefilter"); library("ALL"); data(ALL)
closeto1389_at <- genefinder(ALL, "1389_at", 10, method = "euc")
closeto1389_at[[1]]$indices
round(closeto1389_at[[1]]$dists,1)
featureNames(ALL)[closeto1389_at[[1]]$indices]

########################################
# 7.2 Two types of Cluster Analysis
########################################
library(pvclust)

########################################
# 7.2.1 Single Linkage
########################################

# Example 1. An explanatory example. To illustrate single linkage cluster...
names <- list(c("g1","g2","g3","g4","g5"),c("Coordinate p1 (patient 1 gene expression)","Coordinate p2 (patient 2 gene expression)"))
sl.clus.dat <- matrix(c(1,1,1,1.3,3,2,3,2.4,5,5),ncol = 2, byrow = TRUE,dimnames = names)
plot(sl.clus.dat,
     pch=19,
     col="blue",
     cex=1.4,
     xlim=c(0,6),
     ylim=c(0,6))
text(sl.clus.dat,labels=row.names(sl.clus.dat), pos=4, col="red", cex=1.6)
print(dist(sl.clus.dat,method="euclidian"),digits=3)

sl.out<-hclust(dist(sl.clus.dat,method="euclidian"),method="single")
plot(sl.out,
     lwd=3,
     col="blue",
     ## col.lab = "brown",
     col.axis = "brown",
     ylab="Distance",
     xlab="Clustering of the expression of 5 genes",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 3.5, .5), col = "brown",labels = TRUE, lwd = 4)

# Example 2. Relating data generation processes to cluster trees. It is of...
x<-rnorm(20,0,1)
singlelinkage.out<-hclust(dist(x,method="euclidian"),method="single")
plot(singlelinkage.out,
     lwd=3,
     col="blue",
     col.axis = "brown",
     ylab="Distance",
     xlab="20 genes with normal random distances",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 1.4, .2), col = "brown",labels = TRUE, lwd = 4)

x2 = matrix(rnorm(200,0,1), ncol=20) # 10 samples
fit <- pvclust(x2, method.hclust="average", method.dist="euclidean", nboot=1000)
plot(fit, main="", cex.lab=1.5) # dendogram with p values

x <- c(rnorm(10,0,0.1),rnorm(10,3,0.5),rnorm(10,10,1.0))
sl.out<-hclust(dist(x,method="euclidian"),method="single")
plot(sl.out,
     lwd=3,
     col.axis = "brown",
     col="magenta",
     ylab="Distance",
     xlab="3 clusters of 10 genes each",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 5, 1), col = "brown",labels = TRUE, lwd = 4)

x2 = matrix(c(rnorm(100,0,0.1),rnorm(100,3,0.5),rnorm(100,10,1.0)), ncol=30)
fit <- pvclust(x2, method.hclust="average", method.dist="euclidean", nboot=1000)
plot(fit, main="", cex.lab=1.5) # dendogram with p values

# Example 3. Application to the Golub (1999) data. Recall that the first...
data(golub, package="multtest")
zyxin <- grep("Zyxin",golub.gnames[,2])
ccnd3 <- grep("CCND3",golub.gnames[,2])
clusdata <- data.frame(golub[ccnd3,],golub[zyxin,])
colnames(clusdata)<-c("CCND3 Cyclin D3","Zyxin")
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))
plot(clusdata,
     pch=as.numeric(golubFactor) + 15,
     col=as.numeric(golubFactor) + 1)
legend("topright",
       legend=c("ALL","AML"),
       pch=16:17,
       col=c(2,3))

plot(hclust(dist(clusdata,method="euclidian"),method="single"),
     lwd=3,
     col="blue",
     col.axis = "brown",
     ylab="Distance",
     xlab="Clustering of patients by gene expression",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 1.2, .2), col = "brown",labels = TRUE, lwd = 4)

########################################
# 7.2.2 k-means
########################################

# Example 1. Relating a data generation process to k-means cluster analysis....
data <- rbind(matrix(rnorm(100,0,0.5), ncol = 2),
              matrix(rnorm(100,2,0.5), ncol = 2))
cl <- kmeans(data, 2)
plot(data, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)
initial <- matrix(c(0,0,2,2), nrow = 2, ncol=2, byrow=TRUE)
cl<- kmeans(data, initial, nstart = 10)
n <- 100; nboot<-1000
boot.cl <- matrix(0,nrow=nboot,ncol = 4)

for (i in 1:nboot) {
    dat.star <- data[sample(1:n,replace=TRUE),]
    cl <- kmeans(dat.star, initial, nstart = 10)
    boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
}

quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))
quantile(boot.cl[,3],c(0.025,0.975))
quantile(boot.cl[,4],c(0.025,0.975))

# Example 2. Application to the Golub (1999) data. In the above we found...
data <- data.frame(golub[ccnd3,],golub[zyxin,])
colnames(data) <- c("CCND3 (Cyclin D3)","Zyxin")
cl <- kmeans(data, 2,nstart = 10)
cl
boot.cl <- matrix(0,nrow=nboot,ncol = 4)

plot(data, col = cl$cluster)
points(cl$centers, col = 1:2, pch = 8, cex=2)

for (i in 1:nboot) {
  dat.star <- data[sample(1:nrow(data),replace=TRUE),]
  cl <- kmeans(dat.star, initial, nstart = 10)
  boot.cl[i,] <- c(cl$centers[1,],cl$centers[2,])
}


quantile(boot.cl[,1],c(0.025,0.975))
quantile(boot.cl[,2],c(0.025,0.975))
quantile(boot.cl[,3],c(0.025,0.975))
quantile(boot.cl[,4],c(0.025,0.975))

########################################
# 7.3 The correlation coefficient
########################################

# Example 1. Teaching demonstration. To develop intuition with respect to...
library(TeachingDemos)
run.cor.examp(1000)

# Example 2. Another teaching demonstration. By the function put.points.demo()
put.points.demo()

# Example 3. Application to the Golub (1999) data. We shall illustrate the...
library(multtest); data(golub)
mcm3 <- grep("MCM3",golub.gnames[,2])
mcm3
golub.gnames[mcm3,2]
x <- golub[mcm3[1],]; y <- golub[mcm3[2],]
cor(x,y)
plot(x,y)
cor.test(x,y)

# Example 4. Confidence interval by the bootstrap. Another method to con-...
nboot <- 1000; boot.cor <- matrix(0,nrow=nboot,ncol = 1)
data <- matrix(c(x,y),ncol=2,byrow=FALSE)

for (i in 1:nboot) {
    dat.star <- data[sample(1:nrow(data),replace=TRUE),]
    boot.cor[i,] <- cor(dat.star)[2,1]}

mean(boot.cor)
quantile(boot.cor[,1],c(0.025,0.975))

# Example 5. Application to the Golub (1999) data. The ALL and AML...
library(multtest); data(golub)

corgol <- apply(golub, 1, function(x) cor(x,golub.cl))
o <- order(corgol)
golub.gnames[o[3041:3051],2]

########################################
# 7.4 Principal Components Analysis
########################################

V <- matrix(c(2,1,2.8,2.6),nrow=2, byrow=TRUE)
plot(V,xlim = c(0,3),ylim = c(0,3), col="black", pch=19, cex=1.5, xlab="Variable 1", ylab="Variable 2")
arrows(0,0,2,1, col="red", lwd=2)
arrows(0,0,2.8,2.6, col="blue", lwd=3)

Z <- matrix(c( 1.63, 1.22, -0.40, 0.79, 0.93, 0.97, -1.38, -1.08, -0.17, -0.96, -0.61, -0.93), nrow=6, byrow=TRUE)
K <- eigen(cor(Z))
X<-matrix(c(1.45, 1.09, -0.44, 0.7, 0.8, 0.86, -1.35, -0.99, -0.23, -0.88, -0.64, -0.85), nrow=6, byrow=TRUE)
for(j in 1:2) Z[,j] = (X[,j] - mean(X[,j]))/(0.912871*sd(X[,j]))
print(Z,digits=3)
pca <- princomp(Z, center = TRUE, cor=TRUE, scores=TRUE)
pca$scores
Z %*% K$vec
C <- array(dim=c(6,2))
C[,1] <- pca$scores[,1] * sin(pi/4)
C[,2] <- pca$scores[,1] * cos(pi/4)
plot(Z,xlim = c(-2,2),ylim = c(-2,2), col="blue", pch=19, cex=1.5, xlab="Variable 1", ylab="Variable 2")
arrows(-2,-2,2,2, col="green", lwd=3)
segments(Z[1,1],Z[1,2],C[1,1],C[1,2], col="red", lwd=2)
segments(Z[2,1],Z[2,2],C[2,1],C[2,2], col="red", lwd=2)
segments(Z[3,1],Z[3,2],C[3,1],C[3,2], col="red", lwd=2)
segments(Z[4,1],Z[4,2],C[4,1],C[4,2], col="red", lwd=2)
segments(Z[5,1],Z[5,2],C[5,1],C[5,2], col="red", lwd=2)
segments(Z[6,1],Z[6,2],C[6,1],C[6,2], col="red", lwd=2)

# Example 1. Using R on the above data. It is convenient to store the data of...
V1 <- c( 1.63, -0.40, 0.93, -1.38, -0.17, -0.61)
V2 <- c( 1.22, 0.79, 0.97, -1.08, -0.96, -0.93)
cor.test(V1,V2)

plot(V1,V2, xlim=c(-2,2), ylim=c(-2,2), col="red", pch=15)

Z <- matrix(c( V1, V2), nrow=6, byrow=FALSE)
cor(Z)

K <- eigen(cor(Z))
print(K,digits=2)

Z %*% K$vec[,1]

print(Z %*% K$vec, digits=2)

pca <- princomp(Z, center = TRUE, cor=TRUE, scores=TRUE)
pca$scores

# Example 2. Application to the Golub (1999) data. The first five eigenvalues...
eigen(cor(golub))$values[1:5]

data <- golub; p <- ncol(data); n <- nrow(data) ; nboot<-1000
eigenvalues <- array(dim=c(nboot,p))

for (i in 1:nboot) {
    dat.star <- data[sample(1:n,replace=TRUE),]
    eigenvalues[i,] <- eigen(cor(dat.star))$values
}

for (j in 1:p) {
    print(quantile(eigenvalues[,j],c(0.025,0.975)))
}

for (j in 1:5) {
    cat(j,as.numeric(quantile(eigenvalues[,j],c(0.025,0.975))),"\n" )
}

sum(eigen(cor(golub))$values[1:2])/38*100

pca <- princomp(golub, center = TRUE, cor=TRUE, scores=TRUE)
o <- order(pca$scores[,2])
golub.gnames[o[1:10],2]
golub.gnames[o[3041:3051],2]

# Example 3. Biplot. A useful manner to plot both genes (cases) and patients...
stats::biplot(princomp(data,cor=TRUE),cex=0.5,expand=0.8, xlab="Component 1", ylab="Component 2")

# Example 4. Critical for S-phase. Golub et al. (1999) mention that among...
data(golub, package = "multtest")
golubFactor <- factor(golub.cl)
o1 <- grep("CD",golub.gnames[,2])
o2 <- grep("Op",golub.gnames[,2])
o3 <- grep("MCM",golub.gnames[,2])
o <- c(o1,o2,o3)
length(o)

pt <- apply(golub, 1, function(x) t.test(x ~ golubFactor)$p.value)
oo <- o[pt[o]<0.01]

Z <- as.matrix(scale(golub, center = TRUE, scale = TRUE))
K <- eigen(cor(Z))
P <- Z %*% -K$vec[,1:2]
leu <- data.frame(P[oo,], row.names= oo)

plot(leu,xlim=c(-10,15), ylim=c(-10,10), pch=19, cex=1.2, xlab="Principal Component 1", ylab="Principal Component 2", col="darkgreen")
text(x = leu$X1, y=leu$X2, labels=rownames(leu), pos = 1, col="blue")
fac <- as.integer(oo %in% o1) + 2 * as.integer(oo %in% o2) + 3 * as.integer(oo %in% o3)
text(x = leu$X1, y=leu$X2, labels=fac, pos = 3, col="red")

cl <- hclust(dist(leu,method="euclidian"),method="single")
plot(cl,
     lwd=3,
     col="blue",
     col.axis = "brown",
     ylab="Distance",
     xlab="Clustering of the expression of genes",
     hang=-1,
     main=NA,
     sub=NA,
     axes=FALSE)
axis(side = 2, at = seq(0, 5, 1), col = "brown",labels = TRUE, lwd = 4)

a <- as.integer(rownames(leu)[cl$order])
for (i in 1:length(a)) {
    cat(a[i],golub.gnames[a[i],2],"\n")
}


