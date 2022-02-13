# 1. Illustration of mean and standard deviation.
# (a) Compute the mean and the standard deviation for 1, 1.5, 2, 2.5, 3.
vecA <- c(1, 1.5, 2, 2.5, 3.)

meanA <- mean(vecA) 
# 2

stdA <- sd(vecA)
# [1] 0.7905694

# (b) Compute the mean and the standard deviation for 1, 1.5, 2, 2.5, 30.
vecB <- c(1, 1.5, 2, 2.5, 30.)
meanB <- mean(vecB)
#[1] 7.4

stdB <- sd(vecB)
# [1] 12.64615

# (c) Comment on the differences.


############################################################################
############################################################################
############################################################################

# 2. Testing normality of gene expression. Consider the gene expression values 
# in row 790 and 66 of the Golub et al. (1999) data.

# (a) Produce a box plot for the expression values of the ALL patients
# and comment on the differences. Are there outliers?

data(golub, package = "multtest")
golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))


rowIndex <- 790
boxplot(golub[rowIndex, ] ~ golubFactor, # values
           cex.lab=1.5, 
           xlab="Golub Factor",
           ylab=paste("Row", rowIndex),
           col=c("purple","green")
)
rowIndex <- 66
boxplot(golub[rowIndex, ] ~ golubFactor, # values
           cex.lab=1.5, 
           xlab="Golub Factor",
           ylab=paste("Row", rowIndex),
           col=c("purple","green")
)


# (b) Produce a QQ-plot and formulate a hypothesis about the normality of the genes.
rowIndex <- 790
qqnorm(golub[rowIndex, golubFactor=='ALL'], 
       pch=19, # plot solid circles
       cex.lab=1.5, # make axis labels big
       col="red",
       ylab=paste("Row", rowIndex),
       main=NULL);
qqline(golub[rowIndex, golubFactor=='ALL'], col="blue")

rowIndex <- 66
qqnorm(golub[rowIndex, golubFactor=='ALL'], 
       pch=19, # plot solid circles
       cex.lab=1.5, # make axis labels big
       col="red",
       ylab=paste("Row", rowIndex),
       main=NULL);
qqline(golub[rowIndex, golubFactor=='ALL'], col="blue")



# (c) Compute the mean and the median for the expression values of
# the ALL patients and compare these. Do this for both genes
rowIndex <- 790
mean(golub[rowIndex,golubFactor=="ALL"])
median(golub[rowIndex,golubFactor=="ALL"])
rowIndex <- 66
mean(golub[rowIndex,golubFactor=="ALL"])
median(golub[rowIndex,golubFactor=="ALL"])


############################################################################
############################################################################
############################################################################


# 3. Effect size. An important statistic to measure is the effect size which
# is defined for a sample as x/s. It measures the mean relative to the
# standard deviation, so that its value is large when the mean is large
# and the standard deviation small.


# (a) Determine the five genes with the largest effect size of the ALL
# patients from the Golub et al. (1999) data. Comment on their
# size.

effectSize = apply(golub[,golubFactor=="ALL"], 1, function(x) mean(x)/sd(x))
orderedEffectSize



# (b) Invent a robust variant of the effect size and use it to answer the
# previous question.


