
install.packages(c("TeachingDemos"),repo="http://cran.r-project.org", dep=TRUE)

# install.packages("ctv") # do the install only once!
library("ctv") # loads the "ctv" library into the current R session
install.views("Genetics")


source("http://www.bioconductor.org/biocLite.R")

# installing BiocManager as described on  https://www.bioconductor.org/install/
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")


# using BiocManager instead of biocLite()
# BiocManager::install("ALL")
# BiocManager::install("affy")
packageVersion("BiocManager")
library(ALL)
data(ALL)

# BiocManager::install("multtest")
library(multtest)
data(golub)

class(golub)
str(golub)

objects()
ls()

apropos("diff")

example(boxplot)

1:5
seq(0,1,0.1)


?gl
factor <- gl(3,5)


# 1.6
gene1
sum(gene1)
mean(gene1)/3
sd(gene1)
summary(gene1)


# 1.7

gene1 <- c(1.00,1.50,1.25) 
gene2 <- c(1.35,1.55,1.00)
gene3 <- c(-1.10,-1.50,-1.25)
gene4 <- c(-1.20,-1.30,-1.00)


rowColNames <- list(c("gene1", "gene2", "gene3", "gene4"),
                    c("Eric", "Peter", "Anna"))

geneData <- matrix(c(gene1,gene2,gene3,gene4), nrow=4, ncol=3,
                    byrow=TRUE, dimnames = rowColNames)


geneData[1,2]

dir.create("biology664")
getwd()
write.table(geneData,file=paste(getwd(),"/biology664/geneData.Rdata", sep = ""), append = FALSE)
?write.table

# 1.8
apply(geneData,2,mean)
apply(geneData,1,mean)

meanExpressions <- apply(geneData,1,mean)
o <- order(meanExpressions,decreasing=TRUE)
o

geneData[o, ]

geneData[c(1,2), ]

geneData[c("gene1","gene2"), ]



meanExpressions > 0

geneData[meanExpressions > 0, ]

# 1.9
library(multtest); data(golub)

golub.gnames[1042, ]

nrow(golub)
ncol(golub)

golub[1042,2]
golub[1042,1:27]


golubFactor <- factor(golub.cl, levels=0:1, labels = c("ALL","AML"))


golub[1042,golubFactor=="ALL"]

meanALL <- apply(golub[ ,golubFactor=="ALL"], 1, mean)


cd33 = grep("CD33",golub.gnames[ ,2], ignore.case = TRUE)

golub.gnames[cd33,]

# 1.10

patients.df <-
  data.frame(
    # Define the 3 column names when creating the data.frame
    patientID = c("101", "102", "103", "104"),
    treatment = c("drug", "placebo", "drug", "placebo"),
    age = c(20, 30, 24, 22)
  )

patients.df

patients.df <- data.frame( # Define the 3 column names after creating the data.frame
   c("101", "102", "103", "104"),
   c("drug", "placebo", "drug", "placebo"),
   c(20, 30, 24, 22)
   )


colnames(patients.df) <- c("patientID", "treatment", "age")

patients.df


# 1.11

patients.df[["treatment"]]
patients.df$treatment

# 1.12
patients.df[c(1,3),]
patients.df[patients.df$treatment=='drug',]
treatmentIsDrug = patients.df$treatment=='drug'
treatmentIsDrug

patients.df[treatmentIsDrug,]

patients.df["treatment"]
