data(golub, package = "multtest")

golubFactor <- factor(golub.cl,levels=0:1, labels= c("ALL","AML"))

ccnd3 = grep("CCND3",golub.gnames[ ,2], ignore.case = TRUE)
boxplot(golub[ccnd3,] ~ golubFactor, # values
        cex.lab=1.5, # make axis labels big
        main=NULL, # no title
        xlab="Leukemia subtype",
        ylab="CCND3 (Cyclin D3) Expression",
        col=c("purple","green"))

