############################################################################### 
############################################################################### 
############################################################################### 
############################################################################### 
# 1. Some questions to orientate yourself.
# (a) Use the function class() to find the class to which the following objects belong: golub, golub[1,1],golub.cl, golub.gnames,
# apply(), exp(), golubFactor, plot(), ALL.

library(multtest)
data(golub)
library(ALL)

class(golub)
class(golub[1, 1])
class(golub.cl)
class(golub.gnames)
class(apply)
class(exp)
class(golubFactor)
class(plot)
class(ALL)

# > class(golub)
# [1] "matrix" "array" 
# > class(golub[1, 1])
# [1] "numeric"
# > class(golub.cl)
# [1] "numeric"
# > class(golub.gnames)
# [1] "matrix" "array" 
# > class(apply)
# [1] "function"
# > class(exp)
# [1] "function"
# > class(golubFactor)
# [1] "factor"
# > class(plot)
# [1] "standardGeneric"
# attr(,"package")
# [1] "methods"
# > class(ALL)
# [1] "ExpressionSet"
# attr(,"package")
# [1] "Biobase"

############################################################################### 
############################################################################### 

# (b) Define the meaning of the following abbreviations: rm(), sum(),

# 1) rm() - remove object from memory within an R session
# 2) sum() - computes the sum.
# 3) prod() - produces the product of all elements within a vector
# 4) seq() - produces a sequence of numbers from a lower bound to an upper bound
#            given a certain step size
# 5) sd() - produces the standard deviation of a vector
# 6) nrow() - returns the number of rows in a dataframe or a matrix

############################################################################### 
############################################################################### 


# (c) For what purpose are the following functions useful: grep(),
# agrep(), apply(), gl(), library(), source(), setwd(), history(),
# str().

# 1) grep() Finds the index within an object where a certain condition is met.
#    This can help filter an object into a subset to use in later calculations

?agrep
# 2) agrep() This uses the approximate match of a string within an object.
#    this is useful for discovering and exploring a dataset, or managing data 
#    that has user input or otherwise dirty data

# 3) apply() This is useful to perform a certain operation across a matrix or 
#    a dataframe. Using apply, you only have to use 1 line of code instead of 
#    looping over a structure


# 4)  gl() This is useful for plotting, so I am told, but I can't think of how
#     to use this at the moment

# 5) library() This allows us to import code, methods, or data that others 
#    have produced. Why re-invent the wheel when there are stable, correct, 
#    and useful code that others have produced?

# 6)  source() - allows us to specify the source of a module or other libraries.
#     it also allows us to use code that we have written locally

# 7) setwd() - This allows us to set the working directory for the R environment
#    We can then easily write to files or load data from our current directory.

# 8) history() - This allows us to view the recently executed commands and
#    statements

?str()
# 9) str() - this allows us to briefly view an object and review it 


############################################################################### 
############################################################################### 
############################################################################### 
############################################################################### 

# 2. Standard deviations. Consider the data in the matrix geneData
# that is constructed in Section 1.7. Its small size enables you to check
# your computations even with a pocket calculator. 14

gene1 <- c(1.00,1.50,1.25) 
gene2 <- c(1.35,1.55,1.00)
gene3 <- c(-1.10,-1.50,-1.25)
gene4 <- c(-1.20,-1.30,-1.00)


rowColNames <- list(c("gene1", "gene2", "gene3", "gene4"),
                    c("Eric", "Peter", "Anna"))

geneData <- matrix(c(gene1,gene2,gene3,gene4), nrow=4, ncol=3,
                   byrow=TRUE, dimnames = rowColNames)

# (a) Use apply() to compute the standard deviation of the persons.

apply(geneData,2,sd)
# Eric    Peter     Anna 
# 1.350540 1.690845 1.307032 


# (b) Use apply() to compute the standard deviation of the genes.
apply(geneData,1,sd)
# gene1     gene2     gene3     gene4 
# 0.2500000 0.2783882 0.2020726 0.1527525 


# (c) Order the matrix according to the gene standard deviations.
stdVec <- apply(geneData, 1, sd)
ordered <- order(stdVec, decreasing=FALSE)
geneData[ordered,]
#        Eric Peter  Anna
# gene4 -1.20 -1.30 -1.00
# gene3 -1.10 -1.50 -1.25
# gene1  1.00  1.50  1.25
# gene2  1.35  1.55  1.00

# (d) Which gene has the largest standard deviation?
# gene2

############################################################################### 
############################################################################### 
############################################################################### 
############################################################################### 

# 3. Computations of gene expression means in the Golub data.
# (a) Use apply() to compute the mean gene expression value.


library(multtest)
data(golub)

geneMean <- apply(golub, 1, mean)

# (b) Order the data matrix according to the gene means.

o <- order(geneMean, decreasing = TRUE)

orderedGolub <- golub[o,]



# (c) Give the biological names of the three genes with the largest mean
# expression value.

golub.gnames[o[1:3],]

# [,1]   [,2]                                                                                 [,3]              
# [1,] "5729" "37 kD laminin receptor precursor/p40 ribosome associated protein gene"              "U43901_rna1_s_at"
# [2,] "1707" "RPS14 gene (ribosomal protein S14) extracted from Human ribosomal protein S14 gene" "M13934_cds2_at"  
# [3,] "7097" "GAPD Glyceraldehyde-3-phosphate dehydrogenase"                                      "X01677_f_at"     



