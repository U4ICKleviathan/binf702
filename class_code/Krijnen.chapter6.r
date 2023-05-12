###############################################################################
###############################################################################
# Chapter 6
###############################################################################
###############################################################################

########################################
# 6.1 Probe data
########################################

# Example 1. We will start with a built-in data set called MLL.B from the...
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("affyPLM")
biocLite("ALLMLL")
biocLite("genefilter")
biocLite("limma")

biocLite("annaffy")
biocLite("hgu95av2.db")
biocLite("GEOquery")
biocLite("GO.db")
biocLite("GO")
biocLite("made4")


library(affy)
library(ALLMLL)
data(MLL.B, package = "ALLMLL")
MLL.B

class(MLL.B)
?MLL.B
MLL.B
str(MLL.B)

slotNames(MLL.B)
dim(exprs(MLL.B))
annotation(MLL.B)
probeNames(MLL.B)[1:10]
geneNames(MLL.B)[1:10]
pm(MLL.B,"200000_s_at")[1:4,1:3]

# par(mfrow=c(2,2))
matplot(pm(MLL.B,"200000_s_at"),type="l", xlab="Probe No.", ylab="PM Probe intensity")
# matplot(pm(MLL.B,"200000_s_at") - mm(MLL.B,"200000_s_at"),type="l", xlab="Probe No.", ylab="PM Probe intensity")

hist(MLL.B)

# MA plots for the 6th sample relative to the pseudo-median reference
MAplot(MLL.B, which=c(6), cex.lab=1.5)

# smooth scatter MA plots for the first 4 samples relative to the pseudo-median reference
par(mfrow=c(2,2))   # split the window into 4 parts to plot to sequentially
MAplot(MLL.B, which=c(1,2,3,4), plot.method= "smoothScatter")

# smooth scatter MA plots for the first 4 samples relative to sample 5
MAplot(MLL.B, which=c(1,2,3,4), ref=5, plot.method= "smoothScatter")

par(mfrow=c(1,1))  # revert the canvas back to just one plot per window
image(MLL.B)       # hit the escape key to exit out

# Visualizing RNA degradation
degrade <- AffyRNAdeg(MLL.B)
plotAffyRNAdeg(degrade, col=1:20)

########################################
# 6.2 Preprocessing methods
########################################
bgcorrect.methods()
pmcorrect.methods()
normalize.methods(MLL.B)
express.summary.stat.methods()

# Example 1. The three pre-processing steps can be employed one after...
eset <- expresso(MLL.B,
                 bgcorrect.method="rma",
                 normalize.method="constant",
                 pmcorrect.method="pmonly",
                 summary.method="avgdiff")

# Example 2. Another frequently applied preprocessing method is RMA....
library(affy)
library(affyPLM)
data(MLL.B, package = "ALLMLL")
eset2 <- rma(MLL.B)

mybreaks = seq(from=0, to=max(exprs(eset2)), length.out=30)
hists = apply(exprs(eset2), 2, function(x) {h = hist(x, breaks=mybreaks, plot=FALSE); h$density})
hists.matrix = matrix(unlist(hists), ncol = 20, byrow = FALSE)
matplot(hists.matrix, type="l")

library(limma)

limma::plotMA(eset2,
              array=6,    # plot array 5
              main=NULL,  # No title
              cex.lab=1.5 # Make axis labels big
              )
MAplot(eset2, which=c(6), cex.lab=1.5)

# Example 3. In the sequel we shall frequently work with the ALL data...
data(ALL, package = "ALL")
table(ALL$BT)
?ALL::ALL
slotNames(ALL)
row.names(exprs(ALL))[1:10]

ALL1pp <- ALL1 <- ALL[,ALL$mol == "ALL1/AF4"]
mads <- apply(exprs(ALL1), 2, mad)
meds <- apply(exprs(ALL1), 2, median)
dat <- sweep(exprs(ALL1), 2, meds)
exprs(ALL1pp) <- sweep(dat, 2, mads, FUN="/")

class(ALL1)
boxplot(exprs(ALL1))
boxplot(exprs(ALL1pp))

# The R package "made4" contains a command that will create a quick overview and automatically generate multiple plots
# in a single command
library(made4)
overview(MLL.B)

########################################
# 6.3 Gene filtering
########################################

# Example 1. Filtering by the coefficient of variation. A manner to filter...
cvval <- apply(exprs(ALL1pp),1,function(x){sd(x)/abs(mean(x))})
sum(cvval < 0.2)
ALL1pp[cvval<0.2,]

# Example 2. Combining several filters. It is often desired to combine...
library("genefilter")
f1 <- function(x)(IQR(x)>0.5)
f2 <- pOverA(.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x))<0.1)
f6 <- function(x) (sqrt(10)* abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
library("ALL"); data(ALL)
selected <- genefilter(exprs(ALL[,ALL$BT=="B"]), ff)
sum(selected)

# Example 3. Filtering by t-test and normality. One may also want to...
library("genefilter");library("ALL"); data(ALL)
patientB <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"))
f1 <- function(x) (shapiro.test(x)$p.value > 0.05)
f2 <- function(x) (t.test(x ~ patientB)$p.value < 0.05)
sel1 <- genefilter(exprs(ALL[,patientB==TRUE]), filterfun(f1))
sel2 <- genefilter(exprs(ALL[,patientB==FALSE]), filterfun(f1))
preSelected <- sel1 & sel2
preSelectedALLs <- ALL[preSelected,]
sel3 <- genefilter(exprs(preSelectedALLs), filterfun(f2))
selectedALLs <- preSelectedALLs[sel3,]
dim(selectedALLs)

library(limma)
sel3 <- genefilter(exprs(ALL), filterfun(f2))
x <- matrix(as.integer(c(sel1,sel2,sel3)),ncol = 3,byrow=FALSE)
colnames(x) <- c("sel1","sel2","sel3")
vc <- vennCounts(x, include="both")
vennDiagram(vc, circle.col=c("blue","red","green"), lwd=3)

library("gplots")
heatmap.2(exprs(ALL), scale = "row", col=greenred(75), dendrogram="both", key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5)


########################################
# 6.4 Applications of linear models
########################################

# Example 1. Analysis of variance. We select patients with B-cell leukemia in...
library("ALL"); library("limma");
data(ALL, package = "ALL")
allB <- ALL[,which(ALL$BT %in% c("B","B1","B2"))]
design.ma <- model.matrix(~ 0 + factor(allB$BT))
colnames(design.ma) <- c("B","B1","B2")
fit <- lmFit(allB, design.ma)
fit <- eBayes(fit)
toptab <- topTable(fit, coef=NULL,5,adjust.method="fdr")
print(toptab[,1:6],digits=4)

cont.ma <- makeContrasts(B-B1,B1-B2, levels=factor(allB$BT))
cont.ma

fit1 <- contrasts.fit(fit, cont.ma)
fit1 <- eBayes(fit1)
toptabcon <- topTable(fit1, coef=NULL,5,adjust.method="fdr")
print(toptabcon[,1:5],digits=4)

# Example 2. Summarizing output in HTML format. It is often desired to...
library("annaffy"); library("hgu95av2.db")
anntable <- aafTableAnn(as.character(row.names(toptabcon)), "hgu95av2.db", aaf.handler())
saveHTML(anntable, "~/biology664/ALLB123.html", title = "B-cell 012 ALL")

# Example 3. Using basic R functions. It is also possible to summarize...
library("multtest"); library("annaffy"); library("hgu95av2.db")
library("ALL"); data(ALL, package = "ALL")
ALLB <- ALL[,which(ALL$BT %in% c("B","B1","B2"))]
panova <- apply(exprs(ALLB), 1, function(x) anova(lm(x ~ ALLB$BT))$Pr[1])
genenames <- featureNames(ALLB)[panova<0.000001]
atab <- aafTableAnn(genenames, "hgu95av2.db", aaf.handler()[c(1:3,8:9,11:13)])
saveHTML(atab, file="~/biology664/ANOVAonB-cellGroups.html")
getwd()

# Example 4. Analyzing public available data. The GDS1365 data con-...
library(GEOquery); library(limma); library(hgu95av2.db); library(annaffy)
gds <- getGEO("GDS1365")
eset <- GDS2eSet(gds,do.log2=TRUE)
prot <- pData(eset)$protocol
time <- pData(eset)$time
pval <- apply(exprs(eset)[1:12625,], 1, function(x) anova(lm(x ~ prot * time))$Pr[1:3])
pvalt <- data.frame(t(pval))
colnames(pvalt) <- c("meffprot","mefftime","interaction")
genenames <- featureNames(eset)[pvalt$meffprot< 0.01 & pvalt$mefftime < 0.01 & pvalt$interaction < 0.01]
atab <- aafTableAnn(genenames,"hgu95av2.db",aaf.handler()[c(1:3,8:9,11:13)])
saveHTML(atab, file="~/biology664/Two-way ANOVA protocol by time.html")

########################################
# 6.5 Searching an annotation package
########################################
library("ALL"); data(ALL)
annotation(ALL)

library(hgu95av2.db)
ls("package:hgu95av2.db")
ChrNrOfProbe <- as.list(hgu95av2CHR)
ChrNrOfProbe[1]
?hgu95av2CHR

get("1389_at", env = hgu95av2ACCNUM)
get("1389_at", env = hgu95av2ENTREZID)
get("1389_at", env = hgu95av2SYMBOL)
get("1389_at", env = hgu95av2GENENAME)
get("1389_at", env = hgu95av2UNIGENE)

library(annotate)
genbank("J03779",disp="browser")
genbank(179833,disp="data",type="uid")
get("1389_at", env = hgu95av2CHRLOC)
get("1389_at", env = hgu95av2MAP)


########################################
# 6.6 Using annotation to search literature
########################################
library(hgu95av2.db); library(annotate); library(ALL); data(ALL)
pmid <- get("1389_at",env=hgu95av2PMID)
pubmed(pmid,disp="browser")
absts <- pm.getabst("1389_at", "hgu95av2")
pm.titles(absts)
ne <- pm.abstGrep("neutral endopeptidase",absts[[1]])
pmAbst2HTML(absts[[1]],filename="~/biology664/pmon1389_at.html")

########################################
# 6.7 Searching GO numbers and evidence
########################################
go1389 <- get("1389_at", env = hgu95av2GO)
idl <- lapply(go1389,function(x) x$GOID)
idl[[1]]

library(annotate)
getOntology(go1389,"BP")
getEvidence(go1389)

go1389TAS <- subset(go1389,getEvidence(go1389)=="TAS")

sapply(go1389TAS,function(x) x$GOID)
sapply(go1389TAS,function(x) x$Evidence)
sapply(go1389TAS,function(x) x$Ontology)

########################################
# 6.8 GO parents and children
########################################

# Example 1. Collecting GO information. There are functions to obtain...
GOMFPARENTS$"GO:0003700"
GOMFCHILDREN$"GO:0003700"
go1389 <- get("1389_at", env = hgu95av2GO)
gonr <- getOntology(go1389, "BP")
gP <- getGOParents(gonr)
gC <- getGOChildren(gonr)
gPC <- c(gonr,gP,gC)
pa <- sapply(gP,function(x) x$Parents)
ch <- sapply(gC,function(x) x$Children)
gonrc <- c(gonr,unlist(pa),unlist(ch))
length(gonrc)

# Example 2. Probe selection by GO. A research strategy may be to start...
library(GO.db); library(annotate); library("ALL"); data(ALL)
go1389 <- get("1389_at", env = hgu95av2GO)
gonr <- getOntology(go1389, "BP")
gP <- getGOParents(gonr)
pa <- sapply(gP,function(x) x$Parents)
probes <- mget(unlist(pa),hgu95av2GO2ALLPROBES)
probeNames <- unlist(probes)
ALLpr <- ALL[probeNames,]
dim(exprs(ALLpr))

########################################
# 6.9 Gene filtering by a biological term
########################################

# Example 1. Filter gene by a term. From a biological point of view...
library("GO"); library("annotate"); library("hgu95av2.db")
GOTerm2Tag <- function(term) {
    GTL <- eapply(GOTERM, function(x) {grep(term, x@Term, value=TRUE)})
    Gl <- sapply(GTL, length)
    names(GTL[Gl>0])
}
GOTerm2Tag("transcriptional repressor")
GOTerm2Tag("repressor")

tran1 <- hgu95av2GO2ALLPROBES$"GO:0003714"
tran2 <- hgu95av2GO2ALLPROBES$"GO:0008231"
tran3 <- hgu95av2GO2ALLPROBES$"GO:0017053"
tran <- c(tran1,tran2,tran3)
inboth <- tran %in% row.names(exprs(ALL))
ALLtran <- ALL[tran[inboth],]

GOTERM$"GO:0017053"
dim(exprs(ALLtran))

########################################
# 6.10 Significance per chromosome
########################################

# Example 1. On the expression values of the ALL data we perform a two...
library("ALL"); data(ALL); library("hgu95av2.db")
rawp <- apply(exprs(ALL), 1, function(x) t.test(x ~ ALL$remission)$p.value)
xx <- as.list(hgu95av2CHR)
AffimIDChr19 <- names(xx[xx=="19"])
names(rawp) <- featureNames(ALL)
f <- matrix(NA,2,2)
f[1,1] <- sum(rawp[AffimIDChr19]<0.05); f[1,2] <- sum(rawp[AffimIDChr19]>0.05)
f[2,1] <- sum(rawp<0.05) - f[1,1] ; f[2,2] <- sum(rawp>0.05) - f[1,2]
print(f)
fisher.test(f)
chisq.test(f)



