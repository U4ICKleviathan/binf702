###############################################################################
###############################################################################
# Chapter 9
###############################################################################
###############################################################################

########################################
# 9.1 Using a query language
########################################
## ## UNCOMMENT TO INSTALL
## source("http://bioconductor.org/biocLite.R")
## biocLite("Biostrings")
## biocLite("BSgenome.Celegans.UCSC.ce2")
## biocLite("seqLogo")
## biocLite("annotate")
## biocLite("TFBSTools")
## biocLite("JASPAR2014")
## biocLite("muscle")

## install.packages(c("seqinr"),repo="http://cran.r-project.org",dep=TRUE)

library(seqinr)
choosebank()
choosebank("genbank", timeout=30)
query("ccnd3","k=ccnd3",virtual=TRUE)$nelem
query("ccnd3.human","sp=homo sapiens AND k=ccnd3",virtual=TRUE)$nelem

########################################
# 9.2 Getting information on downloaded sequences
########################################

# Example 1. Let's download sequences related to the species homo sapiens...
choosebank("genbank", timeout=30)
ccnd3.human = query("ccnd3.human","sp=homo sapiens AND k=ccnd3@")
ccnd3.human$nelem

sapply(ccnd3.human$req, getKeyword)
sapply(ccnd3.human$req, getName)

sapply(ccnd3.human$req, getLength)
getSequence(ccnd3.human$req[[1]])[1:15]
getTrans(ccnd3.human$req[[1]])[1:15]
getAnnot(ccnd3.human$req[[1]])

########################################
# 9.3 Computations on sequences
########################################

# Example 1. Frequencies of (di)nucleotides. We shall continue with...
mono = table(getSequence(ccnd3.human$req[[1]]))
mono
mono = table(getSequence(ccnd3.human$req[[1]]),1)
mono
dinucs = count(getSequence(ccnd3.human$req[[1]]),2)
dinucs
decreasing = order(dinucs)
dotchart(dinucs[decreasing], xlab="Dinucleotide Counts")

# Example 2. G + C percentage. We are often interested in the fraction...
GC(getSequence(ccnd3.human$req[[1]]))
GC1(getSequence(ccnd3.human$req[[1]])) # Looking at only codon position 1
GC2(getSequence(ccnd3.human$req[[1]])) # Looking at only codon position 2
GC3(getSequence(ccnd3.human$req[[1]])) # Looking at only codon position 3

require(zoo)
seq = getSequence(ccnd3.human$req[[1]])
# find the 50bp-windowed running GC content fraction for the coding region of the gene
gcContent = rollapply(seq, width=50, by=1, FUN=GC, align="left")
plot(gcContent,         # (x,y) data
     type="l",       # line
     xlab="Nucleotide Index",
     ylab="GC Percentage",
     col="purple",
     lwd=2,          # make line thicker
     cex.lab=1.5)    # make axis labels big
## dev.copy2eps(device=x11,file="GCperc.eps")

# Example 3. Rho and z-scores. The coefficient rho and the corresponding...
# rho = f(XY) / (f(x)*f(Y)) ; X,Y element of {A,C,G,T}
rhos = round(rho(getSequence(ccnd3.human$req[[1]])),2)
decreasing = order(rhos)
dotchart(rhos[decreasing], xlab="Dinucleotide Rhos")
rhos
# CpG is the most under-represented dinucleotide given that it's a very GC-rich gene

zScores = round(zscore(getSequence(ccnd3.human$req[[1]]),modele='base'),2)
decreasing = order(zScores)
dotchart(zScores[decreasing], xlab="Dinucleotide Z-scores")
zScores
# CpG is the most under-represented dinucleotide given that it's a very GC-rich gene

# Example 4. Comparing Amino acid frequencies. We continue with the...
# AF517525.CCND3
table <- table(getTrans(ccnd3.human$req[[1]]))
orderedTable <- table[order(table)]
names(orderedTable) <- aaa(names(orderedTable))
dotchart(orderedTable,pch=19,xlab="Stop and amino-acid-counts", col="red", main=NA)
abline(v=1,lty=2)
## dev.copy2eps(device=x11,file="translationdotchartccnd3hs.eps")

# BC011616.CCND3
table <- table(getTrans(ccnd3.human$req[[2]]))
orderedTable <- table[order(table)]
names(orderedTable) <- aaa(names(orderedTable))
dotchart(orderedTable,pch=19,xlab="Stop and amino-acid-counts", col="red", main=NA)
abline(v=1,lty=2)
## dev.copy2eps(device=x11,file="translationdotchartccnd3hsseq2.eps")

# Example 5. Isoelectric point. The function computePI computes the...

# computePI computes the theoretical isoelectric point of a protein, which is the pH at which the protein
# has a neutral charge
computePI(getTrans(ccnd3.human$req[[1]]))

# pmw computes the protein molecular weight
pmw(getTrans(getSequence(ccnd3.human$req[[1]])))


# Example 6. Hydropathy score. The coefficients ff1, * * * , ff20 are available...

# The hydropathy index of an amino acid is a number representing the hydrophobic or hydrophilic properties of its
# sidechain. It was proposed in 1982 by Jack Kyte and Russell F. Doolittle.
# The larger the number is, the more hydrophobic the amino acid. The most hydrophobic amino acids are isoleucine (4.5)
# and valine (4.2). The most hydrophilic ones are arginine (-4.5) and lysine (-3.9). This is very important in protein
# structure; hydrophobic amino acids tend to be internal (with regard to the protein's 3 dimensional shape)
# while hydrophilic amino acids are more commonly found towards the protein surface.
ccnd3.human.DNAseqs <- sapply(ccnd3.human$req, getSequence)
ccnd3.human.AAseqs <- sapply(ccnd3.human.DNAseqs, getTrans)
data(EXP)
names(EXP$KD) <- sapply(words(), function(x) translate(s2c(x)))
kdc <- EXP$KD[unique(names(EXP$KD))]
# kdc <- -kdc[order(names(kdc))]
kdc <- kdc[order(kdc)]
dotchart(kdc)

getValues = function(seq, coefs) {
  values = numeric()
  for (i in 1:length(seq)) {
    values[i] = coefs[seq[i]]
    }
  return(values)
}

# find the running hydropathy scores for the protein
protein = getTrans(ccnd3.human$req[[1]])
hydro = getValues(protein, kdc)
plot(1:length(protein),hydro,type="l")

# find the 20aa-windowed running hydropathy scores for the protein
mean20 = rollapply(hydro, width=20, by=1, FUN=mean, align="left")
plot(1:length(mean20),mean20,type="l")

# Which sections of the protein are likely buried in the core of the protein, and which sections are likely accessible on the surface?

kdAth = sapply(ccnd3.human.AAseqs, function(x) {mean(getValues(x, kdc))})
print(kdAth,digits=3)

# The data set aaindex of the seqinr library contains more than five hundred sets of coefficients for
# computing specific biophysical/chemical quantities with respect to amino acids.
data(aaindex)
aaindex$CASG920101$D # description of CASG920101 data
aaindex$CASG920101$I # indexed CASG920101 data

#####################################################################################
# 9.4 Matching patterns; Codon, DNA-binding, and RNA-binding motif finding; PWMs
#####################################################################################

# Example 1. Pattern match. In the sequence with NCBI accession number...
library(seqinr)
library(Biostrings)
choosebank("genbank", timeout=30)
ccnd3.human = query("ccnd3.human","sp=homo sapiens AND k=ccnd3@")
ccnd3.human.DNAseqs <- sapply(ccnd3.human$req, getSequence)
ccnd3.human.DNAseq1 <- c2s(ccnd3.human.DNAseqs[[1]])
ccnd3.human.DNAseq1
subseq <- "cccggg"

# For as single-stranded motif like a codon motif or an RNA-binding motif only search the positive strand
countPattern(subseq, ccnd3.human.DNAseq1, max.mismatch = 0)
matchPattern(subseq, ccnd3.human.DNAseq1, max.mismatch = 0)
matchPattern(subseq, ccnd3.human.DNAseq1, max.mismatch = 1)

library(BSgenome.Celegans.UCSC.ce2)

chrII <- Celegans[["chrII"]]
dinucleotideFrequency(chrII)
pattern <- DNAString("TGGGTGTATGTA")

############################################################################
## Virtually transforming DNA with bisulfate and searching for a motif before and after
############################################################################

# Bisulfite sequencing is the use of bisulphite treatment of DNA to determine its pattern of methylation.
# DNA methylation was the first discovered epigenetic mark, and remains the most studied. In animals it
# predominantly involves the addition of a methyl group to the carbon-5 position of cytosine residues of
# the dinucleotide CpG, and is implicated in repression of transcriptional activity.
# Treatment of DNA with bisulphite converts cytosine residues to uracil, but leaves 5-methylcytosine residues
# unaffected.

replace = function(fromSeq, toSeq, dnaString) {
  seq = toString(dnaString)
  newSeq = gsub(fromSeq, toSeq, seq)
  return(DNAString(newSeq))
}

## Transforming CpG to TpG and searching the + strand
plus_strand_bisulfite <- replace("CG", "TG", chrII)
dinucleotideFrequency(plus_strand_bisulfite)
matchPattern(pattern, chrII)
matchPattern(pattern, plus_strand_bisulfite)

## For a dual-strand motif like a dsDNA-binding motif search the positive strand and the negative strand
## However, instead of searching the negative strand, search the positive strand with the revComp(motif),
## which is much faster and more efficient

## Transforming CpG to CpA and searching the + strand with the revComp(motif)
plus_strand2_bisulfite <- replace("CG", "CA", chrII)
dinucleotideFrequency(plus_strand2_bisulfite)
matchPattern(reverseComplement(pattern), chrII)
matchPattern(reverseComplement(pattern), plus_strand2_bisulfite)

# Searching for gapped motifs
# biocLite("BSgenome.Dmelanogaster.UCSC.dm3")
# biocLite("BSgenome.Ecoli.NCBI.20080805")
library(BSgenome.Dmelanogaster.UCSC.dm3)
subject <- Dmelanogaster$chr3R

# Searching for a large gapped motif
Lpattern <- DNAString("AGCTCCGAG")
Rpattern <- DNAString("TTGTTCACA")
matchLRPatterns(Lpattern, Rpattern, 500, subject) # 1 match on the positive strand
matchLRPatterns(reverseComplement(Rpattern), reverseComplement(Lpattern), 500, subject) # 0 matches on the negative strand

# Searching for a smaller gapped motif
Lpattern <- DNAString("AATGCC")
Rpattern <- DNAString("TGGATC")
matchLRPatterns(Lpattern, Rpattern, 10, subject) # 26 matches on the positive strand
matchLRPatterns(reverseComplement(Rpattern), reverseComplement(Lpattern), 10, subject) # 19 matches on the negative strand

# Searching for a pattern in a genome containing Ns:
matchPattern("TGGGTGTCTTT", chrII) # no match
matchPattern("TGGGTGTCTTT", chrII, fixed=FALSE) # 1 match

## Using wildcards ("N") in the pattern on a genome containing Ns:
library(BSgenome.Dmelanogaster.UCSC.dm3)
matchPattern("TTTATGNTTGGTA", Dmelanogaster$chrX, fixed=FALSE) # Problem!

# Use masking to not include occurrences of Ns in the genome in our search
chrX <- maskMotif(Dmelanogaster$chrX, "N")
as(chrX, "Views") # 4 non-masked regions
matchPattern("TTTATGNTTGGTA", chrX, fixed=FALSE) # Not a Problem anymore!
## Can also be achieved with no mask by setting fixed="subject"
masks(chrX) <- NULL
matchPattern("TTTATGNTTGGTA", chrX, fixed="subject")

## Allowing indels in the search
subject <- Celegans$chrI
pattern1 <- DNAString("ACGGACCTAATGTTATC")
## Allowing up to 2 mismatching letters doesnt give any match:
matchPattern(pattern1, subject, max.mismatch=2)
## But allowing up to 2 edit operations gives 3 matches:
matchPattern(pattern1, subject, max.mismatch=2, with.indels=TRUE)

# Search all the worm promoter regions for the TATA Box
Tatabox <- DNAString("TATAAA")
subject <- Celegans$upstream5000
mindex <- vmatchPattern(Tatabox, subject, fixed=FALSE)
start_index <- startIndex(mindex)
hist(unlist(start_index))

# Search all the fly promoter regions for the TATA Box
Tatabox <- DNAString("TATAAA")
subject <- Dmelanogaster$upstream5000
mindex <- vmatchPattern(Tatabox, subject, fixed=FALSE)
start_index <- startIndex(mindex)
hist(unlist(start_index))

# Search all the worm promoter regions for the E-Box
Ebox <- DNAString("CANNTG")
subject <- Celegans$upstream5000
mindex <- vmatchPattern(Ebox, subject, fixed=FALSE)
start_index <- startIndex(mindex)
hist(unlist(start_index))

# Search all the fly promoter regions for the E-box
Ebox <- DNAString("CANNTG")
subject <- Dmelanogaster$upstream5000
mindex <- vmatchPattern(Ebox, subject, fixed=FALSE)
start_index <- startIndex(mindex)
hist(unlist(start_index))

##############################################################
# PWMs - Position Weight Matrices
##############################################################
library("seqLogo")

# HNF4A = Hepatocyte Nuclear Factor 4 Alpha
data(HNF4alpha)
HNF4alpha  # DNAStringSet containing 71 aligned HNF4alpha binding sites

chr3R <- Dmelanogaster$chr3R
chr3R

## Create a PWM from a PFM or directly from a rectangular DNAStringSet object
pcm <- consensusMatrix(HNF4alpha) # position counts mattrix
pcm
pwm <- PWM(pcm) # creates a log-odds score position weight matrix, same as PWM(HNF4alpha)
pwm

# make the logo
pfm = consensusMatrix(HNF4alpha, as.prob=TRUE) # position frequency matrix
pfm
seqLogo::seqLogo(makePWM(pfm[1:4,]))

round(pwm, 2)
maxWeights(pwm)
maxScore(pwm)

pwmRC = reverseComplement(pwm) # reverse complement PWM for searching for the motif on the negative strand
pwmRC

## Score the first 5 start positions
PWMscoreStartingAt(pwm, unmasked(chr3R), starting.at=1:7)

## Search the positive strand
hits <- matchPWM(pwm, chr3R, min.score=.85, with.score=TRUE)
nhits <- countPWM(pwm, chr3R, min.score=.85) # same as length(hits)
nhits
hits

## Use with.score=TRUE to get the scores of the hits:
head(mcols(hits)$score)
min(mcols(hits)$score / maxScore(pwm)) # should be >= min.score

## The scores can also easily be post-calculated if you didn't include "with.score=TRUE" when you called matchPWM():
scores <- PWMscoreStartingAt(pwm, subject(hits), start(hits))

## Search the negative strand for the dsDNA-binding motif as well!!
matchPWM(pwmRC, chr3R, min.score=.85, with.score=TRUE)


# Get the Arnt PWM from the JASPAR database
library(JASPAR2014)
library(TFBSTools)
opts = list()
opts[["name"]] = "Arnt"
pfMatrixList = getMatrixSet(JASPAR2014, opts)
pcmObject = pfMatrixList[[1]] # position count matrix
pcmObject
pcm = as.matrix(pcmObject)

# create the frequency matrix
sums = colSums(pcm)
pfm = pcm * (1/sums) # position frequency matrix
pfm
seqLogo::seqLogo(makePWM(pfm))

pwm = Biostrings::PWM(pcm) # log-odds position weight matrix
maxScore(pwm)
pwm

pwmRC = reverseComplement(pwm) # reverse complement PWM for searching for the motif on the negative strand
pwmRC

hits <- matchPWM(pwm, chr3R, min.score=.85, with.score=TRUE) # Search the positive positive strand for the motif
hits

## Search the negative strand for the dsDNA-binding motif as well!!
matchPWM(pwmRC, chr3R, min.score=.85, with.score=TRUE)

########################################
# 9.5 Pairwise alignments
########################################

# Example 1. Basic recursion. The idea of recursion is to generate a sequence...
x<-double(); x[1] <- 1
for (i in 2:10) {
  x[i] <- 2*x[i-1]-10
}
x[10]

########################################
# Global Alignment with Needleman-Wunsch
########################################

# Example 2. Dynamic programming of DNA sequences. Consider again...
library(seqinr)
x <- s2c("GAATTC"); y <- s2c("GATTA"); d <- 2

# create nucleotide substitution matrix
s <- matrix(data=NA,nrow=length(y),ncol=length(x))
for (i in 1:(nrow(s))) {
  for (j in 1:(ncol(s))) {
    if (y[i]==x[j]) {
      s[i,j]<- 2 # match
    }
    else {
      s[i,j]<- -1 # mismatch
    }
  }
}
rownames(s) <- c(y); colnames(s) <- c(x)
s

needleman.wunsch = function(x, y, sub, d) {
  # create the dynamic programming matrices F and G
  F <- matrix(data=NA,nrow=(length(y)+1),ncol=(length(x)+1))
  G <- matrix(data=NA,nrow=(length(y)+1),ncol=(length(x)+1))

  rownames(F) <- c("",y); colnames(F) <- c("",x)
  rownames(G) <- c("",y); colnames(G) <- c("",x)

  F[,1] <- -seq(0,length(y)*d,d); F[1,] <- -seq(0,length(x)*d,d)
  G[,1] <- rep("u",length(y)+1); G[1,] <- rep("l",length(x)+1)
  G[1,1] = " ";

  for (i in 2:(nrow(F)))
    for (j in 2:(ncol(F)))
    {
        diagonal = F[i-1,j-1]+sub[i-1,j-1]
        up = F[i-1,j]-d
        left = F[i,j-1]-d
        max = max(diagonal, up, left)

        F[i,j] <- max

        if (diagonal == max) {G[i,j] = "d"}
        else if (up == max) {G[i,j] = "u"}
        else {G[i,j] = "l"}
    }

  # Now perform traceback and alignment
  i=nrow(G); j=ncol(G);
  max = max(length(x)+length(y)); ii=max; A1=vector(,max); A2=vector(,max);

  repeat{
      if (G[i,j] == "d") {G[i,j] = "D"; A1[ii] = x[j-1]; A2[ii] = y[i-1]; i=i-1; j=j-1}
      else if (G[i,j] == "u") {G[i,j] = "U"; A1[ii] = "-"; A2[ii] = y[i-1]; i=i-1}
      else {G[i,j] = "L"; A1[ii] = x[j-1]; A2[ii] = "-"; j=j-1}
      ii=ii-1
      if ((i==1) & (j==1)){
          break
      }
  }
  return(list(ScoreMatrix=F, TracebackMatrix=G, Alignment=rbind(A1,A2)[,(ii+1):max], Score=F[nrow(F),ncol(F)]))
}

dpMatrix = needleman.wunsch(x,y,s,d)
dpMatrix

# Example 3. Programming Needleman-Wunsch. For the two sequences...
file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50"
BLOSUM50 <- as.matrix(read.table(file, check.names=FALSE))
x <- s2c("HEAGAWGHEE"); y <- s2c("PAWHEAE"); s <- BLOSUM50[y,x]; d <- 8
dpMatrix = needleman.wunsch(x,y,s,d)
dpMatrix

# Example 4. Global Optimal alignment with Needleman-Wunsch. We may also conveniently use the pairwiseAlignment
# from the Biostrings package which uses backtracing to retrieve the optimal global alignment
globalAlign = pairwiseAlignment(AAString("HEAGAWGHEE"), AAString("PAWHEAE"),
                      substitutionMatrix = "BLOSUM50",gapOpening = 0, gapExtension = -8,
                      scoreOnly = FALSE, type="global")
globalAlign

# Different Alignment with the PAM250 Substitution matrix
globalAlign2 = pairwiseAlignment(AAString("HEAGAWGHEE"), AAString("PAWHEAE"),
                                substitutionMatrix = "PAM250",gapOpening = 0, gapExtension = -8,
                                scoreOnly = FALSE, type="global")
globalAlign2

# Example 5. Comparing with random sequences. To illustrate how the...
allRandomScores <- double()
for (i in 1:1000) {
    x <- c2s(sample(rownames(BLOSUM50),7, replace=TRUE))
    y <- c2s(sample(rownames(BLOSUM50),10, replace=TRUE))
    allRandomScores[i] <- pairwiseAlignment(AAString(x), AAString(y),
                                         substitutionMatrix = "BLOSUM50",gapOpening = 0, gapExtension = -8,
                                         scoreOnly = TRUE, type="global")
}
sum(allRandomScores>1)/1000

# Example 6. Sliding window on Needleman-Wunsch scores. We may also...
choosebank("genbank", timeout=30); library(seqinr)
ccnd3.human = query("ccnd3.human","sp=homo sapiens AND k=ccnd3@")
ccnd3.human.DNAseqs <- sapply(ccnd3.human$req, getSequence)
ccnd3.human.AAseqs <- sapply(ccnd3.human.DNAseqs, getTrans)
x <- c2s(ccnd3.human.AAseqs[[1]])
y <- c2s(ccnd3.human.AAseqs[[1]][50:70])
nwscore <- double() ; n <- length(ccnd3.human.AAseqs[[1]])

for (i in 1:(n-21)) {
    nwscore[i] <- pairwiseAlignment(AAString(c2s(ccnd3.human.AAseqs[[1]][i:(i+20)])),
                                    AAString(y),substitutionMatrix = "BLOSUM50",gapOpening = 0,
                                    gapExtension = -8, scoreOnly = TRUE, type="global")
}

pairwiseAlignment(AAString(y), AAString(y), scoreOnly = TRUE)
max(nwscore)
which.max(nwscore)

########################################
# Local Alignment with Smith-Waterman
########################################
smith.waterman = function(x, y, sub, d) {
  # create the dynamic programming matrices F and G
  F <- matrix(data=NA,nrow=(length(y)+1),ncol=(length(x)+1))
  G <- matrix(data=NA,nrow=(length(y)+1),ncol=(length(x)+1))

  rownames(F) <- c("",y); colnames(F) <- c("",x)
  rownames(G) <- c("",y); colnames(G) <- c("",x)

  F[,1] <- rep(0,length(y)+1); F[1,] <- rep(0,length(x)+1)
  G[,1] <- rep(" ",length(y)+1); G[1,] <- rep(" ",length(x)+1)

  for (i in 2:(nrow(F)))
    for (j in 2:(ncol(F)))
    {
        diagonal = F[i-1,j-1]+sub[i-1,j-1]
        up = F[i-1,j]-d
        left = F[i,j-1]-d
        max = max(0, diagonal, up, left)

        F[i,j] <- max

        if (max == 0) {G[i,j] = " "}
        else if (diagonal == max) {G[i,j] = "d"}
        else if (up == max) {G[i,j] = "u"}
        else {G[i,j] = "l"}
    }

  # Now perform traceback and alignment
  maxF = which(F == max(F), arr.ind = TRUE); i=maxF[1]; j=maxF[2];
  maxL = max(length(x)+length(y)); ii=maxL; A1=vector(,maxL); A2=vector(,maxL);
  repeat{
      if (G[i,j] == "d") {G[i,j] = "D"; A1[ii] = x[j-1]; A2[ii] = y[i-1]; i=i-1; j=j-1}
      else if (G[i,j] == "u") {G[i,j] = "U"; A1[ii] = "-"; A2[ii] = y[i-1]; i=i-1}
      else {G[i,j] = "L"; A1[ii] = x[j-1]; A2[ii] = "-"; j=j-1}
      ii=ii-1
      if (F[i,j] == 0){
          break
      }
  }
  return(list(ScoreMatrix=F, TracebackMatrix=G, Alignment=rbind(A1,A2)[,(ii+1):maxL], Score=F[maxF]))
}

x <- s2c("HEAGAWGHEE"); y <- s2c("PAWHEAE"); s <- BLOSUM50[y,x]; d <- 8
dpMatrix = smith.waterman(x,y,s,d)
dpMatrix

library(Biostrings)
localAlign = pairwiseAlignment(AAString("HEAGAWGHEE"), AAString("PAWHEAE"),
                                substitutionMatrix = "BLOSUM50",gapOpening = 0, gapExtension = -8,
                                scoreOnly = FALSE, type="local")
localAlign

# Overlap (Ends-Free) Alignment with modified Needleman-Wunsch
overlapAlign = pairwiseAlignment(AAString("HEAGAWGHEE"), AAString("PAWHEAE"),
                               substitutionMatrix = "BLOSUM50",gapOpening = 0, gapExtension = -8,
                               scoreOnly = FALSE, type="overlap")
overlapAlign


########################################
# Local Alignment with BLAST
########################################

library(annotate)

## Run NCBI BLAST on an entrez gene ID
blastSequences(17702, program="blastn")

## or x can be a sequence
## hitListSize does not promise that you will get the number of matches you
## want.. It will just try to get that many.
blastSequences(x = "GGCCTTCATTTACCCAAAATG", program="blastn", hitListSize="20")

########################################
# Global Alignment with MUSCLE
########################################
library(muscle)

# align DNA sequences of MAX orthologs from 31 mammalian species
aln <- muscle::muscle(umax)
aln

# align sequences in a fasta file
x <- c("RPLWVAPDGHIFLEAFSPVYK")
y <- c("PLWINPIDGRIILEAFSPLAE")
z <- c("RPLWVAPDGHIILEAFSPVYK")
aaStringSet = AAStringSet(c(x, y, z))
names(aaStringSet) = c("Q38861", "O13768", "O00835")
aaStringSet

aln <- muscle(aaStringSet)
aln


