#reset ctrl+shift+f10
rm(list=ls())
library(ape) # Analysis of phylogenetics and evolution
library(hierfstat) # Hierarchical F-statistics
library(corrplot) # Visualization of correlation matrix
sarco.dna<-read.dna(file="sarco_seq.fasta", format = "fasta")
sarco.dna
class(sarco.dna)
annot <- read.csv("sarcoData.csv", header=TRUE, row.names=1) #read in the annotations

## Multiple sequence alignment

library(ips) #MAFFT is available here
# Requires path to MAFFT binary - set it according to your installation
sarco.mafft <- mafft(x=sarco.dna, method="localpair", maxiterate=100, options="--adjustdirection", exec="C:/Users/Mycolab-2017/Desktop/mafft-win/mafft")
sarco.mafft
# Delete all columns containing at least 25% of gaps
sarco.mafft.ng <- deleteGaps(x=sarco.mafft,gap.max=nrow(sarco.mafft)/4)
# Delete every line (sample) containing at least 20% of missing data
sarco.mafft.ng <- del.rowgapsonly(x=sarco.mafft.ng, threshold=0.2, freq.only=FALSE)
# Delete every alignment position having at least 20% of missing data
sarco.mafft.ng <- del.colgapsonly(x=sarco.mafft.ng, threshold=0.2, freq.only=FALSE)
sarco.mafft.ng
class(sarco.mafft.ng)
image.DNAbin(sarco.mafft.ng) # Plot the alignment
# Check the alignment
checkAlignment(x=sarco.mafft.ng, check.gaps=TRUE, plot=TRUE, what=1:4)
library(adegenet)

# Checking SNPs
# Position of polymorphism within alignment - snpposi.plot() requires input data in form of matrix
snpposi.plot(x=as.matrix(sarco.mafft.ng), codon=FALSE)
# Position of polymorphism within alignment - differentiating codons
snpposi.plot(as.matrix(sarco.mafft.ng))
# When converting to genind object, only polymorphic loci are kept - threshold for polymorphism can be arbitrary (polyThres=...)
sarco.genind <- DNAbin2genind(x=sarco.mafft.ng, polyThres=0.01)
sarco.genind # See it
# Check sequences
# Nucleotide diversity
pegas::nuc.div(x=sarco.mafft.ng)
# Base frequencies
ape::base.freq(x=sarco.mafft.ng)
# GC content
ape::GC.content(x=sarco.mafft.ng)
# Number of times any dimer/trimer/etc oligomers occur in a sequence
seqinr::count(seq=as.character.DNAbin(sarco.dna[["MZ227236"]]), wordsize=3)
#Distance based methods
## Model selection

library(phangorn)
# Conversion to phyDat for phangorn
sarco.phydat <- as.phyDat(sarco.mafft.ng) # Prepare starting tree
modelTest(object=as.phyDat(sarco.mafft.ng), tree=nj(dist.dna(x=sarco.mafft.ng,model="raw")))
# Create the distance matrix
sarco.dist <- dist.dna(x=sarco.mafft.ng, model="JC")
# Check the resulting distance matrix
sarco.dist
class(sarco.dist)
dim(as.matrix(sarco.dist))

#Claster Dendrogram

library(ade4) #Analysis of ecological data, multivariate methods
table.paint(df=as.data.frame(as.matrix(sarco.dist)), cleg=0, clabel.row=0.5, clabel.col=0.5)
# Same visualization, colored
# heatmap() reorders values
# because by default it plots
# also dendrograms on the edges
heatmap(x=as.matrix(sarco.dist), Rowv=NA, Colv=NA, symm=TRUE)
#This is very basic function to make dendrogram
par(mfrow=c(1,1))
plot(hclust(d=sarco.dist, method="complete")) #hierarchical clustering

#UPGMA

# Saving as phylo object (and not hclust) gives more possibilities for further plotting and manipulations
sarco.upgma <- as.phylo(hclust(d=sarco.dist, method="average"))
plot.phylo(x=sarco.upgma, cex=0.75)
title("UPGMA tree") #looks ok, but the branch length is questionable
#Test quality - tests correlation of original distance in the matrix and reconstructed distance from hclust object
plot(x=as.vector(sarco.dist), y=as.vector(as.dist(
cophenetic(sarco.upgma))), xlab="Original pairwise distances",
ylab="Pairwise distances on the tree", main="Is UPGMA appropriate?", pch=20, col=transp(col="black",
                                     alpha=0.1), cex=2)
abline(lm(as.vector(as.dist(cophenetic(sarco.upgma)))~as.vector(sarco.dist)), col="red")
cor.test(x=as.vector(sarco.dist), y=as.vector(as.dist(cophenetic(sarco.upgma))), alternative="two.sided") # Testing the correlation

#Neighbor-Joining tree

sarco.nj <- nj(sarco.dist)
class(sarco.nj)
# Plot a basic tree
plot.phylo(x=sarco.nj, type="phylogram")
sarco.tree <- nj(X=sarco.dist)
plot.phylo(x=sarco.tree, type="unrooted")
title("Unrooted NJ tree")
#Rooting and unrooting trees
plot.phylo(sarco.nj)
print.phylo(sarco.nj)
sarco.nj.rooted <- root.phylo(phy=sarco.nj, interactive=TRUE)
plot.phylo(sarco.nj.rooted)
# Test tree quality - plot original vs. reconstructed distance
plot(x=as.vector(sarco.dist), y=as.vector(as.dist(
  cophenetic(sarco.nj))), xlab="Original pairwise distances",
  ylab="Pairwise distances on the tree", main="Is simple NJ appropriate?", pch=20, col=transp(col="black",
                                                                                          alpha=0.1), cex=2)
abline(lm(as.vector(as.dist(cophenetic(sarco.nj)))~as.vector(sarco.dist)), col="red")
cor.test(x=as.vector(sarco.dist), y=as.vector(as.dist(cophenetic(sarco.nj))), alternative="two.sided") # Testing the correlation
# Linear model for above graph
summary(lm(as.vector(sarco.dist) ~
               as.vector(as.dist(cophenetic(sarco.nj))))) # Prints summary text
# Calculate bootstrap
sarco.tree1 <- bionj(X=sarco.dist)
fit = pml(sarco.tree1, data=sarco.phydat)
fit
methods(class="pml")
fitJC  <- optim.pml(fit, TRUE)
logLik(fitJC)
bs = bootstrap.pml(fitJC, bs=1000, optNni=TRUE,
                   control = pml.control(trace = 0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")
title("bioNJ tree + bootstrap values")
#Maximum parsimony

parsimony(sarco.upgma, sarco.phydat)
parsimony(sarco.nj.rooted, sarco.phydat)
#returns the parsimony score, that is the number of changes which are at least necessary to describe the data for a given tree
#The tree rearrangement implemented are nearest-neighbor interchanges (NNI) and subtree pruning and regrafting (SPR).
treePars  <- optim.parsimony(sarco.upgma, sarco.phydat) #performs tree rearrangements to find trees with a lower parsimony score
treeRatchet  <- pratchet(sarco.phydat, trace = 0) #parsimony ratchet (Nixon 1999) implemented
parsimony(c(treePars, treeRatchet), sarco.phydat)
treeRatchet  <- acctran(treeRatchet, sarco.phydat) #assign branch length to the tree. The branch length are proportional to the number of substitutions / site.
plotBS(midpoint(treeRatchet), type="phylogram")

#Maximum parsimony 2

tre.ini <- nj(dist.dna(sarco.mafft.ng,model="raw")) #reconstruct a tree
tre.ini
parsimony(tre.ini, sarco.phydat) #measure this tree’s parsimony
sarco.pars <- optim.parsimony(tre.ini, sarco.phydat) #the most parsimonious possible tree
sarco.pars
parsimony(sarco.pars, sarco.phydat)
library(ape)
plot(sarco.pars, type="unr", edge.width=2)
title("Maximum-parsimony tree")

#Visualize a tree
library(ggtree)
library(TDbook)
p <- ggtree(sarco.upgma) + geom_tiplab(size=3) #upgma
msaplot(p, sarco.mafft.ng, offset=3, width=15)
p <- ggtree(sarco.nj.rooted) + geom_tiplab(size=3) #nj
msaplot(p, sarco.mafft.ng, offset=3, width=15)
p <- ggtree(treeRatchet) + geom_tiplab(size=3) #Maximum parsimony
msaplot(p, sarco.mafft.ng, offset=3, width=1)
plot.new()
polygon(x=c(0,1,1,0),y=c(0,0,0.5,0.5), col="#FFDD00")
polygon(x=c(0,1,1,0),y=c(0.5,0.5,1,1), col="#0057B7")
text(0.5, 0.25, labels="Слава Україні", cex=2)
