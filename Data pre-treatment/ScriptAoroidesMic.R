# R packages needed
# "phyloseq" to merge ASV_table, ASV_counttable, SampleData, Phylotree 
# "xtable" 
# "phangorn" for phylogenetic trees
# "decipher" performs alignement of multiple sequences following a guided tree
# "microbiome"
# "eulerr"
# "ggplot2"
# "vegan" 

## Phyloseq objec creation:

# ASV_Table - Matrix
ASV_16S1_Table <- read.delim("ASV_Table_16S1_tophyseq.txt", header = TRUE, sep = "\t", row.names=1) 
ASV_16S1_Table<-as.matrix(ASV_16S1_Table)

# Taxonomy table - Matrix
ASV_Tax_16S1 <- read.delim("Taxa_Table_161S1_tophyseq.txt", header = TRUE, sep = "\t", row.names = 1)
ASV_Tax_16S1<-as.matrix(ASV_Tax_16S1)

# Upload the Mapping file
Mapfile=read.delim("MappingFileC.txt" , header=TRUE, sep="\t", row.names=NULL)
rownames(sampledata)=sampledata$SampleID

# Creataion des objets OTU et TAX
ASV = otu_table(ASV_16S1_Table, taxa_are_rows = TRUE)
TAX = tax_table(ASV_Tax_16S1)
sampledata = sample_data(Mapfile)

#Creation object type Phyloseq
physeq16s1Ori = phyloseq(ASV, TAX, sampledata)

# Creation od phylogenetic tree

seqs <- getSequences(TAX)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

physeq16s1=merge_phyloseq(physeq16s1Ori, phy_tree(fitGTR$tree))
"> physeq16s1
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 446 taxa and 78 samples ]
sample_data() Sample Data:       [ 78 samples by 6 sample variables ]
tax_table()   Taxonomy Table:    [ 446 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 446 tips and 444 internal nodes ]"



