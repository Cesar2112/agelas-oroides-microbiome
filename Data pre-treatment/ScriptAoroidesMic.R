# R packages needed
# "phyloseq" to merge ASV_table, ASV_counttable, SampleData, Phylotree 
# "knitr" for table creation  
# "phangorn" for phylogenetic trees
# "decipher" performs alignement of multiple sequences following a guided tree
# "microbiome" series of OTU/ASV operations and analysis
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

#Creation object type Phyloseq
physeq16s1Ori = phyloseq(ASV, TAX, sampledata, phy_tree(fitGTR$tree))
"
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 446 taxa and 30 samples ]
sample_data() Sample Data:       [ 78 samples by 6 sample variables ]
tax_table()   Taxonomy Table:    [ 446 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 446 tips and 444 internal nodes ]"

#Prune taxa not present in any sample (if they exist)
physeq16s1NoZeroTaxa <- prune_taxa(taxa_sums(physeq16s1) > 0, physeq16s1)
physeq16s1<-physeq16s1NoZeroTaxa
"phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 365 taxa and 30 samples ]"
#Create a table with reads number per sample:
nbReads<- sample_sums(physeq16s1)
write.table(nbReads, "nbReads_161S1_MS.tsv", sep = "\t", quote=F, col.names=NA)


#Taxa filtering by prevalecence
#Prevalecence Agelas
prevdf_AO = apply(X = otu_table(physeq16s1),
                  MARGIN = ifelse(taxa_are_rows(physeq16s1), yes = 1, no = 2),
                  FUN = function(x){sum(x > 0)})
prevdf_AO = data.frame(Prevalence = prevdf_AO,
                       TotalAbundance = taxa_sums(physeq16s1),
                       tax_table(physeq16s1))
plyr::ddply(prevdf_AO, "Phylum", function(df1){cbind(mean(df1$Prevalence),
                                                     sum(df1$Prevalence))}) -> dfprev_AO
kable(dfprev_AO)
write.table(dfprev_AO, "PhylumReadsPrev_161S1_AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
filterPhyla = c("Bacteroidetes", "Campilobacterota", "Chlamydiae", "Deinococcus-Thermus", "Rhodothermaeota") 
(physeq16FIN = subset_taxa(physeq16s1, !Phylum %in% filterPhyla))
# physeq16FIN [ 351 taxa and 30 samples ]

#Change Taxa names par ASV#
taxa_names(Physeq_MS)
taxa_names(Physeq_MS) <- paste0("ASV", seq(ntaxa(physeq16FIN)))

#Overview of taxonomic assignement
TaxOverviwAO <-table(tax_table(physeq16FIN)[, 2])
write.table(TaxOverviwAO, "TaxOverview_AO_FN.tsv", sep = "\t", quote=F, col.names=NA)

# Create table, number of features for each phyla and class
PhylaTable_AO<-table(tax_table(physeq16FIN)[, "Phylum"], exclude = NULL)
write.table(PhylaTable_AO, "PhylaTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
ClassTable_AO<-table(tax_table(physeq16FIN)[, "Class"], exclude = NULL)
write.table(ClassTable_AO, "ClassTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
OrderTable_AO<-table(tax_table(physeq16FIN)[, "Order"], exclude = NULL)
write.table(OrderTable_AO, "OrderTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
FamTable_AO<-table(tax_table(physeq16FIN)[, "Family"], exclude = NULL)
write.table(FamTable_AO, "FamilyTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)
GenTable_AO<-table(tax_table(physeq16FIN)[, "Genus"], exclude = NULL)
write.table(GenTable_AO, "GenusTable16s1AO_FN.tsv", sep = "\t", quote=F, col.names=NA)

##Find out the number of a taxon, for example:
Ao <- as.vector(tax_table(physeq16FIN)[, "Phylum"])
b_Ao <- unique(Ao)#10 Phyla
Ao <- as.vector(tax_table(physeq16FIN)[, "Class"])
b_Ao <- unique(Ao)#18 Class
Ao <- as.vector(tax_table(physeq16FIN)[, "Order"])
b_Ao <- unique(Ao)#31 Order

# convert to relative abundance (for abundance barplots)
physeq16FIN_AbRel = transform_sample_counts(physeq16FIN, function(x) x / sum(x) )

               **************** DIVERSITE ALPHA*********************
set.seed(12345)



