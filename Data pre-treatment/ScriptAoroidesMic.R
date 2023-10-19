# R packages needed
# "phyloseq" to merge ASV_table, ASV_counttable, SampleData, Phylotree 
# "knitr" for table creation
# "xtable" Functions converting an R object to an "xtable" object, which can then be printed as a LaTeX or HTML table  
# "phangorn" for phylogenetic trees
# "decipher" performs alignement of multiple sequences following a guided tree
# "microbiome" series of OTU/ASV operations and analysis
# "eulerr" for venn's diagrams
# "ggplot2" graphics management
# "vegan" murtivariate statistics
# "ampvis2" visualise and analyse 16S rRNA amplicon data in different ways
# "microbiomeutilities" a supporting tool for extending the functionality of both phyloseq and microbiome R packages

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

#***** Alpha Diversity*****

set.seed(12345)
#Physeq object used for aplpha diversity has 
plot_richness(physeq16s1Ori, title ="Agelas oroides Alpha Diversity", 
              measures=c("observed","Shannon", "PD", "Chao1")) + 
              geom_point(size=2)+ geom_boxplot()

#Creation of a dataframe storing alpha diversity messures
Div_tabAO <- estimate_richness(physeq16s1Ori, measures=c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
alpha_div_AO <- estimate_pd(physeq16s1Ori)
Div_tabAO <- cbind(sample_data(physeq16s1Ori), alpha_div_AO, Div_tabAO) 
write.table(Div_tabAO, "Div_tab_AO_MS_FN.tsv", sep = "\t", quote=F, col.names=NA)

#Normal distribution test of alpha div messures
shapiro.test(Div_tabAO$PD)
"data:  Div_tabAO$PD
W = 0.60427, p-value = 0.002162"
shapiro.test(Div_tabAO$Chao1)
"data:  Div_tabAO$Chao1
W = 0.9813, p-value = 0.49"    #Chao1 and Observed have normal distribution 
shapiro.test(Div_tabAO$Shannon)
"data:  Div_tabAO$Shannon
W = 0.81979, p-value = 0.0001538"
shapiro.test(Div_tabAO$Observed)
"data:  Div_tabAO$Observed
W = 0.96821, p-value = 0.4915"
shapiro.test(Div_tabAO$Simpson)
"data:  Div_tabAO$Simpson
W = 0.56354, p-value = 0.00000002695"

Anova_AOFraction_ADivChao1 <- aov(Chao1 ~ Fraction, Div_tabAO) 
summary(Anova_AOFraction_ADivChao1)
          "  Df Sum Sq Mean Sq F value Pr(>F)  
Fraction     2  252.3  126.15   2.623  0.091 .
Residuals   27 1298.5   48.09  "

Anova_AOSpoTiss_ADivChao1 <- aov(Chao1 ~ SpoTissue, Div_tabAO) 
summary(Anova_AOSpoTiss_ADivChao1)

"            Df Sum Sq Mean Sq F value Pr(>F)
SpoTissue    2  224.8  112.40   2.289  0.121
Residuals   27 1326.0   49.11           "

#***** Beta Diversity*****
#Sponge body part (endosome, ectosome and general fragment) exemple

# Add metadata to the previously created diversity table
Div_tabAO
Physeq16s1.metaAO <- meta(Div_tabAO)

# Create a list with the factors to compare
pares.spps <- combn(seq_along(SpoT_Levels), 2, simplify = TRUE, FUN = function(i)SpoT_Levels[i])
print(pares.spps)#int(endosome); ext (ectosome); gnl (general)

# divergences calcul 
referenceBetaDiv3 <- apply(abundances(physeq16FIN), 1, median)

div_AOint3 <- divergence(subset_samples(physeq16FIN, SpoTissue == "Int"), referenceBetaDiv3, method = "bray")
div_AOext3 <- divergence(subset_samples(physeq16FIN, SpoTissue == "Ext"), referenceBetaDiv3, method = "bray")
div_AOgnl3 <- divergence(subset_samples(physeq16FIN, SpoTissue == "Gnl"), referenceBetaDiv3, method = "bray")

# transform to dataframe
data.frame(div_AOint3) -> df_div_AOint3
data.frame(div_AOext3) -> df_div_AOext3
data.frame(div_AOgnl3) -> df_div_AOgnl3
mutate(df_div_AOint3, SpoTissue = "Int") -> df_div_AOint3
mutate(df_div_AOext3, SpoTissue = "Ext") -> df_div_AOext3
mutate(df_div_AOgnl3, SpoTissue = "Gnl") -> df_div_AOgnl3
# Add cols to the created dataframe
colnames(df_div_AOint3) <- c("divergence", "SpoTissue")
colnames(df_div_AOext3) <- c("divergence", "SpoTissue")
colnames(df_div_AOgnl3) <- c("divergence", "SpoTissue")
# Graphic with statistical test
p3_3 <- ggviolin(title= "Agelas oroides microbiome divergences Physeq_MS", data = div.boxplot3,
                 x = "SpoTissue", y = "divergence", fill = "SpoTissue",
                 add = "boxplot", palette = c("#fdbf6f","#b2df8a","#a6cee3"))
p3_3 + stat_compare_means(comparisons = pares.SpoT)
p3_3 + stat_compare_means()#Kruskal Wallis p = 0.15

#ordination methods
#exemple Sponge body part
physeq.mds.unifrac <- ordinate(physeq16FIN, method = "MDS", distance = "unifrac")
physeq.mds.unifrac
evals <- physeq.mds.unifrac$values$Eigenvalues
pUnif02b <- plot_ordination(physeq16FIN, physeq.mds.unifrac, color = "SpoTissue", title = "Unifrac distance ordination") +
  labs(col = "A oriodes Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

physeq.mds.Weigtunifrac <- ordinate(physeq16FIN, method = "MDS", distance = "wunifrac")
evals <- physeq.mds.Weigtunifrac$values$Eigenvalues
pWeighUnif02b <- plot_ordination(physeq16FIN, physeq.mds.Weigtunifrac, color = "SpoTissue", title = "Weighted-Unifrac distance ordination") +
  labs(col = "A oriodes Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

physeq.mds.dpcoa <- ordinate(physeq16FIN, method = "MDS", distance = "dpcoa")
evals <- physeq.mds.dpcoa$values$Eigenvalues
pdpcoa_b <- plot_ordination(physeq16FIN, physeq.mds.dpcoa, color = "SpoTissue", title = "DPCoA ordination") +
  labs(col = "Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

physeq.mds.bray02 <- ordinate(physeq16FIN, method = "MDS", distance = "bray")
evals <- physeq.mds.bray02$values$Eigenvalues
pmds02 <- plot_ordination(physeq16FIN, physeq.mds.bray02, color = "SpoTissue", title = "Bray distance ordination") +
  labs(col = "A oriodes Sponge Tissue") +
  coord_fixed(sqrt(evals[2] / evals[1]))

#Adonis statistical test
AdoSpoTPhysq_MS<-adonis2(Esp_bray2 ~ SpoTissue, data = Esp_Braydf2, permutations = 999)
AdoSpoTPhysq_MS
"Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = Esp_bray2 ~ SpoTissue, data = Esp_Braydf2, permutations = 999)
          Df SumOfSqs      R2      F Pr(>F)  
SpoTissue  2  0.10461 0.14791 2.3435  0.038 *
Residual  27  0.60264 0.85209                
Total     29  0.70725 1.00000"

#Abundance Barplots
BPPhylaAgelas_MS <- plot_bar(physeq16FIN_AbRel, fill="Phylum", title = "Agelas oroides Prok Phyla AbRel") +
  scale_fill_manual(values = getPalette(colorCount)) +
  guides(fill=guide_legend(ncol=2))+ 
  geom_bar(stat = "identity", position = "stack") +
  theme(axis.text.x=element_blank())
BPPhylaAgelas_MS +  facet_wrap (~ SpoTissue, scales = "free")
BPPhylaAgelas_MS +  facet_wrap (~ Fraction, scales = "free")
BPPhylaAgelas_MS +  facet_wrap (~ Type, scales = "free")

#Venn diagrams
#create a list with factors
SpoTissueListVenn <- unique(as.character(meta(Physeq_MS_AbRel)$SpoTissue))
print(SpoTissueListVenn)#[1] "Ext" "Gnl" "Int"
list_core <- c() #  empty object to store information

for (n in SpoTissueListVenn){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(Physeq_MS_AbRel, SpoTissue == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.5)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  print(list_core)
}

#export phyloseq data to .tsv for instance for MicrbiomeAnalyst

#ASV table
OTUAOok = as(otu_table(physeq16FIN), "matrix")
write.table(OTUAOok, "OtuTable_Aook.tsv", sep = "\t", quote=F, col.names=NA)
#Taxonomy table
TAXAOok = as(tax_table(physeq16FIN), "matrix")
write.table(TAXAOok, "TaxTable_Aook.tsv", sep = "\t", quote=F, col.names=NA)
#SampleData table
SampleData_MS = as(sample_data(physeq16FIN), "matrix")
write.table(SampleData_MS, "SampleDaTa_Aook.tsv", sep = "\t", quote=F, col.names=NA)
#Tree
tree_MS<-phy_tree(physeq16FIN)
ape::write.tree(tree_MS, "PhySeq_Aook_tree.nwk")

---------------------------------------------------------------------------------------------
