## Pipeline for phylophankton analysis based ont SILVA annotation of 18S rRNA
## R v4.0.2
## library
setwd("../")
library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(picante)
library(dplyr)
library(multcomp)
library(ggalluvial)
library(pheatmap)
library(VennDiagram)
library(ape)
library(ggrepel)
library(FactoMineR)
library(corrplot)
library(caret)
library(indicspecies)
library(patchwork)
library(dismo)
library(gbm)
library(psych)
library(randomForest)
library(linkET)
library(geosphere)
library(ggtree)
library(iNEXT)

### Loading Parameters  & Processing
source("Functions/Loading.R")
species.type <- "zooplankton"
cbbPalette <- Loading_Color()
data <- Loading_Table("Input/table_f_even_depth_w_tax.xls")
group <- Loading_Group("Input/group.txt")
group$Group<- factor(group$Group,
                     levels = c("Spring", "Autumn"))
tree <- Loading_Tree("Input/unrooted_tree.tre")


result <- Data_Process.zoo(data,group,tree)
bac1 <- result[[1]]
bac2 <- result[[2]]
bac3 <- result[[3]]
tree <- result[[4]]
Group_numb <- result[[5]]
Sample_numb <- result[[6]]
wid <- result[[7]]
zoo <- result[[8]]
if (file.exists("Input/Env.txt")) {
    ENV <- read.table("Input/Env.txt",header = TRUE,sep = "\t")
    rownames(ENV) <- ENV[,1]
    ENV <- ENV[,-1]
    ENV_numb <- ncol(ENV)
    ENV <- ENV[group$variable,]
}else{
    ENV <- c()
}
if (file.exists("Input/Geo.txt")) {
    GEO <- read.table("Input/Geo.txt",header = TRUE,sep = "\t",row.names = 1)
}else{
    GEO <- c()
}

## Analysis results
dir.create("Results")
dir.create("Results/Zooplankton")

### Taxonomy statistics
dir.create("Results/Zooplankton/01.ASV.tax")

bac3 <- bac3[,-1]
source("Functions/tax.stat.R")
result <- tax.stat(data[,c(group$variable,"taxonomy")],bac3)    
write.table(result[[1]],
            "Results/Zooplankton/01.ASV.tax/Tax.stat.txt",
            row.names = FALSE,sep = "\t")
pdf("Results/Zooplankton/01.ASV.tax/Tax.anno.per.pdf",width = 4,height = 4)
result[[2]]
dev.off()

### Alpha diversity analysis
dir.create("Results/Zooplankton/02.Alpha")

### Alpha diversity indices calculation
dir.create("Results/Zooplankton/02.Alpha/01.Estimators")
source("Functions/alpha_diversity_indices.R")
otu_Flattening = as.data.frame(rrarefy(bac1, 5102))
bac1 <- otu_Flattening
alpha.bac <- Alpha_diversity_index(bac1,tree)
alpha <- alpha.bac[,-c(3,5,9)]
colnames(alpha) <- c("Observed_species","Chao1","ACE","Shannon","Simpson","Pielou_J",
                     "Pd_faith")
alpha$ACE[is.na(alpha$ACE)] <- alpha$Chao1[is.na(alpha$ACE)]
alpha1 <- data.frame(ID = rownames(alpha),alpha)
bac3 <- read.csv("Results/Zooplankton/01.ASV.tax/bac3.csv",check.names = F,row.names = 1)
write.table(alpha1,"Results/Zooplankton/02.Alpha/01.Estimators/alpha_diversity_indices.txt",
            sep = "\t",row.names = FALSE)

### Phylogentic tree
if (nrow(bac3) > 100) {
    sum <- rowSums(bac3[,1:(ncol(bac3)-4)])
    sum <- data.frame(OTU = rownames(bac3),Sum = sum)
    sum <- sum[order(sum$Sum,decreasing = TRUE),]
    sum <- sum[1:100,]
    bac4 <- bac3[rownames(bac3) %in% sum$OTU,]
    bac5 <- bac1[,colnames(bac1) %in% sum$OTU]
}else{
    bac4 <- bac3
    bac5 <- bac1
}
groupInfo <- split(rownames(bac4),bac4$Phylum)
tree$tip.label <- gsub("(^')|('$)","",tree$tip.label)
tree <- prune.sample(bac5,tree)
tree <- groupOTU(tree, groupInfo)
p <- ggtree(tree,layout = "circular",col="deepskyblue4", branch.length="none") + 
    geom_tiplab(size=3, color="black",hjust = -0.1) + 
    geom_tippoint(size = 3,aes(color = group))
pdf("Results/Zooplankton/01.ASV.tax/Phylogenetic_tree.pdf",width = 8,height = 8)
p
dev.off()

### Differences of alpha diversity indices among different groups
#### Boxplot
dir.create("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_boxplot")
source("Functions/Diff_alpha.R")
result <- diff.alpha.boxplot(alpha1,group)
for (i in 2:ncol(alpha1)) {
    pdf(file = paste("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_boxplot/",
                     colnames(alpha1)[i],".pdf",sep = ""),
        width = 0.8 + 0.8*Group_numb,height = 3.6)
    print(result[[i-1]])
    dev.off()
}
write.table(result[[8]],"Results/Zooplankton/02.Alpha/01.Estimators/div_test.txt",
            row.names = FALSE,sep = "\t")

#### Jitter
dir.create("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_jitter")
source("Functions/Diff_alpha2.R")
result <- diff.alpha.jitter(alpha1,group)
for (i in 2:ncol(alpha1)) {
    pdf(file = paste("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_jitter/",
                     colnames(alpha1)[i],".pdf",sep = ""),
        width = 0.8 + 0.8*Group_numb,height = 3.6)
    print(result[[i-1]])
    dev.off()
}

#### violin
dir.create("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_violin")
source("Functions/Diff_alpha3.R")
result <- diff.alpha.violin(alpha1,group)
for (i in 2:ncol(alpha1)) {
    pdf(file = paste("Results/Zooplankton/02.Alpha/01.Estimators/alpha_div_violin/",
                     colnames(alpha1)[i],".pdf",sep = ""),
        width = 0.8 + 0.8*Group_numb,height = 3.6)
    print(result[[i-1]])
    dev.off()
}

#### Rarefaction curve
dir.create("Results/Zooplankton/02.Alpha/03.Rarefaction")
source("Functions/rare.spec.R")
result <- rare.spec(bac1,group)
pdf("Results/Zooplankton/02.Alpha/03.Rarefaction/rarefaction.pdf",height = 4,width = 5)
result[[1]]
dev.off()
pdf("Results/Zooplankton/02.Alpha/03.Rarefaction/speccum.pdf",height = 4,width = 5)
result[[2]]
dev.off()

## Community analysis
dir.create("Results/Zooplankton/03.Community")
tax.name <- c("Phylum","Genus")

### abundance table
dir.create("Results/Zooplankton/03.Community/01.tax_summary")
source("Functions/tax.summary.R")
result <- tax.summary(bac3)
for (i in 1:2) {
    write.table(result[[i]],
                paste("Results/Zooplankton/03.Community/01.tax_summary/",
                      tax.name[i],".xls",sep = ""),row.names = FALSE,sep = "\t")
}
data.1 <- as.data.frame(t(bac3[,1:(ncol(bac3)-4)]))
aa <- match(rownames(data.1),group$variable)
data.1 <- data.1[aa,]
data.1$Group <- group$Group
data.group <- aggregate(data.1[,1:(ncol(data.1)-1)],
                        list(data.1$Group),mean)
rownames(data.group) <- data.group$Group.1
data.group <- t(data.group[,-1])
data.group <- as.data.frame(cbind(data.group,bac3[,(ncol(bac3)-3):ncol(bac3)]))
result <- tax.summary(data.group)
for (i in 1:2) {
    write.table(result[[i]],
                paste("Results/Zooplankton/03.Community/01.tax_summary/",
                      tax.name[i],"_Group.xls",sep = ""),row.names = FALSE,sep = "\t")
}

### Abundance chart
dir.create("Results/Zooplankton/03.Community/02.Abundance_chart")
source("Functions/Community.s.R")
result <- community.s(bac3,group,species.type)
pdf(file = "Results/Zooplankton/03.Community/02.Abundance_chart/Phylum.pdf",
    width = 0.3*Sample_numb + 1 + max(str_length(bac3$Phylum))*0.12,height = 5)
result
dev.off()
result <- community.s(data.group,group,species.type)
pdf(file = "Results/Zooplankton/03.Community/02.Abundance_chart/Phylum_Group.pdf",
    width = wid*Group_numb + 1 + max(str_length(bac3$Phylum))*0.12,height = 5)
result
dev.off()

source("Functions/Bubblechart.s.R")
result <- bubblechart.s(bac3,group,species.type)
pdf(file = "Results/Zooplankton/03.Community/02.Abundance_chart/Genus.pdf",
    width = 0.22*Sample_numb + 0.5 + 
        max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1,
    height = ifelse(length(unique(bac3$Genus)) > 30, 8.5, 
                    0.25*length(unique(bac3$Genus)) + 1))
result
dev.off()

source("Functions/Bubblechart.g.R")
result <- bubblechart.g(data.group,group,species.type)
pdf(file = "Results/Zooplankton/03.Community/02.Abundance_chart/Genus_Group.pdf",
    width = 0.3*Group_numb + 0.5 + 
        max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1,
    height = ifelse(length(unique(bac3$Genus)) > 30, 8.5, 
                    0.25*length(unique(bac3$Genus)) + 1))
result
dev.off()

source("Functions/heatmap.s.R")
result <- heatmap.s(bac3,group)
pdf(file = "Results/Zooplankton/03.Community/02.Abundance_chart/Genus_heatmap_cluster.pdf",
     width = 0.22*Sample_numb + 1.5 + 
        max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1,
     height = ifelse(length(unique(bac3$Genus)) > 30, 10.5, 
                     0.25*length(unique(bac3$Genus)) + 3))
result
dev.off()

### Venn analysis
dir.create("Results/Zooplankton/02.Alpha/02.Venn")
#### Venn diagram
if (length(levels(group$Group)) <= 5) {
    aa <- levels(group$Group)
    species <- list()
    for (i in 1:length(aa)) {
        dd <- bac2[,group[group$Group == aa[i],1]]
        dd <- dd[which(rowSums(dd) > 0),]
        species[[aa[i]]] <- rownames(dd)
    }
    
    venn.diagram(species,
                 filename = "Results/Zooplankton/02.Alpha/02.Venn/venn.png",
                 height = 5400,width = 5400,
                 resolution = 600,imagetype = "png",units = "px",
                 lwd = 2,lty = 1,fill = cbbPalette[1:length(aa)],cex = 1.5,
                 cat.cex = 2,alpha = 0.8,margin = 0.05,fontface = 2,
                 cat.fontface = 2, print.mode = c("raw","percent"))
    
    source("Functions/overLapper.new.r")
    OLlist <- overLapper(setlist=species, sep="", type="vennsets",keepdups=FALSE)
    counts <- sapply(OLlist$Venn_List, length)
    pdf("Results/Zooplankton/02.Alpha/02.Venn/venn.pdf",
        width=9,height=9,pointsize=16)
    vennPlot(counts=counts,mymain="")
    dev.off()
    
    #### Shared species list
    inter <- get.venn.partitions(species)
    data.tax <- bac3[,(ncol(bac3)-3):ncol(bac3)]
    shared.species <- data.tax[unlist(inter$..values..[1]),]
    shared.species$ID <- rownames(shared.species)
    core_species_id_g <- shared.species[,c("ID",colnames(shared.species)[1:4])]
    write.table(core_species_id_g,
                "Results/Zooplankton/02.Alpha/02.Venn/Shared_species.txt",
                sep = "\t",row.names = FALSE)
}else{
    source("Functions/Flower.g.R")
    pdf(file = "Results/Zooplankton/02.Alpha/02.Venn/Venn.pdf",
        width = 8,height = 8)
    flower.g(bac2,group)
    dev.off()
    
    source("Functions/Core.species.id.g.R")
    core_species_id_g <- core.species.id.g(bac2,group)
    core_species_id_g <- bac3[core_species_id_g,(ncol(bac3)-3):ncol(bac3)]
    core_species_id_g$ID <- rownames(core_species_id_g)
    core_species_id_g <- core_species_id_g[,c("ID",colnames(core_species_id_g)[1:4])]
    write.table(core_species_id_g,
                "Results/Zooplankton/02.Alpha/02.Venn/Shared_species.txt",
                sep = "\t",row.names = FALSE)
}

#### Differences in total abundance of shared species among different groups
source("Functions/Core.species.abun.g.R")
result <- core.species.abun.g(bac2,core_species_id_g,group)
pdf(file = "Results/Zooplankton/02.Alpha/02.Venn/Shared_species_abundance_diff_boxplot.pdf",
    width = wid*Group_numb,height = 4.2)
result[[1]]
dev.off()
pdf(file = "Results/Zooplankton/02.Alpha/02.Venn/Shared_species_abundance_diff_jitter.pdf",
    width = wid*Group_numb,height = 4.2)
result[[2]]
dev.off()
pdf(file = "Results/Zooplankton/02.Alpha/02.Venn/Shared_species_abundance_diff_violin.pdf",
    width = wid*Group_numb,height = 4.2)
result[[3]]
dev.off()
sink("Results/Zooplankton/02.Alpha/02.Venn/Shared_species_abundance_diff_test.txt")
result[[4]]
sink()

#### Ratio of different phyla in shared species
source("Functions/Core_species_tax.g.R")
tax.number <- core.species.tax.g(core_species_id_g)
write.table(tax.number,
            "Results/Zooplankton/02.Alpha/02.Venn/Shared_species_tax_percentage.txt",
            row.names = FALSE,sep = "\t")
pdf(file = "Results/Zooplankton/02.Alpha/02.Venn/Shared_species_tax_percentage.pdf",
    width = 10,height = 8)
pie(tax.number$Num,labels=tax.number$Phylum,
    col=rainbow(length(tax.number$Phylum)),
    main= paste("Percentage of the phyla in shared",species.type))
dev.off()

### Beta diversity analysis
dir.create("Results/Zooplankton/04.Beta")
dir.create("Results/Zooplankton/04.Beta/01.beta_diversity")
dir.create("Results/Zooplankton/04.Beta/02.Distance")
dir.create("Results/Zooplankton/04.Beta/03.PCA")
dir.create("Results/Zooplankton/04.Beta/04.PCoA")
dir.create("Results/Zooplankton/04.Beta/05.NMDS")
dir.create("Results/Zooplankton/04.Beta/08.Beta_test")

#### PCA 
source("Functions/PCA.R")
result <- pca.community(t(bac2),group,species.type)
write.table(result[[1]],"Results/Zooplankton/04.Beta/03.PCA/PCA_top2.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Zooplankton/04.Beta/03.PCA/PCA_group.pdf",
    width = 7.5,height = 5.4)
result[[2]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/03.PCA/PCA_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/03.PCA/PCA_label.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
write.table(result[[5]],"Results/Zooplankton/04.Beta/03.PCA/PCA_score.txt",sep = "\t")

#### PCoA of Bray_curtis distance
source("Functions/PCoA.R")
result <- pcoa.community(t(bac2),group,species.type)
write.table(as.matrix(result[[1]]),
            "Results/Zooplankton/04.Beta/01.beta_diversity/Bray_Curtis_distance.txt",
            sep = "\t")
write.table(result[[2]],"Results/Zooplankton/04.Beta/04.PCoA/PCoA_top2.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Zooplankton/04.Beta/04.PCoA/PCoA_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/04.PCoA/PCoA_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/04.PCoA/PCoA_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()
write.table(result[[6]],"Results/Zooplankton/04.Beta/04.PCoA/PCoA_score.txt",sep = "\t")

#### NMDS for Bray_Curtis distance
source("Functions/NMDS.R")
result <- nmds.community(t(bac2),group,species.type)
write.table(result[[1]],"Results/Zooplankton/04.Beta/05.NMDS/Microbe_stress.txt",
            sep = "\t",row.names = FALSE)
write.table(result[[2]],"Results/Zooplankton/04.Beta/05.NMDS/NMDS_score.txt",
            sep = "\t",row.names = FALSE)
pdf(file = "Results/Zooplankton/04.Beta/05.NMDS/NMDS_group.pdf",
    width = 7.5,height = 5.4)
result[[3]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/05.NMDS/NMDS_ellipse.pdf",
    width = 7.5,height = 5.4)
result[[4]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/05.NMDS/NMDS_label.pdf",
    width = 7.5,height = 5.4)
result[[5]]
dev.off()

### ANOSIM, Adonis, and MRPP tests of bary-curtis distance
source("Functions/beta.test.R")
dist_bray <- vegdist(t(bac2), method = "bray")
result <- beta.test(dist_bray,group)
sink("Results/Zooplankton/04.Beta/08.Beta_test/Adonis.txt")
result[[1]]
sink()
sink("Results/Zooplankton/04.Beta/08.Beta_test/ANOSIM.txt")
result[[2]]
sink()
sink("Results/Zooplankton/04.Beta/08.Beta_test/MRPP.txt")
result[[3]]
sink()

#### Distance heatmap for Bray_Curtis distance
dist_dd <- as.matrix(dist_bray)
pdf("Results/Zooplankton/04.Beta/02.Distance/Distance_heatmap.pdf",
    width = Sample_numb*0.5 + 2,height = Sample_numb*0.5 + 1)
corrplot(dist_dd ,method = "color",p.mat = dist_dd,insig = "p-value",
         sig.level = -1,hclust.method = "average",number.cex = 0.1,
         col = colorRampPalette(c("Red","white","navy"))(100),
         tl.col = "black")
dev.off()

source("Functions/Diff_beta.R")
result <- diff.beta(dist_dd,group)
pdf(file = "Results/Zooplankton/04.Beta/02.Distance/Group_distance_diff_boxplot.pdf",
    width = Group_numb*0.5 + 0.5,height = 4.2)
result[[1]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/02.Distance/Group_distance_diff_jitter.pdf",
    width = Group_numb*0.5 + 0.5,height = 4.2)
result[[2]]
dev.off()
pdf(file = "Results/Zooplankton/04.Beta/02.Distance/Group_distance_diff_violin.pdf",
    width = Group_numb*0.5 + 0.5,height = 4.2)
result[[3]]
dev.off()
sink("Results/Zooplankton/04.Beta/02.Distance/Group_distance_diff_test.txt")
summary(result[[4]])
sink()

### Different abundance chart
dir.create("Results/Zooplankton/03.Community/03.Differential")
source("Functions/top30.R")
result <- top30(bac3)
source("Functions/diff.heatmap.R")
Diff_heat <- diff.heatmap(result,group)
pdf("Results/Zooplankton/03.Community/03.Differential/Phylum_heatmap.pdf",
    width = 0.5*Group_numb + 1 + max(str_length(bac3$Phylum))*0.1,
    height = ifelse(length(unique(bac3$Phylum)) > 30, 7, 
                    0.22*length(unique(bac3$Phylum)) + 1))
Diff_heat[[1]]
dev.off()
pdf("Results/Zooplankton/03.Community/03.Differential/Genus_heatmap.pdf",
    width = 0.5*Group_numb + 1 + 
        max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1,
    height = ifelse(length(unique(bac3$Genus)) > 30, 7, 
                    0.22*length(unique(bac3$Genus)) + 1))
Diff_heat[[3]]
dev.off()

source("Functions/diff.bubble.R")
Diff_heat <- diff.bubble(result,group)
pdf("Results/Zooplankton/03.Community/03.Differential/Phylum_bubblechart.pdf",
    width = 0.5*Group_numb + 1 + max(str_length(bac3$Phylum))*0.1,
    height = ifelse(length(unique(bac3$Genus)) > 30, 7, 
                    0.22*length(unique(bac3$Genus)) + 1))
Diff_heat[[1]]
dev.off()
pdf("Results/Zooplankton/03.Community/03.Differential/Genus_bubblechart.pdf",
    width = 0.5*Group_numb + 1 + 
        max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1,
    height = ifelse(length(unique(bac3$Genus)) > 30, 7, 
                    0.22*length(unique(bac3$Genus)) + 1))
Diff_heat[[3]]
dev.off()

### Indicator identification
dir.create("Results/Zooplankton/03.Community/04.Indicator")
source("Functions/Indicator.R")
source("Functions/tax.summary.R")
result <- tax.summary(bac3)
Ind.result <- Indicator.indval(result,group)
write.table(Ind.result[[1]],
            "Results/Zooplankton/03.Community/04.Indicator/Phylum_indicator.txt",
            sep = "\t")
write.table(Ind.result[[2]],
            "Results/Zooplankton/03.Community/04.Indicator/Genus_indicator.txt",
            sep = "\t")

### SIMPER analysis
dir.create("Results/Zooplankton/03.Community/05.SIMPER")
source("Functions/simper.R")
simper.result <- simper.bac(result,group)
temp <- names(simper.result[[1]])
for (i in 1:length(temp)) {
    write.table(simper.result[[1]][[i]],
                "Results/Zooplankton/03.Community/05.SIMPER/Phylum_SIMPER.txt",
                sep = "\t")
    write.table(simper.result[[2]][[i]],
                "Results/Zooplankton/03.Community/05.SIMPER/Genus_SIMPER.txt",
                sep = "\t")
}

### Wilcox test
dir.create("Results/Zooplankton/03.Community/06.Wilcox")
source("Functions/wilcox.abun.R")
source("Functions/wilcox.plot.R")
source("Functions/wilcox.R")
for (i in 1:length(temp)) {
    dir.create(paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],sep = ""))
    dd <- c()
    for (k in 1:length(levels(group$Group))) {
        a <- ifelse(grepl(levels(group$Group)[k],temp[i]),levels(group$Group)[k],NA)
        dd <- c(dd,a)
    }
    dd <- dd[complete.cases(dd)]
    group1 <- group[group$Group %in% dd,]
    otu <- result[[1]]
    otu <- otu[,c("Phylum",group1$variable)]
    dir.create(paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],"/Phylum",
                     sep = ""))
    set.seed(1111)
    wilcox.result <- wilcox.biomarker(otu,group1)
    write.table(wilcox.result[[1]],
                paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                      "/Phylum/wilcox_all.txt",sep = ""),
                sep = "\t",row.names = FALSE)
    write.table(wilcox.result[[2]],
                paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                      "/Phylum/wilcox_biomarker.txt",sep = ""),
                sep = "\t",row.names = FALSE)
    if (nrow(wilcox.result[[2]]) > 0) {
        Phylum <- wilcox.abun(otu,wilcox.result[[2]],group1)
        write.table(Phylum[[1]],
                    paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                          "/Phylum/Group_mean_abun_biomarker.txt",sep = ""),
                    sep = "\t",row.names = FALSE)
        write.table(Phylum[[2]],
                    paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                          "/Phylum/Abun_change_biomarker.txt",sep = ""),
                    sep = "\t",row.names = FALSE)
        p <- wilcox.plot(Phylum[[1]],Phylum[[2]])
        pdf(paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                  "/Phylum/wilcox_biomarker_abun.pdf",sep = ""),
            width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
        print(p)
        dev.off()
    }
    otu <- result[[2]]
    otu <- otu[,c("Genus",group1$variable)]
    dir.create(paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],"/Genus",
                     sep = ""))
    set.seed(1111)
    wilcox.result <- wilcox.biomarker(otu,group1)
    write.table(wilcox.result[[1]],
                paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                      "/Genus/wilcox_all.txt",sep = ""),
                sep = "\t",row.names = FALSE)
    write.table(wilcox.result[[2]],
                paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                      "/Genus/wilcox_biomarker.txt",sep = ""),
                sep = "\t",row.names = FALSE)
    if (nrow(wilcox.result[[2]]) > 0) {
        Genus <- wilcox.abun(otu,wilcox.result[[2]],group1)
        write.table(Genus[[1]],
                    paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                          "/Genus/Group_mean_abun_biomarker.txt",sep = ""),
                    sep = "\t",row.names = FALSE)
        write.table(Genus[[2]],
                    paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                          "/Genus/Abun_change_biomarker.txt",sep = ""),
                    sep = "\t",row.names = FALSE)
        p <- wilcox.plot(Genus[[1]],Genus[[2]])
        pdf(paste("Results/Zooplankton/03.Community/06.Wilcox/",temp[i],
                  "/Genus/wilcox_biomarker_abun.pdf",sep = ""),
            width = 10,height = 1.8 + 0.18*nrow(wilcox.result[[2]]))
        print(p)
        dev.off()
    }
}

if (!is.null(ENV)) {
    ## Environmental correlation
    dir.create("Results/Zooplankton/05.Environment")
    dir.create("Results/Zooplankton/05.Environment/01.Diff_ENV")
    dir.create("Results/Zooplankton/05.Environment/02.CCA_RDA")
    dir.create("Results/Zooplankton/05.Environment/04.Mantel")
    dir.create("Results/Zooplankton/05.Environment/05.Diversity_heatmap")
    dir.create("Results/Zooplankton/05.Environment/06.TopAbun_heatmap")
    dir.create("Results/Zooplankton/05.Environment/03.ABT")
    
    #### CCA/RDA
    source("Functions/CCA(RDA).R")
    result <- cca.rda(t(bac2),ENV,group)
    write.csv(result[[1]], file = "Results/Zooplankton/05.Environment/02.CCA_RDA/dca.csv")
    write.csv(result[[2]], file = "Results/Zooplankton/05.Environment/02.CCA_RDA/cca.sample.csv")
    write.csv(result[[3]], file = "Results/Zooplankton/05.Environment/02.CCA_RDA/cca.env.csv")
    write.csv(result[[4]], file = "Results/Zooplankton/05.Environment/02.CCA_RDA/cca.species.csv")
    pdf("Results/Zooplankton/05.Environment/02.CCA_RDA/CCA.pdf",width = 7.5,height = 5.4)
    print(result[[5]])
    dev.off()
    write.csv(result[[6]], file = "Results/Zooplankton/05.Environment/02.CCA_RDA/ccaenvfit.csv")
    write.csv(result[[7]],file="Results/Zooplankton/05.Environment/02.CCA_RDA/rda.sample.csv")
    write.csv(result[[8]],file="Results/Zooplankton/05.Environment/02.CCA_RDA/rda.env.csv")
    write.csv(result[[9]],file="Results/Zooplankton/05.Environment/02.CCA_RDA/rda.species.csv")
    pdf("Results/Zooplankton/05.Environment/02.CCA_RDA/RDA.pdf",width = 7.5,height = 5.4)
    print(result[[10]])
    dev.off()
    write.csv(result[[11]],file="Results/Zooplankton/05.Environment/02.CCA_RDA/rdaenvfit.csv")
    
    ### Correlation of environmental factors with alpha diversity
    source("Functions/heatmap.ENV.R")
    result <- heatmap.env(alpha,ENV)
    write.table(result[[1]]$r,
                "Results/Zooplankton/05.Environment/05.Diversity_heatmap/correlation.xls",
                sep="\t",quote=FALSE,col.names=NA)
    write.table(result[[1]]$p,
                "Results/Zooplankton/05.Environment/05.Diversity_heatmap/pvalue.xls",
                sep="\t",quote=FALSE,col.names=NA)
    pdf("Results/Zooplankton/05.Environment/05.Diversity_heatmap/Diversity_heatmap.pdf",
        width = 2.2 + ENV_numb*0.6,height = 5)
    print(result[[2]])
    dev.off()
    
    ### Correlation of environmental factors with dominant species
    result <- top30(bac3)
    result.abun <- heatmap.env(t(result[[1]]),ENV)
    write.table(result.abun[[1]]$r,
                "Results/Zooplankton/05.Environment/06.TopAbun_heatmap/correlation_phylum.xls",
                sep="\t",quote=FALSE,col.names=NA)
    write.table(result.abun[[1]]$p,
                "Results/Zooplankton/05.Environment/06.TopAbun_heatmap/pvalue_phylum.xls",
                sep="\t",quote=FALSE,col.names=NA)
    pdf("Results/Zooplankton/05.Environment/06.TopAbun_heatmap/heatmap_phylum.pdf",
        width = max(str_length(bac3$Phylum))*0.1 + ENV_numb*0.6 + 0.5,
        height = ifelse(nrow(result[[1]]) == 30, 9,
                        0.25*nrow(result[[1]]) + 1.5))
    print(result.abun[[2]])
    dev.off()
    
    result.abun <- heatmap.env(t(result[[3]]),ENV)
    write.table(result.abun[[1]]$r,
                "Results/Zooplankton/05.Environment/06.TopAbun_heatmap/correlation_genus.xls",
                sep="\t",quote=FALSE,col.names=NA)
    write.table(result.abun[[1]]$p,
                "Results/Zooplankton/05.Environment/06.TopAbun_heatmap/pvalue_genus.xls",
                sep="\t",quote=FALSE,col.names=NA)
    pdf("Results/Zooplankton/05.Environment/06.TopAbun_heatmap/heatmap_genus.pdf",
        width = max(str_length(bac3$Genus[complete.cases(bac3$Genus)]))*0.1 + ENV_numb*0.6 + 0.5,
        height = ifelse(nrow(result[[3]]) == 30, 9,
                        0.25*nrow(result[[3]]) + 1.5))
    print(result.abun[[2]])
    dev.off()
    
    #### ABT
    source("Functions/ABT.R")
    ABT.result <- ABT(t(bac2),ENV)
    write.table(ABT.result[[1]],
                "Results/Zooplankton/05.Environment/03.ABT/ABT_result.txt",
                row.names = FALSE,sep = "\t")
    pdf("Results/Zooplankton/05.Environment/03.ABT/ABT_result.pdf",
        width = 5,height = ENV_numb*0.3 + 0.5)
    print(ABT.result[[2]])
    dev.off()
    
    #### Different of environmental factors among different groups
    source("Functions/Diff_ENV.R")
    result <- diff.ENV(ENV,group)
    write.table(result[[1]],
                "Results/Zooplankton/05.Environment/01.Diff_ENV/ENV_diff_test.txt",
                row.names = FALSE,sep = "\t")
    for (i in 1:ENV_numb) {
        pdf(paste("Results/Zooplankton/05.Environment/01.Diff_ENV/diff_",
                  colnames(ENV)[i],".pdf",sep = ""),
            height = 3,width = 0.6*Group_numb + 0.5)
        print(result[[i+1]])
        dev.off()
    }
    
    #### Mantel test and Procrustes analysis
    source("Functions/Mantel.proc.R")
    result <- Mantel.proc(bac2,ENV)
    pdf("Results/Zooplankton/05.Environment/04.Mantel/Procrustes.pdf",width = 6,height = 5)
    print(result[[1]])
    dev.off()
    sink("Results/Zooplankton/05.Environment/04.Mantel/Mantel.txt")
    result[[2]]
    sink()
    sink("Results/Zooplankton/05.Environment/04.Mantel/Procrustes.txt")
    result[[3]]
    sink()
    
    #### Random forest for sample distinguishing
    dir.create("Results/Zooplankton/03.Community/07.RandomForest")
    result <- tax.summary(bac3)
    source("Functions/RandomForest_classification.R")
    tax.name <- c("Phylum","Genus")
    for (i in 1:length(tax.name)) {
        dir.create(paste("Results/Zooplankton/03.Community/07.RandomForest/",
                         tax.name[i],sep = ""))
        set.seed(1111)
        RF.result <- Random.forest.biomarker(result[[i]],group,ENV)
        pdf(paste("Results/Zooplankton/03.Community/07.RandomForest/",
                  tax.name[i],"/classification.pdf",sep = ""),
            height = 3.5,width = 3.5)
        print(RF.result[[1]])
        dev.off()
        importance.otu <- RF.result[[2]]
        importance.otu$ID <- rownames(importance.otu)
        importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
        write.table(importance.otu,
                    paste("Results/Zooplankton/03.Community/07.RandomForest/",
                          tax.name[i],"/importance.txt",sep = ""),
                    sep = "\t",row.names = FALSE)
        pdf(paste("Results/Zooplankton/03.Community/07.RandomForest/",
                  tax.name[i],"/RandomForest_biomarker.pdf",sep = ""),
            height = nrow(importance.otu)*0.3 + 1.85,width = 9)
        print(RF.result[[3]])
        dev.off()
        write.table(RF.result[[4]],
                    paste("Results/Zooplankton/03.Community/07.RandomForest/",
                          tax.name[i],"/biomarker_abundance.txt",sep = ""),
                    sep = "\t")
        aa <- max(str_length(rownames(RF.result[[4]])))
        h.data <- correlate(t(RF.result[[4]][,-ncol(RF.result[[4]])]), ENV,cor.test = TRUE)
        pdf(paste("Results/Zooplankton/03.Community/07.RandomForest/",
                  tax.name[i],"/RandomForest_biomarker_corr_ENV.pdf",sep = ""),
            height = 5,width = (aa*0.1 + ENV_numb*0.3)*3.5/2.5)
        print(RF.result[[5]])
        dev.off()
    }
    
    ### Relative contribution of Environmental factors and geographic locations
    dir.create("Results/Zooplankton/06.Mechanism")
    
    if (!is.null(GEO)) {
        #### Distance-decay similarity of GEO
        dir.create("Results/Zooplankton/06.Mechanism/01.GEO_regression")
        source("Functions/GEO.regression.R")
        result <- GEO.regression(GEO,bac2,group)
        pdf("Results/Zooplankton/06.Mechanism/01.GEO_regression/GEO_regression.pdf",
            width = 4.5,height = 3)
        print(result[[1]])
        dev.off()
        write.table(result[[2]],
                    "Results/Zooplankton/06.Mechanism/01.GEO_regression/GEO_regression.txt",
                    sep = "\t")
        
        #### VPA
        source("Functions/VPA.env.spa.R")
        dir.create("Results/Zooplankton/06.Mechanism/03.VPA")
        for (i in 1:length(levels(group$Group))) {
            otu <- t(bac2[,group[group$Group == levels(group$Group)[i],"variable"]])
            ENV1 <- ENV[group[group$Group == levels(group$Group)[i],"variable"],]
            GEO1 <- GEO[group[group$Group == levels(group$Group)[i],"variable"],]
            bac_dis <- vegdist(otu)
            
            vpa <- vpa.arg2(otu,ENV1,GEO1)
            pdf(file = paste("Results/Zooplankton/06.Mechanism/03.VPA/vpa_",
                             levels(group$Group)[i],".pdf",sep = ""),
                width = 8,height = 6)
            plot(vpa, bg = 2:3,digits = 2,Xnames = c("ENV","GEO"))
            dev.off()
            sink(paste("Results/Zooplankton/06.Mechanism/03.VPA/vpa_",
                       levels(group$Group)[i],".txt",sep = ""))
            print(vpa)
            sink()
        }
    }
    
    
    library(rdacca.hp)
    library(UpSetVP)  
    source("Functions/rdacca.hp.plot2.R")
    dir.create("Results/Zooplankton/06.Mechanism")
    dir.create("Results/Zooplankton/06.Mechanism/04.rdaccahp")
    if (file.exists("Input/vpa_env.txt")) {
      ENV <- read.table("Input/vpa_env.txt",header = TRUE,sep = "\t")
      rownames(ENV) <- ENV[,1]
      ENV <- ENV[,-1]
      ENV_numb <- ncol(ENV)
      ENV <- ENV[group$variable,]
    }else{
      ENV <- c()
    }
    
    if (file.exists("Input/vpa_hydro.txt")) {
      hydro <- read.table("Input/vpa_hydro.txt",header = TRUE,sep = "\t")
      rownames(hydro) <- hydro[,1]
      hydro <- hydro[,-1]
      hydro_numb <- ncol(hydro)
      hydro <- hydro[group$variable,]
    }else{
      hydro <- c()
    }
    
    if (file.exists("Input/Geo.txt")) {
      GEO <- read.table("Input/Geo.txt",header = TRUE,sep = "\t",row.names = 1)
    }else{
      GEO <- c()
    }
    for (i in 1:length(levels(group$Group))) {
      otu <- t(bac2[,group[group$Group == levels(group$Group)[i],"variable"]])
      ENV1 <- ENV[group[group$Group == levels(group$Group)[i],"variable"],]
      GEO1 <- GEO[group[group$Group == levels(group$Group)[i],"variable"],]
      hydro1 <- hydro[group[group$Group == levels(group$Group)[i],"variable"],]
      
      hp <- rdacca.hp.plot(otu,ENV1,GEO1,hydro1)
      pdf(file = paste("Results/Zooplankton/06.Mechanism/04.rdaccahp/RDA_",
                       levels(group$Group)[i],".pdf",sep = ""),
          width = 8,height = 6)
      barplot_hp(hp, col.fill = 'var', 
                 col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA'))
      dev.off()
      
    }
    
    
    #### Distance-decay similarity of ENV
    dir.create("Results/Zooplankton/06.Mechanism/02.ENV_regression")
    source("Functions/ENV.regression.R")
    result <- ENV.regression(ENV,bac2,group)
    pdf("Results/Zooplankton/06.Mechanism/02.ENV_regression/ENV_regression.pdf",
        width = 4.5,height = 3)
    result[[1]]
    dev.off()
    write.table(result[[2]],
                "Results/Zooplankton/06.Mechanism/02.ENV_regression/ENV_regression.txt",
                sep = "\t")
    
    
    ### Ecological quality prediction
    dir.create("Results/Zooplankton/07.EcoPredict")
    
    ### Prediction environmental factors via microbial communities
    dir.create("Results/Zooplankton/07.EcoPredict/03.SingleENV")
    source("Functions/RandomForest_regression.R")
    result <- tax.summary(bac3)
    for (i in 1:length(levels(group$Group))) {
        dir.create(paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                         levels(group$Group)[i],sep = ""))
        group1 <- group[group$Group == levels(group$Group)[i],]
        otu <- result[[2]][,c("Genus",group1$variable)]
        ENV1 <- ENV[group1$variable,]
        for (j in 1:ENV_numb) {
            set.seed(1111)
            dir.create(paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                             levels(group$Group)[i],"/",colnames(ENV)[j],sep = ""))
            group11 <- data.frame(variable = rownames(ENV1),Group = ENV1[,j])
            RF.result <- Random.forest.regression(otu,group11,colnames(ENV1)[j])
            pdf(paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                      levels(group$Group)[i],"/",colnames(ENV)[j],"/Prediction.pdf",sep = ""),
                height = 5,width = 5)
            print(RF.result[[1]])
            dev.off()
            importance.otu <- RF.result[[2]]
            importance.otu$ID <- rownames(importance.otu)
            importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
            write.table(importance.otu,
                        paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                              levels(group$Group)[i],"/",colnames(ENV)[j],
                              "/importance.txt",sep = ""),
                        sep = "\t",row.names = FALSE)
            pdf(paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                      levels(group$Group)[i],"/",colnames(ENV)[j],
                      "/RandomForest_biomarker.pdf",sep = ""),
                height = nrow(importance.otu)*0.3 + 1.85,width = 9)
            print(RF.result[[3]])
            dev.off()
            write.table(RF.result[[4]],
                        paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                              levels(group$Group)[i],"/",colnames(ENV)[j],
                              "/biomarker_abundance.txt",sep = ""),
                        sep = "\t") 
            write.table(RF.result[[5]],
                        paste("Results/Zooplankton/07.EcoPredict/03.SingleENV/",
                              levels(group$Group)[i],"/",colnames(ENV)[j],
                              "/Predicted_value.txt",sep = ""),
                        sep = "\t") 
        } 
    }
    
    ### Prediction of total environmental factors by microbial communities
    dir.create("Results/Zooplankton/07.EcoPredict/01.TotalENV")
    source("Functions/RandomForest_regression2.R")
    temp_pca <- PCA(ENV,scale.unit = FALSE,graph = FALSE)
    group11 <- data.frame(variable = rownames(ENV),Group = temp_pca$ind$coord[,1])
    RF.result <- Random.forest.regression2(result[[2]],group11,
                                           "PC1 of environmental factors",group)
    pdf("Results/Zooplankton/07.EcoPredict/01.TotalENV/Prediction.pdf",
        height = 3.5,width = 4.5)
    print(RF.result[[1]])
    dev.off()
    importance.otu <- RF.result[[2]]
    importance.otu$ID <- rownames(importance.otu)
    importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
    write.table(importance.otu,
                "Results/Zooplankton/07.EcoPredict/01.TotalENV/importance.txt",
                sep = "\t",row.names = FALSE)
    pdf("Results/Zooplankton/07.EcoPredict/01.TotalENV/RandomForest_biomarker.pdf",
        height = nrow(importance.otu)*0.3 + 1.85,width = 9)
    print(RF.result[[3]])
    dev.off()
    write.table(RF.result[[4]],
                "Results/Zooplankton/07.EcoPredict/01.TotalENV/biomarker_abundance.txt",
                sep = "\t") 
    write.table(RF.result[[5]],
                "Results/Zooplankton/07.EcoPredict/01.TotalENV/Predicted_value.txt",
                sep = "\t") 
    write.table(RF.result[[6]],
                "Results/Zooplankton/07.EcoPredict/01.TotalENV/Predicted_accuracy.txt",
                sep = "\t") 
    
    
    ### Prediction of environmental groups by microbial communities
    dir.create("Results/Zooplankton/07.EcoPredict/02.ENVGroup")
    ### Prediction of environmental groups by microbial communities
    if (file.exists("Input/VPA.txt")) {
        ENV.part <- read.delim("Input/VPA.txt",header = FALSE,row.names = 1)
    }else{
        ENV.part <- c()
    }
    if (!is.null(ENV.part) > 0) {
        for (j in 1:nrow(ENV.part)) {
            dir.create(paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                             rownames(ENV.part)[j],sep = ""))
            tt <- as.character(ENV.part[j,])
            tt <- tt[tt != ""]
            temp <- ENV[,tt]
            temp_pca <- PCA(temp,scale.unit = FALSE,graph = FALSE)
            group11 <- data.frame(variable = rownames(ENV),Group = temp_pca$ind$coord[,1])
            RF.result <- Random.forest.regression2(result[[2]],group11,
                                                   paste("PC1 of",rownames(ENV.part)[j]),
                                                   group)
            pdf(paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                      rownames(ENV.part)[j],"/Prediction.pdf",sep = ""),
                height = 3.5,width = 4.5)
            print(RF.result[[1]])
            dev.off()
            importance.otu <- RF.result[[2]]
            importance.otu$ID <- rownames(importance.otu)
            importance.otu <- importance.otu[,c("ID",colnames(importance.otu)[1:(ncol(importance.otu)-1)])]
            write.table(importance.otu,
                        paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                              rownames(ENV.part)[j],"/importance.txt",sep = ""),
                        sep = "\t",row.names = FALSE)
            pdf(paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                      rownames(ENV.part)[j],"/RandomForest_biomarker.pdf",sep = ""),
                height = nrow(importance.otu)*0.3 + 1.85,width = 9)
            print(RF.result[[3]])
            dev.off()
            write.table(RF.result[[4]],
                        paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                              rownames(ENV.part)[j],"/biomarker_abundance.txt",sep = ""),
                        sep = "\t") 
            write.table(RF.result[[5]],
                        paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                              rownames(ENV.part)[j],"/Predicted_value.txt",sep = ""),
                        sep = "\t") 
            write.table(RF.result[[6]],
                        paste("Results/Zooplankton/07.EcoPredict/02.ENVGroup/",
                              rownames(ENV.part)[j],"/Predicted_accuracy.txt",sep = ""),
                        sep = "\t")
        }
    }
}
