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
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)
library(linkET)
library(TITAN2)
library(spaa)
library(ggraph)
library(tidygraph)
library(igraph)
library(ggrepel)
library(multcomp)

### Loading Parameters  & Processing
source("Functions/Loading.R")
cbbPalette <- Loading_Color()
data <- Loading_Table("Input/table_f_even_depth_w_tax.xls")
group <- Loading_Group("Input/group.txt")
group$Group<- factor(group$Group, 
                     levels = c("Spring", "Autumn"))
tree <- Loading_Tree("Input/unrooted_tree.tre")



result <- Data_Process.phyto(data,group,tree)
phy1 <- result[[1]]
phy2 <- result[[2]]
phy3 <- result[[3]]
tree.phy <- result[[4]]
Group_numb <- result[[5]]
Sample_numb <- result[[6]]
wid <- result[[7]]
#phy <- result[[8]]

result <- Data_Process.zoo(data,group,tree)
zoo1 <- result[[1]]
zoo2 <- result[[2]]
zoo3 <- result[[3]]
tree.zoo <- result[[4]]
#zoo <- result[[8]]

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

dir.create("Results/")
dir.create("Results/Phy_Zoo_association")
species.type <- c("Phytoplankton","Zooplankton")

## Neutral model
dir.create("Results/Phy_Zoo_association/1.NeutralModel")
source("Functions/netrual.mod.R")
source("Functions/netrual.plot.R")
source("Functions/per.rich.R")
source("Functions/per.abun.R")
source("Functions/phylum.neutral.R")
phy <- read.csv("Input/phy_ncm.csv",row.names = 1,header = T,check.names = FALSE)
zoo <- read.csv("Input/zoo_ncm.csv",row.names = 1,header = T,check.names = FALSE)

for (j in 1:2) {
    dir.create(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                     species.type[j],sep = ""))
    for (i in 1:length(levels(group$Group))) {
        dir.create(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                         species.type[j],"/",levels(group$Group)[i],sep = ""))
        spp <- if (j == 1) {
            phy
        }else{
            zoo
        }
        spp <- spp[,group$variable[group$Group == levels(group$Group)[i]]]
        result <- neutral.mod(spp,group)
        pdf(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                  species.type[j],"/",levels(group$Group)[i],"/NCM.pdf",sep = ""),
            width = 6,height = 5)
        netrual.plot(result[[1]],result[[2]],result[[4]],result[[5]])
        dev.off()
        pdf(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                        species.type[j],"/",levels(group$Group)[i],
                  "/percentage_richness.pdf",sep = ""),width = 4,height = 4)
        per.rich(result[[2]])
        dev.off()
        pdf(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                  species.type[j],"/",levels(group$Group)[i],
                  "/percentage_abun.pdf",sep = ""),width = 4,height = 4)
        per.abun(result[[1]],result[[3]])
        dev.off()
        
        dd <- if (j == 1) {
            phy
        }else{
            zoo
        }
        
        result[[1]]$Group <- ifelse(result[[1]]$freq >= result[[1]]$Upper,"Above",
                                    ifelse(result[[1]]$freq <= result[[1]]$Lower,"Below","Neutral"))
        
        dd <- data.frame(ID = rownames(dd),taxonomy = dd$taxonomy)
        dd$Phylum <- gsub(".*p__","",dd$taxonomy)
        dd$Phylum <- gsub("\\; f__.*","",dd$Phylum)
        dd2 <- data.frame(ID = rownames(result[[1]]),Group = result[[1]]$Group)
        dd <- merge(dd,dd2)
        
        dd3 <- dd %>%
            group_by(Phylum,Group) %>%
            summarise(Num = length(ID))
        dd4 <- dd3 %>% group_by(Phylum) %>% summarise(Num = Num/sum(Num)*100)
        dd4$Group <- dd3$Group
        dd4$Group <- factor(dd4$Group,levels = c("Above","Neutral","Below"))
        
        dd5 <- as.data.frame(rowSums(spp))
        dd5$ID <- rownames(dd5)
        colnames(dd5) <- c("Num","ID")
        dd6 <- merge(dd,dd5)
        
        dd7 <- dd6 %>%
            group_by(Phylum,Group) %>%
            summarise(Num = sum(Num))
        dd8 <- dd7 %>% group_by(Phylum) %>% summarise(Num = Num/sum(Num)*100)
        dd8$Group <- dd7$Group
        dd8$Group <- factor(dd8$Group,levels = c("Above","Neutral","Below"))
        
        phylum <- phy.neurtal(dd,result[[1]],result[[3]])
        pdf(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                  species.type[j],"/",levels(group$Group)[i],"/Phylum_richness.pdf",sep = ""),
            width = 1 + 0.15*length(unique(dd4$Phylum)),height = 3)
        print(phylum[[1]])
        dev.off()
        pdf(paste("Results/Phy_Zoo_association/1.NeutralModel/",
                  species.type[j],"/",levels(group$Group)[i],"/Phylum_Abun.pdf",sep = ""),
            width = 1 + 0.15*length(unique(dd8$Phylum)),height = 3)
        print(phylum[[2]])
        dev.off()
        
        aa <- as.data.frame(t(result[[3]]))
        aa$ID <- rownames(aa)
        aa <- merge(phylum[[3]],aa)
        write.table(aa,paste("Results/Phy_Zoo_association/1.NeutralModel/",
                             species.type[j],"/",levels(group$Group)[i],"/NCM.txt",sep = ""),
                    row.names = FALSE,sep = "\t")
    }
}

if (length(levels(group$Group)) == 2) {
    tt <- group %>% group_by(Group) %>% summarise(Num = length(variable))
    if (min(tt$Num) >= 12) {
        if (!is.null(ENV)) {
            ### Environmental breadth
            dir.create("Results/Phy_Zoo_association/2.EnvBreadth")
            source("Functions/Env.breadth.R")
            result <- ENV.breadth(phy2,zoo2,group,ENV)
            write.table(result[[1]],file = "Results/Phy_Zoo_association/2.EnvBreadth/EB.txt",
                        sep = "\t",row.names = FALSE)
            pdf("Results/Phy_Zoo_association/2.EnvBreadth/EnvBreadth.pdf",
                width = 0.8 + ENV_numb*0.6,height = 8)
            result[[2]]
            dev.off()
            
            ### Phylogenetic single
            source("Functions/phy.singal.R")
            result <- phy.singal(phy2,zoo2,group,ENV,tree.phy,tree.zoo)
            write.table(result[[1]],"Results/Phy_Zoo_association/2.EnvBreadth/PhyloSignal.txt",
                        row.names = FALSE,sep = "\t")
            pdf("Results/Phy_Zoo_association/2.EnvBreadth/PhyloSignal.pdf",
                width = 0.8 + ENV_numb*0.6,height = 8)
            result[[2]]
            dev.off()
        }
    }
    
    ### Habitat niche
    dir.create("Results/Phy_Zoo_association/3.Niche")
    source("Functions/habitat.niche.R")
    result <- habitat.niche(phy2,zoo2,group)
    write.table(result[[1]],"Results/Phy_Zoo_association/3.Niche/Habitat_niche.txt",
                row.names = FALSE,sep = "\t")
    pdf("Results/Phy_Zoo_association/3.Niche/Habitat_niche.pdf",width = 2,height = 3.5)
    result[[2]]
    dev.off()
    
    ### Dispersal ability
    source("Functions/Dis.ability.R")
    result <- dis.ability(phy2,zoo2,group)
    write.table(result[[1]],"Results/Phy_Zoo_association/3.Niche/Dispersal_ability.txt",
                row.names = FALSE,sep = "\t")
    pdf("Results/Phy_Zoo_association/3.Niche/Dispersal_ability.pdf",width = 2,height = 3.5)
    result[[2]]
    dev.off()
}

### Network analysis
dir.create("Results/Phy_Zoo_association/4.Network")
source("Functions/zi_pi.r")
source("Functions/network11.R")
for (i in 1:length(levels(group$Group))) {
    dir.create(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],sep = ""))
    phy.1 <- phy2[,group[group$Group == levels(group$Group)[i],"variable"]]
    zoo.1 <- zoo2[,group[group$Group == levels(group$Group)[i],"variable"]]
    result <- network.11(phy.1,zoo.1,phy3,zoo3)
    pdf(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
              "/network_abun.pdf",sep = ""),
        width = 6,height = 6)
    print(result[[1]])
    dev.off()
    pdf(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
              "/network_abun_nolabel.pdf",sep = ""),
        width = 6,height = 6)
    print(result[[2]])
    dev.off()
    pdf(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
              "/network_degree.pdf",sep = ""),
        width = 6,height = 6)
    print(result[[3]])
    dev.off()
    pdf(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
              "/network_degree_nolabel.pdf",sep = ""),
        width = 6,height = 6)
    print(result[[4]])
    dev.off()
    write.table(result[[5]],
                paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                    "/edge.csv",sep = ""),sep=",",quote=F,row.names = FALSE)
    write.table(result[[6]], 
                paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                "/nodes_list.txt",sep = ""),sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(result[[7]],
                paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                             "/zi_pi_result.txt",sep = ""), 
                sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(result[[8]],
                paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                      "/net.txt",sep = ""), 
                sep = '\t', row.names = FALSE, quote = FALSE)
    write.table(result[[9]],
                paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                      "/net2.txt",sep = ""), 
                sep = '\t', row.names = FALSE, quote = FALSE)
    pdf(file = paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                     "/zi_pi.pdf",sep = ""),width = 4,height = 3.5)
    print(result[[8]])
    dev.off()
}


###新网络分析
dir.create("Results/Phy_Zoo_association/4.Network")
source("Functions/network_new.R")
for (i in 1:length(levels(group$Group))) {
  dir.create(paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],sep = ""))
  phy.1 <- phy2[,group[group$Group == levels(group$Group)[i],"variable"]]
  zoo.1 <- zoo2[,group[group$Group == levels(group$Group)[i],"variable"]]
  result <- network.new(phy.1,zoo.1,phy3,zoo3)
  
  write.csv(result[[1]],
            paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                  "/network.csv",sep = ""))
  write.csv(result[[2]], 
            paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                  "/edge_list.csv",sep = ""),row.names = FALSE)
  write.csv(result[[3]],
            paste("Results/Phy_Zoo_association/4.Network/",levels(group$Group)[i],
                  "/nodes_list.csv",sep = ""),row.names = FALSE)
}



#### Network stability
dir.create("Results/Phy_Zoo_association/4.Network/Stability")
source("Functions/robustness.R")
source("Functions/info.centrality.R")
source("Functions/cohesion.R")
source("Functions/network.stability.R")
set.seed(111)
group$Group<- factor(group$Group,
                     levels = c("Spring", "Autumn"))
result <- network.stability(phy2,zoo2,group)
write.csv(result[[1]],"Results/Phy_Zoo_association/4.Network/Stability/random_removal_result.csv")
pdf("Results/Phy_Zoo_association/4.Network/Stability/random_removal_results_weighted.pdf",
    width = 5,height = 3)
result[[2]]
dev.off()
pdf("Results/Phy_Zoo_association/4.Network/Stability/random_removal_results_unweighted.pdf",
    width = 5,height = 3)
result[[3]]
dev.off()
write.table(result[[4]],
          "Results/Phy_Zoo_association/4.Network/Stability/robustness.txt",
          sep = "\t",row.names = FALSE)
pdf("Results/Phy_Zoo_association/4.Network/Stability/Robustness.pdf",
    width = 0.7 + Group_numb*0.8,height = 3)
result[[5]]
dev.off()
write.table(result[[6]],
            "Results/Phy_Zoo_association/4.Network/Stability/vulnerability.txt",
            sep = "\t",row.names = FALSE)
pdf("Results/Phy_Zoo_association/4.Network/Stability/vulnerability.pdf",
    0.7 + Group_numb*0.8,height = 3)
result[[7]]
dev.off()
write.table(result[[8]],
            "Results/Phy_Zoo_association/4.Network/Stability/cohesion.txt",
            sep = "\t",row.names = FALSE)
pdf("Results/Phy_Zoo_association/4.Network/Stability/cohesion.pdf",
    width = 0.7 + Group_numb*0.8,height = 3)
result[[9]]
dev.off()











