phy.singal <- function(phy2,zoo2,group,ENV,tree.phy,tree.zoo){
    K <- c()
    P <- c()
    for (j in 1: length(levels(group$Group))) {
        otu <- t(phy2)
        otu <- otu[group[group$Group == levels(group$Group)[i],"variable"],]
        f <- function(x) sum(x==0)
        d <- apply(otu,2,f)
        otu <- otu[,d < (nrow(otu)-3)]
        ENV1 <- ENV[group[group$Group == levels(group$Group)[i],"variable"],]
        res.s <- fast_correlate(otu,ENV1,method = "spearman")
        tree.phy$tip.label <- gsub("\'","",tree.phy$tip.label)
        tree1 <- picante::prune.sample(otu,tree.phy)
        tree1$node.label = NULL
        tree1$edge.length <- tree1$edge.length+(10^-7)
        Ps.s <- multiPhylosignal(res.s[["r"]],multi2di(tree1))
        K <- c(K,Ps.s$K)
        P <- c(P,Ps.s$PIC.variance.P)
    }
    
    Group <- c()
    for (i in 1: length(levels(group$Group))) {
        aa <- rep(levels(group$Group)[i],ENV_numb)
        Group <- c(Group,aa)
    }
    
    Ps <- data.frame(Env = rep(colnames(ENV),length(levels(group$Group))),
                     K = K,P = P,Group = Group,
                     Type = rep("Phytoplankton",length(Group)))
    K <- c()
    P <- c()
    for (j in 1: length(levels(group$Group))) {
        otu <- t(zoo2)
        otu <- otu[group[group$Group == levels(group$Group)[i],"variable"],]
        f <- function(x) sum(x==0)
        d <- apply(otu,2,f)
        otu <- otu[,d < (nrow(otu)-3)]
        ENV1 <- ENV[group[group$Group == levels(group$Group)[i],"variable"],]
        res.s <- fast_correlate(otu,ENV1,method = "spearman")
        tree.zoo$tip.label <- gsub("\'","",tree.zoo$tip.label)
        tree1 <- picante::prune.sample(otu,tree.zoo)
        tree1$node.label = NULL
        tree1$edge.length <- tree1$edge.length+(10^-7)
        Ps.s <- multiPhylosignal(res.s[["r"]],multi2di(tree1))
        K <- c(K,Ps.s$K)
        P <- c(P,Ps.s$PIC.variance.P)
    }
    
    Ps1 <- data.frame(Env = rep(colnames(ENV),length(levels(group$Group))),
                     K = K,P = P,Group = Group,
                     Type = rep("Zooplankton",length(Group)))
    Ps <- as.data.frame(rbind(Ps,Ps1))
    
    Ps$Env <- factor(Ps$Env,levels = colnames(ENV))
    p <- ggplot(Ps,aes(Env,K)) + 
        geom_bar(aes(fill = Group),stat = "identity",position = "dodge") +
        scale_fill_manual(values = cbbPalette,guide = guide_legend(nrow = 1)) +
        ylab("Phylogenetic signal (Blomberg's K)") +
        scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
        facet_wrap(Type~.,scales = "free",ncol = 1) + 
        theme_bw()+
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"), 
              axis.title.x=element_blank(),
              axis.title.y=element_text(colour='black', size=18,face = "bold",
                                        vjust = 1.5),
              axis.text.y=element_text(colour='black',size=13,face = "bold"),
              axis.text.x=element_text(colour = "black",size = 14,face = "bold"),
              strip.text = element_text(colour = "black",size = 18,face = "bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=14,face = "bold"),
              legend.key=element_blank(),legend.position = "top",
              legend.background = element_blank())
    result <- list(Ps,p)
    return(result)
}