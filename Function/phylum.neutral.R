phy.neurtal <- function(dd,bacnlsALL,spp){
    dd3 <- dd %>%
        group_by(Phylum,Group) %>%
        summarise(Num = length(ID))
    dd4 <- dd3 %>% group_by(Phylum) %>% summarise(Num = Num/sum(Num)*100)
    dd4$Group <- dd3$Group
    dd4$Group <- factor(dd4$Group,levels = c("Above","Neutral","Below"))
    
    p1 <- ggplot(data = dd4,aes(x = Phylum, y = Num,fill = Group)) + 
        geom_bar(stat = "identity",width = 0.6) +
        ylab(label = "Proportion (%)") + 
        xlab(label = "") + 
        scale_fill_manual(values = c("#29A6A6","black","#A52A2A")) +
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=8,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 8)) + 
        theme(legend.title = element_blank(),
              legend.position = "top")
    
    dd5 <- as.data.frame(colSums(spp))
    dd5$ID <- rownames(dd5)
    colnames(dd5) <- c("Num","ID")
    dd6 <- merge(dd,dd5)
    
    dd7 <- dd6 %>%
        group_by(Phylum,Group) %>%
        summarise(Num = sum(Num))
    dd8 <- dd7 %>% group_by(Phylum) %>% summarise(Num = Num/sum(Num)*100)
    dd8$Group <- dd7$Group
    dd8$Group <- factor(dd8$Group,levels = c("Above","Neutral","Below"))
    
    p2 <- ggplot(data = dd8,aes(x = Phylum, y = Num,fill = Group)) + 
        geom_bar(stat = "identity",width = 0.6) +
        ylab(label = "Proportion (%)") + 
        xlab(label = "") + 
        scale_fill_manual(values = c("#29A6A6","black","#A52A2A")) +
        theme_bw()+ 
        theme(panel.grid=element_blank()) + 
        theme(panel.border = element_blank()) +
        theme(panel.background=element_rect(fill='transparent', color='black'),
              plot.margin = unit(c(3,5,1,1),"mm")) + 
        theme(axis.text.x=element_text(colour="black",size=8,face = "bold",
                                       angle = 90,vjust = 0.5,hjust = 1)) + 
        theme(axis.text.y=element_text(colour = "black",size = 10)) + 
        theme(axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"))+ 
        theme(axis.title.y = element_text(size = 12,face = "bold",
                                          margin = unit(c(0,1,0,1),"lines"))) + 
        scale_y_continuous(limits = c(0,100.001),expand = c(0,0)) + 
        theme(legend.text = element_text(colour = "black",size = 8)) + 
        theme(legend.title = element_blank(),
              legend.position = "top")
    result <- list(p1,p2,dd)
    return(result)
}