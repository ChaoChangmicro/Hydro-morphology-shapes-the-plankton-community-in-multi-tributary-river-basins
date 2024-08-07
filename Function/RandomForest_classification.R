Random.forest.biomarker <- function(otu,group,ENV){
    rownames(otu) <- otu[,1]
    otu2 <- otu[,-1]
    otu2 <- otu2[rownames(otu2) != "Unclassified",]
    data1 <- t(otu2)
    aa <- match(rownames(data1),group$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group$Group)
    colnames(data1) <- c(rownames(otu2),"Group")
    data1$Group <- factor(data1$Group)
    
    phy.rf <- randomForest(data1[,1:(ncol(data1)-1)],data1$Group,ntree = 1000,
                           importance = TRUE,proximity = TRUE)
    
    p1 <- pheatmap(phy.rf$confusion[,1:length(unique(data1$Group))],cluster_rows = FALSE,
             cluster_cols = FALSE,legend = FALSE,display_numbers = TRUE,
             number_format = "%i",fontsize = 10,number_color = "black",
             main = paste("RF classifier (",(1-round(phy.rf$err.rate[1000,1],4))*100,"% correct)",sep = ""),
             fontsize_number = 26)
    
    imp <- as.data.frame(round(importance(phy.rf), 4))
    imp <- imp[order(imp$MeanDecreaseAccuracy,decreasing = TRUE),]
    if (phy.rf$mtry > 10) {
        imp1 <- imp[1:10,]
    }else{
        imp1 <- imp[1:phy.rf$mtry,]
    }
    if (nrow(imp1) > 0) {
        imp1$ID <- rownames(imp1)
        
        p2<- imp1 %>%
            mutate(ID = fct_reorder(ID,MeanDecreaseAccuracy)) %>%
            ggplot(aes(ID,MeanDecreaseAccuracy)) +
            geom_errorbar(aes(ymin = MeanDecreaseAccuracy - MeanDecreaseGini,
                              ymax = MeanDecreaseAccuracy + MeanDecreaseGini),width = 0.3) + 
            geom_vline(xintercept = 1:nrow(imp),linetype = 3) +
            geom_point(shape = 21,fill = "grey80",size = 5) +
            coord_flip() +
            theme_bw()+ theme(panel.grid=element_blank()) + 
            theme(panel.border = element_blank()) +
            theme(panel.background=element_rect(fill='transparent', color='black'),
                  plot.margin = unit(c(3,5,1,1),"mm")) +
            scale_y_continuous(limits = c(ifelse(min(imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini) > 0,
                                                 min((imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini)*0.9),
                                                 min((imp1$MeanDecreaseAccuracy - imp1$MeanDecreaseGini)*1.1)),
                                          max((imp1$MeanDecreaseAccuracy + imp1$MeanDecreaseGini)*1.1)),
                               expand = c(0,0)) +
            theme(axis.text.y = element_text(size = 10,face = "bold.italic",
                                             colour = "black"),
                  axis.text.x = element_text(size = 12,face = "bold",colour = "black"),
                  axis.title.x = element_text(size = 16,face = "bold",
                                              colour = "black")) +
            xlab("") +ylab("Importance")
        
        abundance <- otu2[rownames(otu2) %in% rownames(imp1),]
        abundance$ID <- rownames(abundance)
        abundance1 <- melt(abundance)
        colnames(group)[1] <- "variable"
        abundance1 <- merge(abundance1,group)
        abundance1$ID <- factor(abundance1$ID,levels = rev(imp1$ID))
        p3 <- ggplot(abundance1) +
            geom_tile(aes(variable,ID,fill = value)) +
            theme_classic() +
            theme(axis.ticks = element_blank(),axis.line = element_blank()) +
            xlab("") + ylab("") +
            facet_grid(.~Group,scales = "free",space = "free") + 
            scale_fill_gradient2('Relative abundance (%)',low = 'black',high = 'red') +
            theme(panel.background = element_rect(color = "black"),
                  legend.position = "bottom",
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  strip.text.x = element_text(colour = "black",face = "bold",size = 12))
        p3 <- p2 + p3 + plot_layout(widths = c(1,2))
        
        reg <- calc_relimp(t(abundance[,-ncol(abundance)]), ENV)
        p4 <- correlate(t(abundance[,-ncol(abundance)]), ENV) %>% 
            qcorrplot(extra_mat = list(importance = reg$importance,
                                       pvalue = reg$p),
                      fixed = FALSE) +
            geom_tile(colour = "grey80", size = 0.25,
                      data = function(data) data[data$p < 0.05, , drop = FALSE]) +
            geom_point(aes(size = importance), fill = NA, shape = 21,
                       data = function(data) data[data$pvalue < 0.05, , drop = FALSE]) +
            scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))
        
        p5 <- ggplot(reg$explained, aes(explained, name)) +
            geom_col(fill = "steelblue") +
            scale_y_discrete(limits = rev(names(t(abundance[,-ncol(abundance)])))) +
            theme_bw() +
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title.y = element_blank())
        
        p6 <- p4 + p5 + plot_layout(nrow = 1,widths = c(2.5, 1),guides = "collect")
    }
    
    if (nrow(imp1) > 0) {
        result <- list(p1,imp1,p3,abundance,p6)
    }else{
        result <- list(p1,imp1)
    }
    
    return(result)
}