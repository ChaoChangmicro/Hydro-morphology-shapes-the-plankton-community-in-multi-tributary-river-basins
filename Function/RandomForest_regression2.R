Random.forest.regression2 <- function(otu,group11,target,group){
    rownames(otu) <- otu[,1]
    otu2 <- otu[,-1]
    otu2 <- otu2[rownames(otu2) != "Unclassified",]
    rownames(otu2) <- gsub("\\W+\\(.*","",rownames(otu2))
    rownames(otu2) <- gsub("\\(.*","",rownames(otu2))
    rownames(otu2) <- gsub("\\[","",rownames(otu2))
    rownames(otu2) <- gsub("\\]","",rownames(otu2))
    rownames(otu2) <- gsub(" ","_",rownames(otu2))
    otu2 <- otu2[!grepl("^(\\d+)",rownames(otu2)),]
    data1 <- t(otu2)
    aa <- match(rownames(data1),group11$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group11$Group)
    colnames(data1) <- c(rownames(otu2),"Group")
    colnames(data1) <- gsub("-","\\.",colnames(data1))
    
    phy.rf <- randomForest(Group~.,data = data1,ntree = 1000,
                           importance = TRUE,proximity = TRUE)
    
    plotdata <- data.frame(x = data1$Group,y = phy.rf$predicted)
    plotdata$variable <- rownames(plotdata)
    plotdata <- merge(plotdata,group)
    
    p1 <- ggplot(plotdata,aes(x,y,fill = Group)) +
       # geom_point(aes(color = Group),size = 3,alpha = 0.7) + 
        geom_smooth(aes(color = Group,fill = Group),
                    method = 'lm', formula = y ~ x,se=FALSE) +
        scale_fill_manual(values = cbbPalette) + 
        scale_color_manual(values = cbbPalette) +
        theme_bw() + 
        labs(x = "Actural value",y = "Predicted value", title = target) +
        theme(axis.title.x = element_text(size = 12,face = "bold",colour = "black"),
              axis.title.y = element_text(size = 12,face = "bold",colour = "black"),
              axis.text.x = element_text(size = 10,colour = "black"),
              axis.text.y = element_text(size = 10,colour = "black"),
              legend.title = element_text(size = 12,face = "bold",colour = "black"),
              legend.text = element_text(size = 10,colour = "black"),
              plot.title = element_text(face = "bold",size = 14,colour = "black",hjust = 0.5))
    
    test <- c()
    for (i in 1:length(unique(plotdata$Group))) {
        m <- lm(y~x,plotdata[plotdata$Group == unique(plotdata$Group)[i],])
        a <- summary(m)
        aa <- c(m$coefficients[2],m$coefficients[1],a$adj.r.squared,a$coefficients[2,4])
        test <- rbind(test,aa)
    }
    test <- as.data.frame(test)
    colnames(test) <- c("slope","intercept","r2","pvalue")
    rownames(test) <- unique(plotdata$Group)
    
    imp <- as.data.frame(round(importance(phy.rf), 4))
    imp <- imp[order(imp$`%IncMSE`,decreasing = TRUE),]
    if (phy.rf$mtry > 10) {
        imp1 <- imp[1:10,]
    }else{
        imp1 <- imp[1:phy.rf$mtry,]
    }
    
    if (nrow(imp1) > 0) {
        imp1$ID <- rownames(imp1)
        
        p2<- imp1 %>%
            mutate(ID = fct_reorder(ID,`%IncMSE`)) %>%
            ggplot(aes(ID,`%IncMSE`)) +
            geom_errorbar(aes(ymin = `%IncMSE` - IncNodePurity,
                              ymax = `%IncMSE` + IncNodePurity),width = 0.3) + 
            geom_vline(xintercept = 1:nrow(imp),linetype = 3) +
            geom_point(shape = 21,fill = "grey80",size = 5) +
            coord_flip() +
            theme_bw()+ theme(panel.grid=element_blank()) + 
            theme(panel.border = element_blank()) +
            theme(panel.background=element_rect(fill='transparent', color='black'),
                  plot.margin = unit(c(3,5,1,1),"mm")) +
            scale_y_continuous(limits = c(ifelse(min(imp1$`%IncMSE` - imp1$IncNodePurity) > 0,
                                                 min((imp1$`%IncMSE` - imp1$IncNodePurity)*0.9),
                                                 min((imp1$`%IncMSE` - imp1$IncNodePurity)*1.1)),
                                          max((imp1$`%IncMSE` + imp1$IncNodePurity)*1.1)),
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
            scale_fill_gradient2("Relative abundance (%)",low = "black",high = "red") +
            theme(panel.background = element_rect(color = "black"),
                  legend.position = "bottom",
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  strip.text.x = element_text(colour = "black",face = "bold",size = 12))
        p3 <- p2 + p3 + plot_layout(widths = c(1,2))
    }
    
    if (nrow(imp1) > 0) {
        result <- list(p1,imp1,p3,abundance,plotdata,test)
    }else{
        result <- list(p1,imp1,test)
    }
    
    return(result)
}