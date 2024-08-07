Random.forest.regression <- function(otu,group,target){
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
    aa <- match(rownames(data1),group$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group$Group)
    colnames(data1) <- c(rownames(otu2),"Group")
    colnames(data1) <- gsub("-","\\.",colnames(data1))
    
    phy.rf <- randomForest(Group~.,data = data1,ntree = 1000,
                           importance = TRUE,proximity = TRUE)
    
    lm_eqn <- function(df){
        m <- lm(y ~ x, df);
        eq <- substitute(italic(r)^2~"="~r2*","~italic(p)~"="~p1,
                         list(r2 = format(summary(m)$r.squared, digits = 3),
                              p1 = format(summary(m)$coefficients[8], digits = 3)))
        as.character(as.expression(eq));
    }
    
    plotdata <- data.frame(x = data1$Group,y = phy.rf$predicted)
    
    m <- lm(y~x,plotdata)
    tem <- summary(m)$r.squared
    
    p1 <- ggplot(plotdata,aes(x,y)) +
        geom_point(shape = 21,color = "black",size = 6) + 
        geom_smooth(method = 'lm', formula = y ~ x) +
        geom_text(x= min(plotdata$x), y = max(plotdata$y),
                  label = lm_eqn(plotdata),parse = TRUE,size = 5,hjust = 0) +
        scale_x_continuous(limits = c(min(plotdata)*0.95,max(plotdata)*1.05)) +
        scale_y_continuous(limits = c(min(plotdata)*0.95,max(plotdata)*1.05)) +
        geom_abline(slope = 1,intercept = 0,linetype = "dashed")+
        theme_bw() + 
        labs(x = "Actural value",y = "Predicted value", title = target) +
        theme(axis.title.x = element_text(size = 16,colour = "black",face = "bold"),
              axis.title.y = element_text(size = 16,colour = "black",face = "bold"),
              axis.text.x = element_text(size = 12,colour = "black"),
              axis.text.y = element_text(size = 12,colour = "black"),
              plot.title = element_text(face = "bold",size = 16,colour = "black",hjust = 0.5))
    
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
        result <- list(p1,imp1,p3,abundance,plotdata,tem)
    }else{
        result <- list(p1,imp1,tem)
    }
    
    return(result)
}