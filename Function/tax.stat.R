tax.stat <- function(bac,bac3){
    tax.stat <- c()
    for (i in 1:(ncol(bac)-1)) {
        dd1 <- c()
        for (j in 1:4) {
            dd <- ifelse(bac3[,i] > 0,bac3[,j + Sample_numb],NA)
            dd <- gsub("Unclassified",NA,dd)
            dd <- dd[complete.cases(dd)]
            dd1 <- c(dd1,length(unique(dd)))
        }
        tax.stat <- cbind(tax.stat,dd1)
    }
    tax.stat <- as.data.frame(tax.stat)
    rownames(tax.stat) <- c("Phylum","Family","Genus","Species")
    colnames(tax.stat) <- colnames(bac)[1:(ncol(bac)-1)]
    dd1 <- c()
    for (j in 1:4) {
        dd <- gsub("Unclassified",NA,bac3[,j + Sample_numb])
        dd <- dd[complete.cases(dd)]
        dd1 <- c(dd1,length(unique(dd)))
    }
    tax.stat$Total <- dd1
    tax.stat$Taxonomy <- rownames(tax.stat)
    tax.stat <- tax.stat[,c("Taxonomy",colnames(tax.stat)[1:(ncol(tax.stat)-1)])]
    
    dd1 <- c()
    for (j in 1:4) {
        dd <- gsub("Unclassified",NA,bac3[,j + Sample_numb])
        dd <- dd[complete.cases(dd)]
        dd1 <- c(dd1,length(dd)/nrow(bac3))
    }
    tax.per <- data.frame(Taxonomy = c("Phylum","Family","Genus","Species"),
                          Per = dd1)
    tax.per$Taxonomy <- factor(tax.per$Taxonomy,levels = c("Phylum","Family","Genus","Species"))
    tax.per$lab <- paste(tax.per$Taxonomy,"-",round(tax.per$Per*100,2),"%",sep = "")
    p <- ggplot(tax.per,aes(Taxonomy,Per,fill = Taxonomy)) +
        geom_bar(stat = "identity",width = 0.9,show.legend = FALSE) +
        coord_polar(theta = "y",start = 0) + 
        geom_text(aes(Taxonomy,y = 0,label = lab),hjust = 1.05,size = 2) + 
        ylim(0,2) + 
        labs(x = NULL,y = NULL) + 
        scale_fill_manual(values = cbbPalette) +
        theme_bw()+
        theme(panel.grid=element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(), 
              axis.title=element_blank(),
              axis.text=element_blank())
    result <- list(tax.stat,p)
    return(result)
}