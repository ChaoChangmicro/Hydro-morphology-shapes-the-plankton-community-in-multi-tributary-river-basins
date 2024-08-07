wilcox.biomarker <- function(otu,group){
    otu2 <- otu
    rownames(otu2) <-otu2[,1]
    otu2 <- otu2[,-1]

    n <- c(1:nrow(otu2))
    p.value <- rep(NA,nrow(otu2))
    for(i in n){
        p.value[i] <- wilcox.test(as.numeric(otu2[i,group[group$Group == unique(group$Group)[1],"variable"]]),
                                  as.numeric(otu2[i,group[group$Group == unique(group$Group)[2],"variable"]]))$p.value
    }
    wilcox <- as.data.frame(cbind(rownames(otu2),p.value))
    wilcox$p.value <- as.numeric(wilcox$p.value)
    wilcox[is.na(wilcox)] <- 0
    wilcox.biomarker <- wilcox[wilcox$p.value < 0.05,]
    wilcox.biomarker <- wilcox.biomarker[order(wilcox.biomarker$p.value),]
    result <- list(wilcox,wilcox.biomarker)
    return(result)
}