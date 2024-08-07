simper.bac <- function(data,group){
    tax.name <- c("Phylum","Genus")
    result <- list()
    for (i in 1:2) {
        otu <- data[[i]]
        rownames(otu) <- otu[,1]
        otu <- otu[,-1]
        otu <- t(otu)
        aa <- match(rownames(otu),group$variable)
        group <- group[aa,]
        sim <- with(group,simper(otu,Group))
        con <- summary(sim)
        result[[i]] <- con
    }
    return(result)
}