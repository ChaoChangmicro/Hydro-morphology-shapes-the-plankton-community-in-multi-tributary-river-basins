tax.summary <- function(data){
    Phylum <- aggregate(data[,1:(ncol(data)-4)],
                             list(data$Phylum),sum)
    colnames(Phylum)[1] <- "Phylum"
    Genus <- aggregate(data[,1:(ncol(data)-4)],
                        list(data$Genus),sum)
    colnames(Genus)[1] <- "Genus"
    result <- list(Phylum,Genus)
    return(result)
}