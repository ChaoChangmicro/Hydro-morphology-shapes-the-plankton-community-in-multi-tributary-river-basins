top30 <- function(abundance){
    result <- list()
    tax.name <- c("Phylum","Family","Genus","Species")
    tax.names <- c("phyla","families","genera","species")
    for (i in 1:4) {
        abundance.1 <- aggregate(abundance[,1:(ncol(abundance)-4)],
                                 list(abundance[,(ncol(abundance)-4+i)]),sum)
        colnames(abundance.1)[1] <- tax.name[i]
        abundance.2 <- as.data.frame(t(abundance.1))
        colnames(abundance.2) <- abundance.2[1,]
        abundance.2 <- abundance.2[-1,]
        abundance.2 <- as.matrix(abundance.2)
        f.abundance <- matrix(as.numeric(abundance.2),nrow = nrow(abundance.2))
        rownames(f.abundance) <- rownames(abundance.2)
        colnames(f.abundance) <- colnames(abundance.2)
        f.abundance <- t(f.abundance)
        sum <- apply(f.abundance,1,sum) 
        f.abundance <- cbind(f.abundance,sum)
        f.abundance <- as.data.frame(f.abundance)
        f.abundance <- f.abundance[order(f.abundance[,"sum"],decreasing = T),]
        f.abundance <- subset(f.abundance, select = -sum)
        f.abundance <- f.abundance[rownames(f.abundance) != "Unclassified",]
        f.abundance.1 <- f.abundance
        if (nrow(f.abundance.1) > 20) {
            f.abundance.1 <- f.abundance.1[1:20,]
        }
        result[[i]] <- f.abundance.1
    }
    return(result)
}