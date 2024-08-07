rare.spec <- function(bac1,group){
    rare <- list()
    for (i in 1:length(levels(group$Group))) {
        tem <- bac1[rownames(bac1) %in% group[group$Group == levels(group$Group)[i],"variable"],]
        rare[[levels(group$Group)[i]]] <- floor(colMeans(tem))
    }
    rare <- iNEXT(rare,q = 0,datatype = "abundance")
    p1 <- ggiNEXT(rare,type = 1) + 
        scale_color_manual(values = cbbPalette) + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "right")
    
    f <- function(x)sum(x==0)
    rare <- list()
    for (i in 1:length(levels(group$Group))) {
        tem <- bac1[rownames(bac1) %in% group[group$Group == levels(group$Group)[i],"variable"],]
        rare[[levels(group$Group)[i]]] <- apply(tem,2,f)
    }
    rare <- iNEXT(rare,q = 0,datatype = "abundance")
    p2 <- ggiNEXT(rare,type = 1) + 
        scale_color_manual(values = cbbPalette) + 
        scale_x_continuous(breaks = c(0,500,1000,1500,2000),labels = c(0,10,20,30,40)) + 
        xlab("No. sampling sites") + 
        theme_bw() +
        theme(panel.grid = element_blank(),
              legend.position = "right",
              axis.text = element_text(colour = "black"))
    result <- list(p1,p2)
    return(result)
}