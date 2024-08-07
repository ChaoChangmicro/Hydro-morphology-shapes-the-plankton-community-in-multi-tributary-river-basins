wilcox.abun <- function(otu,otu.wilcox.biomarker,group){
    colnames(group) <- c("variable","Group")
    data1 <- t(otu[,2:ncol(otu)])
    rownames(data1) <- colnames(otu)[2:ncol(otu)]
    aa <- match(rownames(data1),group$variable)
    group <- group[aa,]
    data1 <- data.frame(data1,Group = group$Group)
    colnames(data1) <- c(otu[,1],"Group")
    data1$Group <- factor(data1$Group)
    abun.bar <- data1[,c(otu.wilcox.biomarker$V1,"Group")] %>% 
        gather(variable,value,-Group) %>% 
        group_by(variable,Group) %>% 
        summarise(Mean = mean(value))
    
    diff1 <- data1[,c(otu.wilcox.biomarker$V1,"Group")] %>% 
        select_if(is.numeric) %>%
        map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')

    diff.mean <- diff1[,c("var","estimate","conf.low","conf.high")]
    aa <- match(diff.mean$var,otu.wilcox.biomarker$V1)
    otu.wilcox.biomarker <- otu.wilcox.biomarker[aa,]
    diff.mean <- data.frame(diff.mean,p.value = otu.wilcox.biomarker$p.value)
    diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                                levels(data1$Group)[2]))
    diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]
    result <- list(abun.bar,diff.mean)
    return(result)
}