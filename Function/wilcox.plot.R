wilcox.plot <- function(abun.bar,diff.mean){
    abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))
    p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
        scale_x_discrete(limits = levels(diff.mean$var)) +
        coord_flip() +
        xlab("") +
        ylab("Mean proportion") +
        theme(panel.background = element_rect(fill = 'transparent'),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=12,face = "bold"),
              axis.text=element_text(colour='black',size=10,face = "bold"),
              legend.title=element_blank(),
              legend.text=element_text(size=12,face = "bold",colour = "black"),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.key.width = unit(0.8,"cm"),
              legend.key.height = unit(0.5,"cm"))
    
    
    for (k in 1:(nrow(diff.mean) - 1)){ 
        p1 <- p1 + annotate('rect', xmin = k+0.5, xmax = k+1.5, ymin = -Inf, ymax = Inf, 
                            fill = ifelse(k %% 2 == 0, 'white', 'gray95'))
    }
    p1 <- p1 + 
        geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
        scale_fill_manual(values=cbbPalette)
    
    
    diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
    diff.mean$p.value <- as.character(diff.mean$p.value)
    p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
        theme(panel.background = element_rect(fill = 'transparent'),
              panel.grid = element_blank(),
              axis.ticks.length = unit(0.4,"lines"), 
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_text(colour='black', size=12,face = "bold"),
              axis.text=element_text(colour='black',size=10,face = "bold"),
              axis.text.y = element_blank(),
              legend.position = "none",
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
        scale_x_discrete(limits = levels(diff.mean$var)) +
        coord_flip() +
        xlab("") +
        ylab("Difference in mean proportions") +
        labs(title="95% confidence intervals") 
    
    for (k in 1:(nrow(diff.mean) - 1)){ 
        p2 <- p2 + annotate('rect', xmin = k+0.5, xmax = k+1.5, ymin = -Inf, ymax = Inf, 
                            fill = ifelse(k %% 2 == 0, 'white', 'gray95'))
    }
    p2 <- p2 +
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                      position = position_dodge(0.8), width = 0.5, size = 0.5) +
        geom_point(shape = 21,size = 3) +
        scale_fill_manual(values = if(length(unique(diff.mean$Group)) == 1){
            if (unique(diff.mean$Group) == levels(group$Group)[2]) {
                cbbPalette[2]
            }else{
                cbbPalette[1] 
            }
        }else{
            c(cbbPalette[1],cbbPalette[2])
        }) +
        geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')
    
    diff.mean$p.value <- as.numeric(diff.mean$p.value)
    diff.mean$p.value <- ifelse(diff.mean$p.value >= 0.0001,round(diff.mean$p.value,6),
                                formatC(diff.mean$p.value,format = "e",digits = 0))
    p3 <- ggplot(diff.mean,aes(var,estimate)) +
        geom_text(aes(y = 0,x = var),label = substr(diff.mean$p.value,1,6),
                  hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
        geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value",
                  srt = 90,fontface = "bold",size = 3) +
        coord_flip() +
        ylim(c(0,1)) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank())
    
    ## 图像拼接
    p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
    return(p)
}