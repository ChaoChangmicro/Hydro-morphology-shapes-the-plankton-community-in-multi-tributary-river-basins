rdacca.hp.plot <- function(ARG_beta,ENV,GEO,hydro){
  
    ARG_hel <- decostand(ARG_beta, method = 'hellinger')
    
    aa <- match(rownames(ARG_beta),rownames(ENV))
    ENV <- ENV[aa,]
    
    
	bb <- match(rownames(ARG_beta),rownames(hydro))
    hydro <- hydro[bb,]
    
	
	
    rownames(GEO) <- rownames(ENV)
    aa <- match(rownames(ARG_beta),rownames(GEO))
    GEO <- GEO[aa,]
    distsite <- dist(GEO[,1:2])
    PCNM.site <- pcnm(distsite)
    
	env  <- cbind(hydro,ENV)
	exp.list <- list(ENV = as.data.frame(env), GEO = as.data.frame(PCNM.site$vectors))
	
	
	#mod1 <- rdacca.hp(ARG_hel, exp.list, method = 'CCA', var.part = TRUE, type = 'R2', scale = FALSE)
	#PP1<- barplot_hp(mod1, col.fill = 'var', 
    #col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA'))
	
	
	mod2 <- rdacca.hp(ARG_hel, exp.list, method = 'RDA', var.part = TRUE, type = 'R2', scale = FALSE)
	PP2<- barplot_hp(mod2, col.fill = 'var', 
    col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA'))
	
	
   # mod3 <- rdacca.hp(ARG_hel, exp.list, method = 'CCA', var.part = TRUE, type = 'adjR2', scale = FALSE)
	#PP3<- barplot_hp(mod3, col.fill = 'var', 
   # col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA'))
	
	
	#mod4 <- rdacca.hp(ARG_hel, exp.list, method = 'RDA', var.part = TRUE, type = 'adjR2', scale = FALSE)
	#PP4<- barplot_hp(mod4, col.fill = 'var', 
    #col.color = c('#8DD3C7', '#FFFFB3', '#BEBADA'))
	
	
	
    #return(PP1,PP2,PP3,PP4)
	return(PP2)
	
}